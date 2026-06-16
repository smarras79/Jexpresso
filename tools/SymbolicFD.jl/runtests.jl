#=---------------------------------------------------------------------------------
# runtests.jl  --  analytic-solution and operator checks for SymbolicFD.jl
#
# Run inside the Jexpresso project environment (so Plots resolves):
#
#       julia --project=. tools/SymbolicFD.jl/runtests.jl
#
# Covers:
#   1. differential operators (∇, ∇⋅, ∇²) vs. exact derivatives + 2nd-order
#      convergence;
#   2. parsing / rank bookkeeping (dimensional and unknown-symbol errors);
#   3. the heat equation  ∂q/∂t = μ∇²q          vs. its analytic decay;
#   4. periodic advection ∂q/∂t + ∇⋅(uq) = 0    returning to itself after one
#      period, with machine-precision mass conservation;
#   5. an end-to-end `solve` smoke test (ascii output, no display).
#---------------------------------------------------------------------------------=#

using Test
include(joinpath(@__DIR__, "src", "SymbolicFD.jl"))
using .SymbolicFD
const S = SymbolicFD

# ---- helpers ---------------------------------------------------------------------
make_mesh(N; xmin = -1.0, xmax = 1.0) =
    S.FDMesh1D(Dict(:npoin => N, :xmin => xmin, :xmax => xmax, :periodic => true))

field(mesh, f) = Field(mesh, 0, [Float64[f(xi) for xi in mesh.x]])

l2(a) = sqrt(sum(abs2, a) / length(a))
rel_l2(num, exact) = l2(num .- exact) / l2(exact)

# integrate an equation WITHOUT any file/plot output, returning (mesh, q0, q)
function integrate_only(eqn, inputs, q0fun; tend, cfl = 0.4, Δt = nothing)
    var, mode, dqdt = parse_equation(eqn, inputs)
    mesh = S.FDMesh1D(inputs)
    rhs! = S.build_rhs(dqdt, mesh)
    q0   = Float64[q0fun(xi) for xi in mesh.x]
    q    = copy(q0)
    dt     = Δt === nothing ? S.stable_dt(inputs, mesh, cfl) : Δt
    nsteps = max(1, ceil(Int, tend / dt)); dt = tend / nsteps
    w = S.RK4Work(mesh.npoin)
    for _ in 1:nsteps
        S.rk4_step!(q, rhs!, dt, w)
    end
    return mesh, q0, q
end

@testset "SymbolicFD" begin

    # ----------------------------------------------------------------------------
    @testset "operators vs. exact derivatives" begin
        m  = make_mesh(256)
        f  = field(m, x -> sinpi(x))                  # sin(πx)
        fx_exact  = [Float64(π) * cospi(x)   for x in m.x]
        fxx_exact = [-Float64(π)^2 * sinpi(x) for x in m.x]

        g = S.gradient(f)
        @test g.rank == 1
        @test rel_l2(g.comp[1], fx_exact) < 1e-3

        l = S.laplacian(f)
        @test l.rank == 0
        @test rel_l2(l.comp[1], fxx_exact) < 1e-3

        # divergence of a vector flux  ∇⋅(u f) = d(u f)/dx  for constant u
        u  = 2.0
        vf = Field(m, 1, [Float64[u * sinpi(x) for x in m.x]])
        d  = S.divergence(vf)
        @test d.rank == 0
        @test rel_l2(d.comp[1], [u * Float64(π) * cospi(x) for x in m.x]) < 1e-3

        # operators annihilate constants exactly
        c = field(m, x -> 3.0)
        @test maximum(abs, S.gradient(c).comp[1])  < 1e-12
        @test maximum(abs, S.laplacian(c).comp[1]) < 1e-12
    end

    # ----------------------------------------------------------------------------
    @testset "2nd-order spatial convergence (∇²)" begin
        err(N) = begin
            m = make_mesh(N)
            num = S.laplacian(field(m, x -> sinpi(x))).comp[1]
            ex  = [-Float64(π)^2 * sinpi(x) for x in m.x]
            l2(num .- ex)
        end
        ratio = err(32) / err(64)
        @test 3.5 < ratio < 4.5          # halving Δx cuts the error by ~4
    end

    # ----------------------------------------------------------------------------
    @testset "parsing & rank bookkeeping" begin
        var, mode, dqdt = parse_equation("∂q/∂t + ∇⋅(uq) = μ∇²q", Dict(:u => [1.0], :μ => 1e-3))
        @test var == "q"
        @test mode == :transient

        # no time derivative -> steady mode, unknown inferred (q not a parameter)
        v2, m2, _ = parse_equation("μ∇²q = f", Dict(:μ => 1.0, :f => x -> 0.0))
        @test v2 == "q"
        @test m2 == :steady

        # a scalar = vector balance must be rejected (∇q is rank 1)
        m = make_mesh(16)
        _, _, bad = parse_equation("∂q/∂t = ∇q", Dict())
        rhs! = S.build_rhs(bad, m)
        @test_throws ErrorException rhs!(zeros(16), ones(16))

        # unknown symbol on a transient eqn must be reported
        @test_throws ErrorException parse_equation("∂q/∂t = a q", Dict())
    end

    # ----------------------------------------------------------------------------
    @testset "heat equation analytic decay" begin
        μ, tend = 0.05, 0.5
        m, q0, q = integrate_only("∂q/∂t = μ∇²q", Dict(:μ => μ, :npoin => 256),
                                  x -> sinpi(x); tend = tend)
        # exact: q(x,t) = exp(-μ π² t) sin(πx)
        exact = [exp(-μ * Float64(π)^2 * tend) * sinpi(x) for x in m.x]
        @test rel_l2(q, exact) < 2e-3
    end

    # ----------------------------------------------------------------------------
    @testset "periodic advection (one period) + mass conservation" begin
        u, tend = 1.0, 2.0          # domain length 2, speed 1 -> one full period
        # nonzero-mean profile so the mass-conservation check is meaningful
        m, q0, q = integrate_only("∂q/∂t + ∇⋅(uq) = 0", Dict(:u => [u], :npoin => 256),
                                   x -> 1.0 + 0.5 * sinpi(x); tend = tend)
        @test rel_l2(q, q0) < 5e-3                       # returns to itself
        Δx = m.Δx
        @test sum(q0) * Δx ≈ 2.0                          # mean 1 over length-2 domain
        @test abs(sum(q) * Δx - sum(q0) * Δx) < 1e-10    # conservative form
    end

    # ----------------------------------------------------------------------------
    @testset "steady (time-independent) problems" begin
        # Laplace  ∇²q = 0,  q(-1)=0, q(1)=1  ->  linear  q = (x+1)/2
        _, mode, node = parse_equation("∇²q = 0", Dict())
        @test mode == :steady
        m  = S.FDMesh1D(Dict(:npoin => 101, :xmin => -1.0, :xmax => 1.0, :periodic => false))
        q, rmax = S.solve_steady(node, m, Dict(:bc_left => 0.0, :bc_right => 1.0))
        @test rel_l2(q, [(xi + 1) / 2 for xi in m.x]) < 1e-10
        @test rmax < 1e-8

        # Poisson  μ∇²q = f,  μ=1, f=-π²sin(πx)  ->  q = sin(πx),  q(±1)=0
        inp = Dict(:μ => 1.0, :f => x -> -Float64(π)^2 * sinpi(x))
        _, mode2, node2 = parse_equation("μ∇²q = f", inp)
        @test mode2 == :steady
        m2 = S.FDMesh1D(Dict(:npoin => 257, :xmin => -1.0, :xmax => 1.0, :periodic => false))
        q2, r2 = S.solve_steady(node2, m2, merge(inp, Dict(:bc_left => 0.0, :bc_right => 0.0)))
        @test rel_l2(q2, [sinpi(xi) for xi in m2.x]) < 2e-3
        @test r2 < 1e-8

        # a steady solve refuses a periodic grid (no Dirichlet data)
        mp = S.FDMesh1D(Dict(:npoin => 32, :periodic => true))
        @test_throws ErrorException S.solve_steady(node, mp, Dict())
    end

    # ----------------------------------------------------------------------------
    @testset "end-to-end solve (ascii, no display)" begin
        dir = mktempdir()
        inputs = Dict(:nsd => 1, :xmin => -1.0, :xmax => 1.0, :npoin => 256,
                      :periodic => true, :u => [1.0], :μ => 1e-3,
                      :q0 => x -> exp(-(x^2) / (2 * 0.1^2)),
                      :tend => 0.5, :CFL => 0.4,
                      :outformat => "ascii", :plot_live => false,
                      :output_dir => dir)
        mesh, q0, q = SymbolicFD.solve("∂q/∂t + ∇⋅(\\mathbf{u}q) = \\mu∇⋅∇(q)", inputs)
        @test length(q) == 256
        @test all(isfinite, q)
        @test isfile(joinpath(dir, "solution.csv"))
        # conservative advection + diffusion both conserve the discrete mass
        @test abs(sum(q) - sum(q0)) * mesh.Δx < 1e-6
        @test maximum(q) ≤ maximum(q0) + 1e-3            # diffusion ⇒ peak decays
    end

end
