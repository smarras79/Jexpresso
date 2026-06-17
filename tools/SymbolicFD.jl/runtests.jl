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
using LinearOperators
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
    @testset "discretization backend (:method)" begin
        # default backend is finite differences
        @test make_mesh(16).disc isa S.FDMethod
        # unknown method is rejected
        @test_throws ErrorException S.FDMesh1D(Dict(:method => :nope))

        # SEM backend reuses Jexpresso's basis (dψ, ξ, ω); only run the numerical
        # check if the Jexpresso package can be loaded in this environment.
        local msem
        sem_ok = true
        try
            msem = S.FDMesh1D(Dict(:method => :sem, :nop => 4, :nelx => 8,
                                   :xmin => -1.0, :xmax => 1.0, :periodic => true))
        catch err
            sem_ok = false
            @warn "SEM backend test skipped (could not load Jexpresso basis here)" err
        end
        if sem_ok
            @test msem.disc isa S.SEMMethod
            @test msem.sem !== nothing
            # spectral first derivative of sin(πx) on the LGL element mesh
            f  = Float64[sinpi(x) for x in msem.x]
            d  = S.deriv1(f, msem, 1)
            ex = [Float64(π) * cospi(x) for x in msem.x]
            @test rel_l2(d, ex) < 1e-3
            # spectral Laplacian (∇⋅∇) via the operator layer
            lap = S.laplacian(Field(msem, 0, f)).comp[1]
            ex2 = [-Float64(π)^2 * sinpi(x) for x in msem.x]
            @test rel_l2(lap, ex2) < 1e-2
        else
            @test_skip sem_ok
        end
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
    @testset "2D structured operators (FD + tensor SEM)" begin
        mesh2(N; method = :fd) = S.FDMesh2D(Dict(:nsd => 2, :method => method,
            :xmin => -1.0, :xmax => 1.0, :ymin => -1.0, :ymax => 1.0,
            :npoinx => N, :npoiny => N, :nelx => N, :nely => N, :nop => 4,
            :periodic => true))
        on(m, f) = Float64[f(m.x[ip], m.y[ip]) for ip in 1:m.npoin]
        exact_x(m)   = on(m, (x, y) -> Float64(π) * cospi(x) * sinpi(y))
        exact_y(m)   = on(m, (x, y) -> Float64(π) * sinpi(x) * cospi(y))
        exact_lap(m) = on(m, (x, y) -> -2 * Float64(π)^2 * sinpi(x) * sinpi(y))

        # --- finite differences (2nd-order central) ---
        m = mesh2(64)
        @test S.nsd(m) == 2
        f = on(m, (x, y) -> sinpi(x) * sinpi(y))
        g = S.gradient(Field(m, 0, f))
        @test g.rank == 1 && length(g.comp) == 2
        @test rel_l2(g.comp[1], exact_x(m)) < 5e-3
        @test rel_l2(g.comp[2], exact_y(m)) < 5e-3
        @test rel_l2(S.laplacian(Field(m, 0, f)).comp[1], exact_lap(m)) < 5e-3
        # divergence of u·f with u=[0.5,1.0] is 0.5 fx + 1.0 fy
        vf = S.fmul(S.vec_const_field(m, [0.5, 1.0]), Field(m, 0, f))
        @test rel_l2(S.divergence(vf).comp[1], 0.5 .* exact_x(m) .+ exact_y(m)) < 5e-3
        # 2nd-order convergence of ∇²
        e1 = l2(S.laplacian(Field(mesh2(32), 0, on(mesh2(32), (x,y)->sinpi(x)*sinpi(y)))).comp[1] .- exact_lap(mesh2(32)))
        e2 = l2(S.laplacian(Field(mesh2(64), 0, on(mesh2(64), (x,y)->sinpi(x)*sinpi(y)))).comp[1] .- exact_lap(mesh2(64)))
        @test 3.5 < e1 / e2 < 4.5

        # --- tensor-product SEM (reuses the SAME 1D Jexpresso basis) ---
        local msem
        sem_ok = true
        try
            msem = mesh2(4; method = :sem)        # 4×4 elements, nop=4
        catch err
            sem_ok = false
            @warn "2D SEM test skipped (could not load Jexpresso basis here)" err
        end
        if sem_ok
            @test msem.sem !== nothing
            fs = on(msem, (x, y) -> sinpi(x) * sinpi(y))
            gs = S.gradient(Field(msem, 0, fs))
            @test rel_l2(gs.comp[1], exact_x(msem)) < 1e-2     # spectral, well under FD
            @test rel_l2(gs.comp[2], exact_y(msem)) < 1e-2
            @test rel_l2(S.laplacian(Field(msem, 0, fs)).comp[1], exact_lap(msem)) < 2e-2
        else
            @test_skip sem_ok
        end
    end

    # ----------------------------------------------------------------------------
    @testset "2D advection-diffusion (structured FD)" begin
        # gaussian blob advected one full period on a periodic box; conservative
        # form ⇒ discrete mass is conserved, diffusion ⇒ the peak decays.
        inputs = Dict(:nsd => 2, :method => :fd, :periodic => true,
                      :xmin => -1.0, :xmax => 1.0, :ymin => -1.0, :ymax => 1.0,
                      :npoinx => 64, :npoiny => 64,
                      :u => [1.0, 1.0], :μ => 1e-3,
                      :q0 => (x, y) -> 1.0 + 0.5 * sinpi(x) * sinpi(y),
                      :tend => 0.3, :CFL => 0.3,
                      :outformat => "ascii", :plot_live => false,
                      :output_dir => mktempdir())
        @vars q u μ
        mesh, q0, q = SymbolicFD.solve(∂t(q) + ∇⋅(u*q) - μ*∇⋅∇(q), inputs)
        @test length(q) == 64 * 64
        @test all(isfinite, q)
        m0 = S.integrate(mesh, q0); m1 = S.integrate(mesh, q)
        @test abs(m1 - m0) / abs(m0) < 1e-6          # conservative
        @test maximum(q) ≤ maximum(q0) + 1e-6        # diffusion ⇒ peak decays
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
    @testset "symbolic equation DSL (no string)" begin
        @vars q u μ
        eq = ∂t(q) + ∇⋅(u*q) - μ*∇⋅∇(q)        # residual form, = 0 implied
        @test eq isa S.Node

        inp = Dict(:u => [1.0], :μ => 1e-3)
        var, mode, node, _ = S.build_from_residual(eq, inp)
        @test var == "q"
        @test mode == :transient

        # symbolic and string forms must build the SAME discrete RHS
        _, _, node_str = parse_equation("∂q/∂t + ∇⋅(u q) = μ∇²q", inp)
        m  = make_mesh(128)
        q0 = Float64[exp(-(x^2) / (2 * 0.1^2)) for x in m.x]
        dq_sym = zeros(m.npoin); S.build_rhs(node,     m)(dq_sym, q0)
        dq_str = zeros(m.npoin); S.build_rhs(node_str, m)(dq_str, q0)
        @test maximum(abs, dq_sym .- dq_str) < 1e-12

        # no ∂t ⇒ steady; ∇⋅∇ collapses to the compact Laplacian
        @vars q f
        _, mode2, node2, _ = S.build_from_residual(∇⋅∇(q) - f, Dict(:f => x -> 0.0))
        @test mode2 == :steady
        @test node2.args[1].op == :lap
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

        # the steady solve is matrix-free: a LinearOperator, no stored matrix
        L, b, _ = S.steady_operator(node, m, Dict(:bc_left => 0.0, :bc_right => 1.0))
        @test L isa AbstractLinearOperator
        @test size(L) == (m.npoin, m.npoin)

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
