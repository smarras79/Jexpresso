#---------------------------------------------------------------------------------
# test/test_sphere_dsgs.jl — the residual-based DynSGS viscosity on the
# spherical shell (src/kernel/operators/sphere_dsgs.jl).
#
#   julia --project=. test/test_sphere_dsgs.jl
#
# What is checked:
#
#   1. the deck switch: :visc_model => DSGS() builds the model, AV() does not,
#      and :μ then holds the dimensionless multipliers [0, 1, 1, 1];
#   2. a resting layer with a zero residual gets ν_e = 0 everywhere — the run
#      starts inviscid — and ν_e ≤ ν_max = C₂Δ(|u|+√φ) always;
#   3. on a state that is NOT a solution (a random field, so the residual is
#      large) ν_e is positive, finite and at the cap: the floor keeps a
#      resting-layer normalisation from blowing up;
#   4. the two-pass RHS of a DynSGS run equals the inviscid RHS plus the
#      viscous term with the same ν_e assembled separately — i.e. splitting
#      the assembly changed nothing — and on a state with ν_e ≡ 0 it equals
#      the plain inviscid RHS;
#   5. the viscous term with DynSGS coefficients can only remove energy:
#      Σ_ip (φu)·RHS_visc ≤ 0;
#   6. the constant-ν path is untouched: the per-element scale defaults to 1.
#
# The fixture is the .msh that ships with the case.
#---------------------------------------------------------------------------------

using Test
using Jexpresso
using Jexpresso: mod_mesh_mesh_driver, build_sphere_metrics, build_sphere_params,
                 build_sphere_dsgs, sphere_dsgs_viscosity!, sphere_dsgs_roll!,
                 sphere_dsgs_cap_estimate, sphere_dsgs_nodal!, sphere_rhs!,
                 DSGS, AV, CL, TOTAL
using PartitionedArrays, MPI
using Printf

const CASE_DIR = joinpath(@__DIR__, "..", "problems", "ShallowWater", "SWsphere_ScottPolvani_DynSGS")
const CASE_MSH = joinpath(CASE_DIR, "cubed_sphere.msh")

@eval Jexpresso begin
    parsed_equations           = "ShallowWater"
    parsed_equations_case_name = "SWsphere_ScottPolvani_DynSGS"
    user_input_file            = "test/test_sphere_dsgs.jl"
end

# the case's equations (fluxes, sources, primitives) and its Ω, g globals
for f in ("user_flux.jl", "user_source.jl", "user_primitives.jl", "user_bc.jl")
    Jexpresso.include(joinpath(CASE_DIR, f))
end

function shell_inputs(nop; kwargs...)
    inputs = Dict{Symbol,Any}(
        :lread_gmsh           => true,
        :gmsh_filename        => CASE_MSH,
        :nop                  => nop,
        :interpolation_nodes  => "lgl",
        :backend              => Jexpresso.CPU(),
        :lspherical_shell     => true,
        :sphere_radius        => 7.1492e7,
        :lproject_to_sphere   => true,
        :lgrid_only           => true,
    )
    Jexpresso.mod_inputs_user_inputs!(inputs, 1)
    inputs[:gmsh_filename] = CASE_MSH
    inputs[:nop]           = nop
    for (k, v) in kwargs
        inputs[k] = v
    end
    return inputs
end

function rest_state(npoin, φref)
    q = zeros(Float64, npoin, 5)
    q[:, 1] .= φref
    return q
end

# a tangential random momentum on top of the resting layer
function noisy_state(mesh, φref, amp)
    npoin = Int(mesh.npoin)
    q = rest_state(npoin, φref)
    for ip = 1:npoin
        x, y, z = mesh.coords[1,ip], mesh.coords[2,ip], mesh.coords[3,ip]
        r = sqrt(x*x + y*y + z*z)
        v = (amp*(2rand() - 1), amp*(2rand() - 1), amp*(2rand() - 1))
        # remove the radial part so the state is on the shell
        vr = (v[1]*x + v[2]*y + v[3]*z)/r^2
        q[ip, 2] = φref*(v[1] - vr*x)
        q[ip, 3] = φref*(v[2] - vr*y)
        q[ip, 4] = φref*(v[3] - vr*z)
        q[ip, 1] = φref*(1 + 0.01*(2rand() - 1))
    end
    return q
end


@testset "spherical shell: DynSGS viscosity" begin

    isfile(CASE_MSH) || error("missing $CASE_MSH — it ships with the case")

    with_mpi() do distribute

        nop  = 4
        a    = 7.1492e7
        Ω    = 1.7585e-4
        φref = (2Ω*0.025*a)^2
        Δt   = 500.0

        base = shell_inputs(nop; :lvisc => true, :ivisc_equations => [2, 3, 4],
                                 :μ => [0.0, 1.0, 1.0, 1.0], :sphere_Omega => Ω)
        mesh, _ = mod_mesh_mesh_driver(base, 1, distribute)
        metrics = build_sphere_metrics(mesh, base; verbose = false)
        npoin   = Int(mesh.npoin)
        nelem   = Int(mesh.nelem)
        neqs    = 4

        q0 = rest_state(npoin, φref)
        qe = rest_state(npoin, φref)

        #--- (1) the switch
        inp_d = copy(base); inp_d[:visc_model] = DSGS()
        sp_d  = build_sphere_params(mesh, metrics, inp_d; neqs = neqs, verbose = false, q0 = q0)
        @test sp_d.dsgs !== nothing
        @test sp_d.μ == [0.0, 1.0, 1.0, 1.0]
        @test sp_d.lvisc
        @test length(sp_d.dsgs.Δ) == nelem && all(>(0), sp_d.dsgs.Δ)

        inp_c = copy(base); inp_c[:visc_model] = AV(); inp_c[:μ] = 1.0e5
        sp_c  = build_sphere_params(mesh, metrics, inp_c; neqs = neqs, verbose = false, q0 = q0)
        @test sp_c.dsgs === nothing
        @test sp_c.μ == [0.0, 1.0e5, 1.0e5, 1.0e5]

        sp_d.Δt = Δt; sp_c.Δt = Δt

        #--- (2) rest: zero residual, zero viscosity; the cap bounds everything
        RHS = zeros(Float64, npoin, neqs)
        sphere_rhs!(RHS, q0, qe, mesh, metrics, sp_d, TOTAL())
        @test maximum(abs, RHS) < 1.0e-12*φref        # the rest state is a steady solution
        cap = sphere_dsgs_cap_estimate(sp_d.dsgs, q0, npoin)
        @test cap ≈ 0.5*maximum(sp_d.dsgs.Δ)*sqrt(φref)
        # the residual is the round-off of that RHS over the floored
        # normalisation: ν_e is not exactly zero, it is nothing against the cap
        @test maximum(sp_d.dsgs.νel) < 1.0e-8*cap

        #--- (3) a random state: the residual is large, ν_e finite and positive
        qn = noisy_state(mesh, φref, 5.0)
        # the history says "the state was at rest": a large BDF2 derivative
        sphere_rhs!(RHS, qn, qe, mesh, metrics, sp_d, TOTAL())
        ν = sp_d.dsgs.νel
        @test all(isfinite, ν)
        @test maximum(ν) > 0
        # never above the element cap
        ngl = Int(mesh.ngl)
        for iel = 1:nelem
            umax = 0.0
            for j = 1:ngl, i = 1:ngl
                ip = mesh.connijk[iel,i,j]
                φ  = qn[ip,1]
                umax = max(umax, sqrt(qn[ip,2]^2 + qn[ip,3]^2 + qn[ip,4]^2)/φ + sqrt(φ))
            end
            @test ν[iel] <= 0.5*sp_d.dsgs.Δ[iel]*umax*(1 + 1e-12)
        end
        @info @sprintf("  random state: ν_e ∈ [%.3e, %.3e] m²/s, %d of %d elements at the cap",
                       minimum(ν), maximum(ν), sp_d.dsgs.ncap, nelem)

        # nodal field: the max over the elements at a node
        νn = sphere_dsgs_nodal!(sp_d.dsgs, mesh)
        @test maximum(νn) ≈ maximum(ν)
        @test minimum(νn) >= minimum(ν)

        #--- (4) the two-pass assembly: inviscid + separately assembled viscous
        RHS_inv = zeros(Float64, npoin, neqs)
        inp_i = copy(base); inp_i[:lvisc] = false
        sp_i  = build_sphere_params(mesh, metrics, inp_i; neqs = neqs, verbose = false, q0 = q0)
        sphere_rhs!(RHS_inv, qn, qe, mesh, metrics, sp_i, TOTAL())

        acc = zeros(Float64, npoin, neqs)
        Jexpresso._sphere_visc_assemble!(acc, qn, mesh.connijk, nelem, ngl,
                                         metrics.Je,
                                         metrics.dξdx, metrics.dξdy, metrics.dξdz,
                                         metrics.dηdx, metrics.dηdy, metrics.dηdz,
                                         metrics.dψ, metrics.ω, sp_d.μ, sp_d.gξ, sp_d.gη, neqs, ν)
        Jexpresso._sphere_dss_scale!(acc, metrics.Minv, npoin, neqs, sp_d.dss)
        # sp_d.dsgs.νel was computed from RHS_inv (the same state, the same
        # history), so the DynSGS RHS must be RHS_inv + acc
        sphere_rhs!(RHS, qn, qe, mesh, metrics, sp_d, TOTAL())
        @test maximum(abs, RHS .- (RHS_inv .+ acc)) < 1.0e-9*maximum(abs, RHS)
        @test maximum(abs, acc[:, 1]) == 0.0          # no mass diffusion

        # after a rolled history that matches the state, the residual is only
        # the spatial part; with q = history the BDF2 derivative vanishes
        sphere_dsgs_roll!(sp_d.dsgs, q0); sphere_dsgs_roll!(sp_d.dsgs, q0)
        sphere_rhs!(RHS, q0, qe, mesh, metrics, sp_d, TOTAL())
        @test maximum(sp_d.dsgs.νel) < 1.0e-8*cap

        #--- (5) the DynSGS viscous term is dissipative
        @test sum(qn[ip, k]*acc[ip, k]*metrics.M[ip] for ip = 1:npoin, k = 2:4) < 0.0

        #--- (6) the constant-ν path: scale defaults to one, same as before
        RHS_c  = zeros(Float64, npoin, neqs)
        sphere_rhs!(RHS_c, qn, qe, mesh, metrics, sp_c, TOTAL())
        acc_c = zeros(Float64, npoin, neqs)
        Jexpresso._sphere_visc_assemble!(acc_c, qn, mesh.connijk, nelem, ngl,
                                         metrics.Je,
                                         metrics.dξdx, metrics.dξdy, metrics.dξdz,
                                         metrics.dηdx, metrics.dηdy, metrics.dηdz,
                                         metrics.dψ, metrics.ω, sp_c.μ, sp_c.gξ, sp_c.gη, neqs,
                                         ones(Float64, nelem))
        Jexpresso._sphere_dss_scale!(acc_c, metrics.Minv, npoin, neqs, sp_c.dss)
        @test maximum(abs, RHS_c .- (RHS_inv .+ acc_c)) < 1.0e-9*maximum(abs, RHS_c)
    end
end
