#---------------------------------------------------------------------------------
# test/test_sphere_forcing.jl — the Scott & Polvani (2007) random small-scale
# forcing and large-scale dissipation on the spherical shell
# (src/kernel/operators/sphere_forcing.jl).
#
#   julia --project=. test/test_sphere_forcing.jl
#
# What is checked, in order of how much it would hurt to get wrong:
#
#   1. the real spherical-harmonic basis is ORTHONORMAL on the grid: YᵀMY/a² = I
#      for the whole forced band — which tests the Legendre recurrence, its
#      normalisation, the (λ, μ) of every node and the SEM quadrature at once;
#   2. u_F = n̂ × ∇ₛY has the energy the amplitude formula assumes,
#      ∫|u_F|² dA = n(n+1) for a unit orthonormal harmonic, and is tangential;
#   3. the ENERGY INJECTION RATE is the ε₀ the deck asked for: exactly (to the
#      gradient's discretisation error) for white-in-time forcing, and
#      statistically for the Markovian forcing, by integrating the linear
#      problem ∂u/∂t = u_F the way the time loop does (one draw per step,
#      forcing held fixed over the step) and reading the kinetic energy slope;
#   4. the two linear dissipations enter the RHS with the right sign and
#      coefficient, and the deck switch really turns everything off;
#   5. the draw is reproducible: same seed, same forcing.
#
# The fixture is the .msh that ships with the case, so the test needs no gmsh.
#---------------------------------------------------------------------------------

using Test
using Jexpresso
using Jexpresso: mod_mesh_mesh_driver, build_sphere_metrics, build_sphere_params,
                 build_sphere_forcing, sphere_forcing_update!, sphere_forcing_apply!,
                 sphere_harmonic_basis, sphere_forcing_diagnostics, sphere_zonal_mean
using PartitionedArrays, MPI
using LinearAlgebra, Printf

const CASE_MSH = joinpath(@__DIR__, "..", "problems", "ShallowWater", "SWsphere_ScottPolvani",
                          "cubed_sphere.msh")

# mod_inputs_user_inputs! reads three module globals that run.jl normally sets
# from the command line before it includes a case.
@eval Jexpresso begin
    parsed_equations           = "ShallowWater"
    parsed_equations_case_name = "SWsphere_ScottPolvani"
    user_input_file            = "test/test_sphere_forcing.jl"
end

function shell_inputs(nop; kwargs...)
    inputs = Dict{Symbol,Any}(
        :lread_gmsh           => true,
        :gmsh_filename        => CASE_MSH,
        :nop                  => nop,
        :interpolation_nodes  => "lgl",
        :backend              => Jexpresso.CPU(),
        :lspherical_shell     => true,
        :sphere_radius        => 7.1492e7,     # Jupiter
        :lproject_to_sphere   => true,
        :lgrid_only           => true,
    )
    Jexpresso.mod_inputs_user_inputs!(inputs, 1)   # rank 1 => no banner spam
    inputs[:gmsh_filename] = CASE_MSH
    inputs[:nop]           = nop
    for (k, v) in kwargs
        inputs[k] = v
    end
    return inputs
end

# a resting layer, q = [φ_ref, 0, 0, 0] with the spare fifth column define_q carries
function rest_state(npoin, φref)
    q = zeros(Float64, npoin, 5)
    q[:, 1] .= φref
    return q
end

# ∫ f dA over the shell, with the assembled diagonal mass matrix
integral(f, M) = sum(M[ip]*f[ip] for ip = 1:length(f))


@testset "spherical shell: Scott & Polvani forcing" begin

    isfile(CASE_MSH) || error("missing $CASE_MSH — it ships with the SWsphere_ScottPolvani case")

    with_mpi() do distribute

        nop  = 5
        a    = 7.1492e7
        Ω    = 1.7585e-4
        φref = (2Ω*0.025*a)^2               # L_D = 0.025 a
        nf   = 8                            # low enough that the grid integrates Y² exactly
        ε0   = 1.0e-5                       # m²/s³

        inputs = shell_inputs(nop;
                              :lsphere_forcing   => true,
                              :forcing_epsilon   => ε0,
                              :forcing_nf        => nf,
                              :forcing_dn        => 4,
                              :forcing_tau       => 0.0,
                              :forcing_seed      => 42,
                              :rayleigh_friction => 0.0,
                              :sphere_Omega      => Ω)

        mesh, _  = mod_mesh_mesh_driver(inputs, 1, distribute)
        metrics  = build_sphere_metrics(mesh, inputs; verbose = false)
        npoin    = Int(mesh.npoin)
        M        = metrics.M

        @test abs(sum(M) - 4π*a^2) < 1e-8*4π*a^2      # the quadrature integrates the sphere

        #-------------------------------------------------------------------------
        # (1) orthonormality of the basis on the grid: (1/a²) YᵀMY = I
        #-------------------------------------------------------------------------
        ns = collect(nf-2:nf+2)
        Y, pairs = sphere_harmonic_basis(mesh, ns)
        nmodes = size(Y, 2)
        @test nmodes == sum(2n+1 for n in ns)
        @test length(pairs) == sum(n+1 for n in ns)

        G = (Y' * (M .* Y)) ./ a^2
        errG = maximum(abs, G - I)
        @info @sprintf("  harmonic basis n ∈ [%d, %d], %d modes: max|YᵀMY/a² - I| = %.3e", ns[1], ns[end], nmodes, errG)
        @test errG < 1.0e-6

        #-------------------------------------------------------------------------
        # (2) u_F of one unit harmonic: ∫|u_F|² dA = n(n+1), tangential
        #-------------------------------------------------------------------------
        sp = build_sphere_params(mesh, metrics, inputs; neqs = 4, verbose = false)
        fc = sp.forcing
        @test fc !== nothing
        @test fc.nband == 5

        for (n, k) in ((nf, pairs[findfirst(p -> p.n == nf && p.m == 3, pairs)].kc),
                       (nf-2, pairs[findfirst(p -> p.n == nf-2 && p.m == 0, pairs)].kc),
                       (nf+2, pairs[findfirst(p -> p.n == nf+2 && p.m == nf+2, pairs)].ks))
            ψ = Y[:, k]
            uF = zeros(Float64, npoin, 4)
            Jexpresso._forcing_velocity_kernel!(uF, ψ, mesh.connijk,
                                                Int(mesh.nelem), Int(mesh.ngl),
                                                metrics.dξdx, metrics.dξdy, metrics.dξdz,
                                                metrics.dηdx, metrics.dηdy, metrics.dηdz,
                                                metrics.nx, metrics.ny, metrics.nz,
                                                metrics.Je, metrics.dψ, metrics.ω, 4)
            Jexpresso._sphere_dss_scale!(uF, metrics.Minv, npoin, 4, sp.dss)

            e2   = integral(uF[:,2].^2 .+ uF[:,3].^2 .+ uF[:,4].^2, M)
            err  = abs(e2 - n*(n+1))/(n*(n+1))
            @info @sprintf("  n = %2d : ∫|n̂×∇Y|² dA = %.6f (exact %d), rel. error %.3e", n, e2, n*(n+1), err)
            @test err < 1.0e-4

            # tangential: u_F = n̂ × ∇ₛψ lies in the DISCRETE tangent plane, and
            # the discrete normal n̂ differs from x̂ by the interpolation error of
            # the shell (~3e-7 relative at nop = 5), so (u_F · x̂) is small to
            # that order, not to round-off. The Lagrange projection after every
            # stage removes what is left.
            tmax = maximum(abs(uF[ip,2]*mesh.coords[1,ip] + uF[ip,3]*mesh.coords[2,ip] +
                               uF[ip,4]*mesh.coords[3,ip])/a for ip = 1:npoin)
            umax = maximum(sqrt(uF[ip,2]^2 + uF[ip,3]^2 + uF[ip,4]^2) for ip = 1:npoin)
            @test tmax < 1.0e-5*umax
        end

        #-------------------------------------------------------------------------
        # (3a) white-in-time: every draw injects exactly ε₀ (given the gradient)
        #-------------------------------------------------------------------------
        q  = rest_state(npoin, φref)
        Δt = 500.0
        for _ = 1:5
            sphere_forcing_update!(fc, q, Δt, mesh, metrics, sp)
            @test abs(fc.Pinj - ε0)/ε0 < 1.0e-4
        end
        # the vorticity forcing written to VTK is ∇ₛ²ψ_F: check on the ψ_F just drawn
        #   ∫ F ψ dA = -∫|∇ψ|² dA = -∫|u_F|² dA
        @test abs(integral(fc.F .* fc.ψ, M) + integral(fc.uF[:,2].^2 .+ fc.uF[:,3].^2 .+ fc.uF[:,4].^2, M)) <
              1.0e-4*integral(fc.uF[:,2].^2 .+ fc.uF[:,3].^2 .+ fc.uF[:,4].^2, M)

        # the linear problem ∂(φu)/∂t = φu_F, one draw per step: KE grows at ε₀
        nsteps = 400
        q  = rest_state(npoin, φref)
        for _ = 1:nsteps
            sphere_forcing_update!(fc, q, Δt, mesh, metrics, sp)
            @. q[:, 2:4] += Δt*φref*fc.uF[:, 2:4]
        end
        KE = sphere_forcing_diagnostics(q, fc, mesh, metrics).KE
        rate = KE/(nsteps*Δt)
        @info @sprintf("  white in time: KE/(NΔt) = %.4e against ε₀ = %.4e (ratio %.4f)", rate, ε0, rate/ε0)
        @test abs(rate/ε0 - 1) < 0.15

        #-------------------------------------------------------------------------
        # (3b) Markovian, τ = 20Δt, FIXED amplitude (:forcing_normalize => false):
        #      the linear-response formula gives ε₀ statistically on the linear
        #      problem, and the gain stays 1
        #-------------------------------------------------------------------------
        inputs_m = copy(inputs); inputs_m[:forcing_tau] = 20Δt; inputs_m[:forcing_normalize] = false
        sp_m = build_sphere_params(mesh, metrics, inputs_m; neqs = 4, verbose = false)
        fcm  = sp_m.forcing
        @test !fcm.lnorm
        nsteps = 3000
        q  = rest_state(npoin, φref)
        for _ = 1:nsteps
            sphere_forcing_update!(fcm, q, Δt, mesh, metrics, sp_m)
            @. q[:, 2:4] += Δt*φref*fcm.uF[:, 2:4]
        end
        @test abs(fcm.r - exp(-1/20)) < 1e-12
        @test fcm.gain == 1.0
        d = sphere_forcing_diagnostics(q, fcm, mesh, metrics)
        rate = d.KE/(nsteps*Δt)
        @info @sprintf("  Markovian (τ = 20Δt), fixed amplitude: KE/(NΔt) = %.4e, mean measured injection = %.4e, ε₀ = %.4e", rate, d.Pmean, ε0)
        @test abs(rate/ε0 - 1) < 0.2
        @test abs(d.Pmean_rel - 1) < 0.2
        @test d.Ro ≈ d.U/(2a*Ω)

        #-------------------------------------------------------------------------
        # (3c) Markovian, RENORMALISED every step (the default, paper Eq. 14):
        #      every step injects exactly ε₀ on the linear problem, whatever the
        #      correlation between the flow and the forcing, and the gain moves
        #      away from 1 to do it
        #-------------------------------------------------------------------------
        inputs_n = copy(inputs); inputs_n[:forcing_tau] = 20Δt
        sp_n = build_sphere_params(mesh, metrics, inputs_n; neqs = 4, verbose = false)
        fcn  = sp_n.forcing
        @test fcn.lnorm
        nsteps = 300
        q  = rest_state(npoin, φref)
        gains = Float64[]
        for _ = 1:nsteps
            sphere_forcing_update!(fcn, q, Δt, mesh, metrics, sp_n)
            @test abs(fcn.Pinj - ε0)/ε0 < 1.0e-10
            push!(gains, fcn.gain)
            @. q[:, 2:4] += Δt*φref*fcn.uF[:, 2:4]
        end
        d = sphere_forcing_diagnostics(q, fcn, mesh, metrics)
        rate = d.KE/(nsteps*Δt)
        @info @sprintf("  Markovian (τ = 20Δt), renormalised: KE/(NΔt) = %.6e, ε₀ = %.6e, gain ∈ [%.3f, %.3f]", rate, ε0, minimum(gains), maximum(gains))
        @test abs(rate/ε0 - 1) < 1.0e-8        # exact: the linear problem IS the energy budget used
        @test abs(d.Pmean_rel - 1) < 1.0e-10
        @test minimum(gains) < 0.999 || maximum(gains) > 1.001   # it did have to act
        @test all(gains .> 0)

        #-------------------------------------------------------------------------
        # (4) the dissipation terms and the deck switch
        #-------------------------------------------------------------------------
        inputs_d = copy(inputs)
        inputs_d[:forcing_epsilon]      = 0.0        # no forcing, only damping
        inputs_d[:rayleigh_friction]    = 1.0e-3
        inputs_d[:radiative_relaxation] = 2.0e-3
        sp_d = build_sphere_params(mesh, metrics, inputs_d; neqs = 4, verbose = false)
        fcd  = sp_d.forcing
        qd   = rest_state(npoin, φref)
        qd[:, 1] .*= 1.1
        qd[:, 2] .= 3.0; qd[:, 3] .= -2.0; qd[:, 4] .= 1.0
        qe   = rest_state(npoin, φref)
        sphere_forcing_update!(fcd, qd, Δt, mesh, metrics, sp_d)
        @test maximum(abs, fcd.uF) == 0.0
        RHS = zeros(Float64, npoin, 4)
        sphere_forcing_apply!(RHS, qd, qe, fcd, npoin)
        @test RHS[:, 1] ≈ -2.0e-3 .* (qd[:, 1] .- φref)
        @test RHS[:, 2] ≈ -1.0e-3 .* qd[:, 2]
        @test RHS[:, 3] ≈ -1.0e-3 .* qd[:, 3]
        @test RHS[:, 4] ≈ -1.0e-3 .* qd[:, 4]

        inputs_off = copy(inputs); inputs_off[:lsphere_forcing] = false
        sp_off = build_sphere_params(mesh, metrics, inputs_off; neqs = 4, verbose = false)
        @test sp_off.forcing === nothing
        RHS0 = copy(RHS)
        sphere_forcing_apply!(RHS, qd, qe, sp_off.forcing, npoin)
        @test RHS == RHS0

        @test_throws ErrorException build_sphere_forcing(mesh, metrics,
            Dict{Symbol,Any}(:lsphere_forcing => true, :forcing_epsilon => 1.0, :forcing_nf => 0))
        @test_throws ErrorException build_sphere_forcing(mesh, metrics,
            Dict{Symbol,Any}(:lsphere_forcing => true, :forcing_epsilon => -1.0, :forcing_nf => 4))

        #-------------------------------------------------------------------------
        # (5) reproducible: two forcings with the same seed draw the same field
        #-------------------------------------------------------------------------
        sp1 = build_sphere_params(mesh, metrics, inputs; neqs = 4, verbose = false)
        sp2 = build_sphere_params(mesh, metrics, inputs; neqs = 4, verbose = false)
        q   = rest_state(npoin, φref)
        for _ = 1:3
            sphere_forcing_update!(sp1.forcing, q, Δt, mesh, metrics, sp1)
            sphere_forcing_update!(sp2.forcing, q, Δt, mesh, metrics, sp2)
        end
        @test sp1.forcing.ψ == sp2.forcing.ψ
        @test sp1.forcing.uF == sp2.forcing.uF
        inputs_s = copy(inputs); inputs_s[:forcing_seed] = 43
        sp3 = build_sphere_params(mesh, metrics, inputs_s; neqs = 4, verbose = false)
        sphere_forcing_update!(sp3.forcing, q, Δt, mesh, metrics, sp3)
        @test sp3.forcing.ψ != sp1.forcing.ψ

        #-------------------------------------------------------------------------
        # (6) the zonal mean of a solid-body rotation u = U cos φ e_λ is U cos φ
        #-------------------------------------------------------------------------
        U  = 10.0
        qz = rest_state(npoin, φref)
        for ip = 1:npoin
            x, y, z = mesh.coords[1,ip], mesh.coords[2,ip], mesh.coords[3,ip]
            λ  = atan(y, x); φl = asin(z/a)
            qz[ip, 2] = φref*(-U*cos(φl)*sin(λ))
            qz[ip, 3] = φref*( U*cos(φl)*cos(λ))
        end
        lat, ū, v̄, φ̄, cnt = sphere_zonal_mean(qz, mesh, metrics; nbins = 18)
        @test all(cnt .> 0)
        @test maximum(abs.(ū .- U .* cosd.(lat))) < 0.05*U     # 10° bands: cos φ varies inside
        @test maximum(abs, v̄) < 1e-10*U
        @test all(φ̄ .≈ φref)
    end
end
