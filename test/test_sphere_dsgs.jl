#---------------------------------------------------------------------------------
# test/test_sphere_dsgs.jl — the Marras-Nazarov DynSGS model on the spherical
# shell (compute_dsgs_viscosity!(::DSGS_SW, ::NSD_2D) in
# src/kernel/physics/SGS.jl, wired in src/kernel/operators/sphere_rhs.jl).
#
#   julia --project=. test/test_sphere_dsgs.jl
#
# The Laplace-Beltrami operator the model feeds is already tested against its
# spherical-harmonic eigenvalues in test_sphere_visc.jl. What is tested HERE is
# the thing DynSGS adds on top of it: WHERE the viscosity goes and HOW BIG it
# is. Those are the two claims the model makes, and each is checkable exactly.
#
#   1  a state that satisfies the discrete equations gets NO viscosity
#   2  μ never exceeds its own C2·Δ·(|u|+c) cap
#   3  μ_res scales as C1·Δ²·‖R‖/‖q−⟨q⟩‖ — the formula, term by term
#   4  the per-equation multipliers of the deck are applied
#   5  SELECTIVITY: an under-resolved feature gets orders of magnitude more
#      viscosity than a smooth resolved one. This is the whole point of using a
#      residual sensor instead of Smagorinsky, and it is the test that a model
#      which merely "damps things" fails.
#   6  a purely ZONAL state does not saturate the cap. This is a regression
#      test: the Cartesian component φw of a zonal jet is identically zero, so
#      normalizing its (large, Coriolis-driven) residual by its own zero spread
#      pinned μ at the cap on every element and killed the Galewsky run at
#      t ≈ 1100 s. The three momentum slots share one denominator for exactly
#      this reason — see the header of compute_dsgs_viscosity!(::DSGS_SW).
#   7  the deck switches, including the guard that rejects a constant-ν :μ
#
# The fixture is the .msh that ships with the SWsphere cases, so the test needs
# no gmsh binary.
#---------------------------------------------------------------------------------

using Test
using Jexpresso
using Jexpresso: mod_mesh_mesh_driver, build_sphere_metrics, build_sphere_params,
                 build_sphere_viscosity, build_sphere_element_size,
                 sphere_dsgs_requested, compute_dsgs_viscosity!,
                 DSGS_SW, NSD_2D
using PartitionedArrays, MPI
using Random, Printf

const CASE_DIR = joinpath(@__DIR__, "..", "problems", "ShallowWater", "SWsphereDSGS")
const CASE_MSH = joinpath(CASE_DIR, "cubed_sphere.msh")

# mod_inputs_user_inputs! reads three module globals that run.jl normally sets
# from the command line before it includes a case.
@eval Jexpresso begin
    parsed_equations           = "ShallowWater"
    parsed_equations_case_name = "SWsphereDSGS"
    user_input_file            = "test/test_sphere_dsgs.jl"
end

#
# Test 5 drives the model with a REAL right-hand side rather than a hand-written
# one, and sphere_rhs! reaches back into the case's own callbacks — so the case
# files have to be loaded into the module, which is what run.jl does before it
# calls the driver. Nothing else in this file needs them.
#
for _f in ("user_flux.jl", "user_source.jl", "user_primitives.jl")
    Base.include(Jexpresso, joinpath(CASE_DIR, _f))
end

function shell_inputs(nop; kwargs...)
    inputs = Dict{Symbol,Any}(
        :lread_gmsh           => true,
        :gmsh_filename        => CASE_MSH,
        :nop                  => nop,
        :interpolation_nodes  => "lgl",
        :backend              => Jexpresso.CPU(),
        :lspherical_shell     => true,
        :sphere_radius        => 6.37122e6,
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

#
# The model, called the way sphere_rhs.jl calls it. Returns μ_dsgs[nelem, neqs].
#
function dsgs(mesh, metrics, sp, q, q1, q2, rhs, Δt)
    fill!(sp.μ_dsgs, 0.0)
    compute_dsgs_viscosity!(sp.μ_dsgs, DSGS_SW(), NSD_2D(),
                            q, q1, q2, rhs, metrics.Minv, sp.μ,
                            sp.davg, sp.ddenom, Δt,
                            mesh.connijk, sp.Δelem, sp.C1, sp.C2,
                            Int(mesh.nelem), Int(mesh.ngl))
    return sp.μ_dsgs
end

#
# A geostrophically trivial but structurally realistic state: a ZONAL flow of
# speed U on top of a uniform geopotential, in the conservative Cartesian form
# the case integrates, q = (φ, φu, φv, φw) with u = U·e_λ.
#
# e_λ = (−sinλ, cosλ, 0), so φw ≡ 0 — which is the degeneracy test 6 is about.
#
function zonal_state(mesh, neqs; φ0 = 1.0e5, U = 80.0)
    npoin = Int(mesh.npoin)
    q = zeros(Float64, npoin, neqs)
    for ip = 1:npoin
        λ = mesh.lon[ip]
        φlat = mesh.lat[ip]
        sλ, cλ = sincos(λ)
        q[ip,1] = φ0
        q[ip,2] = φ0*U*(-sλ)*cos(φlat)
        q[ip,3] = φ0*U*( cλ)*cos(φlat)
        q[ip,4] = 0.0
    end
    return q
end


@testset "spherical shell: DynSGS (DSGS_SW)" begin

    isfile(CASE_MSH) || error("missing $CASE_MSH — it ships with the SWsphereDSGS case")

    with_mpi() do distribute

        nop    = 5
        inputs = shell_inputs(nop; :lvisc => true, :visc_model => DSGS_SW(),
                                   :μ => 1.0, :ivisc_equations => [1,2,3,4],
                                   :dsgs_C1 => 1.0, :dsgs_C2 => 0.5)

        mesh, _ = mod_mesh_mesh_driver(inputs, 1, distribute)
        metrics = build_sphere_metrics(mesh, inputs; verbose = false)
        sp      = build_sphere_params(mesh, metrics, inputs; neqs = 4)

        npoin = Int(mesh.npoin)
        nelem = Int(mesh.nelem)
        ngl   = Int(mesh.ngl)
        neqs  = 4
        R     = Float64(mesh.radius)
        Δt    = 74.0

        @test sp.ldsgs
        @test sp.lvisc
        @test sp.μ == [1.0, 1.0, 1.0, 1.0]

        #-------------------------------------------------------------------------
        # (0) the element length scale. On a cubed sphere with 10 elements per
        #     panel edge, an element spans ~1/40 of a great circle, so its chord
        #     is ~2πR/40 ≈ 1.0e6 m. The cubed-sphere map makes corner elements
        #     smaller than face-centre ones, hence a spread rather than a value.
        #-------------------------------------------------------------------------
        Δe = build_sphere_element_size(mesh)
        @test length(Δe) == nelem
        @test all(>(0.0), Δe)
        @test 0.3e6 < minimum(Δe) && maximum(Δe) < 2.0e6
        # no element is wildly out of scale with its neighbours
        @test maximum(Δe)/minimum(Δe) < 3.0

        q = zonal_state(mesh, neqs)

        # μ_max|e = C2·(Δ_elem/ngl)·(|u| + √φ)∞,e, the model's own cap, rebuilt
        # here independently of the kernel. Used as the yardstick throughout.
        cap = zeros(Float64, nelem)
        for iel = 1:nelem
            Δ = Δe[iel]/ngl
            w = 0.0
            for j = 1:ngl, i = 1:ngl
                ip = mesh.connijk[iel,i,j]
                φ  = q[ip,1]
                u  = sqrt(q[ip,2]^2 + q[ip,3]^2 + q[ip,4]^2)/φ
                w  = max(w, u + sqrt(φ))
            end
            cap[iel] = sp.C2*Δ*w
        end
        capmax = maximum(cap)

        #-------------------------------------------------------------------------
        # (1) A DISCRETELY EXACT STATE GETS NO VISCOSITY.
        #
        #     R = (3q − 4q + q)/(2Δt) − M⁻¹·0 = 0, so μ_res = 0 and
        #     μ = max(0, min(μ_max, 0)) = 0 — however large the cap is. A model
        #     that leaked a floor, or that took the max the wrong way round,
        #     fails here.
        #
        #     "= 0" up to round-off, not bitwise: (3q − 4q + q) is evaluated in
        #     floating point, and 3q carries a rounding error of order eps·q
        #     that the cancellation does not remove. What survives is ~1e-19 of
        #     the cap, which the bound below is 10⁷ times looser than.
        #-------------------------------------------------------------------------
        rhs = zeros(Float64, npoin, neqs)

        μ0 = dsgs(mesh, metrics, sp, q, q, q, rhs, Δt)
        @test maximum(μ0) < 1.0e-12*capmax

        #-------------------------------------------------------------------------
        # (2) THE CAP IS NEVER EXCEEDED, whatever the residual does.
        #
        #     Driven with an absurd RHS — 10¹⁵ times any physical value — μ must
        #     land exactly on C2·Δ·(|u|+√φ) per element and not one bit above.
        #     This is the "degrade to first-order upwind, never worse" property.
        #-------------------------------------------------------------------------
        rng  = MersenneTwister(20250815)
        rhsB = 1.0e15 .* randn(rng, npoin, neqs)

        μB = dsgs(mesh, metrics, sp, q, q, q, rhsB, Δt)

        @test all(μB[iel,ieq] <= cap[iel]*(1 + 1.0e-12)
                  for iel = 1:nelem, ieq = 1:neqs)
        # and it really is AT the cap — the min() branch is live, not vacuous
        @test maximum(abs, μB[:,2] .- cap) <= 1.0e-9*capmax

        #-------------------------------------------------------------------------
        # (3) THE FORMULA, term by term: below the cap, μ = C1·Δ²·‖R‖∞,e/denom.
        #
        #     Built by hand for one element: put a residual on the φ equation
        #     alone, of a size chosen to stay well under the cap, and check the
        #     element coefficient against C1·Δ²·‖R‖/‖φ−⟨φ⟩‖. Also check that
        #     doubling C1 doubles μ and that the neighbouring elements, which
        #     carry no residual, stay at zero — i.e. the model is LOCAL.
        #-------------------------------------------------------------------------
        # a non-degenerate φ field, so ‖φ−⟨φ⟩‖ is a real number and not a floor
        qv = copy(q)
        for ip = 1:npoin
            qv[ip,1] = 1.0e5*(1.0 + 0.05*sin(mesh.lat[ip]))
        end
        #
        # ⟨φ⟩ and ‖φ−⟨φ⟩‖∞ as the KERNEL forms them: a sum over (element, node)
        # pairs, so a node shared by k elements is counted k times. That is the
        # convention every DSGS path in Jexpresso uses — the mean is only there
        # to centre an L∞ norm, not to be a conserved integral — and rebuilding
        # it here as a plain average over unique nodes would leave this check
        # disagreeing with the kernel in the 9th digit for the wrong reason.
        #
        φsum = 0.0
        for iel = 1:nelem, j = 1:ngl, i = 1:ngl
            φsum += qv[mesh.connijk[iel,i,j], 1]
        end
        φavg   = φsum/(nelem*ngl*ngl)
        denomφ = 0.0
        for iel = 1:nelem, j = 1:ngl, i = 1:ngl
            denomφ = max(denomφ, abs(qv[mesh.connijk[iel,i,j],1] - φavg))
        end
        denomφ = max(denomφ, 1.0e-3*abs(φavg)) + 1.0e-16

        iel_t = nelem ÷ 2
        rhsL  = zeros(Float64, npoin, neqs)
        Rwant = 1.0e-4                      # the residual we want, in φ units/s
        for j = 1:ngl, i = 1:ngl
            ip = mesh.connijk[iel_t,i,j]
            # R = |0 − M⁻¹·rhs| = Rwant  ⇒  rhs = Rwant/Minv = Rwant·M
            rhsL[ip,1] = -Rwant*metrics.M[ip]
        end

        μL   = dsgs(mesh, metrics, sp, qv, qv, qv, rhsL, Δt)
        Δ_t  = Δe[iel_t]/ngl
        want = sp.C1*Δ_t*Δ_t*Rwant/denomφ

        @test want < cap[iel_t]             # the check below is of μ_res, not the cap
        #
        # rtol 1e-7, not machine precision: the kernel's BDF2 term is
        # (3φ − 4φ + φ)/(2Δt) evaluated in floating point on φ ~ 1e5, so it
        # contributes ~eps·φ/(2Δt) ≈ 7e-14 of spurious residual against the
        # Rwant = 1e-4 imposed here — a relative floor of ~1e-9. Everything
        # above that is the formula.
        #
        @test isapprox(μL[iel_t,1], want; rtol = 1.0e-7)

        # LOCALITY: an element that shares no node with iel_t sees only that
        # same round-off floor, five hundred times below the bound here, where
        # the element carrying the residual sees 447 m²/s.
        touched = Set(mesh.connijk[iel_t,i,j] for i = 1:ngl, j = 1:ngl)
        for iel = 1:nelem
            any(in(touched), (mesh.connijk[iel,i,j] for i = 1:ngl, j = 1:ngl)) && continue
            @test maximum(@view μL[iel,:]) < 1.0e-6*μL[iel_t,1]
        end

        # C1 is linear in μ_res
        sp2 = build_sphere_params(mesh, metrics,
                                  shell_inputs(nop; :lvisc => true,
                                               :visc_model => DSGS_SW(), :μ => 1.0,
                                               :ivisc_equations => [1,2,3,4],
                                               :dsgs_C1 => 2.0, :dsgs_C2 => 0.5);
                                  neqs = 4)
        μL2 = dsgs(mesh, metrics, sp2, qv, qv, qv, rhsL, Δt)
        @test isapprox(μL2[iel_t,1], 2*μL[iel_t,1]; rtol = 1.0e-10)

        #-------------------------------------------------------------------------
        # (4) the per-equation multipliers of the deck scale the slots.
        #-------------------------------------------------------------------------
        spm = build_sphere_params(mesh, metrics,
                                  shell_inputs(nop; :lvisc => true,
                                               :visc_model => DSGS_SW(),
                                               :μ => [0.0, 1.0, 0.5, 0.25],
                                               :ivisc_equations => 1:4,
                                               :dsgs_C1 => 1.0, :dsgs_C2 => 0.5);
                                  neqs = 4)
        μm = dsgs(mesh, metrics, spm, qv, qv, qv, rhsL, Δt)
        @test μm[iel_t,1] == 0.0
        @test isapprox(μm[iel_t,3], 0.5 *μm[iel_t,2]; rtol = 1.0e-12)
        @test isapprox(μm[iel_t,4], 0.25*μm[iel_t,2]; rtol = 1.0e-12)

        #-------------------------------------------------------------------------
        # (5) SELECTIVITY — the property the model exists for.
        #
        #     Two states, the same everywhere else: one whose φ field is a smooth
        #     l = 2 harmonic (fully resolved at nop = 5), one carrying a sharp
        #     front the grid cannot resolve. Both are advanced through the SAME
        #     BDF2 stencil, so the residual is the discretization error of each.
        #
        #     The residual is formed here from a genuine RHS — the case's own
        #     surface divergence — rather than a hand-written one, so this is the
        #     model responding to real under-resolution.
        #-------------------------------------------------------------------------
        # a params object with the viscous term OFF, so sphere_rhs! returns the
        # plain INVISCID right-hand side. That is what the residual is defined
        # against — feeding it an already-stabilized RHS would make the sensor
        # measure its own output.
        spinv = build_sphere_params(mesh, metrics,
                                    shell_inputs(nop; :lvisc => false);
                                    neqs = 4)
        qe_zero = zeros(Float64, npoin, neqs+1)

        # ∂q/∂t, i.e. M⁻¹·(assembled RHS): what a time integrator consumes.
        function rhs_of(state)
            RHS = zeros(Float64, npoin, neqs)
            Jexpresso.sphere_rhs!(RHS, state, qe_zero,
                                  mesh, metrics, spinv, Jexpresso.TOTAL())
            return RHS
        end

        # ... and the WEAK form of the same thing, which is what
        # compute_dsgs_viscosity! is handed inside sphere_rhs! (the M⁻¹ there
        # happens only after the viscous term is built, and the model applies
        # Minv itself). Getting this wrong is a factor of the mass matrix.
        rhs_weak_of(state) = rhs_of(state) .* metrics.M

        # smooth: φ = φ₀(1 + 0.05 sin(lat)), the l = 1 harmonic. Resolved.
        q_sm = copy(q)
        for ip = 1:npoin
            q_sm[ip,1] = 1.0e5*(1.0 + 0.05*sin(mesh.lat[ip]))
        end
        # sharp: the same amplitude carried by a tanh front one element wide
        # across the equator. Under-resolved by construction.
        Δmin_lat = maximum(Δe)/R
        q_sh = copy(q)
        for ip = 1:npoin
            q_sh[ip,1] = 1.0e5*(1.0 + 0.05*tanh(mesh.lat[ip]/(0.15*Δmin_lat)))
        end

        #
        # TWO Euler steps, so the BDF2 stencil has a genuine three-level history
        # and the residual is what it is supposed to be. With sⁿ⁺¹ = sⁿ + Δt·fⁿ,
        #
        #   (3s² − 4s¹ + s⁰)/(2Δt) − f² = (3f¹ − f⁰)/2 − f²
        #
        # which vanishes to O(Δt²) whenever f varies smoothly from step to step,
        # and does NOT when the solution is carrying grid-scale oscillation. A
        # single Euler step would leave ½·∂q/∂t behind and the test would be
        # measuring how fast the flow evolves instead.
        #
        function bdf2_triple(state)
            s0 = state
            s1 = s0 .+ Δt .* rhs_of(s0)
            s2 = s1 .+ Δt .* rhs_of(s1)
            return s2, s1, s0
        end

        function mu_of(state)
            s2, s1, s0 = bdf2_triple(state)
            return copy(dsgs(mesh, metrics, sp, s2, s1, s0, rhs_weak_of(s2), Δt))
        end

        μ_sm = mu_of(q_sm)
        μ_sh = mu_of(q_sh)

        @info @sprintf("  smooth: mean ν = %.3e, max ν = %.3e m²/s",
                       sum(@view μ_sm[:,2])/nelem, maximum(@view μ_sm[:,2]))
        @info @sprintf("  sharp : mean ν = %.3e, max ν = %.3e m²/s",
                       sum(@view μ_sh[:,2])/nelem, maximum(@view μ_sh[:,2]))

        # The front draws far more viscosity than the smooth field does. The
        # factor measured on this grid is ~2 orders of magnitude; the bound is
        # loose so the test reports a real regression rather than grid noise.
        @test maximum(@view μ_sh[:,2]) > 10*maximum(@view μ_sm[:,2])

        # ... and it is CONCENTRATED, not spread: the elements the front passes
        # through hold most of it, so the mean over the shell stays far below
        # the max. A model that had simply turned itself up everywhere would
        # have mean ≈ max.
        @test sum(@view μ_sh[:,2])/nelem < 0.5*maximum(@view μ_sh[:,2])

        #-------------------------------------------------------------------------
        # (6) REGRESSION: a purely zonal state must not saturate the cap.
        #
        #     φw ≡ 0 for a zonal flow, so its own ‖φw − ⟨φw⟩‖∞ is exactly zero
        #     while its residual is not — the φw equation carries a Coriolis
        #     source the flux divergence cancels only to discretization error.
        #     Normalizing by a floor instead of by the momentum scale inflated
        #     that ratio by ~260× on the Galewsky jet, pinned μ at the cap on
        #     every element, and blew the run up at t ≈ 1100 s. The three
        #     momentum slots therefore share one denominator.
        #-------------------------------------------------------------------------
        # exactly zero third component, as e_λ gives
        @test maximum(abs, @view q[:,4]) == 0.0

        μ_z = mu_of(q)

        @info @sprintf("  zonal : max ν = %.3e m²/s, cap = %.3e m²/s (%.1f%% of cap)",
                       maximum(@view μ_z[:,2]), maximum(cap),
                       100*maximum(@view μ_z[:,2])/maximum(cap))

        for iel = 1:nelem
            @test μ_z[iel,2] < 0.5*cap[iel]
        end

        # the three momentum slots really do share a denominator, so a state
        # ROTATED into a different Cartesian frame gets the same viscosity. This
        # is the invariance the shared denominator buys, and it is what a
        # component-wise normalization does not have.
        @test sp.ddenom[2] == sp.ddenom[3] == sp.ddenom[4]

        #-------------------------------------------------------------------------
        # (7) the deck switches.
        #-------------------------------------------------------------------------
        # the model is a viscosity: :lvisc => false turns it off
        @test !sphere_dsgs_requested(Dict{Symbol,Any}(:lvisc => false,
                                                      :visc_model => DSGS_SW()))
        @test  sphere_dsgs_requested(Dict{Symbol,Any}(:lvisc => true,
                                                      :visc_model => DSGS_SW()))
        # and a deck without :visc_model keeps the constant-ν path
        @test !sphere_dsgs_requested(Dict{Symbol,Any}(:lvisc => true, :μ => 1.0e5))

        # under DynSGS :μ defaults to 1.0, not to 0.0
        μd = build_sphere_viscosity(
            Dict{Symbol,Any}(:lvisc => true, :visc_model => DSGS_SW()), 4)
        @test μd == [0.0, 1.0, 1.0, 1.0]      # :ivisc_equations defaults to 2:neqs

        # a constant-ν value left in :μ is an error, not a 10⁵x overdose
        @test_throws ErrorException build_sphere_viscosity(
            Dict{Symbol,Any}(:lvisc => true, :visc_model => DSGS_SW(),
                             :μ => 1.0e5), 4)
        # ... and the same value without :visc_model is perfectly fine
        @test build_sphere_viscosity(
            Dict{Symbol,Any}(:lvisc => true, :μ => 1.0e5), 4) ==
              [0.0, 1.0e5, 1.0e5, 1.0e5]

        # all-zero multipliers => the residual pass is skipped entirely
        spz = build_sphere_params(mesh, metrics,
                                  shell_inputs(nop; :lvisc => true,
                                               :visc_model => DSGS_SW(),
                                               :μ => 0.0, :ivisc_equations => 1:4);
                                  neqs = 4)
        @test !spz.ldsgs
    end
end
