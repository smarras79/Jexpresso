"""
# Example: Integrating Smagorinsky Model with Spectral Element RHS

This example demonstrates how to integrate the Smagorinsky turbulence model
and viscous fluxes into the Jexpresso RHS computation framework.

Author: Claude Code
Date: 2025-10-24
"""

# Include the new modules
include("../src/kernel/Turbulence/Smagorinsky.jl")
include("../problems/equations/CompEuler/theta/user_viscous_flux.jl")

"""
    build_rhs_with_smagorinsky_2d!(RHS, q, params, time)

Extended RHS computation including Smagorinsky viscous terms for 2D.

This function augments the existing inviscid RHS with viscous contributions.

# Workflow:
1. Compute inviscid RHS (existing code)
2. Compute Smagorinsky eddy viscosity for each element
3. Compute viscous fluxes
4. Add viscous contribution to RHS using weak form
"""
function build_rhs_with_smagorinsky_2d!(RHS, q, params, time)

    # -------------------------------------------------------------------
    # Step 1: Compute inviscid RHS (existing Jexpresso functionality)
    # -------------------------------------------------------------------
    # This is already handled by _build_rhs! in rhs.jl
    # For this example, assume it's already computed and stored in RHS

    # -------------------------------------------------------------------
    # Step 2: Initialize viscous parameters
    # -------------------------------------------------------------------
    ngl = params.mesh.ngl
    neqs = params.neqs
    nelem = params.mesh.nelem

    # Smagorinsky parameters
    C_s = 0.18    # Smagorinsky constant (typical for atmospheric flows)
    Pr_t = 0.85   # Turbulent Prandtl number

    PhysConst = PhysicalConst{Float64}()

    # -------------------------------------------------------------------
    # Step 3: Loop over elements and compute viscous contributions
    # -------------------------------------------------------------------
    for iel = 1:nelem

        # Allocate element-wise arrays
        ν_t = zeros(ngl, ngl)      # Turbulent kinematic viscosity
        μ_t = zeros(ngl, ngl)      # Turbulent dynamic viscosity
        κ_t = zeros(ngl, ngl)      # Turbulent thermal diffusivity

        # Compute Smagorinsky eddy viscosity for this element
        compute_smagorinsky_viscosity_2d!(ν_t, μ_t, κ_t,
                                         q, params.basis.dψ,
                                         params.metrics, params.mesh, iel,
                                         C_s, Pr_t, PhysConst)

        # Add molecular viscosity
        θ_el = zeros(ngl, ngl)
        for j = 1:ngl
            for i = 1:ngl
                ip = params.mesh.connijk[i, j, iel]
                θ_el[i,j] = q[ip, 4] / q[ip, 1]  # θ = ρθ / ρ
            end
        end
        add_molecular_viscosity!(μ_t, κ_t, θ_el, PhysConst)

        # Compute viscous fluxes
        F_visc = zeros(ngl, ngl, neqs)
        G_visc = zeros(ngl, ngl, neqs)

        user_viscous_flux_full_2d!(F_visc, G_visc, q, μ_t, κ_t,
                                  params.basis.dψ, params.metrics,
                                  params.mesh, iel)

        # -------------------------------------------------------------------
        # Step 4: Add viscous contribution to RHS using weak form
        # -------------------------------------------------------------------
        # Weak form: ∫_Ω ∇φ · F_visc dΩ
        #
        # Discrete form:
        # R_visc[kl] = ∑_{ij} ω_i ω_j J_{ij} [dψ_k/dξ|_i ψ_l(η_j) F_ξ_{ij} +
        #                                       ψ_k(ξ_i) dψ_l/dη|_j F_η_{ij}]

        for eq = 1:neqs
            for l = 1:ngl
                for k = 1:ngl
                    ip = params.mesh.connijk[k, l, iel]

                    visc_contrib = 0.0

                    for j = 1:ngl
                        for i = 1:ngl
                            # Jacobian and weights
                            ωJac = params.ω[i] * params.ω[j] * params.metrics.Je[iel,i,j]

                            # Transform viscous flux to reference space
                            F_ξ = (F_visc[i,j,eq] * params.metrics.dξdx[iel,i,j] +
                                   G_visc[i,j,eq] * params.metrics.dξdy[iel,i,j])

                            F_η = (F_visc[i,j,eq] * params.metrics.dηdx[iel,i,j] +
                                   G_visc[i,j,eq] * params.metrics.dηdy[iel,i,j])

                            # Weak form integral: ∫ ∇φ · F dΩ
                            # Uses derivative of test function
                            visc_contrib += ωJac * (params.basis.dψ[i,k] * params.basis.ψ[j,l] * F_ξ +
                                                   params.basis.ψ[i,k] * params.basis.dψ[j,l] * F_η)
                        end
                    end

                    # Add to RHS (note: sign convention depends on formulation)
                    # Typically: dq/dt = -∇·F_inv + ∇·F_visc
                    # In weak form: -∫∇φ·F_inv + ∫∇φ·F_visc
                    RHS[ip, eq] += visc_contrib * params.Minv[k,l]

                end
            end
        end

    end  # End loop over elements

    return nothing
end

"""
    build_rhs_with_smagorinsky_3d!(RHS, q, params, time)

Extended RHS computation including Smagorinsky viscous terms for 3D.
"""
function build_rhs_with_smagorinsky_3d!(RHS, q, params, time)

    ngl = params.mesh.ngl
    neqs = params.neqs
    nelem = params.mesh.nelem

    C_s = 0.18
    Pr_t = 0.85
    PhysConst = PhysicalConst{Float64}()

    for iel = 1:nelem

        # Allocate element-wise arrays (3D)
        ν_t = zeros(ngl, ngl, ngl)
        μ_t = zeros(ngl, ngl, ngl)
        κ_t = zeros(ngl, ngl, ngl)

        # Compute Smagorinsky viscosity
        compute_smagorinsky_viscosity_3d!(ν_t, μ_t, κ_t,
                                         q, params.basis.dψ,
                                         params.metrics, params.mesh, iel,
                                         C_s, Pr_t, PhysConst)

        # Compute viscous fluxes (3D version)
        F_visc = zeros(ngl, ngl, ngl, neqs)
        G_visc = zeros(ngl, ngl, ngl, neqs)
        H_visc = zeros(ngl, ngl, ngl, neqs)

        # TODO: Implement user_viscous_flux_full_3d!
        # user_viscous_flux_full_3d!(F_visc, G_visc, H_visc, ...)

        # Add to RHS using 3D weak form
        for eq = 1:neqs
            for m = 1:ngl
                for l = 1:ngl
                    for k = 1:ngl
                        ip = params.mesh.connijk[k, l, m, iel]

                        visc_contrib = 0.0

                        for kk = 1:ngl
                            for j = 1:ngl
                                for i = 1:ngl
                                    ωJac = (params.ω[i] * params.ω[j] * params.ω[kk] *
                                           params.metrics.Je[iel,i,j,kk])

                                    # Transform to reference space
                                    F_ξ = (F_visc[i,j,kk,eq] * params.metrics.dξdx[iel,i,j,kk] +
                                           G_visc[i,j,kk,eq] * params.metrics.dξdy[iel,i,j,kk] +
                                           H_visc[i,j,kk,eq] * params.metrics.dξdz[iel,i,j,kk])

                                    F_η = (F_visc[i,j,kk,eq] * params.metrics.dηdx[iel,i,j,kk] +
                                           G_visc[i,j,kk,eq] * params.metrics.dηdy[iel,i,j,kk] +
                                           H_visc[i,j,kk,eq] * params.metrics.dηdz[iel,i,j,kk])

                                    F_ζ = (F_visc[i,j,kk,eq] * params.metrics.dζdx[iel,i,j,kk] +
                                           G_visc[i,j,kk,eq] * params.metrics.dζdy[iel,i,j,kk] +
                                           H_visc[i,j,kk,eq] * params.metrics.dζdz[iel,i,j,kk])

                                    visc_contrib += ωJac * (
                                        params.basis.dψ[i,k] * params.basis.ψ[j,l] * params.basis.ψ[kk,m] * F_ξ +
                                        params.basis.ψ[i,k] * params.basis.dψ[j,l] * params.basis.ψ[kk,m] * F_η +
                                        params.basis.ψ[i,k] * params.basis.ψ[j,l] * params.basis.dψ[kk,m] * F_ζ
                                    )
                                end
                            end
                        end

                        RHS[ip, eq] += visc_contrib * params.Minv[k,l,m]
                    end
                end
            end
        end
    end

    return nothing
end

"""
    Example usage in main driver
"""
function example_main()

    println("=" ^ 70)
    println("Spectral Element Compressible Navier-Stokes with Smagorinsky Model")
    println("=" ^ 70)

    # Setup problem (this would use existing Jexpresso initialization)
    # inputs = ...
    # mesh = ...
    # params = ...
    # q = ...  # Initial condition

    # Modified RHS function that includes Smagorinsky
    function rhs_with_smagorinsky!(du, u, params, time)

        # First compute inviscid part (existing code)
        # build_rhs!(params.RHS, u, params, time)

        # Then add Smagorinsky viscous terms
        if params.SD == NSD_2D()
            build_rhs_with_smagorinsky_2d!(params.RHS, params.uaux, params, time)
        elseif params.SD == NSD_3D()
            build_rhs_with_smagorinsky_3d!(params.RHS, params.uaux, params, time)
        end

        # Convert RHS to du format
        # RHStoDU!(du, params.RHS, params.neqs, params.mesh.npoin)
    end

    # Setup ODE problem with modified RHS
    # tspan = (0.0, T_final)
    # prob = ODEProblem(rhs_with_smagorinsky!, q0, tspan, params)

    # Solve using existing time integrators
    # sol = solve(prob, SSPRK43(), ...)

    println("Simulation with Smagorinsky LES complete!")

end

# -------------------------------------------------------------------
# Diagnostic Functions
# -------------------------------------------------------------------

"""
    compute_turbulent_diagnostics(q, ν_t, mesh)

Compute diagnostic quantities for turbulence analysis.
"""
function compute_turbulent_diagnostics(q, ν_t, mesh)

    # Turbulent kinetic energy
    TKE = 0.0

    # Eddy viscosity statistics
    ν_t_mean = mean(ν_t)
    ν_t_max = maximum(ν_t)
    ν_t_min = minimum(ν_t)

    # Effective Reynolds number
    # Re_eff = U L / (ν_mol + ν_t)

    println("Turbulent Diagnostics:")
    println("  Mean eddy viscosity: ", ν_t_mean)
    println("  Max eddy viscosity:  ", ν_t_max)
    println("  Min eddy viscosity:  ", ν_t_min)

    return Dict(
        "nu_t_mean" => ν_t_mean,
        "nu_t_max" => ν_t_max,
        "nu_t_min" => ν_t_min
    )
end

"""
    validate_energy_conservation(q_old, q_new, RHS_inv, RHS_visc, dt, mesh)

Check energy conservation for debugging.
"""
function validate_energy_conservation(q_old, q_new, RHS_inv, RHS_visc, dt, mesh)

    # Total energy at old and new time
    E_old = compute_total_energy(q_old, mesh)
    E_new = compute_total_energy(q_new, mesh)

    # Energy change
    dE = E_new - E_old

    # Expected change from RHS
    dE_expected = dt * (sum(RHS_inv[:, end]) + sum(RHS_visc[:, end]))

    # Relative error
    rel_error = abs(dE - dE_expected) / abs(dE + 1e-16)

    if rel_error > 1e-6
        @warn "Energy conservation check failed!" dE dE_expected rel_error
    end

    return rel_error
end

"""
    compute_total_energy(q, mesh)

Compute total energy in domain.
"""
function compute_total_energy(q, mesh)
    E_total = 0.0
    for ip = 1:mesh.npoin
        ρ = q[ip, 1]
        u = q[ip, 2] / ρ
        v = q[ip, 3] / ρ
        θ = q[ip, 4] / ρ
        # Add kinetic energy and internal energy contributions
        E_total += 0.5 * ρ * (u^2 + v^2) + ρ * θ  # Simplified
    end
    return E_total
end

println("Smagorinsky integration example loaded successfully!")
println("Call example_main() to see usage, or integrate into existing Jexpresso drivers.")
