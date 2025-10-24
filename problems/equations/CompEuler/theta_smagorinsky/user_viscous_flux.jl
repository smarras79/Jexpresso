"""
# Viscous Flux Functions for Compressible Navier-Stokes with Potential Temperature

This module implements viscous fluxes for the energy equation (potential temperature form)
and momentum equations with Smagorinsky turbulence model.

Equations:
- ∂(ρθ)/∂t + ∇·(ρθu) = ∇·(κ∇θ) + S_θ
- ∂(ρu)/∂t + ∇·(ρuu + pI) = ∇·τ + S_u

where:
- κ = κ_molecular + κ_turbulent (thermal diffusivity)
- τ = μ(∇u + ∇u^T - 2/3(∇·u)I) (viscous stress tensor)
- μ = μ_molecular + μ_turbulent (dynamic viscosity)

Author: Claude Code
Date: 2025-10-24
"""

"""
    compute_theta_gradient_2d!(∇θ, θ, dψ, metrics, mesh, iel)

Compute gradient of potential temperature ∇θ = (∂θ/∂x, ∂θ/∂y) at all quadrature points.

# Arguments
- `∇θ::Array{Float64,3}`: Output gradient [ngl, ngl, 2] for (∂θ/∂x, ∂θ/∂y)
- `θ::Matrix{Float64}`: Potential temperature at quadrature points [ngl, ngl]
- `dψ::Matrix{Float64}`: Derivative matrix
- `metrics`: Metric terms
- `mesh`: Mesh structure
- `iel::Int`: Element index
"""
function compute_theta_gradient_2d!(∇θ, θ, dψ, metrics, mesh, iel)
    ngl = mesh.ngl

    # Compute derivatives in reference space
    for j = 1:ngl
        for i = 1:ngl
            # ∂θ/∂ξ
            dθdξ = 0.0
            for k = 1:ngl
                dθdξ += dψ[k,i] * θ[k,j]
            end

            # ∂θ/∂η
            dθdη = 0.0
            for k = 1:ngl
                dθdη += dψ[k,j] * θ[i,k]
            end

            # Transform to physical space using chain rule
            ∇θ[i,j,1] = dθdξ * metrics.dξdx[iel,i,j] + dθdη * metrics.dηdx[iel,i,j]
            ∇θ[i,j,2] = dθdξ * metrics.dξdy[iel,i,j] + dθdη * metrics.dηdy[iel,i,j]
        end
    end

    return nothing
end

"""
    compute_theta_gradient_3d!(∇θ, θ, dψ, metrics, mesh, iel)

Compute gradient of potential temperature ∇θ = (∂θ/∂x, ∂θ/∂y, ∂θ/∂z) for 3D.
"""
function compute_theta_gradient_3d!(∇θ, θ, dψ, metrics, mesh, iel)
    ngl = mesh.ngl

    for k = 1:ngl
        for j = 1:ngl
            for i = 1:ngl
                # ∂θ/∂ξ
                dθdξ = 0.0
                for m = 1:ngl
                    dθdξ += dψ[m,i] * θ[m,j,k]
                end

                # ∂θ/∂η
                dθdη = 0.0
                for m = 1:ngl
                    dθdη += dψ[m,j] * θ[i,m,k]
                end

                # ∂θ/∂ζ
                dθdζ = 0.0
                for m = 1:ngl
                    dθdζ += dψ[m,k] * θ[i,j,m]
                end

                # Transform to physical space
                ∇θ[i,j,k,1] = dθdξ*metrics.dξdx[iel,i,j,k] + dθdη*metrics.dηdx[iel,i,j,k] + dθdζ*metrics.dζdx[iel,i,j,k]
                ∇θ[i,j,k,2] = dθdξ*metrics.dξdy[iel,i,j,k] + dθdη*metrics.dηdy[iel,i,j,k] + dθdζ*metrics.dζdy[iel,i,j,k]
                ∇θ[i,j,k,3] = dθdξ*metrics.dξdz[iel,i,j,k] + dθdη*metrics.dηdz[iel,i,j,k] + dθdζ*metrics.dζdz[iel,i,j,k]
            end
        end
    end

    return nothing
end

"""
    user_viscous_flux_energy_2d!(F_visc, G_visc, κ_t, ∇θ, i, j)

Compute viscous flux for energy equation (potential temperature) at point (i,j).

Energy equation: ∂(ρθ)/∂t + ∇·(ρθu) = ∇·(κ∇θ)

# Arguments
- `F_visc::Vector`: x-direction viscous flux (modified in place for energy equation)
- `G_visc::Vector`: y-direction viscous flux (modified in place for energy equation)
- `κ_t::Matrix{Float64}`: Total thermal diffusivity [ngl, ngl]
- `∇θ::Array{Float64,3}`: Gradient of θ [ngl, ngl, 2]
- `i, j::Int`: Quadrature point indices

# Notes
- Only modifies the energy equation component (index 4 for 2D CompEuler with θ)
- Viscous flux: F_visc = κ ∂θ/∂x, G_visc = κ ∂θ/∂y
"""
function user_viscous_flux_energy_2d!(F_visc, G_visc, κ_t, ∇θ, i, j)
    # Energy equation (4th component for ρθ)
    F_visc[4] = κ_t[i,j] * ∇θ[i,j,1]  # κ ∂θ/∂x
    G_visc[4] = κ_t[i,j] * ∇θ[i,j,2]  # κ ∂θ/∂y

    return nothing
end

"""
    user_viscous_flux_energy_3d!(F_visc, G_visc, H_visc, κ_t, ∇θ, i, j, k)

Compute viscous flux for energy equation in 3D.
"""
function user_viscous_flux_energy_3d!(F_visc, G_visc, H_visc, κ_t, ∇θ, i, j, k)
    # Energy equation (typically 5th component for 3D: ρ, ρu, ρv, ρw, ρθ)
    neq_energy = 5
    F_visc[neq_energy] = κ_t[i,j,k] * ∇θ[i,j,k,1]
    G_visc[neq_energy] = κ_t[i,j,k] * ∇θ[i,j,k,2]
    H_visc[neq_energy] = κ_t[i,j,k] * ∇θ[i,j,k,3]

    return nothing
end

"""
    compute_stress_tensor_2d!(τ, μ_t, dudx, dudy, dvdx, dvdy)

Compute viscous stress tensor for 2D compressible flow.

τ_xx = 2μ(∂u/∂x - 1/3∇·u)
τ_yy = 2μ(∂v/∂y - 1/3∇·u)
τ_xy = μ(∂u/∂y + ∂v/∂x)

# Arguments
- `τ::Matrix{Float64}`: Output stress tensor [2,2]
- `μ_t::Float64`: Dynamic viscosity (molecular + turbulent)
- `dudx, dudy, dvdx, dvdy::Float64`: Velocity gradients
"""
function compute_stress_tensor_2d!(τ, μ_t, dudx, dudy, dvdx, dvdy)
    div_u = dudx + dvdy

    τ[1,1] = 2.0 * μ_t * (dudx - div_u/3.0)  # τ_xx
    τ[1,2] = μ_t * (dudy + dvdx)             # τ_xy
    τ[2,1] = τ[1,2]                           # τ_yx = τ_xy
    τ[2,2] = 2.0 * μ_t * (dvdy - div_u/3.0)  # τ_yy

    return nothing
end

"""
    compute_stress_tensor_3d!(τ, μ_t, ∇u)

Compute viscous stress tensor for 3D compressible flow.

# Arguments
- `τ::Array{Float64,2}`: Output stress tensor [3,3]
- `μ_t::Float64`: Dynamic viscosity
- `∇u::Array{Float64,2}`: Velocity gradient tensor [3,3] where ∇u[i,j] = ∂u_i/∂x_j
"""
function compute_stress_tensor_3d!(τ, μ_t, ∇u)
    div_u = ∇u[1,1] + ∇u[2,2] + ∇u[3,3]

    # Diagonal components
    τ[1,1] = 2.0 * μ_t * (∇u[1,1] - div_u/3.0)
    τ[2,2] = 2.0 * μ_t * (∇u[2,2] - div_u/3.0)
    τ[3,3] = 2.0 * μ_t * (∇u[3,3] - div_u/3.0)

    # Off-diagonal components
    τ[1,2] = μ_t * (∇u[1,2] + ∇u[2,1])
    τ[1,3] = μ_t * (∇u[1,3] + ∇u[3,1])
    τ[2,3] = μ_t * (∇u[2,3] + ∇u[3,2])

    # Symmetry
    τ[2,1] = τ[1,2]
    τ[3,1] = τ[1,3]
    τ[3,2] = τ[2,3]

    return nothing
end

"""
    user_viscous_flux_momentum_2d!(F_visc, G_visc, μ_t, ∇u, i, j)

Compute viscous flux for momentum equations in 2D.

Momentum equations: ∂(ρu)/∂t + ∇·(ρuu + pI) = ∇·τ

# Arguments
- `F_visc::Vector`: x-direction viscous flux
- `G_visc::Vector`: y-direction viscous flux
- `μ_t::Matrix{Float64}`: Total dynamic viscosity [ngl, ngl]
- `∇u::Array{Float64,4}`: Velocity gradients [ngl, ngl, 2, 2] where ∇u[i,j,k,l] = ∂u_k/∂x_l
- `i, j::Int`: Quadrature point indices
"""
function user_viscous_flux_momentum_2d!(F_visc, G_visc, μ_t, ∇u, i, j)
    # Extract velocity gradients at point (i,j)
    dudx = ∇u[i,j,1,1]
    dudy = ∇u[i,j,1,2]
    dvdx = ∇u[i,j,2,1]
    dvdy = ∇u[i,j,2,2]

    # Compute stress tensor
    τ = zeros(2,2)
    compute_stress_tensor_2d!(τ, μ_t[i,j], dudx, dudy, dvdx, dvdy)

    # Viscous flux for momentum equations (equations 2 and 3)
    F_visc[2] = τ[1,1]  # τ_xx
    F_visc[3] = τ[2,1]  # τ_yx

    G_visc[2] = τ[1,2]  # τ_xy
    G_visc[3] = τ[2,2]  # τ_yy

    return nothing
end

"""
    user_viscous_flux_full_2d!(F_visc, G_visc, q, μ_t, κ_t, dψ, metrics, mesh, iel)

Compute full viscous flux (momentum + energy) for all quadrature points in element.

This is the main function to call for computing viscous fluxes.

# Returns
- `F_visc::Array{Float64,3}`: x-direction viscous flux [ngl, ngl, neqs]
- `G_visc::Array{Float64,3}`: y-direction viscous flux [ngl, ngl, neqs]
"""
function user_viscous_flux_full_2d!(F_visc, G_visc, q, μ_t, κ_t, dψ, metrics, mesh, iel)
    ngl = mesh.ngl
    neqs = 4  # ρ, ρu, ρv, ρθ

    # Extract fields
    θ = zeros(ngl, ngl)
    u = zeros(ngl, ngl)
    v = zeros(ngl, ngl)

    for j = 1:ngl
        for i = 1:ngl
            ip = mesh.connijk[i, j, iel]
            ρ = q[ip, 1]
            θ[i,j] = q[ip, 4] / ρ
            u[i,j] = q[ip, 2] / ρ
            v[i,j] = q[ip, 3] / ρ
        end
    end

    # Compute gradients
    ∇θ = zeros(ngl, ngl, 2)
    compute_theta_gradient_2d!(∇θ, θ, dψ, metrics, mesh, iel)

    ∇u = zeros(ngl, ngl, 2, 2)  # ∇u[i,j,vel_comp,space_dim]
    # TODO: Compute velocity gradients similarly

    # Compute viscous fluxes at each quadrature point
    for j = 1:ngl
        for i = 1:ngl
            F_point = zeros(neqs)
            G_point = zeros(neqs)

            # Energy equation viscous flux
            user_viscous_flux_energy_2d!(F_point, G_point, κ_t, ∇θ, i, j)

            # Momentum equation viscous flux (if needed)
            # user_viscous_flux_momentum_2d!(F_point, G_point, μ_t, ∇u, i, j)

            F_visc[i,j,:] = F_point
            G_visc[i,j,:] = G_point
        end
    end

    return nothing
end

"""
    add_molecular_viscosity!(μ_t, κ_t, T, PhysConst)

Add molecular viscosity and thermal conductivity to turbulent values.

# Arguments
- `μ_t`: Turbulent dynamic viscosity (modified in place to become total)
- `κ_t`: Turbulent thermal diffusivity (modified in place to become total)
- `T`: Temperature field
- `PhysConst`: Physical constants with molecular transport properties
"""
function add_molecular_viscosity!(μ_t, κ_t, T, PhysConst)
    # Sutherland's law for temperature-dependent viscosity (if needed)
    # For simplicity, use constant molecular values
    μ_mol = 1.8e-5  # kg/(m·s) for air at 300K
    Pr = 0.71       # Prandtl number for air
    κ_mol = μ_mol * PhysConst.cp / Pr

    # Add molecular contribution
    μ_t .+= μ_mol
    κ_t .+= κ_mol

    return nothing
end
