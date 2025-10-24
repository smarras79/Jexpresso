"""
# Smagorinsky Turbulence Model

Classical Smagorinsky subgrid-scale model for LES.

Eddy viscosity: ν_t = (C_s Δ)² |S|

where:
- C_s = Smagorinsky constant (typically 0.1-0.2)
- Δ = filter width (element size / polynomial order)
- |S| = magnitude of strain rate tensor

Author: Claude Code
Date: 2025-10-24
"""

using LinearAlgebra

"""
    compute_strain_rate_magnitude_2d!(S_mag, u, v, dψ, metrics, mesh, iel)

Compute strain rate magnitude |S| for 2D flow at all quadrature points in element iel.

# Arguments
- `S_mag::Matrix{Float64}`: Output strain rate magnitude [ngl, ngl]
- `u::Matrix{Float64}`: x-velocity at quadrature points [ngl, ngl]
- `v::Matrix{Float64}`: y-velocity at quadrature points [ngl, ngl]
- `dψ::Matrix{Float64}`: Derivative matrix [ngl, ngl]
- `metrics`: Metric terms structure
- `mesh`: Mesh structure
- `iel::Int`: Element index

# Returns
Updates S_mag in place with strain rate magnitude at each point
"""
function compute_strain_rate_magnitude_2d!(S_mag, u, v, dψ, metrics, mesh, iel)
    ngl = mesh.ngl

    # Temporary storage for velocity gradients
    dudξ = zeros(ngl, ngl)
    dudη = zeros(ngl, ngl)
    dvdξ = zeros(ngl, ngl)
    dvdη = zeros(ngl, ngl)

    # Compute derivatives in reference space
    for j = 1:ngl
        for i = 1:ngl
            dudξ[i,j] = 0.0
            dvdξ[i,j] = 0.0
            for k = 1:ngl
                dudξ[i,j] += dψ[k,i] * u[k,j]
                dvdξ[i,j] += dψ[k,i] * v[k,j]
            end

            dudη[i,j] = 0.0
            dvdη[i,j] = 0.0
            for k = 1:ngl
                dudη[i,j] += dψ[k,j] * u[i,k]
                dvdη[i,j] += dψ[k,j] * v[i,k]
            end
        end
    end

    # Transform to physical space and compute strain rate components
    for j = 1:ngl
        for i = 1:ngl
            # Physical space gradients
            dudx = dudξ[i,j] * metrics.dξdx[iel,i,j] + dudη[i,j] * metrics.dηdx[iel,i,j]
            dudy = dudξ[i,j] * metrics.dξdy[iel,i,j] + dudη[i,j] * metrics.dηdy[iel,i,j]
            dvdx = dvdξ[i,j] * metrics.dξdx[iel,i,j] + dvdη[i,j] * metrics.dηdx[iel,i,j]
            dvdy = dvdξ[i,j] * metrics.dξdy[iel,i,j] + dvdη[i,j] * metrics.dηdy[iel,i,j]

            # Divergence
            div_u = dudx + dvdy

            # Strain rate tensor (deviatoric part)
            S_xx = dudx - div_u / 3.0
            S_yy = dvdy - div_u / 3.0
            S_xy = 0.5 * (dudy + dvdx)

            # Magnitude |S| = sqrt(2 S_ij S_ij)
            S_mag[i,j] = sqrt(2.0 * (S_xx^2 + S_yy^2 + 2.0*S_xy^2))
        end
    end

    return nothing
end

"""
    compute_strain_rate_magnitude_3d!(S_mag, u, v, w, dψ, metrics, mesh, iel)

Compute strain rate magnitude |S| for 3D flow at all quadrature points in element iel.
"""
function compute_strain_rate_magnitude_3d!(S_mag, u, v, w, dψ, metrics, mesh, iel)
    ngl = mesh.ngl

    # Temporary storage for velocity gradients
    dudξ = zeros(ngl, ngl, ngl)
    dudη = zeros(ngl, ngl, ngl)
    dudζ = zeros(ngl, ngl, ngl)
    dvdξ = zeros(ngl, ngl, ngl)
    dvdη = zeros(ngl, ngl, ngl)
    dvdζ = zeros(ngl, ngl, ngl)
    dwdξ = zeros(ngl, ngl, ngl)
    dwdη = zeros(ngl, ngl, ngl)
    dwdζ = zeros(ngl, ngl, ngl)

    # Compute derivatives in reference space
    for k = 1:ngl
        for j = 1:ngl
            for i = 1:ngl
                dudξ[i,j,k] = 0.0
                dvdξ[i,j,k] = 0.0
                dwdξ[i,j,k] = 0.0
                for m = 1:ngl
                    dudξ[i,j,k] += dψ[m,i] * u[m,j,k]
                    dvdξ[i,j,k] += dψ[m,i] * v[m,j,k]
                    dwdξ[i,j,k] += dψ[m,i] * w[m,j,k]
                end

                dudη[i,j,k] = 0.0
                dvdη[i,j,k] = 0.0
                dwdη[i,j,k] = 0.0
                for m = 1:ngl
                    dudη[i,j,k] += dψ[m,j] * u[i,m,k]
                    dvdη[i,j,k] += dψ[m,j] * v[i,m,k]
                    dwdη[i,j,k] += dψ[m,j] * w[i,m,k]
                end

                dudζ[i,j,k] = 0.0
                dvdζ[i,j,k] = 0.0
                dwdζ[i,j,k] = 0.0
                for m = 1:ngl
                    dudζ[i,j,k] += dψ[m,k] * u[i,j,m]
                    dvdζ[i,j,k] += dψ[m,k] * v[i,j,m]
                    dwdζ[i,j,k] += dψ[m,k] * w[i,j,m]
                end
            end
        end
    end

    # Transform to physical space and compute strain rate components
    for kk = 1:ngl
        for j = 1:ngl
            for i = 1:ngl
                # Physical space gradients
                dudx = dudξ[i,j,kk]*metrics.dξdx[iel,i,j,kk] + dudη[i,j,kk]*metrics.dηdx[iel,i,j,kk] + dudζ[i,j,kk]*metrics.dζdx[iel,i,j,kk]
                dudy = dudξ[i,j,kk]*metrics.dξdy[iel,i,j,kk] + dudη[i,j,kk]*metrics.dηdy[iel,i,j,kk] + dudζ[i,j,kk]*metrics.dζdy[iel,i,j,kk]
                dudz = dudξ[i,j,kk]*metrics.dξdz[iel,i,j,kk] + dudη[i,j,kk]*metrics.dηdz[iel,i,j,kk] + dudζ[i,j,kk]*metrics.dζdz[iel,i,j,kk]

                dvdx = dvdξ[i,j,kk]*metrics.dξdx[iel,i,j,kk] + dvdη[i,j,kk]*metrics.dηdx[iel,i,j,kk] + dvdζ[i,j,kk]*metrics.dζdx[iel,i,j,kk]
                dvdy = dvdξ[i,j,kk]*metrics.dξdy[iel,i,j,kk] + dvdη[i,j,kk]*metrics.dηdy[iel,i,j,kk] + dvdζ[i,j,kk]*metrics.dζdy[iel,i,j,kk]
                dvdz = dvdξ[i,j,kk]*metrics.dξdz[iel,i,j,kk] + dvdη[i,j,kk]*metrics.dηdz[iel,i,j,kk] + dvdζ[i,j,kk]*metrics.dζdz[iel,i,j,kk]

                dwdx = dwdξ[i,j,kk]*metrics.dξdx[iel,i,j,kk] + dwdη[i,j,kk]*metrics.dηdx[iel,i,j,kk] + dwdζ[i,j,kk]*metrics.dζdx[iel,i,j,kk]
                dwdy = dwdξ[i,j,kk]*metrics.dξdy[iel,i,j,kk] + dwdη[i,j,kk]*metrics.dηdy[iel,i,j,kk] + dwdζ[i,j,kk]*metrics.dζdy[iel,i,j,kk]
                dwdz = dwdξ[i,j,kk]*metrics.dξdz[iel,i,j,kk] + dwdη[i,j,kk]*metrics.dηdz[iel,i,j,kk] + dwdζ[i,j,kk]*metrics.dζdz[iel,i,j,kk]

                # Divergence
                div_u = dudx + dvdy + dwdz

                # Strain rate tensor (deviatoric part)
                S_xx = dudx - div_u / 3.0
                S_yy = dvdy - div_u / 3.0
                S_zz = dwdz - div_u / 3.0
                S_xy = 0.5 * (dudy + dvdx)
                S_xz = 0.5 * (dudz + dwdx)
                S_yz = 0.5 * (dvdz + dwdy)

                # Magnitude |S| = sqrt(2 S_ij S_ij)
                S_mag[i,j,kk] = sqrt(2.0 * (S_xx^2 + S_yy^2 + S_zz^2 +
                                           2.0*(S_xy^2 + S_xz^2 + S_yz^2)))
            end
        end
    end

    return nothing
end

"""
    compute_element_filter_width(mesh, iel, ngl)

Compute filter width Δ for element iel.

For spectral elements: Δ = h_e / ngl
where h_e is characteristic element size.
"""
function compute_element_filter_width(mesh, iel, ngl)
    # Simple estimate: use maximum element dimension
    # For more accurate version, use element volume^(1/d)
    if hasfield(typeof(mesh), :Δx)
        return mesh.Δx[iel] / ngl
    else
        # Fallback: compute from coordinates
        # This is a simplified version - should be improved for general elements
        return 1.0 / ngl  # Placeholder
    end
end

"""
    compute_smagorinsky_viscosity_2d!(ν_t, μ_t, κ_t, q, dψ, metrics, mesh, iel, C_s, Pr_t, PhysConst)

Compute Smagorinsky eddy viscosity and thermal diffusivity for 2D flow.

# Arguments
- `ν_t::Matrix{Float64}`: Output turbulent kinematic viscosity [ngl, ngl]
- `μ_t::Matrix{Float64}`: Output turbulent dynamic viscosity [ngl, ngl]
- `κ_t::Matrix{Float64}`: Output turbulent thermal diffusivity [ngl, ngl]
- `q::Matrix{Float64}`: State vector [ngl*ngl, neqs] for element
- `dψ::Matrix{Float64}`: Derivative matrix
- `metrics`: Metric terms
- `mesh`: Mesh structure
- `iel::Int`: Element index
- `C_s::Float64`: Smagorinsky constant (typically 0.1-0.2)
- `Pr_t::Float64`: Turbulent Prandtl number (typically 0.7-1.0)
- `PhysConst`: Physical constants structure
"""
function compute_smagorinsky_viscosity_2d!(ν_t, μ_t, κ_t, q, dψ, metrics, mesh, iel,
                                          C_s, Pr_t, PhysConst)
    ngl = mesh.ngl

    # Extract velocity field at quadrature points
    u_el = zeros(ngl, ngl)
    v_el = zeros(ngl, ngl)
    ρ_el = zeros(ngl, ngl)

    for j = 1:ngl
        for i = 1:ngl
            ip = mesh.connijk[i, j, iel]
            ρ_el[i,j] = q[ip, 1]
            u_el[i,j] = q[ip, 2] / ρ_el[i,j]
            v_el[i,j] = q[ip, 3] / ρ_el[i,j]
        end
    end

    # Compute strain rate magnitude
    S_mag = zeros(ngl, ngl)
    compute_strain_rate_magnitude_2d!(S_mag, u_el, v_el, dψ, metrics, mesh, iel)

    # Filter width
    Δ = compute_element_filter_width(mesh, iel, ngl)

    # Compute eddy viscosity
    for j = 1:ngl
        for i = 1:ngl
            ν_t[i,j] = (C_s * Δ)^2 * S_mag[i,j]
            μ_t[i,j] = ρ_el[i,j] * ν_t[i,j]

            # Thermal eddy diffusivity (using c_p from PhysConst)
            κ_t[i,j] = μ_t[i,j] * PhysConst.cp / Pr_t
        end
    end

    return nothing
end

"""
    compute_smagorinsky_viscosity_3d!(ν_t, μ_t, κ_t, q, dψ, metrics, mesh, iel, C_s, Pr_t, PhysConst)

Compute Smagorinsky eddy viscosity and thermal diffusivity for 3D flow.
"""
function compute_smagorinsky_viscosity_3d!(ν_t, μ_t, κ_t, q, dψ, metrics, mesh, iel,
                                          C_s, Pr_t, PhysConst)
    ngl = mesh.ngl

    # Extract velocity field at quadrature points
    u_el = zeros(ngl, ngl, ngl)
    v_el = zeros(ngl, ngl, ngl)
    w_el = zeros(ngl, ngl, ngl)
    ρ_el = zeros(ngl, ngl, ngl)

    for k = 1:ngl
        for j = 1:ngl
            for i = 1:ngl
                ip = mesh.connijk[i, j, k, iel]
                ρ_el[i,j,k] = q[ip, 1]
                u_el[i,j,k] = q[ip, 2] / ρ_el[i,j,k]
                v_el[i,j,k] = q[ip, 3] / ρ_el[i,j,k]
                w_el[i,j,k] = q[ip, 4] / ρ_el[i,j,k]
            end
        end
    end

    # Compute strain rate magnitude
    S_mag = zeros(ngl, ngl, ngl)
    compute_strain_rate_magnitude_3d!(S_mag, u_el, v_el, w_el, dψ, metrics, mesh, iel)

    # Filter width
    Δ = compute_element_filter_width(mesh, iel, ngl)

    # Compute eddy viscosity
    for k = 1:ngl
        for j = 1:ngl
            for i = 1:ngl
                ν_t[i,j,k] = (C_s * Δ)^2 * S_mag[i,j,k]
                μ_t[i,j,k] = ρ_el[i,j,k] * ν_t[i,j,k]

                # Thermal eddy diffusivity
                κ_t[i,j,k] = μ_t[i,j,k] * PhysConst.cp / Pr_t
            end
        end
    end

    return nothing
end
