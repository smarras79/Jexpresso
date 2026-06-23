"""
        Kinetic energy preserving two-point flux by Shima et al.
        Adapted for 2D compressible Euler equations
        """
@inline function flux_shima_etal_2d(u_ll, u_rr, uprim_ll, uprim_rr, 
                                    orientation::Integer, γ)
    # Unpack primitive variables [ρ, ρu, ρv, ρE] → [ρ, u, v, p]
    ρ_ll, u_ll_vel, v_ll_vel, p_ll = uprim_ll
    ρ_rr, u_rr_vel, v_rr_vel, p_rr = uprim_rr
    
    # Average each factor of products in flux
    ρ_avg = 0.5 * (ρ_ll + ρ_rr)
    u_avg = 0.5 * (u_ll_vel + u_rr_vel)
    v_avg = 0.5 * (v_ll_vel + v_rr_vel)
    p_avg = 0.5 * (p_ll + p_rr)
    
    # Average of products (for kinetic energy preservation)
    kin_avg = 0.5 * (u_ll_vel * u_rr_vel + v_ll_vel * v_rr_vel)
    
    inv_gamma_minus_one = 1.0 / (γ - 1.0)
    
    if orientation == 1  # x-direction
        pv_avg = 0.5 * (p_ll * u_rr_vel + p_rr * u_ll_vel)
        
        f1 = ρ_avg * u_avg
        f2 = f1 * u_avg + p_avg
        f3 = f1 * v_avg
        f4 = p_avg * u_avg * inv_gamma_minus_one + f1 * kin_avg + pv_avg
    else  # orientation == 2, y-direction
        pv_avg = 0.5 * (p_ll * v_rr_vel + p_rr * v_ll_vel)
        
        f1 = ρ_avg * v_avg
        f2 = f1 * u_avg
        f3 = f1 * v_avg + p_avg
        f4 = p_avg * v_avg * inv_gamma_minus_one + f1 * kin_avg + pv_avg
    end
    
    return f1, f2, f3, f4
end


@inline function cons2prim_2d(u, γ)
    ρ = u[1]
    ρu = u[2]
    ρv = u[3]
    ρE = u[4]
    
    inv_ρ = 1.0 / ρ
    u_vel = ρu * inv_ρ
    v_vel = ρv * inv_ρ
    
    # Pressure from total energy
    p = (γ - 1.0) * (ρE - 0.5 * ρ * (u_vel^2 + v_vel^2))
    
    return ρ, u_vel, v_vel, p
end

function _expansion_inviscid_KEP_twopoint!(u, uprimitive, 
                                           neqs, ngl, dψ, ω,
                                           F, G, S,
                                           Je,
                                           dξdx, dξdy,
                                           dηdx, dηdy,
                                           rhs_el, iel,                                          
                                           ::CL, QT::Inexact, SD::NSD_2D, AD::ContGal)

       
    
    PhysConst = PhysicalConst{Float32}()
    # Loop over quadrature points in the element
    for j = 1:ngl
        for i = 1:ngl
            ωJac = ω[i] * ω[j] * Je[iel, i, j]
            
            # Get conservative and primitive variables at point (i,j)
            u_ij = @SVector [u[i, j, 1], u[i, j, 2], 
                             u[i, j, 3], u[i, j, 4]]
            uprim_ij = (uprimitive[i, j, 1], uprimitive[i, j, 2], 
                        uprimitive[i, j, 3], uprimitive[i, j, 4])
            
            # --- Compute two-point flux divergence in ξ-direction ---
            dFdξ_1, dFdξ_2, dFdξ_3, dFdξ_4 = 0.0, 0.0, 0.0, 0.0
            
            for k = 1:ngl
                # Get state at point (k,j)
                u_kj = @SVector [u[k, j, 1], u[k, j, 2], 
                                 u[k, j, 3], u[k, j, 4]]
                uprim_kj = (uprimitive[k, j, 1], uprimitive[k, j, 2], 
                            uprimitive[k, j, 3], uprimitive[k, j, 4])
                
                # Compute two-point flux in ξ-direction (orientation=1 for x-component)
                f1, f2, f3, f4 = flux_shima_etal_2d(u_ij, u_kj, uprim_ij, uprim_kj, 1, PhysConst.γ)
                
                # Apply derivative matrix
                dFdξ_1 += dψ[k, i] * f1
                dFdξ_2 += dψ[k, i] * f2
                dFdξ_3 += dψ[k, i] * f3
                dFdξ_4 += dψ[k, i] * f4
            end
            
            # --- Compute two-point flux divergence in η-direction ---
            dGdη_1, dGdη_2, dGdη_3, dGdη_4 = 0.0, 0.0, 0.0, 0.0
            
            for k = 1:ngl
                # Get state at point (i,k)
                u_ik = @SVector [u[i, k, 1], u[i, k, 2], 
                                 u[i, k, 3], u[i, k, 4]]
                uprim_ik = (uprimitive[i, k, 1], uprimitive[i, k, 2], 
                            uprimitive[i, k, 3], uprimitive[i, k, 4])
                
                # Compute two-point flux in η-direction (orientation=2 for y-component)
                g1, g2, g3, g4 = flux_shima_etal_2d(u_ij, u_ik, uprim_ij, uprim_ik, 2, PhysConst.γ)
                
                # Apply derivative matrix
                dGdη_1 += dψ[k, j] * g1
                dGdη_2 += dψ[k, j] * g2
                dGdη_3 += dψ[k, j] * g3
                dGdη_4 += dψ[k, j] * g4
            end
            
            # --- Transform to physical space ---
            dξdx_ij = dξdx[iel, i, j]
            dξdy_ij = dξdy[iel, i, j]
            dηdx_ij = dηdx[iel, i, j]
            dηdy_ij = dηdy[iel, i, j]
            
            # For each equation, compute divergence in physical space
            div_flux = MVector{4, Float64}(0.0, 0.0, 0.0, 0.0)
            
            div_flux[1] = dFdξ_1 * dξdx_ij + dGdη_1 * dηdy_ij
            div_flux[2] = dFdξ_2 * dξdx_ij + dGdη_2 * dηdy_ij
            div_flux[3] = dFdξ_3 * dξdx_ij + dGdη_3 * dηdy_ij
            div_flux[4] = dFdξ_4 * dξdx_ij + dGdη_4 * dηdy_ij
            
            # --- Update RHS ---
            rhs_el[iel, i, j, 1] -= ωJac * (div_flux[1] - S[i, j, 1])
            rhs_el[iel, i, j, 2] -= ωJac * (div_flux[2] - S[i, j, 2])
            rhs_el[iel, i, j, 3] -= ωJac * (div_flux[3] - S[i, j, 3])
            rhs_el[iel, i, j, 4] -= ωJac * (div_flux[4] - S[i, j, 4])
        end
    end
end
