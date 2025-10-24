function user_flux!(F, G, SD::NSD_2D, q, mesh::St_mesh; neqs=4, ip=1)

    PhysConst = PhysicalConst{Float64}()
                
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρθ = q[4]
    θ  = ρθ/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    
    Press = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
    F[1] = ρu
    F[2] = ρu*u + Press
    F[3] = ρv*u
    F[4] = ρθ*u
    
    G[1] = ρv
    G[2] = ρu*v
    G[3] = ρv*v + Press
    G[4] = ρθ*v
end

## flux_kennedy_gruber
@inline function user_volume_flux(u_ll, u_rr)
    PhysConst = PhysicalConst{Float64}()
    Temp_ll = (rho_e_ll/rho_ll - 0.5*v1_ll*v1_ll)/PhysConst.cv
    Temp_rr = (rho_e_rr/rho_rr - 0.5*v1_rr*v1_rr)/PhysConst.cv

    p_ll = perfectGasLaw_ρTtoP(PhysConst; ρ=rho_ll, Temp=Temp_ll)
    p_rr = perfectGasLaw_ρTtoP(PhysConst; ρ=rho_rr, Temp=Temp_rr)

    # Unpack left and right state
    rho_ll, rho_v1_ll, rho_v2_ll, rho_e_ll = u_ll
    rho_rr, rho_v1_rr, rho_v2_rr, rho_e_rr = u_rr

    # Average each factor of products in flux
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)
    e_avg = 0.5f0 * (rho_e_ll / rho_ll + rho_e_rr / rho_rr)

    # Calculate fluxes depending on orientation
    
        f1 = rho_avg * v1_avg
        f2 = rho_avg * v1_avg * v1_avg + p_avg
        f3 = rho_avg * v1_avg * v2_avg
        f4 = (rho_avg * e_avg + p_avg) * v1_avg
    
        g1 = rho_avg * v2_avg
        g2 = rho_avg * v2_avg * v1_avg
        g3 = rho_avg * v2_avg * v2_avg + p_avg
        g4 = (rho_avg * e_avg + p_avg) * v2_avg

    return SVector(f1, f2, f3, f4), SVector(g1, g2, g3, g4)
end