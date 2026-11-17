function user_bc_dirichlet!(q,
                            coords,
                            t::AbstractFloat, tag,
                            qbdy::AbstractArray,
                            nx, ny, nz,
                            xmin, xmax,
                            ymin, ymax,
                            zmin, zmax,
                            qe, ::TOTAL)

    qnl = nx*q[2] + ny*q[3] + nz*q[4]
    qbdy[2] = (q[2] - qnl*nx) 
    qbdy[3] = (q[3] - qnl*ny) 
    qbdy[4] = (q[4] - qnl*nz) 
    if (coords[3] < 0.01) 
        qbdy[2] = 0.0
        qbdy[3] = 0.0
        qbdy[4] = 0.0
        # e_tot
        qbdy[5] = 59169*q[1]
        qbdy[6] = 0.0175*q[1]
    end
    
end

function user_bc_dirichlet!(q,
                            coords,
                            t::AbstractFloat, tag,
                            qbdy::AbstractArray,
                            nx, ny, nz,
                            xmin, xmax,
                            ymin, ymax,
                            zmin, zmax,
                            qe, ::PERT)

    # Free-slip wall: remove the normal component of the total momentum. The kernel recognises
    # this projection and constrains only the normal momentum (tangential momentum stays free).
    qnl = nx*(q[2]+qe[2]) + ny*(q[3]+qe[3]) + nz*(q[4]+qe[4])
    qbdy[2] = (q[2]+qe[2] - qnl*nx) - qe[2]
    qbdy[3] = (q[3]+qe[3] - qnl*ny) - qe[3]
    qbdy[4] = (q[4]+qe[4] - qnl*nz) - qe[4]
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords, t::AbstractFloat, tag::String, inputs)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_dirichlet_gpu(q, qe, coords, t, nx, ny, nz, qbdy, lpert)
    T = eltype(q)
    if (lpert)
        qnl = nx*(q[2]+qe[2]) + ny*(q[3]+qe[3]) + nz*(q[4]+qe[4])
        v = (q[3]+qe[3] - qnl*ny) - qe[3]
        w = (q[4]+qe[4] - qnl*nz) - qe[4]
    else
        qnl = nx*(q[2]) + ny*(q[3]) + nz*(q[4])
        
        u = (q[2] - qnl*nx)
        v = (q[3] - qnl*ny)
        w = (q[4] - qnl*nz)
        ρe = qbdy[5]
        ρqt = qbdy[6]
        if (coords[3] < 0.01)
            u = 0.0
            v = 0.0
            w = 0.0
            ρe = 59169*q[1]
            ρqt = 0.0175*q[1]
        end

    end
    return T(qbdy[1]), T(u), T(v), T(w), T(ρe), T(ρqt), T(qbdy[7])
end

#----------------------------------------------------------------------------------------------
# BOMEX prescribed surface fluxes (Siebesma et al. 2003), applied on the "bottom" tag through the
# kernel Neumann path. Active when user_inputs sets :bdy_fluxes => true (and :lwall_model => true,
# which only allocates the per-face τ_f / wθ arrays this call is handed).
# F is the flux INTO the domain per unit area; the kernel adds ∫ F ψ dS to the RHS.
#   hl  (ρ hl, hl = cp T + g z - Lc q_l): sensible heat only, ρ cp w'θl'  (hl carries no vapour latent heat)
#   ρqt: ρ w'qt'
#   momentum: τ = -ρ u*² U_h/|U_h|  (wind at the wall node)
#----------------------------------------------------------------------------------------------
const BOMEX_WTHL  = 8.0e-3     # w'θl' [K m/s]
const BOMEX_WQT   = 5.2e-5     # w'qt' [kg/kg m/s]
const BOMEX_USTAR = 0.28       # u* [m/s]

function user_bc_neumann!(F, q, q1, qe, qe1, tag, coords, τ_f, wθ, SOL_VARS_TYPE, PhysConst; kwargs...)
    tag == "bottom" || return nothing
    lpert = (SOL_VARS_TYPE == PERT())
    ρ  = lpert ? q[1] + qe[1] : q[1]
    ρu = lpert ? q[2] + qe[2] : q[2]
    ρv = lpert ? q[3] + qe[3] : q[3]
    U  = max(sqrt(ρu^2 + ρv^2) / ρ, 0.1)          # floor avoids division by ~0
    τ  = ρ * BOMEX_USTAR^2
    F[2] = -τ * (ρu / ρ) / U
    F[3] = -τ * (ρv / ρ) / U
    F[5] = ρ * PhysConst.cp * BOMEX_WTHL
    F[6] = ρ * BOMEX_WQT
    if length(τ_f) >= 2                            # keep for output/diagnostics
        τ_f[1] = F[2]; τ_f[2] = F[3]
    end
    length(wθ) >= 1 && (wθ[1] = BOMEX_WTHL)
    return nothing
end
