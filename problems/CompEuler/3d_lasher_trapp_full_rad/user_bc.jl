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

    qnl = nx*(q[2]+qe[2]) + ny*(q[3]+qe[3]) + nz*(q[4]+qe[4])
    qbdy[2] = (q[2]+qe[2] - qnl*nx) - qe[2]
    qbdy[3] = (q[3]+qe[3] - qnl*ny) - qe[3]
    qbdy[4] = (q[4]+qe[4] - qnl*nz) - qe[4]
end

#
# Lasher-Trapp et al. (2001): "The cloud is initiated with a Gaussian heat flux (300 W m-2 maximum) at the
# lower boundary." The flux of the energy variable ρhl [J m-3] is a heat flux [W m-2], so no conversion.
# Signature as called by the kernel (BCs.jl). Requires :bdy_fluxes => true and :lwall_model => true
# (the latter only allocates the wall-flux arrays τ_f, wθ that the kernel passes in).
#
const LT_QMAX  = 300.0     # W m-2
const LT_SIGMA = 500.0     # m; keeps the flux < 1% of LT_QMAX at the lateral boundaries (3 km away)

function user_bc_neumann!(F, q, q1, qe, qe1, tag, coords, τ_f, wθ, SOL_VARS_TYPE, PhysConst; kwargs...)
    tag == "bottom" || return nothing
    x, y = coords[1], coords[2]
    xc, yc = 4050.0, 4050.0          # domain centre (8.1 km x 8.1 km)
    F[5] = LT_QMAX * exp(-((x - xc)^2 + (y - yc)^2) / (2 * LT_SIGMA^2))
    return nothing
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
        if (coords[3] < 0.01)
            u = 0.0
            v = 0.0
            w = 0.0
        end

    end
    return T(qbdy[1]), T(u), T(v), T(w), T(qbdy[5]), T(qbdy[6]), T(qbdy[7])
end
