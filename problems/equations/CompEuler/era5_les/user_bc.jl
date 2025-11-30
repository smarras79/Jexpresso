"""
    ERA5 LES - Boundary Conditions

    Author: Jexpresso Development Team
    Date: 2025-11-30
"""

function user_bc_dirichlet!(q,
                            x::AbstractFloat, y::AbstractFloat, z::AbstractFloat,
                            t::AbstractFloat, tag,
                            qbdy::AbstractArray,
                            nx, ny, nz,
                            xmin, xmax,
                            ymin, ymax,
                            zmin, zmax,
                            qe, ::TOTAL)

    # Free-slip lateral boundaries (periodic handled by mesh)
    # No-slip bottom boundary
    qnl = nx*q[2] + ny*q[3] + nz*q[4]
    qbdy[2] = (q[2] - qnl*nx)
    qbdy[3] = (q[3] - qnl*ny)
    qbdy[4] = (q[4] - qnl*nz)

    # Bottom boundary (no-slip)
    if (z < 0.01)
        qbdy[2] = 0.0  # u = 0
        qbdy[3] = 0.0  # v = 0
        qbdy[4] = 0.0  # w = 0
    end

end

function user_bc_dirichlet!(q,
                            x::AbstractFloat, y::AbstractFloat, z::AbstractFloat,
                            t::AbstractFloat, tag,
                            qbdy::AbstractArray,
                            nx, ny, nz,
                            xmin, xmax,
                            ymin, ymax,
                            zmin, zmax,
                            qe, ::PERT)

    # Perturbation formulation
    qnl = nx*(q[2]+qe[2]) + ny*(q[3]+qe[3]) + nz*(q[4]+qe[4])
    qbdy[2] = (q[2]+qe[2] - qnl*nx) - qe[2]
    qbdy[3] = (q[3]+qe[3] - qnl*ny) - qe[3]
    qbdy[4] = (q[4]+qe[4] - qnl*nz) - qe[4]

    # Bottom boundary
    if (z < 0.01)
        qbdy[2] = -qe[2]  # u' = -u_ref
        qbdy[3] = -qe[3]  # v' = -v_ref
        qbdy[4] = -qe[4]  # w' = -w_ref
    end
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray,
                        x::AbstractFloat, y::AbstractFloat, z::AbstractFloat,
                        t::AbstractFloat, tag::String, inputs::Dict)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_dirichlet_gpu(q,qe,x,y,z,t,nx,ny,nz,qbdy,lpert)
    T = eltype(q)

    if (lpert)
        qnl = nx*(q[2]+qe[2]) + ny*(q[3]+qe[3]) + nz*(q[4]+qe[4])
        u = (q[2]+qe[2] - qnl*nx) - qe[2]
        v = (q[3]+qe[3] - qnl*ny) - qe[3]
        w = (q[4]+qe[4] - qnl*nz) - qe[4]

        if (z < 0.01)
            u = -qe[2]
            v = -qe[3]
            w = -qe[4]
        end
    else
        qnl = nx*q[2] + ny*q[3] + nz*q[4]
        u = (q[2] - qnl*nx)
        v = (q[3] - qnl*ny)
        w = (q[4] - qnl*nz)

        if (z < 0.01)
            u = 0.0
            v = 0.0
            w = 0.0
        end
    end

    return T(qbdy[1]), T(u), T(v), T(w), T(qbdy[5]), T(qbdy[6]), T(qbdy[7])
end
