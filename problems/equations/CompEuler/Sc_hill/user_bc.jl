function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny, qe, ::PERT)
    x = coords[1]
    y = coords[2]
    
    # For inflow boundaries (left side typically), maintain reference state + background flow
    if tag == "left" || abs(nx - (-1.0)) < 1e-10
        # Inflow: prescribe reference state with horizontal wind
        qbdy[1] = 0.0  # ρ perturbation = 0
        qbdy[2] = 0.0  # ρu perturbation = 0 (background u from qe)
        qbdy[3] = 0.0  # ρv perturbation = 0
        qbdy[4] = 0.0  # ρhl perturbation = 0
        qbdy[5] = 0.0  # ρqt perturbation = 0
        qbdy[6] = 0.0  # ρqp perturbation = 0
    else
        # For other boundaries (top, right), use slip condition
        qnl = nx*(q[2]+qe[2]) + ny*(q[3]+qe[3])
        qbdy[2] = (q[2]+qe[2] - qnl*nx) - qe[2]
        qbdy[3] = (q[3]+qe[3] - qnl*ny) - qe[3]
    end
end

function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String, qbdy::AbstractArray, nx, ny, qe, ::TOTAL)
    x = coords[1]
    y = coords[2]
    
    # For TOTAL formulation
    if tag == "left" || abs(nx - (-1.0)) < 1e-10
        # Inflow: prescribe reference state
        qbdy[1] = qe[1]  # reference density
        qbdy[2] = qe[2]  # reference ρu
        qbdy[3] = qe[3]  # reference ρv
        qbdy[4] = qe[4]  # reference ρhl
        qbdy[5] = qe[5]  # reference ρqt
        qbdy[6] = qe[6]  # reference ρqp
    else
        # Slip condition
        qnl = nx*q[2] + ny*q[3]
        qbdy[2] = q[2] - qnl*nx
        qbdy[3] = q[3] - qnl*ny
    end
end

function user_bc_dirichlet_gpu(q, qe, x, y, t, nx, ny, qbdy, lpert)
    T = eltype(q)
    u = qbdy[2]
    v = qbdy[3]
    
    # Detect boundary type based on normal vector
    # Left boundary: nx ≈ -1
    # Top boundary: ny ≈ 1
    # Right boundary: nx ≈ 1
    
    if abs(nx + T(1.0)) < T(0.1)  # Left boundary (inflow)
        if lpert
            # Inflow with zero perturbations
            return T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0)
        else
            # Inflow with reference state
            return T(qe[1]), T(qe[2]), T(qe[3]), T(qe[4]), T(qe[5]), T(qe[6])
        end
    else
        # Other boundaries: slip condition
        if lpert
            qnl = nx*(q[2]+qe[2]) + ny*(q[3]+qe[3])
            u = (q[2]+qe[2] - qnl*nx) - qe[2]
            v = (q[3]+qe[3] - qnl*ny) - qe[3]
        else
            qnl = nx*q[2] + ny*q[3]
            u = q[2] - qnl*nx
            v = q[3] - qnl*ny
        end
        return T(qbdy[1]), T(u), T(v), T(qbdy[4]), T(qbdy[5]), T(qbdy[6])
    end
end