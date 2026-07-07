function user_primitives!(u,qe,uprimitive,::TOTAL)
    uprimitive[1] = u[1]
    uprimitive[2] = u[2]/u[1]
    uprimitive[3] = u[3]/u[1]
    uprimitive[4] = u[4]/u[1]
end

function user_primitives!(u,qe,uprimitive,::PERT)
    uprimitive[1] = u[1]+qe[1]
    uprimitive[2] = u[2]/(u[1]+qe[1])
    uprimitive[3] = u[3]/(u[1]+qe[1])
    uprimitive[4] = (u[4]+qe[4])/(u[1]+qe[1])-qe[4]/qe[1]
end

function user_primitives_gpu(u, qe, lpert)
    T = eltype(u)
    if (lpert)
        return T(u[1]+qe[1]), T(u[2]/(u[1]+qe[1])), T(u[3]/(u[1]+qe[1])), T((u[4]+qe[4])/(u[1]+qe[1]) - qe[4]/qe[1])
    else
        return T(u[1]), T(u[2]/u[1]), T(u[3]/u[1]), T(u[4]/u[1])
    end
end

function user_uout!(ip, ET, uout, u, qe; kwargs...)
    # NOTICE: qaux (this function's `u`) only carries `neqs` columns — there
    # is no extra "pressure" state column to read here (unlike q.qn/q.qe,
    # which do get one at initialize() time). Pressure must be recomputed
    # from the equation of state to be written correctly to VTK — needed so
    # read_vtk_restart!/user_read_vtu_point_data! can reconstruct a correct
    # restart pressure, since theta/theta_amr's own user_uout! leaves this
    # "p" slot at 0 (harmless there — those cases never restart from VTK).
    PhysConst = PhysicalConst{Float64}()
    ρ = u[1]
    θ = u[4]/u[1]

    uout[1] = ρ
    uout[2] = u[2]/ρ
    uout[3] = u[3]/ρ
    uout[4] = θ
    uout[5] = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
end

function user_read_vtu_point_data!(q, vars, ip_map, mesh)
    ρ_arr = vars["ρ"]
    u_arr = vars["u"]
    w_arr = vars["w"]   # second velocity component (named "w", not vertical z-velocity)
    θ_arr = vars["θ"]
    p_arr = vars["p"]

    for ip = 1:mesh.npoin
        j = ip_map[ip]
        ρ = ρ_arr[j]

        q.qn[ip, 1]   = ρ
        q.qn[ip, 2]   = ρ * u_arr[j]
        q.qn[ip, 3]   = ρ * w_arr[j]
        q.qn[ip, 4]   = ρ * θ_arr[j]
        q.qn[ip, end] = p_arr[j]
    end
end
