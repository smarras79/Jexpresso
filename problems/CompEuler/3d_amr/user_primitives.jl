function user_primitives!(u,qe,uprimitive,::TOTAL)
    uprimitive[1] = u[1]
    uprimitive[2] = u[2]/u[1]
    uprimitive[3] = u[3]/u[1]
    uprimitive[4] = u[4]/u[1]
    uprimitive[5] = u[5]/u[1]
end

function user_primitives!(u,qe,uprimitive,::PERT)
    uprimitive[1] = u[1]+qe[1]
    uprimitive[2] = u[2]/(u[1]+qe[1])
    uprimitive[3] = u[3]/(u[1]+qe[1])
    uprimitive[4] = u[4]/(u[1]+qe[1])
    uprimitive[5] = (u[5]+qe[5])/(u[1]+qe[1])-qe[5]/qe[1]
end

function user_primitives_gpu(u,qe,lpert)
    T = eltype(u)
    if (lpert)
        return T(u[1]+qe[1]), T(u[2]/(u[1]+qe[1])), T(u[3]/(u[1]+qe[1])), T(u[4]/(u[1]+qe[1])), T((u[5]+qe[5])/(u[1]+qe[1]) - qe[5]/qe[1])
    else
        return T(u[1]), T(u[2]/u[1]), T(u[3]/u[1]), T(u[4]/u[1]), T(u[5]/u[1])
    end
end


function user_uout!(ip, ET, uout, u, qe; kwargs...)

    #
    # IMPORTANT NOTICE:
    # uout has length equivalent to 'qoutvars' defined in your initialize.jl
    # If you have not defined qoutvars in initialize.jl, then uout will only
    # be allocated to the length of u by default.
    #
    # NOTICE: qaux (this function's `u`) only carries `neqs` columns — there
    # is no extra "pressure" state column to read here. uout[6] ("p") must be
    # recomputed from the equation of state (rather than qe[end], which is
    # just the constant background reference pressure) so it's the actual
    # current total pressure — needed so read_vtk_amr_restart!/
    # user_read_vtu_point_data! can reconstruct a correct restart state.
    PhysConst = PhysicalConst{Float64}()
    if ET == TOTAL()
        ρ = u[1]
        θ = u[5]/ρ
        uout[1] = ρ
        uout[2] = u[2]/ρ
        uout[3] = u[3]/ρ
        uout[4] = u[4]/ρ
        uout[5] = θ
        uout[6] = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
    elseif ET == PERT()
        ρ = u[1]+qe[1]
        θ = (u[5]+qe[5])/ρ
        uout[1] = ρ
        uout[2] = u[2]/ρ
        uout[3] = u[3]/ρ
        uout[4] = u[4]/ρ
        uout[5] = θ - qe[5]/qe[1]
        uout[6] = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
    end
end

function user_read_vtu_point_data!(q, vars, ip_map, mesh)
    ρ_arr = vars["ρ"]
    u_arr = vars["u"]
    v_arr = vars["v"]
    w_arr = vars["w"]
    θ_arr = vars["θ"]
    p_arr = vars["p"]

    for ip = 1:mesh.npoin
        j = ip_map[ip]
        ρ = ρ_arr[j]

        q.qn[ip, 1]   = ρ
        q.qn[ip, 2]   = ρ * u_arr[j]
        q.qn[ip, 3]   = ρ * v_arr[j]
        q.qn[ip, 4]   = ρ * w_arr[j]
        q.qn[ip, 5]   = ρ * θ_arr[j]
        q.qn[ip, end] = p_arr[j]
    end
end
