function user_primitives!(u, qe, uprimitive, ::TOTAL)
    PhysConst = PhysicalConst{Float64}()
    uprimitive[1] = u[1]
    uprimitive[2] = u[2]/u[1]
    uprimitive[3] = u[3]/u[1]
    uprimitive[4] = u[4]/u[1]
    uprimitive[5] = u[5]/u[1]
    uprimitive[6] = u[6]/u[1]
    uprimitive[7] = u[7]/u[1]
end

function user_primitives!(u,qe,uprimitive,::PERT)
    uprimitive[1] = u[1]+qe[1]
    uprimitive[2] = (u[2]+qe[2])/(u[1]+qe[1])-qe[2]/qe[1]
    uprimitive[3] = u[3]/(u[1]+qe[1])
    uprimitive[4] = u[4]/(u[1]+qe[1])
    uprimitive[5] = (u[5]+qe[5])/(u[1]+qe[1])-qe[5]/qe[1]
    uprimitive[6] = (u[6]+qe[6])/(u[1]+qe[1])-qe[6]/qe[1]
    uprimitive[7] = (u[7]+qe[7])/(u[1]+qe[1])-qe[7]/qe[1]
end

function user_primitives_gpu(u,qe,lpert)
    T = eltype(u)
    if (lpert)
        return T(u[1]+qe[1]), T(u[2]/(u[1]+qe[1])), T(u[3]/(u[1]+qe[1])), T(u[4]/(u[1]+qe[1])), T((u[5]+qe[5])/(u[1]+qe[1]) - qe[5]/qe[1])
    else
        return T(u[1]), T(u[2]/u[1]), T(u[3]/u[1]), T(u[4]/u[1]), T(u[5]/u[1])
    end
end

function user_read_vtu_point_data!(q, vars, ip_map, mesh)
    ρ_arr  = vars["ρ"]
    u_arr  = vars["ρu"]    # velocity u  (despite the name)
    v_arr  = vars["ρv"]
    w_arr  = vars["ρw"]
    hl_arr = vars["hl"]    # specific moist enthalpy hl/ρ
    qt_arr = vars["ρqt"]   # specific total water qt
    qp_arr = vars["ρqp"]   # specific precipitation qp
    P_arr  = vars["pressure"]

    for ip = 1:mesh.npoin
        j  = ip_map[ip]
        ρ  = ρ_arr[j]
        u  = u_arr[j]
        v  = v_arr[j]
        w  = w_arr[j]
        hl = hl_arr[j]
        qt = qt_arr[j]
        qp = qp_arr[j]

        q.qn[ip, 1]   = ρ
        q.qn[ip, 2]   = ρ * u
        q.qn[ip, 3]   = ρ * v
        q.qn[ip, 4]   = ρ * w
        q.qn[ip, 5]   = ρ * hl
        q.qn[ip, 6]   = ρ * qt
        q.qn[ip, 7]   = ρ * qp
        q.qn[ip, end] = P_arr[j]
    end
end

function user_uout!(ip, ET, uout, u, qe; mp = mp)

    if ET == TOTAL()
        uout[1] = u[1]
        uout[2] = u[2]/u[1]
        uout[3] = u[3]/u[1]
        uout[4] = u[4]/u[1]
        uout[5] = u[5]/u[1]
        uout[6] = u[6]/u[1]
        uout[7] = u[7]/u[1]

        uout[16] = u[2]/u[1]-qe[2]/qe[1]
        uout[17] = u[3]/u[1]-qe[3]/qe[1]
        uout[18] = u[4]/u[1]-qe[4]/qe[1]
        uout[19] = u[5]/u[1]-qe[5]/qe[1]
        uout[20] = u[end]
    elseif ET == PERT()
        ρ_tot = u[1] + qe[1]
        uout[1] = ρ_tot
        uout[2] = (u[2]+qe[2])/ρ_tot
        uout[3] = (u[3]+qe[3])/ρ_tot
        uout[4] = (u[4]+qe[4])/ρ_tot
        uout[5] = (u[5]+qe[5])/ρ_tot - qe[5]/qe[1]
        uout[6] = (u[6]+qe[6])/ρ_tot - qe[6]/qe[1]
        uout[7] = (u[7]+qe[7])/ρ_tot - qe[7]/qe[1]

        uout[16] = (u[2]+qe[2])/ρ_tot - qe[2]/qe[1]
        uout[17] = (u[3]+qe[3])/ρ_tot - qe[3]/qe[1]
        uout[18] = (u[4]+qe[4])/ρ_tot - qe[4]/qe[1]
        uout[19] = (u[5]+qe[5])/ρ_tot - qe[5]/qe[1]
        uout[20] = u[end]
    end
    uout[8] = mp.Tabs[ip]
    uout[9] = mp.qn[ip]
    uout[10] = mp.qc[ip]
    uout[11] = mp.qi[ip]
    uout[12] = mp.qr[ip]
    uout[13] = mp.qs[ip]
    uout[14] = mp.qg[ip]
    uout[15] = mp.qsatt[ip]
        
end

