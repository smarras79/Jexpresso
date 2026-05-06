function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)

    @info " Initialize fields for 2D moist CompEuler (hl formulation) ................. "

    qvars    = ("ρ", "ρu", "ρv", "hl", "ρqt", "ρqp")
    qoutvars = ["ρ", "ρu", "ρv", "hl", "ρqt", "ρqp", "T", "qn", "qc", "qi", "qr", "qs", "qg",
                "u_prime", "v_prime", "hl_prime", "qt_prime"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend];
                 neqs=length(qvars), qoutvars=qoutvars)

    data       = read_sounding(inputs[:sounding_file])
    background = interpolate_sounding(inputs[:backend], mesh.npoin, mesh.y, data)

    if (inputs[:backend] == CPU())

        PhysConst = PhysicalConst{Float64}()

        # ── Step 1: hydrostatic column integration (moist) ───────────────────
        # Uses virtual potential temperature θv = θ*(1 + 0.61*qv)
        x_unique = sort(unique(round.(mesh.x, digits=2)))
        p_hydro  = zeros(mesh.npoin)

        for xi in x_unique
            col_mask = abs.(mesh.x .- xi) .< 0.01
            col_pts  = findall(col_mask)
            col_z    = mesh.y[col_pts]
            sort_idx = sortperm(col_z)
            col_pts  = col_pts[sort_idx]
            col_z    = col_z[sort_idx]

            ip_sfc = col_pts[1]
            p_hydro[ip_sfc] = background[ip_sfc, 5]   # col 6 = pressure

            for k = 1:length(col_pts)-1
                ip_lo = col_pts[k]
                ip_hi = col_pts[k+1]
                dz    = col_z[k+1] - col_z[k]

                θ_lo  = background[ip_lo, 1]
                qv_lo = background[ip_lo, 2] / 1000.0   # g/kg → kg/kg
                θv_lo = θ_lo * (1.0 + 0.61*qv_lo)

                θ_hi  = background[ip_hi, 1]
                qv_hi = background[ip_hi, 2] / 1000.0
                θv_hi = θ_hi * (1.0 + 0.61*qv_hi)

                ρ_lo  = perfectGasLaw_θPtoρ(PhysConst; θ=θv_lo, Press=p_hydro[ip_lo])

                p_hi = p_hydro[ip_lo] - ρ_lo * PhysConst.g * dz
                for _ = 1:20
                    ρ_hi     = perfectGasLaw_θPtoρ(PhysConst; θ=θv_hi, Press=p_hi)
                    p_hi_new = p_hydro[ip_lo] - 0.5*(ρ_lo + ρ_hi) * PhysConst.g * dz
                    if abs(p_hi_new - p_hi) < 1e-10; break; end
                    p_hi = p_hi_new
                end
                p_hydro[ip_hi] = p_hi
            end
        end

        # ── Step 2: initialize q ──────────────────────────────────────────────
        for iel_g = 1:mesh.nelem
            for j = 1:mesh.ngl, i = 1:mesh.ngl

                ip  = mesh.connijk[iel_g, i, j]
                y   = mesh.y[ip]

                θ    = background[ip, 1]
                qv   = background[ip, 2] / 1000.0   # kg/kg
                u    = background[ip, 3]
                p    = p_hydro[ip]

                θv   = θ * (1.0 + 0.61*qv)
                Tv   = θv * (p / PhysConst.pref)^(PhysConst.Rair / PhysConst.cp)
                ρ    = perfectGasLaw_TPtoρ(PhysConst; Temp=Tv, Press=p)
                hl   = PhysConst.cp * Tv + PhysConst.g * y

                ρref  = ρ
                hlref = hl
                uref  = u

                if inputs[:SOL_VARS_TYPE] == PERT()
                    q.qn[ip, 1] = ρ    - ρref
                    q.qn[ip, 2] = ρ*u  - ρref*uref
                    q.qn[ip, 3] = 0.0
                    q.qn[ip, 4] = ρ*hl - ρref*hlref
                    q.qn[ip, 5] = ρ*qv - ρref*qv
                    q.qn[ip, 6] = 0.0
                    q.qn[ip, end] = p
                else
                    q.qn[ip, 1] = ρ
                    q.qn[ip, 2] = ρ*u
                    q.qn[ip, 3] = 0.0
                    q.qn[ip, 4] = ρ*hl
                    q.qn[ip, 5] = ρ*qv
                    q.qn[ip, 6] = 0.0
                    q.qn[ip, end] = p
                end

                q.qe[ip, 1] = ρref
                q.qe[ip, 2] = ρref*uref
                q.qe[ip, 3] = 0.0
                q.qe[ip, 4] = ρref*hlref
                q.qe[ip, 5] = ρref*qv
                q.qe[ip, 6] = 0.0
                q.qe[ip, end] = p
            end
        end

        @info "max hl_prime at t=0: $(maximum(q.qn[:,4]))"
        @info "min hl_prime at t=0: $(minimum(q.qn[:,4]))"
        @info "max qt at t=0:       $(maximum(q.qn[:,5]))"

        # ── Step 3: semi-infinite elements ────────────────────────────────────
        for iel_g = 1:mesh.nelem_semi_inf
            for j = 1:mesh.ngr, i = 1:mesh.ngl

                ip  = mesh.connijk_lag[iel_g, i, j]
                y   = mesh.y[ip]

                θ    = background[ip, 1]
                qv   = background[ip, 2] / 1000.0
                u    = background[ip, 3]
                p    = p_hydro[ip]

                θv   = θ * (1.0 + 0.61*qv)
                Tv   = θv * (p / PhysConst.pref)^(PhysConst.Rair / PhysConst.cp)
                ρ    = perfectGasLaw_TPtoρ(PhysConst; Temp=Tv, Press=p)
                hl   = PhysConst.cp * Tv + PhysConst.g * y

                ρref  = ρ
                hlref = hl
                uref  = u

                if inputs[:SOL_VARS_TYPE] == PERT()
                    q.qn[ip, 1] = ρ    - ρref
                    q.qn[ip, 2] = ρ*u  - ρref*uref
                    q.qn[ip, 3] = 0.0
                    q.qn[ip, 4] = ρ*hl - ρref*hlref
                    q.qn[ip, 5] = ρ*qv - ρref*qv
                    q.qn[ip, 6] = 0.0
                    q.qn[ip, end] = p
                else
                    q.qn[ip, 1] = ρ
                    q.qn[ip, 2] = ρ*u
                    q.qn[ip, 3] = 0.0
                    q.qn[ip, 4] = ρ*hl
                    q.qn[ip, 5] = ρ*qv
                    q.qn[ip, 6] = 0.0
                    q.qn[ip, end] = p
                end

                q.qe[ip, 1] = ρref
                q.qe[ip, 2] = ρref*uref
                q.qe[ip, 3] = 0.0
                q.qe[ip, 4] = ρref*hlref
                q.qe[ip, 5] = ρref*qv
                q.qe[ip, 6] = 0.0
                q.qe[ip, end] = p
            end
        end

    end  # CPU

    @info "Initialize moist fields ........................ DONE"
    @info maximum(q.qn[:,1:6]), minimum(q.qn[:,1:6])
    return q
end