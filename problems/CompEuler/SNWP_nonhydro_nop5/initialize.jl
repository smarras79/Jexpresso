using Base

function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """

                """
    @info " Initialize fields for 2D CompEuler with θ equation ........................ "

    qvars = ["dρ", "dρu", "dρv", "dρθ"]
    qoutvars = ["ρ", "ρu", "ρv", "θ", "θ_prime", "Press", "u_total", "v_total"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)

    data       = read_sounding(inputs[:sounding_file])
    background = interpolate_sounding(inputs[:backend], mesh.npoin, mesh.y, data)

    if (inputs[:backend] == CPU())

        PhysConst = PhysicalConst{Float64}()

        # ── Step 1: get all unique x-columns ──────────────────────────────────
        x_unique = sort(unique(round.(mesh.x, digits=2)))

        # ── Step 2: for each column, integrate hydrostatic balance upward ─────
        p_hydro = zeros(mesh.npoin)

        for xi in x_unique
            col_mask = abs.(mesh.x .- xi) .< 0.01
            col_pts  = findall(col_mask)
            col_z    = mesh.y[col_pts]
            sort_idx = sortperm(col_z)
            col_pts  = col_pts[sort_idx]
            col_z    = col_z[sort_idx]

            ip_sfc = col_pts[1]
            p_hydro[ip_sfc] = background[ip_sfc, 5]

            for k = 1:length(col_pts)-1
                ip_lo = col_pts[k]
                ip_hi = col_pts[k+1]
                dz    = col_z[k+1] - col_z[k]
                θ_lo  = background[ip_lo, 1]
                θ_hi  = background[ip_hi, 1]
                ρ_lo  = perfectGasLaw_θPtoρ(PhysConst; θ=θ_lo, Press=p_hydro[ip_lo])

                # iterative corrector
                p_hi = p_hydro[ip_lo] - ρ_lo * PhysConst.g * dz  # initial guess
                for _ = 1:20
                    ρ_hi     = perfectGasLaw_θPtoρ(PhysConst; θ=θ_hi, Press=p_hi)
                    p_hi_new = p_hydro[ip_lo] - 0.5*(ρ_lo + ρ_hi) * PhysConst.g * dz
                    if abs(p_hi_new - p_hi) < 1e-10
                        break
                    end
                    p_hi = p_hi_new
                end
                p_hydro[ip_hi] = p_hi
            end
        end

        # ── Step 3: initialize q using hydrostatically consistent pressure ────
        for iel_g = 1:mesh.nelem
            for j=1:mesh.ngl, i=1:mesh.ngl

                ip = mesh.connijk[iel_g,i,j]

                θ    = background[ip, 1]
                u    = background[ip, 3]
                v    = 0.0
                p    = p_hydro[ip]

                ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ, Press=p)
                ρref = ρ

                if inputs[:SOL_VARS_TYPE] == PERT()
                    q.qn[ip,1] = ρ - ρref
                    q.qn[ip,2] = ρ*u - ρref*u
                    q.qn[ip,3] = ρ*v - ρref*v
                    q.qn[ip,4] = ρ*θ - ρref*θ
                    q.qn[ip,end] = p
                else
                    q.qn[ip,1] = ρ
                    q.qn[ip,2] = ρ*u
                    q.qn[ip,3] = ρ*v
                    q.qn[ip,4] = ρ*θ
                    q.qn[ip,end] = p
                end
                q.qe[ip,1] = ρref
                q.qe[ip,2] = ρref*u
                q.qe[ip,3] = ρref*v
                q.qe[ip,4] = ρref*θ
                q.qe[ip,end] = p
            end
        end

        @info "max θ_prime at t=0: $(maximum(q.qn[:,4]))"
        @info "min θ_prime at t=0: $(minimum(q.qn[:,4]))"

        # ── Step 4: same for semi-infinite elements if present ────────────────
        for iel_g = 1:mesh.nelem_semi_inf
            for j=1:mesh.ngr, i=1:mesh.ngl

                ip = mesh.connijk_lag[iel_g,i,j]

                θ    = background[ip, 1]
                u    = background[ip, 3]
                v    = 0.0
                p    = p_hydro[ip]

                ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ, Press=p)
                ρref = ρ

                if inputs[:SOL_VARS_TYPE] == PERT()
                    q.qn[ip,1] = ρ - ρref
                    q.qn[ip,2] = ρ*u - ρref*u
                    q.qn[ip,3] = ρ*v - ρref*v
                    q.qn[ip,4] = ρ*θ - ρref*θ
                    q.qn[ip,end] = p
                else
                    q.qn[ip,1] = ρ
                    q.qn[ip,2] = ρ*u
                    q.qn[ip,3] = ρ*v
                    q.qn[ip,4] = ρ*θ
                    q.qn[ip,end] = p
                end
                q.qe[ip,1] = ρref
                q.qe[ip,2] = ρref*u
                q.qe[ip,3] = ρref*v
                q.qe[ip,4] = ρref*θ
                q.qe[ip,end] = p
            end
        end

    else
        if (inputs[:SOL_VARS_TYPE] == PERT())
            lpert = true
        else
            lpert = false
        end
        PhysConst = PhysicalConst{TFloat}()

        k = initialize_gpu!(inputs[:backend])
        k(q.qn, q.qe, mesh.x, mesh.y, background, PhysConst, lpert; ndrange = (mesh.npoin))
    end

    @info "Initialize fields for system of 2D CompEuler with θ equation ........................ DONE"
    @info maximum(q.qn[:,1:4]), minimum(q.qn[:,1:4])
    return q
end

@kernel function initialize_gpu!(qn, qe, x, y, background, PhysConst, lpert)
    ip = @index(Global, Linear)

    T = eltype(x)

    θ    = T(background[ip, 1])
    u    = T(background[ip, 3])
    v    = T(0.0)
    p    = T(background[ip, 5])

    ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ, Press=p)
    ρref = ρ

    if lpert
        qn[ip,1] = ρ - ρref
        qn[ip,2] = ρ*u - ρref*u
        qn[ip,3] = ρ*v - ρref*v
        qn[ip,4] = ρ*θ - ρref*θ
    else
        qn[ip,1] = ρ
        qn[ip,2] = ρ*u
        qn[ip,3] = ρ*v
        qn[ip,4] = ρ*θ
    end
    qn[ip,end] = p

    qe[ip,1] = ρref
    qe[ip,2] = ρref*u
    qe[ip,3] = ρref*v
    qe[ip,4] = ρref*θ
    qe[ip,end] = p
end