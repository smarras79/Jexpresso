function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)

    @info " Initialize fields for 3D CompEuler with θ equation ........................ "

    qvars = ("ρ", "ρu", "ρv", "hl", "ρqt", "ρqp")
    qoutvars = ["ρ", "ρu", "ρv", "hl", "ρqt", "ρqp", "T", "qn", "qc", "qi", "qr", "qs", "qg","u_prime", "v_prime", "hl_prime", "qt_prime"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)

    if (inputs[:backend] == CPU())
        PhysConst = PhysicalConst{Float64}()

        if inputs[:lrestart] == true
            q.qn, q.qe = read_output(mesh.SD, inputs[:restart_input_file_path], inputs, mesh.npoin, HDF5(); nvar=length(qvars))
            PhysConst = PhysicalConst{Float64}()
            for ip=1:mesh.npoin
                ρ  = q.qn[ip,1]
                ρθ = q.qn[ip,5]
                θ  = ρθ/ρ
                P = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
                q.qn[ip,end] = P
                ρe  = q.qe[ip,1]
                ρθe = q.qe[ip,5]
                θe  = ρθe/ρ
                Pe = perfectGasLaw_ρθtoP(PhysConst, ρ=ρe, θ=θe)
                q.qe[ip,end] = Pe
            end

        else
            xc = 0.0
            zc = 2000.0
            θc   =   3.0
            rx = 10000.0
            rz = 1500.0

            data       = read_sounding(inputs[:sounding_file])
            background = interpolate_sounding(inputs[:backend], mesh.npoin, mesh.y, data)

            # ----------------------------------------------------------------
            # HYDROSTATIC REBALANCING
            # Integrate dP/dz = -ρg upward for each unique x-column
            # using iterative predictor-corrector (20 iterations per level)
            # ----------------------------------------------------------------
            p_hydro = zeros(Float64, mesh.npoin)

            # find unique x-columns
            x_unique = sort(unique(round.(mesh.x, digits=2)))

            for xval in x_unique
                col_idx = findall(i -> abs(mesh.x[i] - xval) < 0.01, 1:mesh.npoin)
                # sort by height
                sort!(col_idx, by = i -> mesh.y[i])

                # surface pressure from sounding at lowest point
                ip0    = col_idx[1]
                z0     = mesh.y[ip0]
                qv0    = background[ip0, 2] / 1000.0
                θ0     = background[ip0, 1]
                Tv0    = θ0 * (1.0 + 0.61*qv0)
                p_hydro[ip0] = background[ip0, 5]   # use sounding P at surface

                # integrate upward
                for k = 2:length(col_idx)
                    ip_lo = col_idx[k-1]
                    ip_hi = col_idx[k]
                    dz    = mesh.y[ip_hi] - mesh.y[ip_lo]

                    qv_lo = background[ip_lo, 2] / 1000.0
                    θ_lo  = background[ip_lo, 1]
                    Tv_lo = θ_lo * (1.0 + 0.61*qv_lo)
                    T_lo  = Tv_lo / (1.0 + 0.61*qv_lo) * (1.0 + 0.61*qv_lo)  # = Tv_lo
                    ρ_lo  = perfectGasLaw_TPtoρ(PhysConst; Temp=Tv_lo, Press=p_hydro[ip_lo])

                    # predictor
                    p_hi = p_hydro[ip_lo] - ρ_lo * PhysConst.g * dz

                    # corrector iterations
                    for _ = 1:20
                        qv_hi = background[ip_hi, 2] / 1000.0
                        θ_hi  = background[ip_hi, 1]
                        Tv_hi = θ_hi * (1.0 + 0.61*qv_hi)
                        ρ_hi  = perfectGasLaw_TPtoρ(PhysConst; Temp=Tv_hi, Press=p_hi)
                        ρ_avg = 0.5*(ρ_lo + ρ_hi)
                        p_hi  = p_hydro[ip_lo] - ρ_avg * PhysConst.g * dz
                    end
                    p_hydro[ip_hi] = p_hi
                end
            end
            # ----------------------------------------------------------------

            for ip = 1:mesh.npoin
                x, y, z = mesh.x[ip], mesh.y[ip], mesh.z[ip]

                r = sqrt( (x - xc)^2/(rx^2) + (y - zc)^2/(rz^2) )
                Δθ = 0.0
                if r <= 1
                    Δθ = θc * cospi(r/2)^2
                end

                θ_ref  = background[ip, 1]
                qv_ref = background[ip, 2] / 1000.0
                u_ref  = background[ip, 3]
                v_ref  = background[ip, 4]
                pref   = p_hydro[ip]           # ← use hydrostatic pressure

                θ_ref2 = θ_ref / (1.0 + 0.61*qv_ref)
                θ1     = θ_ref + Δθ
                Tref   = θ_ref / (PhysConst.pref/pref)^(PhysConst.Rair/PhysConst.cp)
                T      = θ1    / (PhysConst.pref/pref)^(PhysConst.Rair/PhysConst.cp)
                Tv     = T    * (1.0 + 0.61*qv_ref)
                Tv_ref = Tref * (1.0 + 0.61*qv_ref)

                ρ    = perfectGasLaw_TPtoρ(PhysConst; Temp=Tv,     Press=pref)
                ρref = perfectGasLaw_TPtoρ(PhysConst; Temp=Tv_ref, Press=pref)

                hl      = PhysConst.cp*Tv     + PhysConst.g*y
                hl_ref  = PhysConst.cp*Tv_ref + PhysConst.g*y
                u = u_ref
                v = 0.0
                pref_m = ρref * Tv_ref * PhysConst.Rair

                if inputs[:SOL_VARS_TYPE] == PERT()
                    q.qn[ip,1] = ρ - ρref
                    q.qn[ip,2] = ρ*u    - ρref*u
                    q.qn[ip,3] = ρ*v    - ρref*v_ref
                    q.qn[ip,4] = ρ*hl   - ρref*hl_ref
                    q.qn[ip,5] = ρ*qv_ref - ρref*qv_ref
                    q.qn[ip,6] = 0.0
                    q.qn[ip,end] = pref_m

                    q.qe[ip,1] = ρref
                    q.qe[ip,2] = ρref*u_ref
                    q.qe[ip,3] = ρref*v_ref
                    q.qe[ip,4] = ρref*hl_ref
                    q.qe[ip,5] = ρref*qv_ref
                    q.qe[ip,6] = 0.0
                    q.qe[ip,end] = pref_m
                else
                    q.qn[ip,1] = ρ
                    q.qn[ip,2] = ρ*u
                    q.qn[ip,3] = ρ*v
                    q.qn[ip,4] = ρ*hl
                    q.qn[ip,5] = ρ*qv_ref
                    q.qn[ip,6] = 0.0
                    q.qn[ip,end] = pref_m

                    q.qe[ip,1] = ρref
                    q.qe[ip,2] = ρref*u_ref
                    q.qe[ip,3] = ρref*v_ref
                    q.qe[ip,4] = ρref*hl_ref
                    q.qe[ip,5] = ρref*qv_ref
                    q.qe[ip,6] = 0.0
                    q.qe[ip,end] = pref_m
                end
            end
        end

        if inputs[:CL] == NCL()
            if inputs[:SOL_VARS_TYPE] == PERT()
                q.qn[:,2] .= q.qn[:,2]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,3] .= q.qn[:,3]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,4] .= q.qn[:,4]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,5] .= q.qn[:,5]./(q.qn[:,1] + q.qe[:,1])
                q.qe[:,5] .= q.qe[:,5]./q.qe[:,1]
            else
                q.qn[:,2] .= q.qn[:,2]./q.qn[:,1]
                q.qn[:,3] .= q.qn[:,3]./q.qn[:,1]
                q.qn[:,4] .= q.qn[:,4]./q.qn[:,1]
                q.qn[:,5] .= q.qn[:,5]./q.qn[:,1]
                q.qe[:,5] .= q.qe[:,5]./q.qe[:,1]
            end
        end

    else
        # GPU branch — unchanged
        if (inputs[:SOL_VARS_TYPE] == PERT())
            lpert = true
        else
            lpert = false
        end
        data = read_sounding(inputs[:sounding_file])
        background = interpolate_sounding(inputs[:backend], mesh.npoin, mesh.z, data)
        PhysConst = PhysicalConst{TFloat}()
        xc = TFloat((maximum(mesh.x) + minimum(mesh.x))/2)
        zc = TFloat(2000.0)
        rz = TFloat(1500.0)
        rx = TFloat(10000.0)
        θref = TFloat(300.0)
        θc   = TFloat(2.0)
        k = initialize_gpu!(inputs[:backend])
        k(q.qn, q.qe, background, mesh.x, mesh.y, mesh.z, xc, rx, rz, zc, θc, PhysConst, lpert; ndrange = (mesh.npoin))
    end

    @info maximum(q.qe[:,end]), minimum(q.qe[:,end])
    @info " Initialize fields for 3D CompEuler with θ equation ........................ DONE "
    return q
end

@kernel function initialize_gpu!(qn, qe, background, x, y, z, xc, rx, rz, zc, θc, PhysConst, lpert)
    ip = @index(Global, Linear)

    T = eltype(x)
    x1 = x[ip]
    y1 = y[ip]
    z1 = z[ip]

    r = T(sqrt( (x1 - xc)^2/(rx^2) + (z1 - zc)^2/(rz^2) ))

    Δθ = T(0.0) #K
    if r <= 1
        Δθ = T(θc*cospi(r/2)^2)
    end
    θ_ref = background[ip,1]
    qv_ref = T(background[ip,2]/T(1000))
    u_ref = background[ip,3]
    v_ref = background[ip,4]
    pref = background[ip,5]
    θv_ref = θ_ref*(T(1) + T(0.608)*qv_ref)
    θ = θ_ref + Δθ
    θv = θv_ref + Δθ
    p    = PhysConst.pref*(T(1.0) - PhysConst.g*T(z1/(PhysConst.cp*θv)))^(PhysConst.cpoverR) #Pa
    Tref = T(θ_ref / (PhysConst.pref/pref)^(PhysConst.Rair/PhysConst.cp))
    Tabs = T(θ / (PhysConst.pref/pref)^(PhysConst.Rair/PhysConst.cp))
    ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θv,    Press=pref)    #kg/m³
    ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θv_ref, Press=pref) #kg/m³
    hl = PhysConst.cp*Tabs + PhysConst.g*z1
    hl_ref = PhysConst.cp*Tvref + PhysConst.g*z1
    u = u_ref
    v = v_ref
    w = T(0.0)
    pref_m = ρref*PhysConst.Rair*Tref + ρref*qv_ref*PhysConst.Rvap*Tref

    if (lpert)
        qn[ip,1] = T(ρ - ρref)
        qn[ip,2] = ρ*u - ρref*u
        qn[ip,3] = ρ*v - ρref*v
        qn[ip,4] = ρ*w - ρref*w
        qn[ip,5] = ρ*hl - ρref*hl_ref
        qn[ip,6] = ρ*qv_ref - ρref*qv_ref
        qn[ip,7] = T(0.0)
        qn[ip,end] = pref_m
    else
        qn[ip,1] = ρ
        qn[ip,2] = ρ*u
        qn[ip,3] = ρ*v
        qn[ip,4] = ρ*w
        qn[ip,5] = ρ*hl
        qn[ip,6] = ρ*qv_ref
        qn[ip,7] = T(0.0)
        qn[ip,end] = pref_m
    end

                    #Store initial background state for plotting and analysis of pertuebations
    qe[ip,1] = ρref
    qe[ip,2] = ρref*u
    qe[ip,3] = ρref*v
    qe[ip,4] = ρref*w
    qe[ip,5] = ρref*hl_ref
    qe[ip,6] = ρref*qv_ref
    qe[ip,7] = 0.0
    qe[ip,end] = pref_m

end
