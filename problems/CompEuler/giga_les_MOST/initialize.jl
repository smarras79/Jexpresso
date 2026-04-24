function initialize(SD::NSD_3D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)
    """

            """
    @info " Initialize fields for 3D CompEuler with θ equation ........................ "
    
    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is the second dimension of q = define_q()
    # 
    #---------------------------------------------------------------------------------
    qvars = ("ρ", "ρu", "ρv", "ρw", "hl", "ρqt", "ρqp")
    qoutvars = ["ρ", "ρu", "ρv", "ρw", "hl", "ρqt", "ρqp", "T", "qn", "qc", "qi", "qr", "qs", "qg", "qsatt","u_prime", "v_prime", "w_prime", "hl_prime", "pressure"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)
    #---------------------------------------------------------------------------------
    
    if (inputs[:backend] == CPU())
        PhysConst = PhysicalConst{Float64}()
        if inputs[:lrestart] == true
            #
            # READ RESTART HDF5:
            #
            q.qn, q.qe = read_output(mesh.SD, inputs[:restart_input_file_path], inputs, mesh.npoin, HDF5(); nvar=length(qvars))
            PhysConst = PhysicalConst{Float64}()
        
            for ip=1:mesh.npoin
                z = mesh.z[ip]
                ρ  = q.qn[ip,1]
                hl = q.qn[ip,5] / ρ
                qv = q.qn[ip,6] / ρ
                T  = (hl - PhysConst.g*z) / PhysConst.cp
                Tv = T*(1+0.61*qv)
                q.qn[ip,end] = ρ*Tv*PhysConst.Rair
            
                ρe  = q.qe[ip,1]
                hle = q.qe[ip,5] / ρe
                qve = q.qe[ip,6] / ρe
                Te  = (hle - PhysConst.g*z) / PhysConst.cp
                Tve = Te*(1+0.61*qve)
                q.qe[ip,end] = ρe*Tve*PhysConst.Rair
            end
        
        else
            #
            # INITIAL STATE from scratch:
            #
            data       = read_sounding(inputs[:sounding_file])
            background = interpolate_sounding(inputs[:backend],mesh.npoin,mesh.z,data) 
            
            # data_u                 = read_sounding("./data_files/GLES_initial_u.dat")
            # data_u_reordered       = zeros(size(data_u))
            # data_u_reordered[:,1] .= data_u[:,2]*1000
            # data_u_reordered[:,2] .= data_u[:,1]
            # background_u           = interpolate_sounding(inputs[:backend],mesh.npoin,mesh.z,data_u_reordered)

            # data_qv                 = read_sounding("./data_files/GLES_initial_qv.dat")
            # data_qv_reordered       = zeros(size(data_qv))
            # data_qv_reordered[:,1] .= data_qv[:,2]*1000
            # data_qv_reordered[:,2] .= data_qv[:,1]
            # background_qv           = interpolate_sounding(inputs[:backend],mesh.npoin,mesh.z,data_qv_reordered)
            
            amp = 0.1 #K

            for ip = 1:mesh.npoin
            
                x, y, z = mesh.x[ip], mesh.y[ip], mesh.z[ip]

                rand_noise = 0.0 #K
                if z < 1000.0 # change to 300m later to be consistent to the ref paper
                    rand_noise = 2*amp*(rand() - 1.0)
                end

                # #else
                # qv_ref = background_qv[ip,1]/1000
                # #end
                # u_ref  = background_u[ip,1]
                # pref   = background[ip,1]
                # T      = T_ref + rand_noise
                # Tv     = T*(1+0.61*qv_ref) 
                # Tv_ref = T_ref*(1+0.61*qv_ref) 
                # ρ      = perfectGasLaw_TPtoρ(PhysConst; Temp=T,    Press=pref)    #kg/m³
                # # ρ      = perfectGasLaw_TPtoρ(PhysConst; Temp=Tv,    Press=pref)    #kg/m³
                # # ρref   = perfectGasLaw_TPtoρ(PhysConst; Temp=Tv_ref, Press=pref) #kg/m³
                # ρref   = perfectGasLaw_TPtoρ(PhysConst; Temp=T_ref, Press=pref) #kg/m³
                # # hl is the liqui/ice static energy used in SAM
                # # use potential is less conservative
                # hl     = PhysConst.cp*T + PhysConst.g*z
                # hl_ref = PhysConst.cp*T_ref + PhysConst.g*z
                # u      = u_ref
                # v_ref  = 0.0
                # v      = v_ref
                # w      = 0.0
                # # pref_m = ρref*Tv_ref*PhysConst.Rair
                # # dry pressure
                # pref_m = ρref*T_ref*PhysConst.Rair#ρref*PhysConst.Rair*Tref + ρref*qv_ref*PhysConst.Rvap*Tref

                

                u_ref  = background[ip,4]
                v_ref  = background[ip,5]
                w_ref  = 0.0
                pref   = background[ip,1]
                θ_ref  = background[ip,2]
                qv_ref = background[ip,3]/1000

                T_ref  = θ_ref*(pref/PhysConst.pref)^(1/PhysConst.cpoverR)
                Tv_ref = T_ref*(1+0.61*qv_ref)

                
                ρref   = perfectGasLaw_TPtoρ(PhysConst; Temp=Tv_ref, Press=pref) #kg/m³
                hl_ref = PhysConst.cp*Tv_ref + PhysConst.g*z
                pref_m = ρref*Tv_ref*PhysConst.Rair

                u      = u_ref
                v      = v_ref
                w      = w_ref
                qv     = qv_ref

                T      = T_ref + rand_noise
                Tv     = T*(1+0.61*qv)

                ρ      = perfectGasLaw_TPtoρ(PhysConst; Temp=Tv, Press=pref) #kg/m³
                hl     = PhysConst.cp*Tv + PhysConst.g*z



                if inputs[:SOL_VARS_TYPE] == PERT()
                    q.qn[ip,1] = ρ - ρref
                    q.qn[ip,2] = ρ*u - ρref*u
                    q.qn[ip,3] = ρ*v - ρref*v
                    q.qn[ip,4] = ρ*w - ρref*w
                    q.qn[ip,5] = ρ*hl - ρref*hl_ref#ρ*θ - ρref*θref
                    q.qn[ip,6] = ρ*qv_ref-ρref*qv_ref
                    q.qn[ip,7] = 0.0
                    q.qn[ip,end] = pref_m #+ ρ*qv_ref*PhysConst.Rvap*T

                
                    #Store initial background state for plotting and analysis of pertuebations
                    q.qe[ip,1] = ρref
                    q.qe[ip,2] = ρref*u_ref
                    q.qe[ip,3] = ρref*v_ref
                    q.qe[ip,4] = ρref*w
                    q.qe[ip,5] = ρref*hl_ref
                    q.qe[ip,6] = ρref*qv_ref
                    q.qe[ip,7] = 0.0
                    q.qe[ip,end] = pref_m #+ ρref*qv*PhysConst.Rvap*Tref
                else
                    q.qn[ip,1] = ρ
                    q.qn[ip,2] = ρ*u
                    q.qn[ip,3] = ρ*v
                    q.qn[ip,4] = ρ*w
                    q.qn[ip,5] = ρ*hl
                    q.qn[ip,6] = 0.0
                    q.qn[ip,6] = ρ*qv#0.0
                    q.qn[ip,7] = 0.0
                    q.qn[ip,end] = pref_m #+ ρ*qv*PhysConst.Rvap*T

                    #Store initial background state for plotting and analysis of pertuebations
                    q.qe[ip,1] = ρref
                    q.qe[ip,2] = ρref*u_ref
                    q.qe[ip,3] = ρref*v_ref
                    q.qe[ip,4] = ρref*w_ref
                    q.qe[ip,5] = ρref*hl_ref
                    q.qe[ip,6] = ρref*qv_ref
                    q.qe[ip,7] = 0.0
                    q.qe[ip,end] = pref_m #+ ρref*qv*PhysConst.Rvap*Tref
                end
                #end
            end
        end


    else
        if (inputs[:SOL_VARS_TYPE] == PERT())
            lpert = true
        else
            lpert = false
        end
        data = read_sounding(inputs[:sounding_file])
        background = interpolate_sounding(inputs[:backend],mesh.npoin,mesh.z,data)
        PhysConst = PhysicalConst{TFloat}()
        xc = TFloat((maximum(mesh.x) + minimum(mesh.x))/2)
        zc = TFloat(2000.0) #m
        rz = TFloat(1500.0) #m
        rx = TFloat(10000.0)
        θref = TFloat(300.0) #K
        θc   =   TFloat(2.0) #K
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
    hl_ref = PhysConst.cp*Tref + PhysConst.g*z1
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
