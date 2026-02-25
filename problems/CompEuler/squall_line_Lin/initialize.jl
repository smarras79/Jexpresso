function initialize(SD::NSD_3D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
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
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars))
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
            #
            # INITIAL STATE from scratch:
            #
            xc = (maximum(mesh.x) + minimum(mesh.x))/2
            zc = 2000.0 #m
        
            θc   =   3.0 #K
            qvc = 2e-3
            rx = 2400.0
            rz = 1500.0
            data = read_sounding(inputs[:sounding_file])
            background = interpolate_sounding(inputs[:backend],mesh.npoin,mesh.z,data) 
            balanced = zeros(mesh.npoin,1)
            @info minimum(background[:,3]), maximum(background[:,3])           
            for ip = 1:mesh.npoin
            
                x, y, z = mesh.x[ip], mesh.y[ip], mesh.z[ip]
                 
                r = sqrt( (x - xc)^2/(rx)^2 )
            
                ΔT = 0.0 #K
                Δqv = 0.0
                if (r <= 1) && (z <= 1199.0) && (z>= 350.0)
                    ΔT = θc*sinpi((r+1)/2)/((z-200)/200)
                    Δqv = qvc*sinpi((r+1)/2)/((z-200)/200) 
                end
                T_ref = background[ip,1]
                qv_ref = background[ip,2]/1000
                u_ref = background[ip,3]
                v_ref = background[ip,4]
                pref = background[ip,5]
                T = T_ref + ΔT
                qv = qv_ref + Δqv
                #=θ_ref = T_ref * (PhysConst.pref/pref)^(PhysConst.Rair/PhysConst.cp)
                θ = T * (PhysConst.pref/pref)^(PhysConst.Rair/PhysConst.cp)
                θv_ref = θ_ref*(1 + 0.61*qv_ref)
                θv = θ*(1 + 0.61*qv)=#
                Tv = T*(1 + 0.61*qv)
                Tv_ref = T_ref*(1 + 0.61*qv_ref)
                ρ    = perfectGasLaw_TPtoρ(PhysConst; Temp=Tv,    Press=pref)    #kg/m³
                ρref = perfectGasLaw_TPtoρ(PhysConst; Temp=Tv_ref, Press=pref) #kg/m³
                hl = PhysConst.cp*T + PhysConst.g*z
                hl_ref = PhysConst.cp*T_ref + PhysConst.g*z
                u = u_ref
                v = v_ref
                w = 0.0
                pref_m = ρref*Tv_ref*PhysConst.Rair#ρref*PhysConst.Rair*T_ref + ρref*qv_ref*PhysConst.Rvap*T_ref
                p_m = ρ*PhysConst.Rair*T + ρ*qv_ref*PhysConst.Rvap*T
                ##attempt to assure balance
                #=for iter=1:10
                    θ_ref = T_ref * (PhysConst.pref/p_m)^(PhysConst.Rair/PhysConst.cp)
                    θ = T * (PhysConst.pref/pref_m)^(PhysConst.Rair/PhysConst.cp)
                    θv_ref = θ_ref*(1 + 0.61*qv_ref)
                    θv = θ*(1 + 0.61*qv)
                    ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θv,    Press=p_m)    #kg/m³
                    ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θv_ref, Press=pref_m)
                    pref_m = ρref*PhysConst.Rair*T_ref + ρref*qv_ref*PhysConst.Rvap*T_ref
                    p_m = ρ*PhysConst.Rair*T + ρ*qv_ref*PhysConst.Rvap*T
                    #@info θ_ref, θv_ref, ρref, pref_m
                end=#
                if inputs[:SOL_VARS_TYPE] == PERT()
                    q.qn[ip,1] = ρ - ρref
                    q.qn[ip,2] = ρ*u - ρref*u
                    q.qn[ip,3] = ρ*v - ρref*v
                    q.qn[ip,4] = ρ*w - ρref*w
                    q.qn[ip,5] = ρ*hl - ρref*hl_ref#ρ*θ - ρref*θref
                    q.qn[ip,6] = ρ*qv-ρref*qv_ref
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
                    q.qn[ip,6] = ρ*qv_ref#0.0
                    q.qn[ip,7] = 0.0
                    q.qn[ip,end] = pref_m #+ ρ*qv*PhysConst.Rvap*T

                    #Store initial background state for plotting and analysis of pertuebations
                    q.qe[ip,1] = ρref
                    q.qe[ip,2] = ρref*u_ref
                    q.qe[ip,3] = ρref*v_ref
                    q.qe[ip,4] = ρref*w
                    q.qe[ip,5] = ρref*hl_ref
                    q.qe[ip,6] = ρref*qv_ref
                    q.qe[ip,7] = 0.0
                    q.qe[ip,end] = pref_m #+ ρref*qv*PhysConst.Rvap*Tref
                end
                #end
            end
        end
        @info maximum(q.qn[:,end]), minimum(q.qn[:,end]), maximum(q.qn[:,2]), minimum(q.qn[:,2]), maximum(q.qn[:,6]), minimum(q.qn[:,6]), maximum(q.qn[:,1]), minimum(q.qn[:,1])    
        if inputs[:CL] == NCL()
            if inputs[:SOL_VARS_TYPE] == PERT()
                q.qn[:,2] .= q.qn[:,2]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,3] .= q.qn[:,3]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,4] .= q.qn[:,4]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,5] .= q.qn[:,5]./(q.qn[:,1] + q.qe[:,1])
            
                #Store initial background state for plotting and analysis of pertuebations
                q.qe[:,5] .= q.qe[:,5]./q.qe[:,1]
            else
                q.qn[:,2] .= q.qn[:,2]./q.qn[:,1]
                q.qn[:,3] .= q.qn[:,3]./q.qn[:,1]
                q.qn[:,4] .= q.qn[:,4]./q.qn[:,1]
                q.qn[:,5] .= q.qn[:,5]./q.qn[:,1]

                #Store initial background state for plotting and analysis of pertuebations
                q.qe[:,5] .= q.qe[:,5]./q.qe[:,1]
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
