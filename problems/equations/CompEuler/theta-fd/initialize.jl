function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """

            """
    @info " Initialize fields for 2D CompEuler with θ equation ........................ "
    
    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is the second dimension of q = define_q()
    # 
    #---------------------------------------------------------------------------------
    qvars = ("ρ", "ρu", "ρv", "ρθ")
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
                ρθ = q.qn[ip,4]
                θ  = ρθ/ρ
                P = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
                q.qn[ip,end] = P
            
                ρe  = q.qe[ip,1]
                ρθe = q.qe[ip,4]
                θe  = ρθe/ρ
                Pe = perfectGasLaw_ρθtoP(PhysConst, ρ=ρe, θ=θe)
                q.qe[ip,end] = Pe
            end
        
        else
            #
            # INITIAL STATE from scratch:
            #
            xc = (maximum(mesh.x) + minimum(mesh.x))/2
            yc = 2500.0 #m
            r0 = 2000.0 #m
        
            θref = 300.0 #K
            θc   =   2.0 #K
            for ip = 1:mesh.npoin
            
                x, y = mesh.x[ip], mesh.y[ip]
                r = sqrt( (x - xc)^2 + (y - yc)^2 )
            
                Δθ = 0.0 #K
                if r < r0
                    Δθ = θc*(1.0 - r/r0)
                end
                θ = θref + Δθ
                p    = PhysConst.pref*(1.0 - PhysConst.g*y/(PhysConst.cp*θ))^(PhysConst.cpoverR) #Pa
                pref = PhysConst.pref*(1.0 - PhysConst.g*y/(PhysConst.cp*θref))^(PhysConst.cpoverR)
                ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ,    Press=p)    #kg/m³
                ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θref, Press=pref) #kg/m³

                u = 0.0
                v = 0.0

                if inputs[:SOL_VARS_TYPE] == PERT()
                    q.qn[ip,1] = ρ - ρref
                    q.qn[ip,2] = ρ*u - ρref*u
                    q.qn[ip,3] = ρ*v - ρref*v
                    q.qn[ip,4] = ρ*θ - ρref*θ
                    q.qn[ip,end] = p
                
                    #Store initial background state for plotting and analysis of pertuebations
                    q.qe[ip,1] = ρref
                    q.qe[ip,2] = u
                    q.qe[ip,3] = v
                    q.qe[ip,4] = ρref*θref
                    q.qe[ip,end] = pref
                else
                    q.qn[ip,1] = ρ
                    q.qn[ip,2] = ρ*u
                    q.qn[ip,3] = ρ*v
                    q.qn[ip,4] = ρ*θ
                    q.qn[ip,end] = p

                    #Store initial background state for plotting and analysis of pertuebations
                    q.qe[ip,1] = ρref
                    q.qe[ip,2] = u
                    q.qe[ip,3] = v
                    q.qe[ip,4] = ρref*θref
                    q.qe[ip,end] = pref
                end
                #end
            end
        end
    
        if inputs[:CL] == NCL()
            if inputs[:SOL_VARS_TYPE] == PERT()
                q.qn[:,2] .= q.qn[:,2]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,3] .= q.qn[:,3]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,4] .= q.qn[:,4]./(q.qn[:,1] + q.qe[:,1])
            
                #Store initial background state for plotting and analysis of pertuebations
                q.qe[:,4] .= q.qe[:,4]./q.qe[:,1]
            else
                q.qn[:,2] .= q.qn[:,2]./q.qn[:,1]
                q.qn[:,3] .= q.qn[:,3]./q.qn[:,1]
                q.qn[:,4] .= q.qn[:,4]./q.qn[:,1]

                #Store initial background state for plotting and analysis of pertuebations
                q.qe[:,4] .= q.qe[:,4]./q.qe[:,1]
            end
        end

        if (inputs[:lwrite_initial] == true)
            outvarsref = ("drho_init", "du_init", "dv_init", "dtheta_init", "dp_init")
            write_vtk_ref(SD, mesh, q.qn.-q.qe, "initial_state", inputs[:output_dir]; nvar=length(q.qn[1,:]), outvarsref=outvarsref)
        
            outvarsref = ("rho_ref", "u_ref", "v_ref", "theta_ref", "p_ref")    
            write_vtk_ref(SD, mesh, q.qe, "REFERENCE_state", inputs[:output_dir]; nvar=length(q.qe[1,:]), outvarsref=outvarsref)
        end
    else
        if (inputs[:SOL_VARS_TYPE] == PERT())
            lpert = true
        else
            lpert = false
        end
        PhysConst = PhysicalConst{TFloat}()
        xc = TFloat((maximum(mesh.x) + minimum(mesh.x))/2)
        yc = TFloat(2500.0) #m
        rθ = TFloat(2000.0) #m

        θref = TFloat(300.0) #K
        θc   =   TFloat(2.0) #K
        k = initialize_gpu!(inputs[:backend])
        k(q.qn, q.qe, mesh.x, mesh.y, xc, rθ, yc, θref, θc, PhysConst,lpert; ndrange = (mesh.npoin))
    end
    @info " Initialize fields for 2D CompEuler with θ equation ........................ DONE "
    
    return q
end

@kernel function initialize_gpu!(qn, qe, x, y, xc, rθ, yc, θref, θc, PhysConst,lpert)
    ip = @index(Global, Linear)

    T = eltype(x)
    x = x[ip]
    y = y[ip]
    r = sqrt( (x - xc)^2 + (y - yc)^2 )
    Δθ = T(0.0) #K
    if r < rθ
        Δθ = T(θc*(T(1.0) - r/rθ))
    end
    θ = θref + Δθ
    p    = PhysConst.pref*(T(1.0) - PhysConst.g*y/(PhysConst.cp*θ))^(PhysConst.cpoverR) #Pa
    pref = PhysConst.pref*(T(1.0) - PhysConst.g*y/(PhysConst.cp*θref))^(PhysConst.cpoverR)
    ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ,    Press=p)    #kg/m³
    ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θref, Press=pref) #kg/m³

    u = T(0.0)
    v = T(0.0)
    if (lpert)
        qn[ip,1] = ρ - ρref
        qn[ip,2] = ρ*u
        qn[ip,3] = ρ*v
        qn[ip,4] = ρ*θ - ρref*θref
        qn[ip,end] = p
    else
        qn[ip,1] = ρ
        qn[ip,2] = ρ*u
        qn[ip,3] = ρ*v
        qn[ip,4] = ρ*θ
        qn[ip,end] = p
    end

                    #Store initial background state for plotting and analysis of pertuebations
    qe[ip,1] = ρref
    qe[ip,2] = u
    qe[ip,3] = v
    qe[ip,4] = ρref*θref
    qe[ip,end] = pref

end
