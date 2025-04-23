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
    qvars    = ("ρ", "ρu", "ρv", "ρθ", "qtr", "qtr2")
    qoutvars = ("A", "B", "C", "ρθ", "YY", "ZZ")
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, qoutvars, TFloat, inputs[:backend]; neqs=length(qvars))
    #---------------------------------------------------------------------------------
    if (inputs[:backend] == CPU())    
        PhysConst = PhysicalConst{Float64}()
        if (inputs[:case] === "rtb")
        
            xc = (maximum(mesh.x) + minimum(mesh.x))/2
            yc = 2500.0 #m
            r0 = 2000.0 #m

            xc2   = -2000.0
            yc2   =  4000.0 #m
            r02   =  1000.0 #m
            qtrc2 =     1.0 #K
            
            θref = 300.0 #K
            θc   =   2.0 #K
            qtrref = 0.0 #K
            qtrc   = 1.0 #K
            for iel_g = 1:mesh.nelem
                for j=1:mesh.ngl, i=1:mesh.ngl
                
                    ip = mesh.connijk[iel_g,i,j]
                    x, y = mesh.x[ip], mesh.y[ip]
                    r = sqrt( (x - xc)^2 + (y - yc)^2 )
                    Δθ   = 0.0 #K
                    Δqtr = 0.0 #K
                    if r < r0
                        Δθ = θc*(1.0 - r/r0)
                        Δqtr = qtrc
                    end

                #Tracer 2
                    r = sqrt( (x - xc2)^2 + (y - yc2)^2 )
                    Δqtr2 = 0.0 #K
                    if r < r02
                        Δqtr2 = qtrc2
                    end
                    qtr  = qtrref + Δqtr
                    qtr2 = qtrref + Δqtr2

                    θ   = θref + Δθ
                    p   = PhysConst.pref*(1.0 - PhysConst.g*y/(PhysConst.cp*θ))^(PhysConst.cpoverR) #Pa
                    pref = PhysConst.pref*(1.0 - PhysConst.g*y/(PhysConst.cp*θref))^(PhysConst.cpoverR)
                    ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ,    Press=p)    #kg/m³
                    ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θref, Press=pref) #kg/m³

                    u = 0.0
                    v = 0.0

                    if inputs[:SOL_VARS_TYPE] == PERT()
                        q.qn[ip,1] = ρ - ρref
                        q.qn[ip,2] = ρ*u
                        q.qn[ip,3] = ρ*v
                        q.qn[ip,4] = ρ*θ - ρref*θref
                        q.qn[ip,5] = qtr - qtrref
                        q.qn[ip,6] = qtr2 - qtrref
                        q.qn[ip,end] = p
                    
                        #Store initial background state for plotting and analysis of pertuebations
                        q.qe[ip,1] = ρref
                        q.qe[ip,2] = u
                        q.qe[ip,3] = v
                        q.qe[ip,4] = ρref*θref
                        q.qe[ip,5] = qtrref
                        q.qe[ip,6] = qtrref
                        q.qe[ip,end] = pref
                    else
                        q.qn[ip,1] = ρ
                        q.qn[ip,2] = ρ*u
                        q.qn[ip,3] = ρ*v
                        q.qn[ip,4] = ρ*θ
                        q.qn[ip,5] = qtr
                        q.qn[ip,6] = qtr2
                        q.qn[ip,end] = p

                        #Store initial background state for plotting and analysis of pertuebations
                        q.qe[ip,1] = ρref
                        q.qe[ip,2] = u
                        q.qe[ip,3] = v
                        q.qe[ip,4] = ρref*θref
                        q.qe[ip,5] = qtrref
                        q.qe[ip,6] = qtrref
                        q.qe[ip,end] = pref
                    end
                end
            end
        
        else
            error(" ERROR: CompEuler: initialize.jl:\n assign value to inputs[:case]")
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

        #
        # Write reference to VTK:
        #  
        if (inputs[:lwrite_initial] == true)
            outvarsref = Array{Union{Nothing, String}}(nothing, q.neqs)
            for i = 1:length(outvarsref)
                outvarsref[i] = string(qvars[i], "_ref")
            end
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
        
        θref = 300.0 #K
        θc   =   2.0 #K
        
        k = initialize_gpu!(inputs[:backend])
        k(q.qn, q.qe, mesh.x, mesh.y, xc, rθ, yc, θref, θc, PhysConst,lpert; ndrange = (mesh.npoin))
    end
    @info " Initialize fields for 2D CompEuler with θ equation ........................ DONE "
    
    return q
end

@kernel function initialize_gpu!(qn, qe, x, y, xc, rθ, yc, θref, θc, PhysConst,lpert)
    ip = @index(Global, Linear)

    T = eltype(x)
    xc2   = T(-2000.0)
    yc2   =  T(4000.0) #m
    r02   =  T(1000.0) #m
    qtrc2 =     T(1.0) #K

    θref = T(300.0) #K
    θc   =   T(2.0) #K
    qtrref = T(0.0) #K
    qtrc   = T(1.0) #K

    x = x[ip]
    y = y[ip]
    r = sqrt( (x - xc)^2 + (y - yc)^2 )
    Δθ = T(0.0) #K
    Δqtr = T(0.0)
    if r < rθ
        Δθ = T(θc*(T(1.0) - r/rθ))
        Δqtr = qtrc
    end

    r = sqrt( (x - xc2)^2 + (y - yc2)^2 )
    Δqtr2 = T(0.0) #K
    if r < r02
        Δqtr2 = qtrc2
    end
    qtr  = qtrref + Δqtr
    qtr2 = qtrref + Δqtr2

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
        qn[ip,5] = qtr - qtrref
        qn[ip,6] = qtr2 - qtrref
        qn[ip,end] = p
    else
        qn[ip,1] = ρ
        qn[ip,2] = ρ*u
        qn[ip,3] = ρ*v
        qn[ip,4] = ρ*θ
        qn[ip,5] = qtr
        qn[ip,6] = qtr2 
        qn[ip,end] = p
    end

                    #Store initial background state for plotting and analysis of pertuebations
    qe[ip,1] = ρref
    qe[ip,2] = u
    qe[ip,3] = v
    qe[ip,4] = ρref*θref
    qe[ip,5] = qtrref
    qe[ip,6] = qtrref
    qe[ip,end] = pref

end
