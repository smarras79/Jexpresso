function initialize(SD::NSD_2D, PT::CompEuler, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
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
    
        xc = minimum(mesh.x)#7000.0#(maximum(mesh.x) + minimum(mesh.x))/2
        yc = 3000.0 #m
        rx = 4000.0 #m
        ry = 2000.0 #m
    
        θref = 300.0 #K
        θc   = -15.0 #K
        for iel_g = 1:mesh.nelem
            for j=1:mesh.ngl, i=1:mesh.ngl
            
                ip = mesh.connijk[iel_g,i,j]
                x, y = mesh.x[ip], mesh.y[ip]
                r = sqrt( (x - xc)^2/rx^2 + (y - yc)^2/ry^2 )
            
                Δθ = 0.0 #K
                if r <= 1.0
                    Δθ = 0.5*θc*(1.0 + cos(π*r))
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
                    q.qn[ip,2] = ρ*u
                    q.qn[ip,3] = ρ*v
                    q.qn[ip,4] = ρ*θ - ρref*θref
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
        
    else
        if (inputs[:SOL_VARS_TYPE] == PERT())
            lpert = true
        else
            lpert = false
        end
        PhysConst = PhysicalConst{TFloat}()
        xc = TFloat(minimum(mesh.x))
        yc = TFloat(3000.0) #m
        rx = TFloat(4000.0) #m
        ry = TFloat(2000.0)
        θref = TFloat(300.0) #K
        θc   =   TFloat(-15.0) #K
        k = initialize_gpu!(inputs[:backend])
        k(q.qn, q.qe, mesh.x, mesh.y, xc, rx, ry, yc, θref, θc, PhysConst,lpert; ndrange = (mesh.npoin))
    end
    @info " Initialize fields for 2D CompEuler with θ equation ........................ DONE "
    
    return q
end

@kernel function initialize_gpu!(qn, qe, x, y, xc, rx, ry, yc, θref, θc, PhysConst,lpert)
    ip = @index(Global, Linear)

    T = eltype(x)
    x = x[ip]
    y = y[ip]
    r = sqrt( (x - xc)^2/rx^2 + (y - yc)^2/ry^2 )

    Δθ = T(0.0) #K
    if r <= T(1.0)
        Δθ = T(0.5)*θc*(T(1.0) + cos(T(π)*r))
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
