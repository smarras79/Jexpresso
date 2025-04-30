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
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat; neqs=length(qvars))
    #---------------------------------------------------------------------------------


    PhysConst = PhysicalConst{Float64}()
    if (inputs[:case] === "rtb")
        
        xc = (maximum(mesh.x) + minimum(mesh.x))/2
        yc = 2500.0 #m
        r0 = 2000.0 #m
        
        θref = 300.0 #K
        θc   =   2.0 #K
        for iel_g = 1:mesh.nelem
            for j=1:mesh.ngl, i=1:mesh.ngl
                
                ip = mesh.connijk[iel_g,i,j]
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
    
    @info " Initialize fields for 2D CompEuler with θ equation ........................ DONE "
    
    return q
end
