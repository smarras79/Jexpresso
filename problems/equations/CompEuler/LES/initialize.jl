function initialize(SD::NSD_3D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0
        @info " Initialize analytic neutral ABL ........................ "
    end
    
    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is the second dimension of q = define_q()
    # 
    #---------------------------------------------------------------------------------
    qvars    = ["ρ", "ρu", "ρv", "ρw", "ρθ"]
    qoutvars = ["ρ", "u", "v", "w", "θ", "p"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)
    #---------------------------------------------------------------------------------
    if (inputs[:backend] == CPU())    
        PhysConst = PhysicalConst{Float64}()
        
        #
        # INITIAL STATE from scratch:
        #    
        zi  = 840.0 #m
        θ1  = 289.0 #K
        θ2  = 297.5 #K
        Ts  = 290.4
        H   = PhysConst.Rair*Ts/PhysConst.g
        amp = 0.25
        for ip = 1:mesh.npoin
            
            z = mesh.z[ip]

            if z <= zi
                θ = θ1
            else
                θ = θ2 + cbrt(z - zi)
            end
            p    = PhysConst.pref*exp(-z/H)
            pref = p

            randnoise = 0.0
            if z < 800.0
                randnoise = 2*amp*(rand() - 1.0)
            end
            θ    = θ + randnoise
            
            ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ, Press=p)    #kg/m³
            ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θ, Press=pref) #kg/m³

            u = 10.0
            v = 0.0
            w = 0.0

            if inputs[:SOL_VARS_TYPE] == PERT()
                q.qn[ip,1] = ρ - ρref
                q.qn[ip,2] = ρ*u - ρref*u
                q.qn[ip,3] = ρ*v - ρref*v
                q.qn[ip,4] = ρ*w - ρref*w
                q.qn[ip,5] = ρ*θ - ρref*θ
                q.qn[ip,end] = p
                
                #Store initial background state for plotting and analysis of pertuebations
                q.qe[ip,1] = 0.0*ρref
                q.qe[ip,2] = 0.0*u
                q.qe[ip,3] = 0.0*v
                q.qe[ip,4] = 0.0*w
                q.qe[ip,5] = 0.0*ρref*θ
                q.qe[ip,end] = 0.0*pref
            else
                q.qn[ip,1] = ρ
                q.qn[ip,2] = ρ*u
                q.qn[ip,3] = ρ*v
                q.qn[ip,4] = ρ*w
                q.qn[ip,5] = ρ*θ
                q.qn[ip,end] = p

                #Store initial background state for plotting and analysis of pertuebations
                q.qe[ip,1] = 0.0*ρref
                q.qe[ip,2] = 0.0*u
                q.qe[ip,3] = 0.0*v
                q.qe[ip,4] = 0.0*w
                q.qe[ip,5] = 0.0*ρref*θ
                q.qe[ip,end] = 0.0*pref
            end
        end
    end
    if rank == 0
        @info " Initialize analytic neutral ABL ........................ DONE "
    end
    
    return q
end


function user_get_adapt_flags(inputs, old_ad_lvl, q, qe, connijk, nelem, ngl)
    adapt_flags = KernelAbstractions.zeros(CPU(), TInt, Int64(nelem))
    ips         = KernelAbstractions.zeros(CPU(), TInt, ngl * ngl)
    tol         = 1.0
    max_level   = inputs[:amr_max_level] 
    
    for iel = 1:nelem
        m = 1
        for i = 1:ngl
            for j = 1:ngl
                ips[m] = connijk[iel, i, j]
                m += 1
            end
        end
        # @info q[ips,4] - qe[ips,4]
        theta      = q[ips, 4] ./ q[ips, 1]
        theta_ref  = qe[ips, 4] ./ qe[ips, 1]
        dtheta     = theta - theta_ref
        # @info dtheta
        if any(dtheta .> tol) && (old_ad_lvl[iel] < max_level)
            adapt_flags[iel] = refine_flag
        end
        if all(dtheta .< tol)
            adapt_flags[iel] = coarsen_flag
        end
    end
    return adapt_flags
end
