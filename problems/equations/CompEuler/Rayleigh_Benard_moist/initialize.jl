function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """

            """
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0
        @info " Initialize fields for 2D CompEuler with θ equation ........................ "
    end
    
    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: while these names can be arbitrary, the length of this tuple
    # defines neqs, which is the second dimension of q = define_q()
    # 
    #---------------------------------------------------------------------------------
    qvars    = ["ρ", "ρu", "ρv", "ρhl", "ρqt", "ρqp"]
    qoutvars = ["ρ", "u", "w", "hl", "qt", "qp", "qn", "T"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)
    #---------------------------------------------------------------------------------
    if (inputs[:backend] == CPU())    
        PhysConst = PhysicalConst{Float64}()
        if inputs[:lrestart] == true
            #
            # READ RESTART HDF5:
            #
            q.qn, q.qe = read_output(mesh.SD, inputs[:restart_input_file_path], inputs, mesh.npoin, HDF5(); nvar=length(qvars), qoutvars=qoutvars)
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
            comm = MPI.COMM_WORLD
            max_x = MPI.Allreduce(maximum(mesh.x), MPI.MAX, comm)
            min_x = MPI.Allreduce(minimum(mesh.x), MPI.MIN, comm)
            xc = (max_x + min_x)/2
            yc = 2500.0 #m
            r0 = 2000.0 #m
        
            θref = 300.0 #K
            θc   =   0.0 #K
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
                T = θ/((PhysConst.pref/p)^(PhysConst.Rair/PhysConst.cp))
                u = 0.0
                v = 0.0
                hl = T*PhysConst.cp + PhysConst.g*y
                if inputs[:SOL_VARS_TYPE] == PERT()
                    q.qn[ip,1] = ρ - ρref
                    q.qn[ip,2] = ρ*u - ρref*u
                    q.qn[ip,3] = ρ*v - ρref*v
                    q.qn[ip,4] = 0.0#ρ*θ - ρref*θ
                    q.qn[ip,5] = 0.0
                    q.qn[ip,6] = 0.0
                    q.qn[ip,end] = p
                
                    #Store initial background state for plotting and analysis of pertuebations
                    q.qe[ip,1] = ρref
                    q.qe[ip,2] = u
                    q.qe[ip,3] = v
                    q.qe[ip,4] = ρref*hl#ρref*θref
                    q.qe[ip,5] = 0.0
                    q.qe[ip,6] = 0.0
                    q.qe[ip,end] = pref
                else
                    q.qn[ip,1] = ρ
                    q.qn[ip,2] = ρ*u
                    q.qn[ip,3] = ρ*v
                    q.qn[ip,4] = ρ*hl
                    q.qn[ip,5] = 0.0
                    q.qn[ip,6] = 0.0
                    q.qn[ip,end] = p

                    #Store initial background state for plotting and analysis of pertuebations
                    q.qe[ip,1] = ρref
                    q.qe[ip,2] = u
                    q.qe[ip,3] = v
                    q.qe[ip,4] = ρref*hl
                    q.qe[ip,5] = 0.0
                    q.qe[ip,6] = 0.0
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
    if rank == 0
        @info " Initialize fields for 2D CompEuler with θ equation ........................ DONE "
    end
    # @mystop("my stop at mesh.jl L135")
    
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
