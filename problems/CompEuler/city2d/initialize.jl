function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0
        @info " Initialize fields for 2D city2d (flow around buildings) ........................ "
    end

    #---------------------------------------------------------------------------------
    # Solution variables
    #---------------------------------------------------------------------------------
    qvars    = ["ρ", "ρu", "ρv", "ρθ"]
    qoutvars = ["ρ", "u", "w", "θ", "p"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)
    #---------------------------------------------------------------------------------

    if (inputs[:backend] == CPU())
        PhysConst = PhysicalConst{Float64}()
        if inputs[:lrestart] == true
            q.qn, q.qe = read_output(mesh.SD, inputs[:restart_input_file_path], inputs, mesh.npoin, HDF5(); nvar=length(qvars))
            for ip=1:mesh.npoin
                ρ  = q.qn[ip,1]
                ρθ = q.qn[ip,4]
                θ  = ρθ/ρ
                P  = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
                q.qn[ip,end] = P

                ρe  = q.qe[ip,1]
                ρθe = q.qe[ip,4]
                θe  = ρθe/ρ
                Pe  = perfectGasLaw_ρθtoP(PhysConst, ρ=ρe, θ=θe)
                q.qe[ip,end] = Pe
            end
        else
            #
            # INITIAL STATE: neutral atmosphere with uniform horizontal inflow.
            #
            θref = 300.0   # K
            u0   = 10.0    # m/s

            for ip = 1:mesh.npoin
                x, y = mesh.x[ip], mesh.y[ip]

                θ    = θref
                p    = PhysConst.pref*(1.0 - PhysConst.g*y/(PhysConst.cp*θ))^(PhysConst.cpoverR)
                pref = p
                ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ,    Press=p)
                ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θref, Press=pref)

                u = u0
                v = 0.0

                if inputs[:SOL_VARS_TYPE] == PERT()
                    q.qn[ip,1] = ρ - ρref
                    q.qn[ip,2] = ρ*u - ρref*u
                    q.qn[ip,3] = ρ*v - ρref*v
                    q.qn[ip,4] = ρ*θ - ρref*θref
                    q.qn[ip,end] = p
                else
                    q.qn[ip,1] = ρ
                    q.qn[ip,2] = ρ*u
                    q.qn[ip,3] = ρ*v
                    q.qn[ip,4] = ρ*θ
                    q.qn[ip,end] = p
                end

                # Background state used by the sponge: store conservative
                # equilibrium values so (q - qe) is the perturbation.
                q.qe[ip,1]   = ρref
                q.qe[ip,2]   = ρref*u
                q.qe[ip,3]   = ρref*v
                q.qe[ip,4]   = ρref*θref
                q.qe[ip,end] = pref
            end
        end

        if inputs[:CL] == NCL()
            if inputs[:SOL_VARS_TYPE] == PERT()
                q.qn[:,2] .= q.qn[:,2]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,3] .= q.qn[:,3]./(q.qn[:,1] + q.qe[:,1])
                q.qn[:,4] .= q.qn[:,4]./(q.qn[:,1] + q.qe[:,1])
                q.qe[:,4] .= q.qe[:,4]./q.qe[:,1]
            else
                q.qn[:,2] .= q.qn[:,2]./q.qn[:,1]
                q.qn[:,3] .= q.qn[:,3]./q.qn[:,1]
                q.qn[:,4] .= q.qn[:,4]./q.qn[:,1]
                q.qe[:,4] .= q.qe[:,4]./q.qe[:,1]
            end
        end
    end

    if rank == 0
        @info " Initialize fields for 2D city2d (flow around buildings) ........................ DONE "
    end

    return q
end
