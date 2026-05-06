function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)

    @info " Initialize fields for 2D CompEuler with θ equation ........................ "

    qvars    = ["ρ", "ρu", "ρv", "ρθ", "qtr", "qtr2"]
    qoutvars = ["ρ", "u", "v", "θ", "qtr", "qtr2"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)

    if (inputs[:backend] == CPU())
        PhysConst = PhysicalConst{Float64}()
        if (inputs[:case] === "rtb")
            comm  = MPI.COMM_WORLD
            max_x = MPI.Allreduce(maximum(mesh.x), MPI.MAX, comm)
            min_x = MPI.Allreduce(minimum(mesh.x), MPI.MIN, comm)

            # ── Jang & Hong bubble parameters ──────────────────────────────
            xc     = (max_x + min_x)/2   # x center = 10km (domain center)
            yc     = 2000.0              # z center = 2km
            R      = 2000.0              # radius = 2km (circular)
            θref   = 300.0               # K isentropic background
            θc     = 2.0                 # K perturbation amplitude
            # tracers not used in Jang & Hong — set to zero
            qtrref = 0.0
            qtrc   = 0.0
            xc2    = 0.0
            yc2    = 0.0
            r02    = 0.0
            qtrc2  = 0.0

            for iel_g = 1:mesh.nelem
                for j=1:mesh.ngl, i=1:mesh.ngl

                    ip = mesh.connijk[iel_g,i,j]
                    x, y = mesh.x[ip], mesh.y[ip]

                    # Jang & Hong: circular bubble, linear profile
                    r    = sqrt((x - xc)^2 + (y - yc)^2)
                    Δθ   = 0.0
                    Δqtr = 0.0
                    if r < R
                        Δθ = θc*(1.0 - r/R)
                    end

                    # tracer 2 — inactive
                    Δqtr2 = 0.0

                    qtr  = qtrref + Δqtr
                    qtr2 = qtrref + Δqtr2

                    θ    = θref + Δθ
                    p    = PhysConst.pref*(1.0 - PhysConst.g*y/(PhysConst.cp*θ))^(PhysConst.cpoverR)
                    pref = PhysConst.pref*(1.0 - PhysConst.g*y/(PhysConst.cp*θref))^(PhysConst.cpoverR)
                    ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ,    Press=p)
                    ρref = perfectGasLaw_θPtoρ(PhysConst; θ=θref, Press=pref)

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
                q.qe[:,4] .= q.qe[:,4]./q.qe[:,1]
            else
                q.qn[:,2] .= q.qn[:,2]./q.qn[:,1]
                q.qn[:,3] .= q.qn[:,3]./q.qn[:,1]
                q.qn[:,4] .= q.qn[:,4]./q.qn[:,1]
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
        xc  = TFloat((maximum(mesh.x) + minimum(mesh.x))/2)
        yc  = TFloat(2000.0)
        rθ  = TFloat(2000.0)
        θref = TFloat(300.0)
        θc   = TFloat(2.0)
        k = initialize_gpu!(inputs[:backend])
        k(q.qn, q.qe, mesh.x, mesh.y, xc, rθ, yc, θref, θc, PhysConst, lpert; ndrange=(mesh.npoin))
    end

    @info " Initialize fields for 2D CompEuler with θ equation ........................ DONE "
    return q
end