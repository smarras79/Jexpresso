using Base

function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)
    """

                """
    @info " Initialize fields for 2D CompEuler with θ equation ........................ "

    qvars = ["dρ", "dρu", "dρv", "dρθ"]
    qoutvars = ["ρ", "ρu", "ρv", "θ", "θ_prime", "Press"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)

    data       = read_sounding(inputs[:sounding_file])
    background = interpolate_sounding(inputs[:backend], mesh.npoin, mesh.y, data)

    if (inputs[:backend] == CPU())

        PhysConst = PhysicalConst{Float64}()

        for iel_g = 1:mesh.nelem
            for j=1:mesh.ngl, i=1:mesh.ngl

                ip = mesh.connijk[iel_g,i,j]

                θ    = background[ip, 1]
                u    = background[ip, 3]
                v    = 0.0
                p    = background[ip, 5]

                ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ, Press=p)
                ρref = ρ

        

                if inputs[:SOL_VARS_TYPE] == PERT()
                    q.qn[ip,1] = ρ - ρref
                    q.qn[ip,2] = ρ*u - ρref*u
                    q.qn[ip,3] = ρ*v - ρref*v
                    q.qn[ip,4] = ρ*θ - ρref*θ
                    q.qn[ip,end] = p
                else
                    q.qn[ip,1] = ρ
                    q.qn[ip,2] = ρ*u
                    q.qn[ip,3] = ρ*v
                    q.qn[ip,4] = ρ*θ
                    q.qn[ip,end] = p
                end
                #Store initial background state for plotting and analysis of perturbations
                q.qe[ip,1] = ρref
                q.qe[ip,2] = ρref*u
                q.qe[ip,3] = ρref*v
                q.qe[ip,4] = ρref*θ
                q.qe[ip,end] = p
            end
        end

        @info "max θ_prime at t=0: $(maximum(q.qn[:,4]))"
        @info "min θ_prime at t=0: $(minimum(q.qn[:,4]))"

        for iel_g = 1:mesh.nelem_semi_inf
            for j=1:mesh.ngr, i=1:mesh.ngl

                ip = mesh.connijk_lag[iel_g,i,j]

                θ    = background[ip, 1]
                u    = background[ip, 3]
                v    = 0.0
                p    = background[ip, 5]

                ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ, Press=p)
                ρref = ρ

                if inputs[:SOL_VARS_TYPE] == PERT()
                    q.qn[ip,1] = ρ - ρref
                    q.qn[ip,2] = ρ*u - ρref*u
                    q.qn[ip,3] = ρ*v - ρref*v
                    q.qn[ip,4] = ρ*θ - ρref*θ
                    q.qn[ip,end] = p
                else
                    q.qn[ip,1] = ρ
                    q.qn[ip,2] = ρ*u
                    q.qn[ip,3] = ρ*v
                    q.qn[ip,4] = ρ*θ
                    q.qn[ip,end] = p
                end
                #Store initial background state for plotting and analysis of perturbations
                q.qe[ip,1] = ρref
                q.qe[ip,2] = ρref*u
                q.qe[ip,3] = ρref*v
                q.qe[ip,4] = ρref*θ
                q.qe[ip,end] = p
            end
        end

    else
        if (inputs[:SOL_VARS_TYPE] == PERT())
            lpert = true
        else
            lpert = false
        end
        PhysConst = PhysicalConst{TFloat}()

        k = initialize_gpu!(inputs[:backend])
        k(q.qn, q.qe, mesh.x, mesh.y, background, PhysConst, lpert; ndrange = (mesh.npoin))
    end

    @info "Initialize fields for system of 2D CompEuler with θ equation ........................ DONE"
    @info maximum(q.qn[:,1:4]), minimum(q.qn[:,1:4])
    return q
end

@kernel function initialize_gpu!(qn, qe, x, y, background, PhysConst, lpert)
    ip = @index(Global, Linear)

    T = eltype(x)

    θ    = T(background[ip, 1])
    u    = T(background[ip, 3])
    v    = T(0.0)
    p    = T(background[ip, 5])

    ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ, Press=p)
    ρref = ρ

    if lpert
        qn[ip,1] = ρ - ρref
        qn[ip,2] = ρ*u - ρref*u
        qn[ip,3] = ρ*v - ρref*v
        qn[ip,4] = ρ*θ - ρref*θ
    else
        qn[ip,1] = ρ
        qn[ip,2] = ρ*u
        qn[ip,3] = ρ*v
        qn[ip,4] = ρ*θ
    end
    qn[ip,end] = p

    qe[ip,1] = ρref
    qe[ip,2] = ρref*u
    qe[ip,3] = ρref*v
    qe[ip,4] = ρref*θ
    qe[ip,end] = p
end