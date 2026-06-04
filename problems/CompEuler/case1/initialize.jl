function initialize(SD, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)
    @info " Initialize fields for 1D CompEuler (case1: DSGS sound pulse) ........ "

    PhysConst = PhysicalConst{Float64}()

    qvars = ["ρ", "ρu", "ρE"]
    q = define_q(SD,
                 mesh.nelem, mesh.npoin, mesh.ngl,
                 qvars,
                 TFloat, inputs[:backend];
                 neqs=length(qvars))

    γ  = PhysConst.γ
    xs = 1.5
    ωsq = 0.125^2

    for iel_g = 1:mesh.nelem
        for i = 1:mesh.ngl
            ip = mesh.connijk[iel_g, i, 1, 1]
            x  = mesh.coords[ip, 1]

            ρ = 1.0
            u = 0.0
            p = exp(-log(2) * ((x - xs)^2)/ωsq) + 1.0

            q.qn[ip, 1] = ρ
            q.qn[ip, 2] = ρ*u
            q.qn[ip, 3] = p/(γ - 1.0) + 0.5*ρ*u*u

            # Background (unperturbed) state used for plotting / Dirichlet bc
            q.qe[ip, 1] = 1.0
            q.qe[ip, 2] = 0.0
            q.qe[ip, 3] = 1.0/(γ - 1.0)
        end
    end

    @info " Initialize fields for 1D CompEuler (case1: DSGS sound pulse) ........ DONE "

    return q
end
