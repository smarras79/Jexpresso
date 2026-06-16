function initialize(SD, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)
    println(" Initialize fields for 1D CompEuler (sod1d: shock tube) ............ ")

    PhysConst = PhysicalConst{Float64}()
    γ = PhysConst.γ

    qvars = ["ρ", "ρu", "ρE"]
    q = define_q(SD,
                 mesh.nelem, mesh.npoin, mesh.ngl,
                 qvars,
                 TFloat, inputs[:backend];
                 neqs=length(qvars))

    # Sod (1978) initial data
    ρL, uL, pL = 1.000, 0.0, 1.0
    ρR, uR, pR = 0.125, 0.0, 0.1
    xshock_initial = 0.5

    for iel_g = 1:mesh.nelem
        for i = 1:mesh.ngl
            ip = mesh.connijk[iel_g, i, 1, 1]
            x  = mesh.coords[ip, 1]

            if x < xshock_initial
                ρ, u, p = ρL, uL, pL
            else
                ρ, u, p = ρR, uR, pR
            end

            q.qn[ip, 1] = ρ
            q.qn[ip, 2] = ρ*u
            q.qn[ip, 3] = p/(γ - 1.0) + 0.5*ρ*u*u

            # Background state used by Dirichlet BCs / plotting (left state on
            # the left half, right state on the right half).
            q.qe[ip, 1] = ρ
            q.qe[ip, 2] = ρ*u
            q.qe[ip, 3] = p/(γ - 1.0) + 0.5*ρ*u*u
        end
    end

    println(" Initialize fields for 1D CompEuler (sod1d: shock tube) ............ DONE ")

    return q
end
