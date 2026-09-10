#---------------------------------------------------------------------------------
# Brio-Wu MHD shock tube (Brio & Wu, JCP 75:400, 1988), the setting of Dao &
# Nazarov (2022, JSC 92:77), Sec. 5.2:
#
#     (ρ, u, p, Bx, By) = (1,     0, 1,   0.75,  1)   x ∈ [0, 0.5)
#                         (0.125, 0, 0.1, 0.75, -1)   x ∈ [0.5, 1]
#
# v = w = Bz = 0, γ = 2, domain (0, 1), final time 0.1 (Brio & Wu's 0.2 on
# (−1, 1), the state of the paper's Fig. 2). Output variables
# (user_uout! in user_primitives.jl): ρ, u, v, p, By — the panels of the
# classical Brio-Wu figure — against the reference solution of
# reference_hll.dat (see user_analytic.jl).
#---------------------------------------------------------------------------------
function initialize(SD, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    println(" Initialize fields for 1D ideal MHD (brioWu1d: Brio-Wu shock tube) ............ ")

    qvars    = ["ρ", "ρu", "ρv", "ρE", "ρw", "Bx", "By", "Bz"]
    qoutvars = ["ρ", "u", "v", "p", "By"]
    q = define_q(SD,
                 mesh.nelem, mesh.npoin, mesh.ngl,
                 qvars,
                 TFloat, inputs[:backend];
                 neqs=length(qvars), qoutvars=qoutvars)

    ρL, uL, pL, ByL = 1.0,   0.0, 1.0,  1.0
    ρR, uR, pR, ByR = 0.125, 0.0, 0.1, -1.0
    Bx = 0.75
    x0 = 0.5

    for iel_g = 1:mesh.nelem
        for i = 1:mesh.ngl
            ip = mesh.connijk[iel_g, i, 1, 1]
            x  = mesh.coords[1, ip]
            if x < x0
                ρ, u, p, By = ρL, uL, pL, ByL
            else
                ρ, u, p, By = ρR, uR, pR, ByR
            end
            E = p/(γ_mhd - 1.0) + 0.5*ρ*u*u + 0.5*(Bx*Bx + By*By)
            q.qn[ip, 1] = ρ
            q.qn[ip, 2] = ρ*u
            q.qn[ip, 3] = 0.0
            q.qn[ip, 4] = E
            q.qn[ip, 5] = 0.0
            q.qn[ip, 6] = Bx
            q.qn[ip, 7] = By
            q.qn[ip, 8] = 0.0
            for ieq = 1:8
                q.qe[ip, ieq] = q.qn[ip, ieq]
            end
        end
    end

    println(" Initialize fields for 1D ideal MHD (brioWu1d: Brio-Wu shock tube) ............ DONE ")
    return q
end
