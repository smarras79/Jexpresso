function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)

    @info " Initialize fields for 2D CompEuler with θ equation ........................ "

    PhysConst = PhysicalConst{Float64}()

    qvars = ("ρ", "ρu", "ρv", "hl", "ρqt", "ρqp")
    qoutvars = ["ρ", "ρu", "ρv", "hl", "ρqt", "ρqp", "T", "qn", "qc", "qi", "qr", "qs", "qg", "u_prime", "v_prime", "hl_prime", "qt_prime"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)

    data       = read_sounding(inputs[:sounding_file])
    background = interpolate_sounding(inputs[:backend], mesh.npoin, mesh.y, data)

    for iel_g = 1:mesh.nelem
        for j = 1:mesh.ngl, i = 1:mesh.ngl

            ip = mesh.connijk[iel_g, i, j]
            y  = mesh.y[ip]

            θ_ref  = background[ip, 1]
            qv_ref = background[ip, 2] / 1000.0
            u_ref  = background[ip, 3]   # this will be -12 m/s from sounding
            pref   = background[ip, 5]

            θv_ref = θ_ref * (1.0 + 0.61 * qv_ref)
            Tref   = θ_ref / (PhysConst.pref / pref)^(PhysConst.Rair / PhysConst.cp)
            Tv_ref = Tref * (1.0 + 0.61 * qv_ref)
            ρref   = perfectGasLaw_TPtoρ(PhysConst; Temp=Tv_ref, Press=pref)
            hl_ref = PhysConst.cp * Tv_ref + PhysConst.g * y
            pref_m = ρref * Tv_ref * PhysConst.Rair

            # No bubble perturbation — pure background state
            ρ  = ρref
            hl = hl_ref
            u  = u_ref
            v  = 0.0

            if inputs[:SOL_VARS_TYPE] == PERT()
                q.qn[ip, 1] = ρ - ρref
                q.qn[ip, 2] = ρ * u - ρref * u_ref
                q.qn[ip, 3] = ρ * v - ρref * 0.0
                q.qn[ip, 4] = ρ * hl - ρref * hl_ref
                q.qn[ip, 5] = ρ * qv_ref - ρref * qv_ref
                q.qn[ip, 6] = 0.0
                q.qn[ip, end] = pref_m
            else
                q.qn[ip, 1] = ρ
                q.qn[ip, 2] = ρ * u
                q.qn[ip, 3] = ρ * v
                q.qn[ip, 4] = ρ * hl
                q.qn[ip, 5] = ρ * qv_ref
                q.qn[ip, 6] = 0.0
                q.qn[ip, end] = pref_m
            end

            q.qe[ip, 1] = ρref
            q.qe[ip, 2] = ρref * u_ref
            q.qe[ip, 3] = ρref * 0.0
            q.qe[ip, 4] = ρref * hl_ref
            q.qe[ip, 5] = ρref * qv_ref
            q.qe[ip, 6] = 0.0
            q.qe[ip, end] = pref_m
        end
    end

    @info " Initialize fields for 2D CompEuler with θ equation ........................ DONE"

    return q
end