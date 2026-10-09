#---------------------------------------------------------------------------------
# MHD rotor of Balsara & Spicer (1999, JCP 149:270) as set up by Dao & Nazarov (2022, JSC 92:77, §5.5):
# ambient (ρ, u, p, B) = (1, 0, 1, (5/√(4π), 0)); disc r < r₀: ρ = 10, u = (u₀/r₀)(0.5 − y, x − 0.5);
# taper r₀ ≤ r < r₁: ρ = 1 + 9f, u = (f u₀/r)(0.5 − y, x − 0.5), f = (r₁ − r)/(r₁ − r₀);
# r₀ = 0.1, r₁ = 0.115, u₀ = 2 (Tóth 2000's "first rotor"; not printed by Dao & Nazarov), γ = 1.4.
#---------------------------------------------------------------------------------
function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank == 0 && @info " Initialize fields for 2D ideal GLM-MHD (rotor) ........... "

    qvars    = ["ρ", "ρu", "ρv", "ρE", "ρw", "Bx", "By", "Bz", "ψ"]
    qoutvars = ["ρ", "u", "v", "w", "p", "Bx", "By", "Bz", "ψ", "T", "pmag", "Mach"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)

    inputs[:backend] != CPU()        && error(" problems/MHD/rotorDaoNazarov2022: only the CPU backend is supported.")
    inputs[:SOL_VARS_TYPE] != TOTAL() && error(" problems/MHD/rotorDaoNazarov2022: only SOL_VARS_TYPE = TOTAL() is supported.")

    γm1        = γ_mhd - 1.0
    r0, r1, u0 = 0.1, 0.115, 2.0
    p, Bx, By  = 1.0, 5.0/sqrt(4.0*π), 0.0
    ch_local   = 0.0
    for ip = 1:mesh.npoin
        dx, dy = mesh.coords[1,ip] - 0.5, mesh.coords[2,ip] - 0.5
        r = sqrt(dx*dx + dy*dy)
        if r < r0
            ρ, ω = 10.0, u0/r0
        elseif r < r1
            f    = (r1 - r)/(r1 - r0)
            ρ, ω = 1.0 + 9.0*f, f*u0/r
        else
            ρ, ω = 1.0, 0.0
        end
        u, v = -ω*dy, ω*dx

        q.qn[ip,1] = ρ
        q.qn[ip,2] = ρ*u
        q.qn[ip,3] = ρ*v
        q.qn[ip,4] = p/γm1 + 0.5*ρ*(u*u + v*v) + 0.5*(Bx*Bx + By*By)
        q.qn[ip,5] = 0.0
        q.qn[ip,6] = Bx
        q.qn[ip,7] = By
        q.qn[ip,8] = 0.0
        q.qn[ip,9] = 0.0
        q.qn[ip,end] = p
        for ieq = 1:length(qvars)
            q.qe[ip,ieq] = q.qn[ip,ieq]
        end
        q.qe[ip,end] = p

        ch_local = max(ch_local, sqrt(u*u + v*v) + sqrt((γ_mhd*p + Bx*Bx + By*By)/ρ))  # |v| + c_f bound
    end
    c_h_mhd[] = MPI.Allreduce(ch_local, MPI.MAX, comm)   # GLM cleaning speed, constant in time

    if rank == 0
        (abs(mesh.xmin) > 1e-8 || abs(mesh.xmax - 1.0) > 1e-8 || abs(mesh.ymin) > 1e-8 || abs(mesh.ymax - 1.0) > 1e-8) &&
            @warn " problems/MHD/rotorDaoNazarov2022: the rotor is centered in the unit square, but this mesh spans [$(mesh.xmin), $(mesh.xmax)] × [$(mesh.ymin), $(mesh.ymax)]."
        @info " GLM divergence-cleaning speed c_h = $(c_h_mhd[])"
        @info " Initialize fields for 2D ideal GLM-MHD (rotor) ........... DONE"
    end
    return q
end
