#---------------------------------------------------------------------------------
# MHD blast wave (Balsara & Spicer 1999b; Balsara 2004 §7.3), centered in the unit square (0, 1)²:
# ρ = 1, v = 0, p = 1000 for r < 0.1 and 0.1 outside, B = (100/√(4π), 0, 0), ψ = 0, γ = 1.4.
# Heaviside-Lorentz B (magnetic pressure ½|B|²): ambient β = 2p/|B|² = 2.51e-4.
#---------------------------------------------------------------------------------
function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank == 0 && @info " Initialize fields for 2D ideal GLM-MHD (Balsara-Spicer blast) ........... "

    qvars    = ["ρ", "ρu", "ρv", "ρE", "ρw", "Bx", "By", "Bz", "ψ"]
    qoutvars = ["ρ", "u", "v", "w", "p", "Bx", "By", "Bz", "ψ", "T", "pmag", "Mach", "β"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)

    inputs[:backend] != CPU()        && error(" problems/MHD/blastBalsaraSpicer1999: only the CPU backend is supported.")
    inputs[:SOL_VARS_TYPE] != TOTAL() && error(" problems/MHD/blastBalsaraSpicer1999: only SOL_VARS_TYPE = TOTAL() is supported.")

    γm1      = γ_mhd - 1.0
    xc, yc   = 0.5, 0.5
    r0       = 0.1
    ρ0       = 1.0
    pin, pout = 1000.0, 0.1
    Bx, By   = 100.0/sqrt(4.0*π), 0.0
    ch_local = 0.0
    for ip = 1:mesh.npoin
        r = sqrt((mesh.coords[1,ip] - xc)^2 + (mesh.coords[2,ip] - yc)^2)
        p = r < r0 ? pin : pout

        q.qn[ip,1] = ρ0
        q.qn[ip,2] = 0.0
        q.qn[ip,3] = 0.0
        q.qn[ip,4] = p/γm1 + 0.5*(Bx*Bx + By*By)
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

        ch_local = max(ch_local, sqrt((γ_mhd*p + Bx*Bx + By*By)/ρ0))   # fast-speed bound, fluid at rest
    end
    c_h_mhd[] = MPI.Allreduce(ch_local, MPI.MAX, comm)

    if rank == 0
        (abs(mesh.xmin) > 1e-8 || abs(mesh.xmax - 1.0) > 1e-8 || abs(mesh.ymin) > 1e-8 || abs(mesh.ymax - 1.0) > 1e-8) &&
            @warn " problems/MHD/blastBalsaraSpicer1999: the blast is centered at (0.5, 0.5) of (0, 1)², but this mesh spans [$(mesh.xmin), $(mesh.xmax)] × [$(mesh.ymin), $(mesh.ymax)]."
        @info " GLM divergence-cleaning speed c_h = $(c_h_mhd[])"
        @info " Initialize fields for 2D ideal GLM-MHD (Balsara-Spicer blast) ........... DONE"
    end
    return q
end
