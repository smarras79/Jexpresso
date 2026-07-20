#---------------------------------------------------------------------------------
# Two-dimensional magnetized Kelvin-Helmholtz instability for the ideal
# GLM-MHD equations.
#
# Chan et al., Frontiers in Physics 10:898028 (2022), Section 3.2.1, Eq. (10):
#
#   domain [-1, 1]², doubly periodic, run to T_final = 15, with
#
#     B̃(x,y) = tanh(15y + 7.5) - tanh(15y - 7.5)     (smoothed step, Eq. (5))
#
#     ρ  = 1/2 + 3/4·B̃        p  = 1                 ψ  = 0
#     u  = (B̃ - 1)/2          v  = sin(2πx)/10       w  = 0
#     Bx = 0                  By = 0.125              Bz = 0
#
# The divergence-cleaning speed c_h is set to the maximum wave speed of the
# initial condition and kept constant throughout the simulation (the paper's
# standard choice, after [Derigs et al. / Rueda-Ramírez et al.]).
#---------------------------------------------------------------------------------
function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0
        @info " Initialize fields for 2D ideal GLM-MHD (magnetized KHI) ........... "
    end

    #---------------------------------------------------------------------------------
    # Solution variables:
    #
    # NOTICE: the length of qvars defines neqs. Slot 4 MUST carry the total
    # energy ρE (not ρw) because Jexpresso's shared 2D kernels assume the
    # energy lives in slot 4 — see the header of user_flux.jl.
    #---------------------------------------------------------------------------------
    qvars    = ["ρ", "ρu", "ρv", "ρE", "ρw", "Bx", "By", "Bz", "ψ"]
    qoutvars = ["ρ", "u", "v", "w", "p", "Bx", "By", "Bz", "ψ", "T"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)
    #---------------------------------------------------------------------------------

    if (inputs[:backend] != CPU())
        error(" problems/MHD/kelvinHelmholtzChan2022: only the CPU backend is supported for now.")
    end
    if (inputs[:SOL_VARS_TYPE] != TOTAL())
        error(" problems/MHD/kelvinHelmholtzChan2022: only SOL_VARS_TYPE = TOTAL() is supported.")
    end

    γ   = γ_mhd
    γm1 = γ - 1.0

    ch_local = 0.0

    for ip = 1:mesh.npoin

        x, y = mesh.x[ip], mesh.y[ip]

        B̃ = tanh(15.0*y + 7.5) - tanh(15.0*y - 7.5)

        ρ  = 0.5 + 0.75*B̃
        p  = 1.0
        u  = 0.5*(B̃ - 1.0)
        v  = sinpi(2.0*x)/10.0
        w  = 0.0
        Bx = 0.0
        By = 0.125
        Bz = 0.0
        ψ  = 0.0

        ρE = p/γm1 + 0.5*ρ*(u*u + v*v + w*w) + 0.5*(Bx*Bx + By*By + Bz*Bz) + 0.5*ψ*ψ

        q.qn[ip,1] = ρ
        q.qn[ip,2] = ρ*u
        q.qn[ip,3] = ρ*v
        q.qn[ip,4] = ρE
        q.qn[ip,5] = ρ*w
        q.qn[ip,6] = Bx
        q.qn[ip,7] = By
        q.qn[ip,8] = Bz
        q.qn[ip,9] = ψ
        q.qn[ip,end] = p

        # Store the initial state as the background state (used only for
        # perturbation output/analysis).
        for ieq = 1:length(qvars)
            q.qe[ip,ieq] = q.qn[ip,ieq]
        end
        q.qe[ip,end] = p

        # Local maximum wave speed |v| + c_f, with the fast magnetosonic
        # speed maximized over propagation directions: c_f ≤ sqrt(a² + b²),
        # a² = γp/ρ (sound), b² = |B|²/ρ (Alfvén).
        vmag = sqrt(u*u + v*v + w*w)
        cf   = sqrt(γ*p/ρ + (Bx*Bx + By*By + Bz*Bz)/ρ)
        ch_local = max(ch_local, vmag + cf)
    end

    #
    # GLM divergence-cleaning speed: max wave speed of the IC over the whole
    # (global) domain, constant in time.
    #
    c_h_mhd[] = MPI.Allreduce(ch_local, MPI.MAX, comm)

    if rank == 0
        @info " GLM divergence-cleaning speed c_h = $(c_h_mhd[])"
        @info " Initialize fields for 2D ideal GLM-MHD (magnetized KHI) ........... DONE"
    end

    return q
end
