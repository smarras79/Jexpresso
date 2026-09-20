#---------------------------------------------------------------------------------
# Two-dimensional Orszag-Tang vortex for the ideal GLM-MHD equations.
#
# A. Bormanis, C. A. Leon, A. Scheinker,
# "Solving the Orszag-Tang vortex magnetohydrodynamics problem with
#  physics-constrained convolutional neural networks",
# Phys. Plasmas 31, 012101 (2024), Section I A, Eqs. (6)-(9).
#
#   domain [0,1]², doubly periodic, γ = 5/3, run to T_final = 1, with
#
#     ρ(r,0) = γ²/(4π)                                   (Eq. 6)
#     p(r,0) = γ/(4π)                                    (Eq. 7)
#     v(r,0) = (-sin(2πy), sin(2πx))                     (Eq. 8)
#     B(r,0) = 1/(2π√(4π)) ∇⊥(cos(4πx)/2 + cos(2πy))     (Eq. 9)
#
#   with ∇⊥ = (∂/∂y, -∂/∂x). Carrying out the curl of Eq. (9) gives the
#   familiar closed form of the Orszag-Tang magnetic field,
#
#     Bx = -B₀ sin(2πy),   By = B₀ sin(4πx),   B₀ = 1/√(4π),
#
#   i.e. B = ∇×(A_z ẑ) with the vector potential
#   A_z = B₀ (cos(4πx)/(4π) + cos(2πy)/(2π)), so the discrete initial
#   field is divergence free to the accuracy of the SEM derivative.
#
#   Numerically: ρ = 25/(36π) ≈ 0.2210485, p = 5/(12π) ≈ 0.1326291,
#   B₀ ≈ 0.2820948, so that the initial sound speed is exactly
#   a = sqrt(γp/ρ) = 1, max|v| = √2 and the minimum plasma beta is
#   β = p/(½|B|²_max) = 5/3.
#
# The out-of-plane components w and Bz are identically zero for this
# problem (as is ψ), but they are carried along because the equation set
# implemented here is the full nine-field GLM-MHD system.
#
# The divergence-cleaning speed c_h is set to the maximum wave speed of
# the initial condition and kept constant throughout the simulation (the
# standard GLM choice, after Derigs et al. / Rueda-Ramírez et al.). For
# this initial condition it evaluates to c_h ≈ 2.603.
#---------------------------------------------------------------------------------
function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0
        @info " Initialize fields for 2D ideal GLM-MHD (Orszag-Tang vortex) ........... "
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
        error(" problems/MHD/orszagTangBormanis2024: only the CPU backend is supported for now.")
    end
    if (inputs[:SOL_VARS_TYPE] != TOTAL())
        error(" problems/MHD/orszagTangBormanis2024: only SOL_VARS_TYPE = TOTAL() is supported.")
    end

    γ   = γ_mhd
    γm1 = γ - 1.0

    #
    # Paper Eqs. (6)-(7): ρ and p are uniform at t = 0.
    #
    ρ  = γ*γ/(4.0*π)     # = 25/(36π) ≈ 0.2210485
    p  = γ/(4.0*π)       # =  5/(12π) ≈ 0.1326291
    B₀ = 1.0/sqrt(4.0*π) #            ≈ 0.2820948

    ψ  = 0.0
    w  = 0.0
    Bz = 0.0

    ch_local = 0.0

    for ip = 1:mesh.npoin

        x, y = mesh.x[ip], mesh.y[ip]

        # Eq. (8)
        u  = -sinpi(2.0*y)
        v  =  sinpi(2.0*x)

        # Eq. (9), evaluated in closed form (see the header)
        Bx = -B₀*sinpi(2.0*y)
        By =  B₀*sinpi(4.0*x)

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

    #
    # The initial condition above is written for the paper's unit square:
    # every field has period 1 in both x and y, so a mesh spanning anything
    # other than [0,1]² either tiles the vortex or truncates it (and then
    # the periodic BCs are inconsistent with the data). Warn loudly rather
    # than silently solving a different problem.
    #
    xmin, xmax = mesh.xmin, mesh.xmax   # global (MPI-reduced) extents, set in sem_setup
    ymin, ymax = mesh.ymin, mesh.ymax

    if rank == 0
        if (abs(xmin) > 1.0e-8 || abs(xmax - 1.0) > 1.0e-8 ||
            abs(ymin) > 1.0e-8 || abs(ymax - 1.0) > 1.0e-8)
            @warn string(" problems/MHD/orszagTangBormanis2024: the Orszag-Tang initial condition ",
                         "is defined on the unit square [0,1]², but this mesh spans ",
                         "[$(xmin), $(xmax)] × [$(ymin), $(ymax)]. ",
                         "Point :gmsh_filename at the mesh that ships with this case, ",
                         "problems/MHD/orszagTangBormanis2024/OT_32x32_periodic.msh, ",
                         "or regenerate it from the companion .geo file.")
        end
        @info " GLM divergence-cleaning speed c_h = $(c_h_mhd[])"
        @info " Initialize fields for 2D ideal GLM-MHD (Orszag-Tang vortex) ........... DONE"
    end

    return q
end
