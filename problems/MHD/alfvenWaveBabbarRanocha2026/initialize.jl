#---------------------------------------------------------------------------------
# The smooth ALFVÉN WAVE: the GLM-MHD convergence test of
# Babbar & Ranocha, J. Comput. Phys. 563 (2026) 115060, §5.5.1, after
# Tóth, J. Comput. Phys. 161 (2000) 605. A circularly polarized Alfvén wave
# travels at the angle α to the x axis through a doubly periodic box. It is an
# exact nonlinear solution of ideal MHD: |B|² = 1.01 and p = 0.1 are uniform, so
# the total pressure is too, and the error can be measured against it at any
# time.
#
#   x⊥ = x cos α + y sin α,     B⊥ = 0.1 sin(2π(x⊥ + t)),     α = π/4
#
#   ρ  = 1,    p = 0.1,    γ = 5/3
#   v  = (−B⊥ sin α,  B⊥ cos α,  0.1 cos(2π(x⊥ + t)))
#   B  = (cos α + v₁,  sin α + v₂,  v₃)      (B∥ = 1 along (cos α, sin α))
#   ψ  = 0
#
# on [0, 1/cos α] × [0, 1/sin α] = [0, √2]², periodic, to t = 2.
#
# This is the solution exactly as the paper's runs compute it: the
# reproducibility repository (github.com/Arpit-Babbar/paper_crk_nonconservative,
# run/convergence_alfven_wave) takes it from Tenkai.jl's ExactSolutionAlfvenWave,
# which is Trixi.jl's initial_condition_convergence_test for
# IdealGlmMhdEquations2D. Note that the third velocity component is
# 0.1 cos(2π(x⊥ + t)) itself, not B⊥ times it.
#
# The wave moves along −(cos α, sin α) at the Alfvén speed B∥/√ρ = 1, and the
# box holds one wavelength along each side, so the solution is back to the
# initial condition at every integer time: t = 2 is two periods.
#---------------------------------------------------------------------------------
const AW_ALPHA = π/4     # propagation angle
const AW_AMP   = 0.1     # wave amplitude, in v and in B
const AW_RHO0  = 1.0
const AW_P0    = 0.1

# The exact state at (x, y, t), as the conserved vector of this case
# (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ) followed by the pressure.
@inline function aw_state(x, y, t, γ)
    sα, cα = sincos(AW_ALPHA)
    xp = x*cα + y*sα
    Bp = AW_AMP*sinpi(2.0*(xp + t))

    ρ  = AW_RHO0
    u  = -Bp*sα
    v  =  Bp*cα
    w  = AW_AMP*cospi(2.0*(xp + t))
    Bx = cα + u
    By = sα + v
    Bz = w
    ψ  = 0.0
    p  = AW_P0

    ρE = p/(γ - 1.0) + 0.5*ρ*(u*u + v*v + w*w) + 0.5*(Bx*Bx + By*By + Bz*Bz) + 0.5*ψ*ψ
    return (ρ, ρ*u, ρ*v, ρE, ρ*w, Bx, By, Bz, ψ, p)
end

function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank == 0 && @info " Initialize fields for 2D ideal GLM-MHD (smooth Alfvén wave) ........... "

    #---------------------------------------------------------------------------------
    # Slot 4 MUST carry the total energy ρE — see the header of user_flux.jl.
    #---------------------------------------------------------------------------------
    qvars    = ["ρ", "ρu", "ρv", "ρE", "ρw", "Bx", "By", "Bz", "ψ"]
    qoutvars = ["ρ", "u", "v", "w", "p", "Bx", "By", "Bz", "ψ", "T"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend];
                 neqs = length(qvars), qoutvars = qoutvars)

    inputs[:backend] == CPU() ||
        error(" problems/MHD/alfvenWaveBabbarRanocha2026: only the CPU backend is supported for now.")
    inputs[:SOL_VARS_TYPE] == TOTAL() ||
        error(" problems/MHD/alfvenWaveBabbarRanocha2026: only SOL_VARS_TYPE = TOTAL() is supported.")

    γ = γ_mhd

    ch_local = 0.0
    for ip = 1:mesh.npoin
        s = aw_state(mesh.coords[1,ip], mesh.coords[2,ip], 0.0, γ)
        for ieq = 1:length(qvars)
            q.qn[ip,ieq] = s[ieq]
            q.qe[ip,ieq] = s[ieq]      # background state, for perturbation output only
        end
        q.qn[ip,end] = s[10]
        q.qe[ip,end] = s[10]

        ρ  = s[1]
        u  = s[2]/ρ; v = s[3]/ρ; w = s[5]/ρ
        B2 = s[6]^2 + s[7]^2 + s[8]^2
        cf = sqrt(γ*s[10]/ρ + B2/ρ)
        ch_local = max(ch_local, sqrt(u*u + v*v + w*w) + cf)
    end

    # GLM divergence-cleaning speed: the maximum wave speed of the initial
    # condition over the whole (global) domain, constant in time, as in
    # smoothVortex: 0.1 + sqrt(γ p + |B|²) = 1.185 here. ψ = 0 and ∇·B = 0 in
    # the exact solution, so c_h only shapes the discretization error.
    c_h_mhd[] = MPI.Allreduce(ch_local, MPI.MAX, comm)

    if rank == 0
        # The wave is periodic on [0, 1/cos α] × [0, 1/sin α] only; on any other
        # box the initial condition jumps across the periodic seam and no
        # scheme converges below that jump.
        Lx = mesh.xmax - mesh.xmin
        Ly = mesh.ymax - mesh.ymin
        sα, cα = sincos(AW_ALPHA)
        if abs(Lx*cα - round(Lx*cα)) > 1.0e-8 || abs(Ly*sα - round(Ly*sα)) > 1.0e-8 ||
           round(Lx*cα) < 1 || round(Ly*sα) < 1
            @warn string(" problems/MHD/alfvenWaveBabbarRanocha2026: the box [", mesh.xmin, ", ", mesh.xmax,
                         "] x [", mesh.ymin, ", ", mesh.ymax, "] is not a whole number of wavelengths ",
                         "along x and y (1/cos α = ", 1/cα, ", 1/sin α = ", 1/sα, "): the wave is not periodic on it.")
        end
        @info " GLM divergence-cleaning speed c_h = $(c_h_mhd[])"
        @info " Initialize fields for 2D ideal GLM-MHD (smooth Alfvén wave) ........... DONE"
    end

    return q
end
