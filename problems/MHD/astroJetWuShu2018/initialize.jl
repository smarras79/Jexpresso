#---------------------------------------------------------------------------------
# Two-dimensional magnetized astrophysical jet for the ideal GLM-MHD equations.
#
#   K. Wu, C.-W. Shu, "A provably positive discontinuous Galerkin method for
#   multidimensional ideal magnetohydrodynamics", SIAM J. Sci. Comput. 40(5)
#   (2018) B1302-B1329, Example 5.6 ("Astrophysical jets"),
#
# which magnetizes the Mach 800 gas dynamical jet of
#
#   D. S. Balsara, "Self-adjusting, positivity preserving high order schemes
#   for hydrodynamics and magnetohydrodynamics", JCP 231 (2012) 7504-7517.
#
# INITIAL CONDITION — the whole domain [-0.5,0.5] x [0,1.5] is a uniform,
# static, magnetized ambient medium:
#
#     ρ = 0.1γ = 0.14        (γ = 1.4)
#     p = 1
#     v = (0, 0, 0)
#     B = (0, B_a, 0)        B_a = √200 by default (plasma β_a = 2p/B_a² = 10⁻²)
#
# There is nothing else at t = 0: everything that happens is driven through the
# nozzle {y = 0, |x| ≤ 0.05} by the boundary condition of user_bc.jl, which
# injects
#
#     ρ = γ = 1.4,  p = 1,  v = (0, 800, 0),  B = (0, B_a, 0)
#
# i.e. a beam ten times denser than the ambient gas at exactly Mach 800 (the
# beam sound speed is sqrt(γ p/ρ) = sqrt(γ·1/γ) = 1).
#
# The out-of-plane components w and Bz stay identically zero for this problem
# (nothing in the 2D fluxes couples them once they vanish), as does ψ up to the
# divergence error the discretization generates; they are carried because the
# implemented system is the full nine-field GLM-MHD one.
#
# THE DIVERGENCE-CLEANING SPEED. c_h is the standard GLM choice — the maximum
# wave speed, held constant — but it must be taken over the initial condition
# AND the prescribed inflow state. The domain is at rest at t = 0, so the
# initial condition alone gives c_h = c_f,ambient ≈ 37.9 (B_a = √200) while the
# beam entering from the very first step carries 800 + c_f,beam ≈ 812. Using
# the smaller value would make the ψ waves slower than the flow that generates
# the divergence error, i.e. cleaning that cannot keep up.
#---------------------------------------------------------------------------------
function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs, OUTPUT_DIR::String, TFloat)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    if rank == 0
        @info " Initialize fields for 2D ideal GLM-MHD (magnetized astrophysical jet) ........... "
    end

    #---------------------------------------------------------------------------------
    # Solution variables.
    #
    # NOTICE: the length of qvars defines neqs. Slot 4 MUST carry the total
    # energy ρE (not ρw) — see the header of user_flux.jl.
    #---------------------------------------------------------------------------------
    qvars    = ["ρ", "ρu", "ρv", "ρE", "ρw", "Bx", "By", "Bz", "ψ"]
    qoutvars = ["ρ", "u", "v", "w", "p", "Bx", "By", "Bz", "ψ", "T",
                "log10rho", "log10p", "beta", "Mach"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars), qoutvars=qoutvars)
    #---------------------------------------------------------------------------------

    if (inputs[:backend] != CPU())
        error(" problems/MHD/astroJetWuShu2018: only the CPU backend is supported for now.")
    end
    if (inputs[:SOL_VARS_TYPE] != TOTAL())
        error(" problems/MHD/astroJetWuShu2018: only SOL_VARS_TYPE = TOTAL() is supported.")
    end

    amb = aj_ambient_state()     # (ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ)
    jet = aj_jet_state()

    for ip = 1:mesh.npoin
        for ieq = 1:length(qvars)
            q.qn[ip,ieq] = amb[ieq]
            # The reference/background state: the ambient medium. It is used
            # for output of perturbations and as the SCALE of the
            # "did the boundary routine prescribe this component?" test in
            # kernel/boundaryconditions/BCs.jl. It is NOT subtracted from the
            # DynSGS residual (:dsgs_reference is left at its default false)
            # and it is NOT a weight in the viscous operator
            # (:dsgs_ref_weight is false): a uniform reference makes both a
            # no-op, since ∇q_e ≡ 0.
            q.qe[ip,ieq] = amb[ieq]
        end
        q.qn[ip,end] = AJ_P_AMB
        q.qe[ip,end] = AJ_P_AMB
    end

    #
    # GLM divergence-cleaning speed: the larger of the ambient and the injected
    # beam wave speed. Both states are uniform, so no reduction over the mesh
    # is needed and every rank computes the same number — but the Allreduce is
    # kept so that the value is provably identical across ranks even if a
    # future variant makes the initial condition non-uniform.
    #
    ch_local  = max(aj_wave_speed(amb), aj_wave_speed(jet))
    c_h_mhd[] = MPI.Allreduce(ch_local, MPI.MAX, comm)

    #
    # The initial and boundary conditions above are written for the paper's
    # domain [-0.5,0.5] x [0,1.5]: the nozzle half-width 0.05 and the outflow
    # boundaries are absolute positions, not fractions, so a mesh spanning
    # anything else solves a different problem. Warn loudly rather than
    # silently doing so.
    #
    xmin, xmax = mesh.xmin, mesh.xmax   # global (MPI-reduced) extents, set in sem_setup
    ymin, ymax = mesh.ymin, mesh.ymax

    if rank == 0
        if (abs(xmin + 0.5) > 1.0e-8 || abs(xmax - 0.5) > 1.0e-8 ||
            abs(ymin)       > 1.0e-8 || abs(ymax - 1.5) > 1.0e-8)
            @warn string(" problems/MHD/astroJetWuShu2018: this case is defined on ",
                         "[-0.5, 0.5] x [0, 1.5], but the mesh spans ",
                         "[$(xmin), $(xmax)] x [$(ymin), $(ymax)]. ",
                         "Point :gmsh_filename at one of the meshes that ship with the case ",
                         "(AJ_40x60.msh, AJ_100x150.msh) or regenerate one from AJ.geo.")
        end

        Ba   = aj_Ba[]
        uj   = aj_ujet[]
        βa   = 2.0*AJ_P_AMB/(Ba*Ba)
        c_j  = sqrt(γ_mhd*AJ_P_JET/AJ_RHO_JET)
        ρE_j = jet[4]
        @info @sprintf(" γ = %.4f   ambient (ρ, p) = (%.6g, %.6g)   c_amb = %.6g",
                       γ_mhd, AJ_RHO_AMB, AJ_P_AMB, sqrt(γ_mhd*AJ_P_AMB/AJ_RHO_AMB))
        @info @sprintf(" nozzle |x| <= %.4g : beam (ρ, p, v) = (%.6g, %.6g, %.6g), c_jet = %.6g, Mach = %.6g",
                       AJ_XNOZZLE, AJ_RHO_JET, AJ_P_JET, uj, c_j, uj/c_j)
        @info @sprintf(" B = (0, B_a, 0), B_a = %.8g (B_a² = %.6g), plasma β_a = %.3e", Ba, Ba*Ba, βa)
        @info @sprintf(" v_A: ambient %.6g, beam %.6g   |   c_f: ambient %.6g, beam %.6g",
                       Ba/sqrt(AJ_RHO_AMB), Ba/sqrt(AJ_RHO_JET),
                       aj_wave_speed(amb), aj_wave_speed(jet) - uj)
        @info @sprintf(" GLM divergence-cleaning speed c_h = %.6g (max of the ambient and the beam wave speed)", c_h_mhd[])
        # The number that says how hard this test is: the fraction of the
        # beam's total energy that IS the pressure. p = (γ-1)(ρE - KE - ME),
        # so a relative error of this size in ρE zeroes the pressure.
        @info @sprintf(" beam ρE = %.8g, of which p/(γ-1) = %.4g (%.3e of the total): a relative error of that order in ρE gives p < 0",
                       ρE_j, AJ_P_JET/(γ_mhd - 1.0), (AJ_P_JET/(γ_mhd - 1.0))/ρE_j)
        if aj_smooth[] > 0.0
            @info @sprintf(" nozzle lip SMOOTHED over w = %.4g (JEXPRESSO_AJ_SMOOTH); the paper's condition is the sharp top hat w = 0", aj_smooth[])
        end
        @info " Initialize fields for 2D ideal GLM-MHD (magnetized astrophysical jet) ........... DONE"
    end

    return q
end
