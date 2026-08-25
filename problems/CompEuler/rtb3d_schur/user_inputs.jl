function user_inputs()

    #---------------------------------------------------------------------------
    # CompEuler/rtb3d_schur -- a 3D rising thermal bubble in a 10 km cube, built
    # to MEASURE the scalar Schur stage solve (:imex_schur) against the
    # five-field one on a mesh where it should pay.
    #
    # WHY THIS MESH. The column preconditioner is exact for the vertical
    # acoustic operator and does nothing for the horizontal one, so what it
    # leaves for the Krylov iteration scales with the grid's acoustic
    # anisotropy h_x/h_z. 20 x 20 x 80 over 10 km cubed is h_x = h_y = 500 m
    # against h_z = 125 m, i.e. 4:1, and at nop = 4 the smallest LGL gaps are
    # 86.4 m and 21.6 m.
    #
    # FIRST REAL MEASUREMENT, and it does not match the projection this case was
    # built on. On the 10 x 10 x 40 variant (same 4:1, 270,681 points), one rank,
    # Delta t = 0.6, CFL_h = 0.526, rtol 1e-8:
    #
    #                    s/step (steady)   cold iterations   total wall
    #     five-field          19.46              20            225.8 s
    #     scalar Schur        16.15              61            196.9 s
    #                        -> 1.21x
    #
    # The reduction IS faster, and for the opposite reason to the one predicted.
    # test/hevi/test_schur_precond.jl swept the mock at production stiffness and
    # measured H needing FEWER iterations than the five-field operator -- 2.67x
    # fewer at this very anisotropy. On a real mesh it needs 3x MORE (61 against
    # 20). The whole gain comes from per-iteration cost, where one implicit field
    # instead of five cuts the matvec, the gather/scatter and the
    # orthogonalisation, and it is large enough to absorb tripling the count.
    #
    # So the mock's ITERATION RATIO did not transfer -- it inverted. Two likely
    # reasons, neither yet tested: this run sits at 18% of the stability limit
    # where the five-field system is already easy at 20 iterations, while the
    # mock sweep sat at lam*rate = 10.7; and production preconditions [1, 4, 5]
    # rather than all five, because hevi_choose_vars drops rho_u and rho_v as
    # exactly the identity on a z-aligned mesh, so the five-field preconditioner
    # may capture its operator more completely than the scalar one captures H.
    #
    # TREAT THE 1.21x AS A LOWER BOUND ON ONE RANK AND NOTHING MORE. Two effects
    # pull opposite ways and neither is small: one rank has no MPI reduce at all,
    # which is the one part of the stage solve the reduction cannot shrink (6.8%
    # on the 64x64x60 profile), so its absence FLATTERS Schur; and CFL_h = 0.526
    # is half what the full mesh runs at, which keeps the stage solve a smaller
    # share of the step and so UNDERSTATES it. The number that settles it is the
    # full mesh on 25 ranks.
    #
    # HOW TO RUN THE A/B
    #
    #   DBG_SCHUR=0 mpiexec -n 25 julia --project=. src/Jexpresso.jl CompEuler rtb3d_schur
    #   DBG_SCHUR=1 mpiexec -n 25 julia --project=. src/Jexpresso.jl CompEuler rtb3d_schur
    #
    # 25 ranks is what tools/pick_nranks.jl recommends for this mesh (a 5 x 5
    # rank grid over the 20 x 20 element columns: 16 columns and 84k points
    # each). 16 ranks is marginally better per rank; 4, 5, 8, 10 and 20 also
    # work. DO NOT pick a count that tool does not list -- with :lxy_partition
    # the usable parallelism is bounded by nelemx*nelemy = 400 and the rest
    # leave ranks owning zero elements.
    #
    # MEMORY. Measured here on the full mesh: a ONE-RANK run reaches 13.4 GB
    # during IMEX3D setup and was OOM-killed on a 16 GB machine, having got as
    # far as the 3D operator. The mesh read and the high-order node population
    # alone take 253 s and several GB before the solver is even built. This case
    # is meant for >= 4 ranks; on 25 each rank holds 1/25 of the points. Use
    # DBG_MESH=10x10x40 for a single-rank smoke test.
    #
    # READ THE HEARTBEAT, NOT THE TOTAL. :lstep_heartbeat reports s/step
    # measured since the previous heartbeat, so after the first few steps it is
    # a steady-state rate with the JIT excluded, which total wall time is not.
    # DBG_TEND=60 is plenty for a timing run; the default 1000 s is there for
    # the physics.
    #
    # WHAT THE COMPARISON IS, AND IS NOT. :imex_schur forces the ADVECTIVE
    # Theta row -- the reduction does not close on one scalar with the flux
    # form -- so DBG_SCHUR=0 and DBG_SCHUR=1 are two different SPLITTINGS, not
    # two ways of solving the same one. They differ by 0.06% of the flux form
    # (test/hevi/test_theta_advective.jl). That is the right comparison for
    # wall clock and the wrong one for reading a state difference as an error.
    # What IS pinned: given the same operator, the two stage solves return the
    # same five fields to ~1e-12 (test/imex3d/test_schur_stage.jl).
    #---------------------------------------------------------------------------

    _vdiff = parse(Bool, get(ENV, "DBG_VDIFF", "false"))

    scheme = Symbol(get(ENV, "DBG_SCHEME", "imex"))
    scheme in (:imex, :hevi, :explicit) ||
        error("rtb3d_schur: DBG_SCHEME must be imex, hevi or explicit; got $scheme")

    #---------------------------------------------------------------------------
    # WHERE Δt = 0.6 COMES FROM.
    #
    # With every acoustic term implicit, Δt is set by the ADVECTIVE and viscous
    # rates on the smallest LGL gap, which here is 21.6 m in z:
    #
    #   advection at |v| = 20 m/s   0.93 1/s  x 1.44 (SEM correction) = 1.33
    #   SGS diffusion, mu = 125     0.27 1/s  x 1.44^2                = 0.56
    #   ---------------------------------------------------------------------
    #   explicit half                                                   1.9 1/s
    #
    # ARS343 is neutral on the reachable wedge out to zE = 2.83, so Δt_max is
    # about 1.5 s and 0.6 sits at ~40% of it. The margin is for |v|: it is the
    # one input to that estimate that is a forecast rather than a measurement,
    # and on a bubble at rest it is exactly zero at t = 0.
    #
    # The setup report recomputes all of this on the real mesh and prints the
    # limit, so a run that is over it is visible before it diverges rather than
    # after. Sweep with DBG_DT and a short DBG_TEND; divergence shows up as a
    # DomainError out of perfectGasLaw_rhothetatoP, i.e. pressure going
    # negative.
    #
    # FOR A SCHUR MEASUREMENT, Δt IS NOT NEUTRAL. The Krylov iteration count
    # grows linearly with CFL_h = γΔt·c/h_x, which is 1.05 here at Δt = 0.6. A
    # larger Δt puts more of the step inside the stage solve, which is exactly
    # what the Schur reduction addresses -- so it FLATTERS the comparison, and
    # a smaller one understates it. Comparing at one Δt, as here, is the honest
    # default; sweeping DBG_DT tells you how the saving moves with it.
    #---------------------------------------------------------------------------
    Δt_default = scheme === :imex ? "0.6" : (scheme === :hevi ? "0.2" : "0.1")
    Δt   = parse(Float64, get(ENV, "DBG_DT",   Δt_default))
    tend = parse(Float64, get(ENV, "DBG_TEND", "1000.0"))

    ode_solver = scheme === :imex     ? IMEX_ARK(:ARS343)      :
                 scheme === :hevi     ? HEVI_ARK(:ARS232)      :
                                        CarpenterKennedy2N54()

    inputs = Dict(
        :ode_solver           => ode_solver,
        :Δt                   => Δt,
        :tinit                => 0.0,
        :tend                 => tend,
        :diagnostics_at_times => (0:100:tend),
        :lstep_heartbeat      => true,
        :lsource              => true,
        :restart_time         => 1.0e7,
        :lrestart             => false,

        #-----------------------------------------------------------------------
        # THE SWITCH THIS CASE EXISTS FOR.
        #
        # :imex_schur eliminates rho, rho_u, rho_v and rho_w exactly and solves
        # one Helmholtz equation in P = beta*Theta -- Np unknowns instead of
        # 5*Np -- then rebuilds the five fields pointwise. The saving has two
        # parts: per iteration the matvec, the gather/scatter and the
        # orthogonalisation all scale with the implicit field count (73% of the
        # stage solve on the 64x64x60 profile, against 6.8% for the MPI reduce,
        # which does not scale), and per solve the iteration count falls by the
        # factor in the header.
        #
        # Measured 1.21x per step on the small variant (see the header), on the
        # per-iteration saving alone -- the iteration count went the WRONG way,
        # 61 against 20. Not yet measured on the full mesh or on many ranks.
        #-----------------------------------------------------------------------
        :imex_schur           => parse(Bool, get(ENV, "DBG_SCHUR", "false")),

        :implicit_vdiff       => _vdiff,
        :lcfl_report          => true,
        :lcfl_report_every    => parse(Int, get(ENV, "DBG_CFL_EVERY", "0")),
        :imex_verify          => parse(Bool, get(ENV, "DBG_VERIFY", "true")),
        :imex_warm_start      => parse(Bool, get(ENV, "DBG_WARM", "true")),
        :imex_rtol            => parse(Float64, get(ENV, "DBG_RTOL", "1.0e-8")),
        # 20 is right while CFL_h stays below ~3; it is 1.05 here. The setup
        # report prints CFL_h, the grid anisotropy behind it, and an advised m.
        :imex_restart         => parse(Int,     get(ENV, "DBG_RESTART", "20")),
        :imex_maxiter         => parse(Int,     get(ENV, "DBG_MAXITER", "200")),
        # With :imex_schur on, this selects the column preconditioner for the
        # SCALAR system rather than the five-field one; :none measures what it
        # buys, on either path.
        :imex_precond         => Symbol(get(ENV, "DBG_PRECOND", "column")),
        # The largest flow speed the run is expected to reach. It cannot be
        # measured at setup: this scheme's explicit half is dominated by
        # ADVECTION, which on a bubble at rest is exactly zero at t = 0, so an
        # estimate from the initial state would call any Δt safe. A 2 K / 2000 m
        # bubble reaches |w| of roughly 12 m/s; 20 is an upper bound.
        :imex_umax            => 20.0,
        # All six faces of this box are no-flux, and the implicit operator has
        # to know: it carries the horizontal mass flux, and letting sound run
        # out through the side walls would leave a stiff residual in the
        # explicit half exactly where rhs! projects the normal momentum out.
        :imex_lateral_walls   => Symbol(get(ENV, "DBG_LATWALL", "auto")),
        :imex_wall_flux       => true,
        # :RS is right here -- beta and thetabar barely move for a 2 K bubble on
        # a 300 K background, so :PS would pay for a refactorisation and recover
        # a coefficient that was never stale.
        :imex_linearization   => Symbol(get(ENV, "DBG_LIN", _vdiff ? "PS" : "RS")),
        :imex_update_freq     => parse(Int, get(ENV, "DBG_UPDFREQ", "5")),
        :imex_stability_guard => Symbol(get(ENV, "DBG_GUARD", "warn")),
        :imex_monitor         => parse(Bool, get(ENV, "DBG_IMEXMON", "false")),

        :hevi_verify          => parse(Bool, get(ENV, "DBG_VERIFY", "true")),
        :hevi_linearization   => Symbol(get(ENV, "DBG_LIN", _vdiff ? "PS" : "RS")),
        :hevi_update_freq     => parse(Int, get(ENV, "DBG_UPDFREQ", "5")),
        :hevi_wall_flux       => true,

        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,

        #---------------------------------------------------------------------------
        # Physical parameters -- CompEuler/rtb_imex's, unchanged.
        #---------------------------------------------------------------------------
        :lvisc                => parse(Bool, get(ENV, "DBG_VISC", "true")),
        :visc_model           => AV(),
        :μ                    => [0.0, 125.0, 125.0, 125.0, 125.0],
        :energy_equation      => "theta",

        #---------------------------------------------------------------------------
        # Mesh. Regenerate with either of
        #
        #   gmsh -3 problems/CompEuler/rtb3d_schur/rtb3d_20x20x80.geo
        #   julia --project=. problems/CompEuler/rtb3d_schur/generate_mesh.jl
        #
        # The second needs no gmsh on PATH -- it drives the SDK GridapGmsh
        # already ships.
        #
        # :lxy_partition is stated rather than left to its default because this
        # case is meant to be run on many ranks and it is the flag that makes
        # the decomposition COLUMNAR. On the p4est space-filling-curve
        # partition the assembled RHS and the mass matrix pick up different
        # ghost multiplicities at some rank-shared nodes, the implicit operator
        # stops being skew, and the run grows at a few 1/s with every
        # self-check still passing. IMEX3D refuses that combination outright.
        #---------------------------------------------------------------------------
        # DBG_MESH picks the size. 20x20x80 is the case; 10x10x40 is the SAME
        # 4:1 aspect ratio at 270,681 points instead of 2,106,081, for a quick
        # check that a build runs before spending cluster time on it -- and it
        # fits on one rank, which the full mesh does not (see the memory note
        # in the header).
        :lread_gmsh           => true,
        :gmsh_filename        => string("./problems/CompEuler/rtb3d_schur/rtb3d_",
                                        get(ENV, "DBG_MESH", "20x20x80"), ".msh"),
        :lxy_partition        => true,
        :lwarp                => false,
        :lstretch             => false,

        #---------------------------------------------------------------------------
        # Output -- one tree per stage solve, so an A/B leaves both solutions on
        # disk rather than the second overwriting the first.
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :output_dir           => scheme !== :imex ?
                                   (scheme === :hevi ? "./output_hevi/" : "./output_explicit/") :
                                 (parse(Bool, get(ENV, "DBG_SCHUR", "false")) ?
                                   "./output_imex_schur/" : "./output_imex_full/"),
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        :loutput_pert         => true,

        #---------------------------------------------------------------------------
        # AMR: off. IMEX3D refuses :ladapt -- adaptation invalidates the column
        # topology the preconditioner is built on, its factorised matrices and
        # its gather/scatter plan. :linitial_refine must be false as well:
        # mesh.jl consults :lxy_partition ONLY when both are off.
        #---------------------------------------------------------------------------
        :linitial_refine      => false,
        :init_refine_lvl      => 1,
        :ladapt               => false,
    ) #Dict
    return inputs
end
