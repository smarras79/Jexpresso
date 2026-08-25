function user_inputs()

    #---------------------------------------------------------------------------
    # CompEuler/rtb2d_schur -- THE DESKTOP ANALOGUE of CompEuler/rtb3d_schur:
    # the same rising thermal bubble, the same mesh spacing, the same Delta t,
    # the same Schur A/B, on 1/16 of the points. Built to be re-runnable on a
    # laptop or a workstation instead of a queue slot.
    #
    # IT IS A SLAB, NOT AN nsd == 2 MESH, AND THAT IS NOT A SHORTCUT
    # --------------------------------------------------------------
    # HEVI and IMEX3D refuse a 2D mesh outright (params_setup.jl: "IMEX3D is 3D
    # only"). The refusal is protective. build_column_topology reads coords[3,:]
    # and builds its column catalogue from distinct (x, y) pairs, which a 2D
    # mesh -- coords shaped (2, npoin), vertical along y -- cannot supply; every
    # element kernel is a triple i,j,k loop over ngl^3 nodes of
    # connijk[iel,i,j,k], which in 2D is shaped (nelem, ngl, ngl, 1); and the
    # operator is contracted with dzeta/d{x,y,z}, which a 2D mesh ALLOCATES and
    # never fills. That last one is the reason the guard is an error rather than
    # a fallback: without it the acoustic operator would assemble to zero and
    # the run would be silently wrong rather than dead.
    #
    # So this is the idiom rtb_hevi and rtb_imex already use: a 3D mesh ONE
    # ELEMENT THICK in y with free slip front and back, and a CYLINDRICAL bubble
    # in initialize.jl. v is zero at t = 0, free slip holds it at zero on both y
    # faces, and nothing in a y-invariant problem generates it -- so the answer
    # is exactly the 2D one while the solver runs its ordinary 3D path. Nothing
    # about the Schur reduction, the preconditioner or the profile is special-
    # cased for this case.
    #
    # WHAT IS HELD FIXED FROM THE 3D CASE, so the two are comparable
    # --------------------------------------------------------------
    # Same 10 km x 10 km in x-z, same 20 x 80 elements there, so the SAME
    # spacing (h_x = 500 m, h_z = 125 m) and the SAME 4:1 acoustic anisotropy --
    # the ratio that decides how much the Schur reduction can save, because the
    # column preconditioner is exact vertically and does nothing horizontally.
    # Same nop = 4, same Delta t = 0.6, same CFL_h ~ 1.05, same 2 K bubble, same
    # viscosity. The only change is that 20 elements of y collapse to 1.
    #
    # y is one element 500 m wide, ISOTROPIC WITH h_x on purpose. A thinner slab
    # would raise h_x/h_y and hand the horizontal operator a stiffness the 3D
    # case does not have -- which is the one thing this benchmark exists to hold
    # fixed.
    #
    # SIZE. 81 x 5 x 321 = 130,005 gridpoints against the 3D case's 2,106,081.
    # DBG_MESH=10x1x40 gives 33,005 for a smoke test.
    #
    # HOW TO RUN THE A/B
    #
    #   ./problems/CompEuler/rtb2d_schur/run_rtb2d.sh full
    #   ./problems/CompEuler/rtb2d_schur/run_rtb2d.sh schur
    #   ./problems/CompEuler/rtb2d_schur/run_rtb2d.sh schur-ref
    #
    # or by hand, where NR is 4 (see the :lxy_partition note below -- the slab
    # has only 20 element columns, so more ranks than that is not an option and
    # even 10 is already communication-bound):
    #
    #   DBG_SCHUR=0 mpiexec -n 4 julia --project=. src/Jexpresso.jl CompEuler rtb2d_schur
    #   DBG_SCHUR=1 mpiexec -n 4 julia --project=. src/Jexpresso.jl CompEuler rtb2d_schur
    #
    # WHAT TO EXPECT, AND WHAT DOES NOT CARRY OVER FROM 3D
    # ----------------------------------------------------
    # The per-iteration savings should carry: they are a property of solving for
    # one field instead of five, and the anisotropy that sets them is identical.
    # Two things do NOT:
    #
    #   * the MPI reduce is the one part of the stage solve the reduction cannot
    #     shrink, and on 4 ranks over 20 columns it is a bigger share of a
    #     smaller problem than it was on 25 ranks over 400. That works AGAINST
    #     Schur here, so a smaller speedup than 3D's 2.82x is the expected
    #     result, not a regression.
    #   * absolute s/step means nothing across the two -- 130,005 points against
    #     2,106,081.
    #
    # Read the profile SPLIT, not just the step time. That is what transfers.
    #
    # READ THE HEARTBEAT, NOT THE TOTAL. :lstep_heartbeat reports s/step since
    # the previous heartbeat, so after a few steps it is a steady-state rate
    # with the JIT excluded, which total wall time is not. DBG_TEND=60 is plenty
    # for timing; the default 1000 s is there for the physics.
    #
    # WHAT THE COMPARISON IS, AND IS NOT. :imex_schur forces the ADVECTIVE Theta
    # row -- the reduction does not close on one scalar with the flux form -- so
    # DBG_SCHUR=0 and DBG_SCHUR=1 are two different SPLITTINGS, not two ways of
    # solving the same one. They differ by 0.06% of the flux form
    # (test/hevi/test_theta_advective.jl). That is the right comparison for wall
    # clock and the wrong one for reading a state difference between
    # output_imex_full/ and output_imex_schur/ as an error. What IS pinned:
    # given the same operator, the two stage solves return the same five fields
    # to ~1e-12 (test/imex3d/test_schur_stage.jl).
    #
    # PLOT dθ, NOT θ. θ is the total -- 300 K of background with a 2 K bubble on
    # it -- so it autoscales to the stratification and the bubble vanishes.
    #---------------------------------------------------------------------------

    _vdiff = parse(Bool, get(ENV, "DBG_VDIFF", "false"))

    scheme = Symbol(get(ENV, "DBG_SCHEME", "imex"))
    scheme in (:imex, :hevi, :explicit) ||
        error("rtb2d_schur: DBG_SCHEME must be imex, hevi or explicit; got $scheme")

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
        # Measured on the 3D case this is the analogue of, 25 ranks: with the
        # scalar matvec fixed the step went 9.581 -> 3.398 s, 2.82x. See the
        # README here for what to expect on this mesh and what does not carry
        # over from that one.
        #-----------------------------------------------------------------------
        :imex_schur           => parse(Bool, get(ENV, "DBG_SCHUR", "false")),

        # WHICH FORM OF H THE SCALAR MATVEC USES. Default on, and there is no
        # accuracy reason to turn it off: the two forms agree to 1.9e-16
        # (test/hevi/test_schur_kernel.jl).
        #
        # Off, `schur_H!` is two full five-field operator applications -- the
        # reference form, correct by construction because it reuses the verified
        # operator, and 2.05x the cost of ONE application. That is what the
        # first cluster profile of this path measured, and it is why the matvec
        # came out 46% slower than the five-field solve it replaced even though
        # the orthogonalisation fell 12.3x, the banded solve 6.9x and the MPI
        # reduce 8.9x. On, the same H costs 0.36x one application.
        #
        # DBG_SCHUR_KERN=0 with DBG_SCHUR=1 is therefore the third leg of the
        # A/B: it separates "the reduction is sound" from "the matvec is fast".
        :imex_schur_kernel    => parse(Bool, get(ENV, "DBG_SCHUR_KERN", "true")),

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
        #   gmsh -3 problems/CompEuler/rtb2d_schur/rtb2d_20x1x80.geo
        #   julia --project=. problems/CompEuler/rtb2d_schur/generate_mesh.jl
        #
        # The second needs no gmsh on PATH -- it drives the SDK GridapGmsh
        # already ships.
        #
        # :lxy_partition is stated rather than left to its default because it
        # is the flag that makes
        # the decomposition COLUMNAR. On the p4est space-filling-curve
        # partition the assembled RHS and the mass matrix pick up different
        # ghost multiplicities at some rank-shared nodes, the implicit operator
        # stops being skew, and the run grows at a few 1/s with every
        # self-check still passing. IMEX3D refuses that combination outright.
        #
        # ON A SLAB IT ALSO CAPS THE PARALLELISM, and low: z is never
        # partitioned and nelemy is 1, so there are 20 element columns and that
        # is the ceiling. tools/pick_nranks.jl lists 1, 2, 4, 5 and 10 as the
        # counts that leave no rank empty; 4 is the default here. Ten cores does
        # NOT mean ten ranks on this mesh -- at 10 each rank owns 2 columns and
        # 13k points and the halo per element rises to 3.0, so the run becomes
        # communication-bound and the profile stops being about the stage solve.
        # That ceiling is a property of a 2D problem, not of this deck.
        #---------------------------------------------------------------------------
        # DBG_MESH picks the size. 20x1x80 is the case, 130,005 points;
        # 10x1x40 is the SAME 4:1 aspect ratio at 33,005 points, for a smoke
        # test on a single core. Both are committed. KEEP nelemy = 1: more than
        # one element in y is a thin 3D box, not the 2D analogue, and since the
        # bubble is a cylinder nothing new happens physically for the cost.
        :lread_gmsh           => true,
        :gmsh_filename        => string("./problems/CompEuler/rtb2d_schur/rtb2d_",
                                        get(ENV, "DBG_MESH", "20x1x80"), ".msh"),
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
        # DEAD KEY, kept for parity with the other decks. mod_inputs.jl gives
        # :loutput_pert a default and nothing in src/ ever reads it -- what
        # actually decides whether the perturbation reaches the VTU is
        # user_uout! in user_primitives.jl, which writes it into slot 6, and
        # qoutvars in initialize.jl, which names that slot. Setting this true
        # does nothing either way.
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
