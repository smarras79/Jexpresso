function user_inputs()

    # Implicit vertical diffusion is opt-in per run: DBG_VDIFF=1. Read into a
    # local here rather than inline because two keys below have to agree about
    # it -- switching it on also switches the linearisation to :PS, without
    # which the implicit operator would carry no diffusion at all.
    _vdiff = parse(Bool, get(ENV, "DBG_VDIFF", "false"))

    #---------------------------------------------------------------------------
    # The IMEX analogue of CompEuler/rtb_hevi: same domain, same bubble, same
    # mesh, same viscosity, same output cadence, same tend. Only the time
    # integrator differs, which is the point -- a wall-clock comparison between
    # schemes means nothing unless everything else is bit-identical.
    #
    #   :imex      IMEX_ARK(:ARS343)       ALL acoustics implicit, 3D  -> ./output_imex/
    #   :hevi      HEVI_ARK(:ARS232)       vertical acoustics implicit -> ./output_hevi/
    #   :explicit  CarpenterKennedy2N54()  nothing implicit            -> ./output_explicit/
    #
    # ONE switch rather than three sets of commented-out lines, because swapping
    # comments by hand is how the mesh or the viscosity quietly ends up
    # different between the two runs being compared.
    #
    # EACH SCHEME RUNS AT ITS OWN LARGEST SAFE STEP. Timing them at a common Δt
    # would measure cost per step, which is not the question. Every implicit
    # scheme here costs MORE per step and pays for it by taking a bigger one:
    #
    #   explicit  5 rhs!, no solve
    #   hevi      3 rhs!, 3 operator applications, 2 banded column solves
    #   imex      4 rhs!, 4 operator applications, 3 preconditioned Krylov
    #             solves -- each of which is itself several operator
    #             applications plus several column solves
    #
    # so only wall-clock to a fixed physical time answers it.
    #---------------------------------------------------------------------------
    scheme = Symbol(get(ENV, "DBG_SCHEME", "imex"))
    scheme in (:imex, :hevi, :explicit) ||
        error("rtb_imex: DBG_SCHEME must be imex, hevi or explicit; got $scheme")

    #---------------------------------------------------------------------------
    # WHERE Δt = 1.0 COMES FROM, AND WHY IT IS NOT LARGER.
    #
    # With every acoustic term implicit, Δt is set by the ADVECTIVE and viscous
    # rates. On this mesh (Δx = 1000 m, Δy = 5000 m, Δz = 200 m elements at
    # nop = 4, so the smallest LGL gaps are 172.7 / 863.5 / 34.5 m):
    #
    #   advection at |v| = 20 m/s   0.72 1/s   x 1.44 (SEM correction) = 1.04
    #   SGS diffusion, mu = 125     0.22 1/s   x 1.44^2                = 0.46
    #   ----------------------------------------------------------------------
    #   explicit half                                                    1.5 1/s
    #   acoustics (implicit)       12.5 1/s   x 1.44                  = 18   1/s
    #
    # ARS343 is neutral on the reachable wedge out to zE = 2.83, which at
    # these rates is Δt = 1.89 s (ark_wedge_dt_max, measured). At 1.0 s the run
    # sits at 53% of that, which is the margin the vertical advective rate
    # needs: it is proportional to |v|, and |v| is the one quantity here that is
    # a forecast rather than a measurement.
    #
    # BUT STABILITY IS NOT THE ONLY REASON THIS IS 1.0 AND NOT 1.89. The Krylov
    # iteration count grows LINEARLY with Δt -- measured on this exact mesh:
    #
    #     Δt     γΔt     CFL_h   iterations/solve (rtol 1e-8)
    #     0.35   0.153   0.31     11
    #     0.50   0.218   0.44     13
    #     1.00   0.436   0.88     22
    #     2.00   0.872   1.75     42
    #     4.00   1.744   3.50     86
    #
    # where CFL_h = γΔt·c/h_x is the HORIZONTAL acoustic Courant number the
    # vertical preconditioner leaves behind. So cost per step is linear in Δt
    # and cost per unit SIMULATED time approaches a floor: past the point where
    # the four rhs! evaluations stop dominating the step, a bigger Δt buys
    # nothing at all. On this mesh that crossover is right about here.
    # See src/kernel/solvers/hevi/README_IMEX3D.md.
    #
    # The setup report recomputes all of this for the real mesh and prints the
    # limit, so a run that is over it is visible before it diverges rather than
    # after. Sweep with DBG_DT and a short DBG_TEND (100 s is enough) to find
    # the real limit; divergence shows up as a DomainError out of
    # perfectGasLaw_rhothetatoP, which is the pressure going negative.
    #---------------------------------------------------------------------------
    Δt_default = scheme === :imex ? "1.0" : (scheme === :hevi ? "0.35" : "0.2")
    Δt   = parse(Float64, get(ENV, "DBG_DT",   Δt_default))
    tend = parse(Float64, get(ENV, "DBG_TEND", "1000.0"))

    ode_solver = scheme === :imex     ? IMEX_ARK(:ARS343)      :
                 scheme === :hevi     ? HEVI_ARK(:ARS232)      :
                                        CarpenterKennedy2N54()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # Time integration. :lcfl_report below prints the stability table for
        # whichever scheme is running, including what each split would buy on
        # this mesh, so a run that is silently over its limit is visible first.
        #---------------------------------------------------------------------------
        :ode_solver           => ode_solver,
        :Δt                   => Δt,

        :tinit                => 0.0,
        :tend                 => tend,
        :diagnostics_at_times => (0:100:tend),
        # The heartbeat is the timing instrument: it reports s/step measured
        # over the interval since the PREVIOUS heartbeat, so after the first few
        # steps it is a steady-state rate with the JIT excluded, which total
        # wall time is not. The figure of merit is
        #
        #     wall to t = tend  =  (s/step) x tend / Δt
        #
        # It also removes the "has it hung?" question: the integrator warm-up
        # prints PATIENCE and then nothing when it finishes, so without a
        # heartbeat the run is silent for minutes while it is working normally.
        :lstep_heartbeat      => true,
        :lsource              => true,
        :restart_time         => 1.0e7,
        :lrestart             => false,

        #---------------------------------------------------------------------------
        # IMEX3D options. All optional; these are the defaults, spelled out so
        # the case documents them. See src/kernel/solvers/hevi/README.md.
        #---------------------------------------------------------------------------
        # :implicit_vdiff  MAKE THE VERTICAL SGS DIFFUSION IMPLICIT AS WELL.
        #
        # HEVI and the 3D IMEX remove the acoustic limit. What is left is
        # whatever was second, and on a mesh refined in z that is the SGS
        # diffusion of the vertical derivative, nu/dz^2. It is invisible at
        # startup -- the CFL report is taken on a laminar sounding where nu_t
        # is ~0 -- and it arrives on its own schedule as the boundary layer
        # spins up, which is what a blow-up at a fixed MODEL time looks like.
        #
        # This puts d/dz(mu d/dz) on u, v, w and theta into the same column
        # operator that already carries the vertical acoustics, so it costs a
        # wider band and no new solve. The horizontal and cross terms stay
        # explicit; they are dz/dx of what is removed.
        #
        # It REQUIRES :PS under a dynamic closure -- the coefficient is read
        # from the closure, which returns its molecular floor at t = 0, so :RS
        # would freeze it there and the operator would carry no diffusion at
        # all. DBG_VDIFF=1 therefore switches the linearisation with it (and
        # DBG_LIN still overrides).
        #---------------------------------------------------------------------------
        :implicit_vdiff       => _vdiff,
        #---------------------------------------------------------------------------
        :lcfl_report          => true,
        # :lcfl_report_every  re-print that table every N steps. :lcfl_report shows it
        #                   ONCE, at startup, where nu_t is ~0 on a laminar sounding --
        #                   so its viscous row says diffusion never limits dt, and says
        #                   it truthfully, at t = 0. Watch that row cross the line as the
        #                   boundary layer spins up. Collective and cheap; 0 = off.
        :lcfl_report_every    => parse(Int, get(ENV, "DBG_CFL_EVERY", "0")),
        :imex_verify          => parse(Bool, get(ENV, "DBG_VERIFY", "true")),
        #-----------------------------------------------------------------------
        # The Krylov solve. These four decide the cost of the whole scheme.
        #
        # :imex_rtol is the one to think about. It has to be tight enough that
        # the stage solve is not the leading error term -- the tableau is third
        # order, so a stage residual of 1e-8 relative is comfortably below the
        # O(Δt^4) local error and does not touch the observed order -- and loose
        # enough not to spend iterations on digits the step is going to throw
        # away. Loosening it to 1e-6 typically saves one or two iterations per
        # stage; loosening it to 1e-4 starts showing up as a loss of order.
        #
        # :imex_warm_start  START EACH STAGE SOLVE FROM THE PREVIOUS ONE'S ANSWER
        #                   rather than from zero. On by default, and the single
        #                   largest lever on this scheme's cost.
        #
        #                   Consecutive ARK stages differ by O(dt*f) while the
        #                   right-hand side is the WHOLE deviation u - qe, so the
        #                   guess arrives several orders down. This operator is
        #                   skew and GMRES on it converges LINEARLY -- iterations
        #                   go as log(residual/tol) -- so those orders come
        #                   straight off the count. Measured end to end through
        #                   the ARK stepper: 5.0 iterations/solve against 21.7
        #                   cold, with the two states agreeing to 1.3e-10.
        #
        #                   It cannot change what the solve returns: the tolerance
        #                   is unchanged and is measured against the true residual.
        #                   DBG_WARM=0 turns it off to measure the difference.
        :imex_warm_start      => parse(Bool, get(ENV, "DBG_WARM", "true")),
        # :imex_restart is NOT the knob it looks like. The usual advice -- raise
        # the restart length, because GMRES(m) stagnates on restart -- is for
        # elliptic operators. This one is skew: its spectrum is a line segment,
        # there is no clustering for a longer Krylov polynomial to exploit, and
        # restarting costs almost nothing. Measured: 71 iterations at m = 20
        # against 68 at m = 240. Leave it at 20 and spend the memory elsewhere.
        #-----------------------------------------------------------------------
        :imex_rtol            => parse(Float64, get(ENV, "DBG_RTOL", "1.0e-8")),
        :imex_restart         => parse(Int,     get(ENV, "DBG_RESTART", "20")),
        :imex_maxiter         => parse(Int,     get(ENV, "DBG_MAXITER", "200")),
        # :column is HEVI's own banded column solve over all five equations --
        # a line-relaxation smoother in the stiffest direction, already
        # factorised. :none is there to MEASURE what it buys: run both and
        # compare the iteration counts the setup report prints. Without it this
        # case takes tens of iterations per stage instead of a handful.
        :imex_precond         => Symbol(get(ENV, "DBG_PRECOND", "column")),
        #-----------------------------------------------------------------------
        # The largest flow speed the run is expected to reach, in m/s.
        #
        # This is the one input to the stability estimate that cannot be
        # measured at setup, and it matters more here than it ever did for
        # HEVI. HEVI's explicit half is dominated by horizontal SOUND, which is
        # a property of the mesh and the base state and is therefore known at
        # t = 0. This scheme's explicit half is dominated by ADVECTION, which at
        # t = 0 on a bubble at rest is exactly zero -- so an estimate made from
        # the initial state would report that any Δt whatsoever is safe.
        #
        # 20 m/s is an upper bound for this benchmark: the 2 K / 2000 m bubble
        # reaches |w| of roughly 12 m/s by t = 700 s.
        #-----------------------------------------------------------------------
        :imex_umax            => 20.0,
        # Free-slip lateral walls, taken from the deck's OWN boundary faces --
        # left/right (x) and front/back (y) here, all six faces of this box
        # being no-flux. The implicit operator has to know about them: it now
        # carries the horizontal mass flux, and letting sound run out through
        # the side walls would leave a stiff residual in the explicit half at
        # exactly the nodes where rhs! projects the normal momentum out.
        # :none is right for a laterally periodic case (and is what :auto picks
        # there anyway, since a periodic face is not a boundary face).
        :imex_lateral_walls   => Symbol(get(ENV, "DBG_LATWALL", "auto")),
        :imex_wall_flux       => true,
        # :RS freezes the operator's coefficients at qe; :PS refreshes them
        # from the solution every :imex_update_freq steps and refactorises the
        # preconditioner. As on rtb_hevi, :RS is right for THIS case -- beta and
        # thetabar barely move for a 2 K bubble on a 300 K background, so :PS
        # pays for a refactorisation and recovers a coefficient that was never
        # stale. :PS earns its cost when the solution departs far from any fixed
        # reference state.
        :imex_linearization   => Symbol(get(ENV, "DBG_LIN", _vdiff ? "PS" : "RS")),
        :imex_update_freq     => parse(Int, get(ENV, "DBG_UPDFREQ", "5")),
        # :warn prints the wedge amplification and lets the run proceed; :error
        # makes it fatal. Left at :warn because the wedge bound is itself
        # conservative -- it treats the viscous part of the explicit spectrum,
        # which is real and damping, as if it were imaginary.
        :imex_stability_guard => Symbol(get(ENV, "DBG_GUARD", "warn")),
        # Per-N-solves line reporting the average Krylov iteration count. Off by
        # default; the setup report already gives one measurement, and this is
        # for watching it drift as the flow develops.
        :imex_monitor         => parse(Bool, get(ENV, "DBG_IMEXMON", "false")),

        #---------------------------------------------------------------------------
        # HEVI options, for the :hevi arm of the comparison.
        #---------------------------------------------------------------------------
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
        # Physical parameters -- CompEuler/theta's, unchanged, as rtb_hevi uses.
        #---------------------------------------------------------------------------
        :lvisc                => parse(Bool, get(ENV, "DBG_VISC", "true")),
        :visc_model           => AV(),
        :μ                    => [0.0, 125.0, 125.0, 125.0, 125.0],
        :energy_equation      => "theta",

        #---------------------------------------------------------------------------
        # Mesh -- byte-identical to rtb_hevi's, deliberately: this case exists to
        # be compared against that one, and a comparison across two meshes is
        # not a comparison. Regenerate with
        #
        #   gmsh -3 problems/CompEuler/rtb_imex/rtb_10x1x10.geo
        #
        # (the .geo is named for an older element count and now writes the
        # 10 x 1 x 50 mesh: x in [-5, 5] km, y in [0, 5] km, z in [0, 10] km).
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/CompEuler/rtb_imex/hevi_10x1x50.msh",
        :lwarp                => false,
        :lstretch             => false,

        #---------------------------------------------------------------------------
        # Output -- one tree per scheme, so an A/B leaves both solutions on disk
        # to compare rather than the second overwriting the first.
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :output_dir           => scheme === :imex ? "./output_imex/" :
                                 scheme === :hevi ? "./output_hevi/" : "./output_explicit/",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        :loutput_pert         => true,

        #---------------------------------------------------------------------------
        # AMR: off. IMEX3D refuses :ladapt -- adaptation invalidates the column
        # topology the preconditioner is built on, its factorised matrices and
        # its gather/scatter plan.
        #---------------------------------------------------------------------------
        :linitial_refine      => false,
        :init_refine_lvl      => 1,
        :ladapt               => false,
    ) #Dict
    return inputs
end
