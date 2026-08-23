function user_inputs()

    # Implicit vertical diffusion is opt-in per run: DBG_VDIFF=1. Read into a
    # local here rather than inline because two keys below have to agree about
    # it -- switching it on also switches the linearisation to :PS, without
    # which the implicit operator would carry no diffusion at all.
    _vdiff = parse(Bool, get(ENV, "DBG_VDIFF", "false"))

    #---------------------------------------------------------------------------
    # HEVI or fully explicit -- ONE switch, right here.
    #
    #   lexplicit = false  ->  HEVI_ARK(:ARS232)        -> ./output_hevi/
    #   lexplicit = true   ->  CarpenterKennedy2N54()   -> ./output_explicit/
    #
    # One switch rather than two commented-out lines, because a wall-clock
    # comparison is only worth something if mesh, viscosity, output cadence and
    # tend are identical between the two runs, and swapping comments by hand is
    # how one of those quietly ends up different.
    #
    # Each scheme runs at ITS OWN largest safe step. Timing them at a common Δt
    # would measure cost per step, which is not the question: HEVI costs more
    # per step (ARS232: 3 rhs! evaluations, 3 applications of the implicit
    # operator and 2 column solves, against CarpenterKennedy2N54's 5 rhs! and no
    # solve) and pays for it by taking a bigger one. Only wall-clock to a fixed
    # physical time answers it.
    #---------------------------------------------------------------------------
    lexplicit = false
    #lexplicit = true

    Δt   = parse(Float64, get(ENV, "DBG_DT", lexplicit ? "0.2" : "0.35"))
    tend = parse(Float64, get(ENV, "DBG_TEND", "1000.0"))

    inputs = Dict(
        #---------------------------------------------------------------------------
        # Time integration.  :lcfl_report below prints the stability table for
        # whichever scheme is running, so a run that is silently over its limit
        # is visible before it diverges.
        #---------------------------------------------------------------------------
        :ode_solver           => get(ENV,"DBG_SPLIT","0")=="1" ? SPLIT_EXPLICIT() : (lexplicit ? CarpenterKennedy2N54() : HEVI_ARK(:ARS232)),
        :Δt                   => Δt,

        :tinit                => 0.0,
        :tend                 => tend,
        # Diagnostics every 50 s, and a per-step heartbeat.
        #
        # Both matter for a case like this: at Δt = 0.05 the FIRST diagnostic
        # output is 1000 steps in, and the integrator warm-up above it does not
        # print when it finishes (that line is commented out in
        # TimeIntegrators.jl). So without a heartbeat the run is silent for
        # minutes after "PATIENCE: ONLY DONE ON 1st RUN", which is
        # indistinguishable from a hang -- and it is not one.
        #
        # For a quick smoke run, set tend = 50.0 above.
        #
        # The heartbeat is also the timing instrument: it reports s/step
        # measured over the interval since the PREVIOUS heartbeat, so after the
        # first few steps it is a steady-state rate with the JIT excluded.
        :diagnostics_at_times => (0:100:tend),
        :lstep_heartbeat      => true,
        :lsource              => true,
        :restart_time         => 1.0e7,
        :lrestart             => false,
        #---------------------------------------------------------------------------
        # HEVI options. All optional; these are the defaults, spelled out so the
        # case documents them. See src/kernel/solvers/hevi/README.md.
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
        :lcfl_report          => true,   # print the stability table at startup
        # :lcfl_report_every  re-print that table every N steps. :lcfl_report shows it
        #                   ONCE, at startup, where nu_t is ~0 on a laminar sounding --
        #                   so its viscous row says diffusion never limits dt, and says
        #                   it truthfully, at t = 0. Watch that row cross the line as the
        #                   boundary layer spins up. Collective and cheap; 0 = off.
        :lcfl_report_every    => parse(Int, get(ENV, "DBG_CFL_EVERY", "0")),
        :hevi_verify          => parse(Bool, get(ENV,"DBG_VERIFY","true")),   # setup self-check; cheap, leave it on
        #-----------------------------------------------------------------------
        # LINEARISATION OF THE IMPLICIT OPERATOR -- the switch, and why this
        # deck picks :PS.
        #
        #   :RS  coefficients frozen at the reference state qe for the whole run
        #   :PS  coefficients refreshed from the solution every
        #        :hevi_update_freq steps, then refactorised
        #
        # Giraldo, de Braganca Alves, Kelly, Kang & Reinecke (arXiv:2311.11425)
        # time three HEVI variants in NUMA -- a spectral element code, so the
        # comparison transfers -- and find LINEAR HEVI fastest by a wide margin:
        # 5x NHEVI-LU and 10x NHEVI-GMRES, with the second- and third-order ARK
        # pairs the best of the five they test. This code is already linear
        # HEVI, so that result is a result FOR this scheme, not an argument to
        # change it. :PS is the sub-variant they recommend, refreshed every five
        # steps.
        #
        # THIS DECK USES :RS ANYWAY, because on THIS case :PS is measurably
        # worse. Measured, 1 rank, tend = 100 s:
        #
        #   :RS   0.207 s/step    max|u| = 2.6760e+00
        #   :PS   0.233 s/step    max|u| = 2.6760e+00     +12.6% for nothing
        #
        # and with the stability guard off both are still stable at Δt = 0.45,
        # so :PS buys no step size either. That is what should be expected here:
        # beta = gamma p/(rho theta) and thetabar barely move for a 2 K bubble
        # on a 300 K background, so refreshing them costs a refactorisation and
        # recovers a coefficient that was never stale.
        #
        # :PS earns its cost when the solution departs far from any fixed
        # reference state -- their 100-day baroclinic instability, real data, or
        # a whole-atmosphere run with large day/night temperature swings. Switch
        # with DBG_LIN=PS or by editing the line below.
        #-----------------------------------------------------------------------
        :hevi_linearization   => Symbol(get(ENV, "DBG_LIN", _vdiff ? "PS" : "RS")),
        :hevi_update_freq     => parse(Int, get(ENV, "DBG_UPDFREQ", "5")),
        :hevi_wall_flux       => true,   # zero implicit vertical mass flux at floor/lid
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # Physical parameters
        #
        # CompEuler/theta's viscosity, unscaled, because this is now theta's
        # domain, theta's resolution and theta's bubble -- see initialize.jl.
        #
        # It was previously 15, theta's 125 scaled by the ratio of minimum LGL
        # gaps (21.6/172.7) for the old 1 km domain. Scaling μ that way keeps
        # the CELL Reynolds number fixed, which is the right thing for keeping
        # a bubble edge smooth, but it does NOT keep the solution looking the
        # same: how diffused the bubble looks is set by sqrt(μ t)/r0, and
        # shrinking r0 by 8x while scaling μ by 8x moves that ratio from
        # theta's 0.18 to 0.41. Same-domain, same-μ removes the question.
        #---------------------------------------------------------------------------
        :lvisc                => parse(Bool, get(ENV,"DBG_VISC","true")),
        :visc_model           => AV(),
        :μ                    => [0.0, 125.0, 125.0, 125.0, 125.0],
        :energy_equation      => "theta",
        #---------------------------------------------------------------------------
        # Mesh -- theta's 10 km box at theta's resolution, one element thick.
        #
        # Regenerate with:  gmsh -3 problems/CompEuler/rtb_hevi/rtb_10x1x10.geo
        #
        # rtb_8x1x24.msh (1 km box, Δz 3x finer than Δx) is kept alongside for
        # the HEVI PERFORMANCE claim, which needs an anisotropic mesh. This one
        # is isotropic and exists for CORRECTNESS: explicit here must reproduce
        # the 2D theta answer.
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        #:gmsh_filename        => "./problems/CompEuler/rtb_hevi/rtb_10x1x10.msh",
        :gmsh_filename        => "./problems/CompEuler/rtb_hevi/hevi_10x1x50.msh",
        :lwarp                => false,
        :lstretch             => false,
        #---------------------------------------------------------------------------
        # Output
        #---------------------------------------------------------------------------
        # Separate output trees, so an A/B leaves both solutions on disk to
        # compare rather than the second overwriting the first.
        :outformat            => "vtk",
        :output_dir           => lexplicit ? "./output_explicit/" : "./output_hevi/",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        :loutput_pert         => true,
        #---------------------------------------------------------------------------
        # AMR: off. HEVI refuses :ladapt -- adaptation invalidates the column
        # topology, the factorised column matrices and the gather/scatter plan.
        #---------------------------------------------------------------------------
        :linitial_refine      => false,
        :init_refine_lvl      => 1,
        :ladapt               => false,
    ) #Dict
    return inputs
end
