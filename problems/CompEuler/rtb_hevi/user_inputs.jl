function user_inputs()

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

    Δt   = parse(Float64, get(ENV, "DBG_DT", lexplicit ? "0.5" : "0.19"))
    tend = parse(Float64, get(ENV, "DBG_TEND", "1000.0"))

    inputs = Dict(
        #---------------------------------------------------------------------------
        # Time integration.  :lcfl_report below prints the stability table for
        # whichever scheme is running, so a run that is silently over its limit
        # is visible before it diverges.
        #---------------------------------------------------------------------------
        :ode_solver           => lexplicit ? CarpenterKennedy2N54() : HEVI_ARK(:ARS232),
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
        :lcfl_report          => true,   # print the stability table at startup
        :hevi_verify          => true,   # setup self-check; cheap, leave it on
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
        :lvisc                => true,
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
        :linitial_refine      => true,
        :init_refine_lvl      => 1,
        :ladapt               => false,
    ) #Dict
    return inputs
end
