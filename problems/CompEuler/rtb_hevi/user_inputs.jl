#==============================================================================
 rtb_hevi -- rising thermal bubble, run with the HEVI time integrator.

 The dry nonhydrostatic benchmark (Robert 1993; Giraldo & Restelli, JCP 227,
 2008): a warm cosine bubble, θ' = 0.5 K over a 250 m radius centred at
 (500, 350) m, released into a neutrally stratified atmosphere at rest inside
 a 1000 m box. It accelerates upward under its own buoyancy and rolls into a
 mushroom by t ≈ 700 s.

 WHY THIS CASE, FOR THIS INTEGRATOR
 ----------------------------------
 It exercises exactly what the vertically-implicit split touches:

   * the vertical acoustic terms HEVI moves to the implicit side;
   * the buoyancy coupling -g·dρ inside the implicit operator;
   * a base state at rest, so any imbalance the split introduced would show up
     as the whole domain drifting rather than as a subtle error in the bubble.

 The mesh is DELIBERATELY anisotropic -- 8 x 1 x 24 elements, so Δz = 41.7 m
 against Δx = 125 m. On an isotropic bubble mesh HEVI buys nothing, because
 there is no vertical acoustic term to remove that the horizontal one does not
 match. Here the vertical spacing is 3x finer and the split is worth a measured
 1.75x in Δt and ~2.5x in wall-clock -- 0.185 s/step at Δt = 0.02 explicit
 against 0.132 s/step at Δt = 0.035 under ARS232. (At Δt = 0.05 it is 3.50x,
 but that is 97% of the joint stability limit and is not a step to run at; see
 the note on Δt below.) That is the claim this case is here to make checkable.

 One element in y: the bubble is a cylinder along y, so the solution is
 y-invariant and the third dimension exists only because HEVI is 3D-only.
 |v| staying at round-off is a free correctness check.

 TO COMPARE AGAINST THE EXPLICIT SCHEME
 --------------------------------------
 Set `lexplicit = true` at the top of user_inputs() below. The explicit run needs 2.5x the number of
 steps to reach t = 700 s, and the two solutions should be indistinguishable --
 the standalone version of this case checks exactly that.

 A version of this case that needs no mesh, no Gridap and no MPI stack --
 and that asserts all of the above rather than leaving it to the eye -- is
 test/hevi/test_rtb.jl.
==============================================================================#
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

    # Step sizes.
    #
    # HEVI at 0.035 s with ARS232 -- 70% of the joint IMEX limit, NOT at it.
    #
    # ark_joint_dt_max puts the neutral limit at 0.0515 s on this mesh, where the
    # explicit half of the split sees 26.7 1/s and the implicit half 69.4 1/s.
    # This case shipped at 0.05, which is 97% of that, and it was stable on 1 and
    # 3 ranks and divergent at t ~ 5 s on 10. At the limit the scheme is EXACTLY
    # neutral, and three things move it by a few percent: the explicit rate is a
    # node-spacing estimate calibrated against one measured spectrum, it assumes
    # the flow speeds at t = 0, and the assembly at rank-shared nodes is mildly
    # partition-dependent so the effective rates shift with the rank count. None
    # of those matters at 70%; all of them decide the run at 97%.
    #
    # The setup report now prints Δt as a percentage of the limit and warns above
    # 85%, and still refuses outright if max|R| > 1.
    #
    # This case ran :ARS343 at 0.03 until the joint region was measured. ARS343
    # has the largest EXPLICIT imaginary radius of the ARS family (2.83 against
    # ARS232's 1.73) and a joint Δt_max of 0.0004 s: it amplifies by 2.4% per
    # step here, which took ~15 s of model time to become visible and moved
    # with the rank count. Judge a tableau by ark_joint_amplification.
    #
    # The explicit limit measured between 0.02 and 0.04, so 0.02 has margin.
    #
    # To find either limit, sweep :Δt with a short :tend (say 30.0) and watch
    # for the pressure going negative. The CFL table printed at startup reports
    # the CFL each run is actually sitting at.
    Δt   = lexplicit ? 0.02 : 0.035
    tend = 700.0

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
        :diagnostics_at_times => (0:50:tend),
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
        # A small constant artificial viscosity keeps the bubble edge smooth on
        # a mesh this coarse (16 points across the bubble diameter). The value
        # is CompEuler/theta's 125 scaled by the ratio of minimum LGL gaps,
        # 21.6/172.7 -- that case is the same benchmark on a 10 km domain.
        # It is small enough not to bind Δt: the parabolic limit here is ~0.5 s
        # against the 0.1 s the acoustic terms allow.
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :visc_model           => AV(),
        :μ                    => [0.0, 15.0, 15.0, 15.0, 15.0],
        :energy_equation      => "theta",
        #---------------------------------------------------------------------------
        # Mesh
        #
        # Regenerate with:  gmsh -3 problems/CompEuler/rtb_hevi/rtb_8x1x24.geo
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/CompEuler/rtb_hevi/rtb_8x1x24.msh",
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
