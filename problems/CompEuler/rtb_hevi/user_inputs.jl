#==============================================================================
 rtb_hevi -- rising thermal bubble, run with the HEVI time integrator.

 THIS CASE IS CompEuler/theta AS A 3D SLAB. Same 10 km x 10 km domain, same
 10 x 10 element resolution (Δx = Δz = 1000 m), same bubble (θ' = 2 K, linear
 taper, r0 = 2000 m, centred at 2500 m), same μ = 125, same tend = 1000 s --
 extruded one element in y because HEVI is 3D-only. Run it explicitly and it
 must reproduce theta; run it under HEVI and it must reproduce theta too. Any
 difference is the time integrator, and nothing else.

 VERIFIED, max|u| at 100 s intervals to t = 1000 s, 1 rank:

   explicit 3D slab vs 2D theta   within 1.6% at t = 100, under 1% thereafter
   HEVI     vs explicit 3D slab   within 0.07% at every sample

 (The slab is not bit-identical to theta and should not be: theta runs
 :SOL_VARS_TYPE => PERT() and this case runs the TOTAL default, and the 3D
 operators are not the 2D ones. Agreement at this level is the check.)

 That is the whole design. An earlier version of this case ran a 250 m bubble
 in a 1 km box at μ = 15 and looked far more diffused than theta at μ = 125,
 which reads as a code problem and is not one: how smeared a bubble looks is
 set by the diffusion length relative to the bubble, sqrt(μ t)/r0. theta sits
 at sqrt(125*1000)/2000 = 0.18; that case sat at sqrt(15*700)/250 = 0.41, more
 than twice as diffused on an 8x smaller μ. Matching the geometry is what
 makes the two comparable, so the geometry is matched here.

 WHAT THIS CASE CAN AND CANNOT SHOW
 ----------------------------------
 It exercises what the vertically-implicit split touches:

   * the vertical acoustic terms HEVI moves to the implicit side;
   * the buoyancy coupling -g·dρ inside the implicit operator;
   * a base state at rest, so any imbalance the split introduced would show up
     as the whole domain drifting rather than as a subtle error in the bubble.

 It CANNOT show HEVI to be faster, and it is not meant to. This mesh is
 isotropic, and on an isotropic mesh a vertically-implicit split buys nothing:
 there is no vertical acoustic term to remove that the horizontal one does not
 match. The setup report measures the implicit half at 2.8 1/s against the
 explicit half at 5.7 1/s, and ARS232's joint limit (0.284 s) lands BELOW
 CarpenterKennedy2N54's (~0.69 s).

 MEASURED here, 1 rank, to t = 1000 s:

   explicit  Δt = 0.5    2000 steps   0.040 s/step    81.8 s
   HEVI      Δt = 0.19   5263 steps   0.022 s/step   148.1 s     1.81x SLOWER

 HEVI is cheaper PER STEP even so (ARS232's 3 rhs! and 2 column solves against
 CarpenterKennedy2N54's 5 rhs!, and the column solves are trivial on 10 vertical
 elements) -- it just needs 2.6x as many of them. Cost per step was never the
 question; wall-clock to a fixed physical time is, and on this geometry the
 explicit scheme wins it.

 For the performance claim, point :gmsh_filename at rtb_8x1x24.msh -- 1 km box,
 8 x 1 x 24 elements, Δz = 41.7 m against Δx = 125 m -- and restore Δt = 0.02
 explicit / 0.035 HEVI with tend = 700. On that mesh the split is worth a
 measured 1.75x in Δt and ~2.5x in wall-clock: 0.185 s/step explicit against
 0.132 s/step under ARS232. Both meshes are kept in this directory precisely
 because one case cannot answer both questions.

 One element in y: the bubble is a cylinder along y, so the solution is
 y-invariant and the third dimension exists only because HEVI is 3D-only.
 |v| staying at round-off is a free correctness check.

 TO COMPARE AGAINST THE EXPLICIT SCHEME
 --------------------------------------
 Set `lexplicit = true` at the top of user_inputs() below. Each scheme runs at
 its own largest safe step, and the two solutions should be indistinguishable.

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
    # MEASURED on rtb_10x1x10.msh by the setup report (it prints these):
    #
    #   explicit half of the split   5.7 1/s
    #   implicit half of the split   2.8 1/s
    #   ARS232 joint neutral limit   0.284 s   -> 0.19 here, 67% of it
    #   CarpenterKennedy2N54 limit  ~0.69 s    -> 0.50 here, theta's step
    #
    # Read the first two numbers before reading anything else. The implicit
    # half is SMALLER than the explicit half, because this mesh is isotropic
    # (Δx = Δz = 1000 m): there is no vertical acoustic term for the split to
    # remove that the horizontal one does not match. HEVI's limit (0.284 s) is
    # consequently BELOW the explicit scheme's (~0.69 s) -- ARS232's explicit
    # imaginary radius is 1.73 against CarpenterKennedy2N54's 3.95, and on an
    # isotropic mesh that is the whole story. Expect HEVI to be SLOWER here.
    #
    # That is not a defect of the case; it is the point of it. This geometry
    # is theta's, so an explicit run must reproduce the 2D theta answer and a
    # HEVI run must reproduce it too. Correctness is what an isotropic mesh can
    # check. For the PERFORMANCE claim -- 1.75x in Δt and ~2.5x in wall-clock --
    # switch :gmsh_filename to rtb_8x1x24.msh, whose Δz is 3x finer than its Δx,
    # and restore Δt = 0.02 explicit / 0.035 HEVI with tend = 700.
    #
    # The setup report prints Δt as a percentage of the joint limit, warns above
    # 85%, and refuses outright if max|R| > 1. It refuses on the JOINT region,
    # not the explicit half: ARS343 has the largest explicit imaginary radius of
    # the ARS family (2.83 vs ARS232's 1.73) and a joint Δt_max of 0.0004 s on
    # the anisotropic mesh. Judge a tableau by ark_joint_amplification.
    #
    # Both steps are overridable for a sweep: DBG_DT=... DBG_TEND=...
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
        :gmsh_filename        => "./problems/CompEuler/rtb_hevi/rtb_10x1x10.msh",
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
