#---------------------------------------------------------------------------------
# SWEEP SWITCHES. Every knob this case is being debugged on is an environment
# variable, so a variant is one line on the command line and never a deck edit:
#
#   JEXPRESSO_M7_SENSOR   "legacy" (default) | "residual"
#   JEXPRESSO_M7_NORMS    "domain" (default) | "rank"
#   JEXPRESSO_M7_MU1      slot-1 (β∇ρ) multiplier              default 1.0
#   JEXPRESSO_M7_MU       slots 2-4 (μ, κ) multiplier          default 4.0
#   JEXPRESSO_M7_CMAX     :dsgs_Cmax, the cap constant         default 0.5
#   JEXPRESSO_M7_CMIN     :dsgs_Cmin, background floor         default 0.0
#   JEXPRESSO_M7_FILTER   modal-filter blend μ_x, 0 = off      default 0.0
#   JEXPRESSO_M7_DT       :Δt                                  default 5.0e-8
#   JEXPRESSO_M7_TEND     :tend                                default 3.5e-3
#   JEXPRESSO_M7_REF      initial refinement level, 0 = off    default 0
#
# Every default reproduces the deck as it stands, so an unset environment is
# the baseline run and nothing below changes behaviour on its own.
#---------------------------------------------------------------------------------
_m7_s(k, d)  = get(ENV, k, d)
_m7_f(k, d)  = parse(Float64, get(ENV, k, string(d)))
_m7_i(k, d)  = parse(Int,     get(ENV, k, string(d)))

function user_inputs()

    m7_tend   = _m7_f("JEXPRESSO_M7_TEND", 3.5e-3)
    m7_filter = _m7_f("JEXPRESSO_M7_FILTER", 0.005)
    m7_ref    = _m7_i("JEXPRESSO_M7_REF", 0)

    inputs = Dict(
        #---------------------------------------------------------------------------
        #
        # 2D CompEuler: Mach-7 supersonic flow over a forward-facing step.
        #
        # THIS IS CompEuler/ffs_step WITH M∞ = 7 INSTEAD OF 3, and nothing
        # else changed that is not forced by that one number. Same tunnel,
        # same mesh, same boundary conditions, same fluxes, same DynSGS
        # settings. The two decks differ in exactly four values:
        #
        #     M∞     3      ->  7        (initialize.jl, ffs_freestream)
        #     Δt     1.0e-7 ->  5.0e-8   scaled with |u|+c, see below
        #     tend   8.0e-3 ->  3.5e-3   the same tunnel flow-through count
        #     output 5.0e-5 ->  2.5e-5   the same cadence in flow-throughs
        #
        # WHY IT EXISTS. rampCaoEtAl2021 on an unstretched grid goes
        # non-finite at t = 4.8e-7 (step 481 at its Δt = 1e-9), in slot 1,
        # on every rank at once. Two things are new in that deck at the same
        # time — the Mach number (7.7) and the ramp geometry/grid — so the
        # failure names neither. ffs_step runs, so raising ONLY the Mach
        # number on ffs_step splits the question in two:
        #
        #   * this deck reaches tend -> Mach 7 is not by itself the problem,
        #     and the ramp's grid, its BCs and its leading-edge/corner
        #     treatment are what to look at;
        #   * this deck dies -> the same failure is reproduced on a
        #     configuration whose mesh, BCs, fluxes and DynSGS settings are
        #     all validated at Mach 3, which is a much smaller thing to
        #     debug — and the levers are then tried here, not on the ramp.
        #
        # Mach 7 is the rung below the ramp's 7.7, close enough that a deck
        # surviving here makes the Mach number an unlikely sole culprit.
        #
        # Wind tunnel 3 m x 1 m with a 0.2 m step 0.6 m from the inflow; the
        # tunnel is filled with, and fed from the left by, a uniform Mach-7
        # stream of air at p = 101325 Pa, T = 293 K (|u| = 2400.4 m/s,
        # against 1028.7 m/s at Mach 3). The thermodynamic state is that of
        # the Loci/STREAM "2D Supersonic Forward Step" tutorial, i.e. the
        # dimensional form of Emery (1968) / Woodward & Colella (1984), also
        # Section 5.1 of Nazarov & Hoffman (IJNMF 71:339-357, 2013); only M∞
        # is raised. The gas stays calorically perfect — ideal-gas Euler, no
        # dissociation, no vibrational excitation — so this is a NUMERICAL
        # test at Mach 7, not a physical model of Mach-7 air.
        #
        # A bow shock stands off the step, reflects from the roof, and the
        # reflections merge into a Mach stem; at Mach 7 the shock stands
        # closer to the step and the shock layer is thinner. Everything
        # interesting in this case is a discontinuity, which drives the two
        # decisions below — both inherited unchanged from ffs_step.
        #
        # (1) :energy_equation => "energy".  Slot 4 is ρE, not ρθ. ρθ is an
        #     entropy variable: it is conserved across a contact but NOT
        #     across a shock, so the Euler-θ system carries the wrong shock
        #     speed no matter how it is stabilized. The entropy jump across
        #     the bow shock is larger at Mach 7, so this matters more here.
        #
        # (2) :visc_model => DSGS().  Residual-based artificial viscosity as
        #     shock capturing, Nazarov & Hoffman eq. (3.4)-(3.7): the
        #     viscosity is proportional to the local residual of the
        #     conservation laws, so it appears at the shocks and stays near
        #     zero in the smooth 90% of the field. A constant AV() coefficient
        #     large enough to hold the Mach-7 shocks would smear the whole
        #     domain. The total-energy branch of compute_dsgs_viscosity! (2D)
        #     is selected by :energy_equation above.
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(),
        :tinit                => 0.0,
        # One tunnel flow-through is 3 m / 2400.4 m/s = 1.25e-3 s at Mach 7
        # (2.92e-3 s at Mach 3). ffs_step integrates 2.74 of them; 3.5e-3
        # here is 2.80, i.e. the same picture at the same stage.
        :tend                 => m7_tend,          # ≈ 2.8 tunnel flow-throughs
        :lrestart             => false,
        :restart_time         => 0.0,
        # CFL. The grid is h = 0.025 m with :nop => 4, so the tightest LGL
        # node spacing is ≈ 0.0043 m; against the free-stream wave speed
        # |u| + c ≈ 1372 m/s that puts the ADVECTIVE limit near 1.4e-6 s.
        # Those are the numbers for the RAW gmsh grid. :init_refine_lvl
        # below halves the element size once per level, so every length
        # above shrinks by 2^lvl and the viscous limit by 4^lvl. All the
        # measured constants in this deck — the 14.95 cap, the 0.22, and
        # the whole sweep below — are at :init_refine_lvl => 1, i.e.
        # h = 0.0125 m, smallest LGL gap 0.00216 m, advective limit
        # 1.6e-6 s. Re-measure them if you change the level.
        #
        # The binding constraint is NOT advective, it is VISCOUS. DynSGS
        # saturates its own μ_max bound at the step corner (measured: μ =
        # 14.4 Pa·s against a local cap of 14.95), and μΔt/(ρΔx²) at
        # Δt = 5e-7 is already ≈ 0.22 with :μ => 1.0. Scaling :μ up without
        # scaling Δt down therefore blows the viscous limit — see the note
        # on :μ below. The pair (Δt, :μ) has to move together.
        #
        # WHAT MACH 7 DOES TO THAT. The free-stream characteristic speed is
        # |u| + c, and it is the one number this case changes:
        #
        #     Mach 3:  1028.7 + 342.9 = 1371.6 m/s
        #     Mach 7:  2400.4 + 342.9 = 2743.3 m/s      exactly 2x
        #
        # so the advective limit halves and Δt halves with it, 1.0e-7 ->
        # 5.0e-8. The VISCOUS number comes out UNCHANGED by that pair of
        # moves: μ_max = C_max·Δ·(|u|+c) doubles with the wave speed while
        # Δt halves, so μΔt/(ρΔx²) lands where ffs_step measured it. That is
        # why :μ below is not touched — ffs_step's warning that dissipation
        # helps only when Δt is cut to match applies here in reverse, Δt
        # having already been cut to pay for the dissipation Mach 7 brings
        # on its own.
        :Δt                   => _m7_f("JEXPRESSO_M7_DT", 5.0e-8),
        :diagnostics_at_times => (0:(m7_tend/140):m7_tend),
        # Wall-clock note, not a setting: at Δt = 5.0e-8 the diagnostics
        # above are 500 steps apart and the whole run is 70000 steps — the
        # same step count and the same per-step cost as ffs_step, so the
        # same order of wall clock. A long silence after "Integrator warm-up
        # with real callbacks" is the run working, not a hang.
        # `JEXPRESSO_STEP_HEARTBEAT=1` turns on a per-step trace without
        # editing this deck.
        #
        # FOR A FIRST LOOK, the question this deck exists to answer does not
        # need the full run: the step corner is where ffs_step failed in
        # every row of the sweep below, and it fails early or not at all, so
        # :tend => 2.0e-4 (4000 steps) already answers it.
        :lsource              => false,
        :SOL_VARS_TYPE        => TOTAL(),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 3,               # polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters / constants
        #---------------------------------------------------------------------------
        :energy_equation      => "energy",        # slot 4 is ρE — see note (1)
        :lvisc                => true,
        # DynSGS sensor: "legacy" = the sensor this case was validated with
        # (the assembled RHS against a fixed BDF2 of the stage state, in
        # effect a |∂ₜq| sensor); "residual" (the default) = the element-wise
        # strong residual with the stage-consistent stencil, DSGS.md §1.2.
        :visc_model           => DSGS(),          # residual-based shock capturing
        # JEXPRESSO_M7_SENSOR. "legacy" is R ≈ |∂ₜq| (rhs.jl:205), a RATE
        # sensor: it fires on a moving front and is blind to an oscillation
        # that is merely standing there. "residual" is the stage-consistent
        # element residual, which is O(1/h) at a discontinuity whether or not
        # it is moving. That difference is the first hypothesis for the
        # per-element beads along the oblique shock.
        :dsgs_sensor          => _m7_s("JEXPRESSO_M7_SENSOR", "legacy"),
        # Startup hold OFF — the single difference that this case cannot
        # absorb. 7dd6f0c holds the coefficient at zero until the BDF2
        # history is a time derivative (default 2 steps, 3 rotations),
        # because the sensor was reading a SMOOTH initial condition as
        # unresolved and pinning ν at its cap on step one. This initial
        # condition is not smooth: a Mach-3 stream is started impulsively
        # against the step, and the whole transient is at the step face and
        # the convex corner. Holding ν at zero there integrates the most
        # violent steps of the run with no dissipation at all, and the
        # oscillation it plants at the corner is what the rest of the run
        # has to carry. sm/newmaster predates the hold, never holds, and
        # runs this case to t = 8e-3; with the hold on it dies at 1.46e-3,
        # in that corner. 0 restores the older behaviour: the sensor fires
        # from the first call on a history seeded from the initial
        # condition, which reads 1.5x the forward difference of q — an
        # over-estimate of the rate, which at an impulsive start is the
        # side to err on.
        # It matters more at Mach 7, not less: the impulsive start is
        # twice as violent, so the steps the hold would integrate without
        # any dissipation are exactly the ones this deck can least afford.
        :dsgs_hold_steps      => 0,
        # Per-equation multiplier on the DynSGS coefficient, INHERITED FROM
        # ffs_step UNCHANGED — see the Δt note above for why Mach 7 does not
        # by itself call for more. The method is parameter-free, so 1.0 is
        # the paper's own setting; the ×4 on the momentum and energy slots
        # is the Mach-3 case's, and it is measured, not guessed — see the
        # sweep below. EVERY ROW OF THAT SWEEP IS AT MACH 3: read its Δt
        # column in Mach-3 units and halve it for this deck.
        #
        # NOTE slot 1 is NOT zero: the total-energy DynSGS carries the
        # density diffusion β∇ρ of eq. (3.3), and it is what keeps the
        # density jump across the bow shock from ringing (the Euler-θ path
        # drops it, following Marras eq. 10). It is load-bearing here —
        # setting it to 0.0 is the single most destabilising change
        # measured on this case — and at Mach 7 the density jump it has to
        # hold is 5.5x across a normal shock rather than 3.9x.
        #
        # WHAT THE STEP CORNER COSTS. (0.6, 0.2) is a convex corner, i.e. a
        # geometric singularity sitting in an expansion fan, and it governs
        # how long this case survives. Runs to t = 2e-3 s, all with the
        # corner BC of user_bc.jl in place, failing on p < 0 in soundSpeed:
        #
        #   :μ [1,1,1,1]  Δt 5.0e-7   fails 4-6e-4     (viscous CFL ≈ 0.22)
        #   :μ [1,1,1,1]  Δt 1.25e-7  fails   4e-4     Δt is NOT the cause
        #   :μ [1,4,4,4]  Δt 5.0e-7   fails  <2e-4     viscous limit breached
        #   :μ [1,8,8,8]  Δt 5.0e-7   fails  <2e-4     ditto, worse
        #   :μ [0,1,1,1]  Δt 5.0e-7   fails  <2e-4     β∇ρ off — worst of all
        #   :μ [1,1,1,1]  Δt 5.0e-7   fails 6-8e-4     with :nop => 3
        #   :μ [1,4,4,4]  Δt 1.25e-7  past 8e-4        <- the only row that
        #                                                     survives; set here
        #
        # WHAT THE SWEEP WAS MEASURED UNDER, because two DynSGS defaults
        # moved after it and neither is a setting of this deck:
        #
        #   * :dsgs_norms was rank-local; it defaults to "domain" since
        #     10b177a (2026-09-12). The domain denominator is the spread
        #     over the WHOLE field, which the bow shock sets, so it is
        #     larger than the local spread on the rank holding the step
        #     corner — the same :μ therefore buys LESS viscosity there
        #     than it did in this table. Pinned explicitly below.
        #   * the normalization floor rose from 1e-3 of each variable's
        #     physical scale to the scale itself, 7dd6f0c (2026-09-11),
        #     :dsgs_rel. It binds only where a variable is nearly
        #     uniform, which here is the free stream at startup, not the
        #     developed field.
        #
        # Both move ν DOWN relative to the rows above, so treat the
        # survival times as optimistic and :μ [1,4,4,4] as the floor of
        # what this case needs, not the ceiling.
        #
        # The pattern: dissipation helps only when Δt is cut to match, and
        # cutting Δt alone does nothing. If THIS deck fails at the corner,
        # the levers in order are (a) Δt down another factor 2 WITH :μ up to
        # [1,8,8,8], (b) :nop => 3, (c) :init_refine_lvl => 2 or the ref = 2
        # grid, (d) a Woodward & Colella corner entropy fix. Raising :μ
        # alone will only blow the viscous limit sooner — that is measured,
        # not a guess.
        # JEXPRESSO_M7_MU1 (slot 1) and JEXPRESSO_M7_MU (slots 2-4). Both
        # multiply the coefficient AFTER the cap (SGS.jl:1481), so they still
        # bite when μ is pinned at μ_cap — but they cannot help where the
        # sensor has switched the coefficient off, which is what makes the
        # sensor the first thing to test.
        :μ                    => [_m7_f("JEXPRESSO_M7_MU1", 1.0),
                                  _m7_f("JEXPRESSO_M7_MU",  4.0),
                                  _m7_f("JEXPRESSO_M7_MU",  4.0),
                                  _m7_f("JEXPRESSO_M7_MU",  4.0)],
        # JEXPRESSO_M7_CMAX. μ_cap = Cmax·Δ·ρ_max·(|u|+c) (SGS.jl:1477). The
        # coarse grid survived 3x longer than the refined one at the SAME Δt,
        # which is what a binding cap looks like: Δ doubled, so did the cap.
        :dsgs_Cmax            => _m7_f("JEXPRESSO_M7_CMAX", 0.5),
        # JEXPRESSO_M7_CMIN. An UNCONDITIONAL floor μ_fl = Cmin·Δ·ρ_max·(|u|+c)
        # (SGS.jl:1478), i.e. a bounded cell Reynolds number of 1/Cmin, applied
        # whatever the sensor says. This is the knob for the structures that
        # convect downstream of the corner: an expansion fan is a SMOOTH
        # solution, so a residual sensor correctly returns ~0 in it, and the
        # entropy layer the corner sheds is a contact-type feature that never
        # self-heals — nothing else in this deck damps it. Unlike the modal
        # filter it acts on every mode and scales with Δ and the wave speed.
        # brioWu1d runs 0.06; the flux-emergence case 0.03. Try 0.02-0.05.
        :dsgs_Cmin            => _m7_f("JEXPRESSO_M7_CMIN", 0.0),
        # Artificial Prandtl number P of eq. (3.7): κ = P/(γ-1)·μ. Nazarov &
        # Hoffman use P ≈ 0.1.
        :Pr                   => 0.1,
        # Scope of the DynSGS normalising scales ⟨q⟩ and ‖q−⟨q⟩‖ — the
        # paper's Ω, i.e. the whole domain, which is also the default since
        # 10b177a. Pinned rather than left implicit: this case is run on 64
        # ranks and the sweep above was measured under the OLD rank-local
        # default, so the scope has to be visible in the deck to be
        # comparable. Costs 2 Allreduce of three doubles per RHS call, 10
        # per step here; no effect on a serial run. "rank" reverts to the
        # pre-September-2026 behaviour and makes the answer depend on the
        # partition. See ENVIRONMENT_VARIABLES.md.
        # JEXPRESSO_M7_NORMS. NOTE only "domain" and "rank" exist for the
        # CompEuler total-energy kernel — it takes a lglobal_norms::Bool
        # (SGS.jl:1333). "element" is DSGS_MHD only (SGS.jl:633), so the
        # per-element normalization that would remove the weak-feature
        # starvation outright is NOT available here; "rank" only shrinks Ω,
        # and makes the answer depend on the partition.
        :dsgs_norms           => _m7_s("JEXPRESSO_M7_NORMS", "domain"),
        #---------------------------------------------------------------------------
        # Mesh
        #
        # ffs_step_transfinite.msh is a three-block transfinite quad mesh of
        # the fluid L-shape, h = 0.025 m uniform (4032 elements) — a
        # byte-for-byte copy of ffs_step's, kept in this directory so the
        # case is self-contained. Its physical
        # curve groups — "inflow", "outflow", "wall" — are the tags that
        # reach user_bc_dirichlet!. Regenerate at twice the resolution by
        # setting ref = 2 in ffs_step_transfinite.geo.
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/CompEuler/ffs_step_M7/ffs_step_transfinite.msh",
        #---------------------------------------------------------------------------
        # Plotting
        #---------------------------------------------------------------------------
        # JEXPRESSO_M7_FILTER: the Boyd-Vandeven modal filter, blended as
        # F = μ_x·(L W L⁻¹) + (1-μ_x)·I (filter.jl:791-800).
        #
        # AT :nop => 4 THIS IS A TOP-MODE KILLER AND NOTHING ELSE. The
        # Boyd-Vandeven transfer function only acts on k > 2n/3, which at
        # n = 4 leaves the weights [1, 1, 1, 0.9957, 0] — modes 0-2 exactly
        # untouched, mode 3 cut by 0.4% at FULL strength, mode 4 annihilated.
        # ("exp" and "quad" reduce to the same thing at this order.) So μ_x
        # is not an amplitude, it is a RATE: the top mode decays as
        # (1-μ_x) per RHS call, and filter! runs inside rhs! (rhs.jl:799),
        # i.e. once per RK stage, five times a step. The top mode therefore
        # e-folds in 1/(5·μ_x) steps:
        #
        #   μ_x = 0.01   -> 20 steps      μ_x = 1e-4 -> 2000 steps
        #   μ_x = 0.005  -> 40 steps      μ_x = 8e-6 -> one flow-through
        #
        # Anything fast enough to catch a Gibbs mode is a complete P4 -> P3
        # truncation over the run; anything gentle enough to leave the mode
        # alive is too slow to matter. 0.005 is the value to try.
        :lfilter              => (m7_filter > 0.0),
        :mu_x                 => m7_filter,
        :mu_y                 => m7_filter,
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        #:output_dir           => "/scratch/smarras/smarras/output/shock/",
        :output_dir           => "./output",
        :loutput_pert         => false,           # plot the total state
        # Numerical schlieren from ρ, computed at output times only
        # (kernel/physics/schlieren.jl). Adds two point-data fields to the
        # VTU on top of :outvars —
        #   schlieren_grad_rho  |∇ρ| [kg/m⁴], quantitative
        #   schlieren           exp(-k|∇ρ|/max|∇ρ|), the picture
        # For the familiar dark-shock look, colour "schlieren" with a
        # REVERSED greyscale in ParaView. This is the field to look at on
        # this case: the bow shock, its reflection off the roof, the Mach
        # stem and the slip line downstream of the triple point are all
        # density features, and the exponential map keeps the weak ones
        # visible next to the strong bow shock — which at Mach 7 is stronger
        # still, so the exponential map earns its place here even more.
        :lschlieren           => true,
        :schlieren_k          => 20.0,            # contrast; Hadjadj uses 10-100
        #---------------------------------------------------------------------------
        # AMR off: the mesh already resolves the shocks at h/nop = 1/80, and
        # DynSGS is what handles what is left under-resolved.
        #---------------------------------------------------------------------------
        :linitial_refine      => (m7_ref > 0),
        :init_refine_lvl      => max(m7_ref, 1),
        :ladapt               => false,
    ) #Dict

    return inputs

end
