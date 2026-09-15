function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        #
        # 2D CompEuler: Mach-3 supersonic flow over a forward-facing step.
        #
        # Wind tunnel 3 m x 1 m with a 0.2 m step 0.6 m from the inflow; the
        # tunnel is filled with, and fed from the left by, a uniform Mach-3
        # stream of air at p = 101325 Pa, T = 293 K (|u| ≈ 1029 m/s), per the
        # Loci/STREAM "2D Supersonic Forward Step" tutorial. This is the
        # dimensional form of Emery (1968) / Woodward & Colella (1984), also
        # Section 5.1 of Nazarov & Hoffman (IJNMF 71:339-357, 2013).
        #
        # A bow shock stands off the step, reflects from the roof, and the
        # reflections merge into a Mach stem. Everything interesting in this
        # case is a discontinuity, which drives the two decisions below.
        #
        # (1) :energy_equation => "energy".  Slot 4 is ρE, not ρθ. ρθ is an
        #     entropy variable: it is conserved across a contact but NOT
        #     across a shock, so the Euler-θ system carries the wrong shock
        #     speed no matter how it is stabilized.
        #
        # (2) :visc_model => DSGS().  Residual-based artificial viscosity as
        #     shock capturing, Nazarov & Hoffman eq. (3.4)-(3.7): the
        #     viscosity is proportional to the local residual of the
        #     conservation laws, so it appears at the shocks and stays near
        #     zero in the smooth 90% of the field. A constant AV() coefficient
        #     large enough to hold the Mach-3 shocks would smear the whole
        #     domain. The total-energy branch of compute_dsgs_viscosity! (2D)
        #     is selected by :energy_equation above.
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(),
        :tinit                => 0.0,
        :tend                 => 8.0e-3,          # ≈ 2.7 tunnel flow-throughs
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
        :Δt                   => 1.25e-7,
        :diagnostics_at_times => (0:5.0e-5:8.0e-3),
        # Wall-clock note, not a setting: at Δt = 1.25e-7 the diagnostics
        # above are 3200 steps apart, so the CFL/VTK lines are ~35-40 min
        # apart and the whole run is 64000 steps, order 12 h on one core.
        # A long silence after "Integrator warm-up with real callbacks" is
        # the run working, not a hang. If you ever want to watch it step,
        # `JEXPRESSO_STEP_HEARTBEAT=1` turns on a per-step trace without
        # editing this deck.
        :lsource              => false,
        :SOL_VARS_TYPE        => TOTAL(),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,               # polynomial order
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
        :dsgs_sensor          => "legacy",
        # Per-equation multiplier on the DynSGS coefficient. The method is
        # parameter-free, so 1.0 is the paper's own setting; the ×4 on the
        # momentum and energy slots is this case's, and it is measured, not
        # guessed — see the sweep below.
        #
        # NOTE slot 1 is NOT zero: the total-energy DynSGS carries the
        # density diffusion β∇ρ of eq. (3.3), and it is what keeps the
        # density jump across the bow shock from ringing (the Euler-θ path
        # drops it, following Marras eq. 10). It is load-bearing here —
        # setting it to 0.0 is the single most destabilising change
        # measured on this case.
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
        # cutting Δt alone does nothing. If this case still fails downstream
        # of the corner for you, the next lever is NOT more :μ — it is the
        # ref = 2 grid, :nop => 3, or a Woodward & Colella corner entropy
        # fix. Raising :μ further without lowering Δt will only blow the
        # viscous limit sooner.
        :μ                    => [1.0, 4.0, 4.0, 4.0],
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
        :dsgs_norms           => "domain",
        #---------------------------------------------------------------------------
        # Mesh
        #
        # ffs_step_transfinite.msh is a three-block transfinite quad mesh of
        # the fluid L-shape, h = 0.025 m uniform (4032 elements). Its physical
        # curve groups — "inflow", "outflow", "wall" — are the tags that
        # reach user_bc_dirichlet!. Regenerate at twice the resolution by
        # setting ref = 2 in ffs_step_transfinite.geo.
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        #:gmsh_filename        => "./meshes/gmsh_grids/ffs_step_transfinite.msh",
        :gmsh_filename        => "./problems/CompEuler/ffs_step/ffs_step_transfinite.msh",
        #---------------------------------------------------------------------------
        # Plotting
        #---------------------------------------------------------------------------
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
        # visible next to the strong bow shock.
        :lschlieren           => true,
        :schlieren_k          => 20.0,            # contrast; Hadjadj uses 10-100
        #---------------------------------------------------------------------------
        # AMR off: the mesh already resolves the shocks at h/nop = 1/80, and
        # DynSGS is what handles what is left under-resolved.
        #---------------------------------------------------------------------------
        :linitial_refine      => true,
        :init_refine_lvl      => 1,
        :ladapt               => false,
    ) #Dict

    return inputs

end
