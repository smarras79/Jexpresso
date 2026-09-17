function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        #
        # 2D CompEuler: Mach-7 LAMINAR flow over a circular cylinder.
        #
        # CompEuler/shock_circle raised to Mach 7 and made viscous, with the
        # grid clustered at the wall so the boundary layer and the surface
        # heat flux are resolved. It is the curved-wall rung on the way to
        # rampCaoEtAl2021, and it is the first deck in this series written
        # around what the ffs_step_M7 experiments actually showed.
        #
        # WHAT THOSE EXPERIMENTS SHOWED, because it is why this deck is
        # configured the way it is:
        #
        #   * filleting the step corner changed NOTHING (same death time to
        #     the digit), so a geometric singularity was not the cause;
        #   * a COMPLETE removal of the top mode by the modal filter bought
        #     1.47x, so the defect does not live in the top mode either;
        #   * the finer grid died EARLIER, and the element-scale checkerboard
        #     was present in the UNDISTURBED FREE STREAM upstream of the bow
        #     shock, where no physical feature exists and no shock-capturing
        #     sensor can legitimately fire.
        #
        # What is left is the discretization itself. A collocation CG
        # integrates the nonlinear flux with the same LGL rule it interpolates
        # on, so the flux is ALIASED: energy is transferred into the
        # grid-scale modes and CG has nothing to take it back out. That is
        # true at any Mach number; what Mach 7 changes is the PRICE, because
        #
        #     p = (γ-1)(ρE - ½ρ|u|²)
        #
        # is a difference of nearly equal numbers and a relative error in ρE
        # or ρu leaves that subtraction amplified by
        #
        #     (γ-1)·ρE/p = 1 + γ(γ-1)M²/2
        #
        # which is 3.5 at Mach 3 and 14.65 at Mach 7 — it grows like M². The
        # same aliasing error buys four times the pressure error, the pressure
        # drives the momentum flux, and the loop closes everywhere at once.
        #
        # THE STABILISING CHOICE IS THEREFORE :lkep, NOT MORE VISCOSITY, and
        # it is set in this deck rather than left to an environment variable
        # because it is the deck's answer, not a sweep.
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(),
        :tinit                => 0.0,
        # The bow shock and the stagnation boundary layer establish on the
        # body timescale R/u∞ = 1.28e-4 s. 1.0e-3 s is about eight of those —
        # a settled shock layer and a settled wall heat flux — and about half
        # a box flow-through. 66700 steps at the Δt below.
        :tend                 => 1.0e-3,
        :lrestart             => false,
        :restart_time         => 0.0,
        # The grid's smallest cell is 2.23 mm (at the wall), so at :nop => 4
        # the tightest LGL node gap is 3.84e-4 m. Against max(|u|+c), which is
        # ~1792 m/s in the free stream and cannot exceed √(2h₀) + c₀ ≈ 2380
        # anywhere, Δt = 1.5e-8 s is an advective CFL of about 0.08.
        #
        # DO NOT TRUST THE PRINTED "Viscous CFL" ON THIS MESH. computeCFL
        # forms it as max(ν) over the whole mesh times Δt over min(Δx)² over
        # the whole mesh (soundSpeed.jl), and on a graded grid those two are
        # at OPPOSITE ENDS: max ν is in the 77 mm far-field cells, where
        # μ_cap ∝ Δ is largest, and min Δx is in the 2.2 mm wall cells. The
        # first run printed 0.312 that way, while the true per-cell parabolic
        # number ν_i Δt/Δx_i² is at most 8.1e-3 anywhere in the domain —
        # a factor of 38 of pure diagnostic artifact. The number is correct
        # on the near-uniform meshes it was written against (ffs_step) and
        # meaningless here.
        #
        # The real viscous margins, per cell, at this Δt: 8.1e-3 from DynSGS
        # in the wall cells, 2.4e-4 in the far field, and 0.032 from the
        # MOLECULAR viscosity at the wall, where the low density and the
        # 300 K wall give ν = μ/ρ = 0.32 m²/s. All comfortable.
        # 3.75e-9: halved AGAIN, and the only variable that has ever moved
        # this case. The scaling is the reason, and it is the one measurement
        # in this series that says the failure may be removable at all:
        #
        #   run 1  Δt 1.50e-8   887 steps   t_fail 1.33e-5
        #   run 3  Δt 7.50e-9  3969 steps   t_fail 2.98e-5
        #
        # Halving Δt multiplied the steps by 4.47 and the PHYSICAL survival
        # time by 2.24. Compare what the two null hypotheses predict:
        #
        #   a hard physical limit (the shock simply cannot form)  -> x1.00
        #   fixed damage per step (purely numerical, unbounded)   -> x0.50
        #   MEASURED                                              -> x2.24
        #
        # Better than both. The damage per unit PHYSICAL time FELL when Δt
        # fell, which is what a Δt-dependent instability looks like and what a
        # hard limit does not. So there may be a Δt at which this runs.
        #
        # This run is the test of that, and it is worth something either way:
        # if the trend holds, t_fail lands near 6.7e-5; if it saturates near
        # 3e-5, Δt is exhausted and the answer is structural, not a step size.
        #
        # Why Δt and not more dissipation: the two runs so far differ only in
        # how much dissipation they had, and the one with LESS died sooner
        # (398 steps against 887). So the failure responds to the numerics
        # rather than sitting at a fixed physical time, and Δt is the lever
        # that buys margin through a violent transient without touching the
        # boundary layer this case exists to resolve. Both runs died with the
        # shock layer only half formed — 9.4 mm and 20.9 mm of flow travel
        # against a 42 mm standoff — so the whole difficulty is the FORMATION
        # of the normal shock, not any developed state.
        :Δt                   => 3.75e-9,
        :diagnostics_at_times => (0:2.0e-6:1.0e-3),   # dense: the first µs is the hard part
        :lsource              => false,
        :SOL_VARS_TYPE        => TOTAL(),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # Physical parameters / constants
        #---------------------------------------------------------------------------
        :energy_equation      => "energy",        # slot 4 is ρE, not ρθ
        :lvisc                => true,
        #---------------------------------------------------------------------------
        # THE STABILISATION. Flux differencing with an entropy-conservative
        # two-point volume flux.
        #
        # With :lkep the inviscid RHS is not assembled from the pointwise flux
        # but from symmetric two-point volume fluxes
        # (_expansion_inviscid_KEP! in rhs.jl), fed by the user_fluxaux!
        # methods in user_flux.jl. That construction removes the
        # aliasing-driven transfer into the grid-scale modes by construction —
        # it is aimed at the exact term the header argues is killing these
        # runs, and unlike a filter or a viscosity floor it costs nothing in
        # the boundary layer, which matters because the boundary layer is what
        # this case exists to resolve.
        #
        # ranocha() is entropy conservative. kennedy_gruber() is
        # kinetic-energy preserving only; central_euler() is the plain central
        # flux written in flux-differencing form and is the honest control — if
        # this case behaves identically under central_euler(), the two-point
        # machinery is not what is helping.
        #---------------------------------------------------------------------------
        :lkep                 => true,
        :volume_flux          => ranocha(),
        #---------------------------------------------------------------------------
        # Molecular viscosity: Sutherland, standard air constants. μ(125 K) =
        # 8.656e-6 Pa·s, which with the free stream of initialize.jl gives
        # Re_D = 1.01e4 and a laminar boundary layer δ ~ 2.8 mm. This is REAL
        # viscosity, not a numerical device: without it there is no boundary
        # layer and no surface heat flux to be stable about.
        #---------------------------------------------------------------------------
        :lsutherland          => true,
        :sutherland_muref     => 1.716e-5,        # Pa·s
        :sutherland_Tref      => 273.15,          # K
        :sutherland_S         => 110.4,           # K
        :Pr_lam               => 0.71,            # molecular Prandtl number
        #---------------------------------------------------------------------------
        # Shock capturing, on top of the molecular viscosity.
        #
        # :dsgs_sensor => "residual" — the DEFAULT, and the right one here for
        # two independent reasons, both of which the other decks in this
        # series get wrong for their own historical reasons:
        #
        #   1. A cylinder's bow shock becomes STATIONARY. "legacy" is
        #      R ≈ |∂ₜq| (rhs.jl:205), a rate sensor: it fires on a moving
        #      front and goes quiet on a standing one, which is precisely the
        #      steady state this case is integrating towards. The element
        #      residual is O(1/h) at a discontinuity whether or not it moves.
        #   2. A no-slip isothermal wall is a strong Dirichlet constraint on
        #      three of the four slots. The residual path explicitly zeroes
        #      the residual at Dirichlet nodes (_dsgs_boundary_pairs! in
        #      rhs.jl) because otherwise the constraint force reads as
        #      under-resolution and pins ν at its cap along the whole wall —
        #      measured on the rising bubble, and it blew that run up. The
        #      legacy branch returns before that correction.
        #
        # :dsgs_Cmax => 0.1 rather than the default 0.5 is rampCaoEtAl2021's
        # setting and for its reason: the cap has to be low enough that the
        # residual viscosity cannot smear the laminar boundary layer it is
        # sitting next to. :μ = [1,1,1,1] is the parameter-free setting; the
        # x4 that ffs_step carries is that case's, measured on an inviscid
        # run with no boundary layer to protect.
        #
        # NO :dsgs_Cmin. A background floor is the obvious lever for the
        # free-stream checkerboard, but it is a viscosity applied everywhere
        # including inside the boundary layer, and it would corrupt the wall
        # heat flux this case exists to measure. :lkep is meant to make it
        # unnecessary. If the free stream still quilts, :dsgs_Cmin => 0.01 is
        # the first thing to try, and the heat flux must then be re-checked.
        #---------------------------------------------------------------------------
        :visc_model           => DSGS(),
        #---------------------------------------------------------------------------
        # :ldsgs_nodal IS OFF, AND MUST STAY OFF WHILE :dsgs_sensor IS
        # "residual". Tried as run 4 and it was a 17x regression — 231 steps
        # against 3969, and the failure went GLOBAL (reported nodes scattered
        # to the domain corner at (0,-1)) instead of staying on the stagnation
        # streamline. The reason is in the kernel's own header (SGS.jl:2156):
        #
        #   "the residual is the assembled (lumped-mass) nodal residual
        #    R_i = |BDF2(q)_i - M^-1_i rhs_i|"
        #
        # and DSGS.md §1.2 says what that quantity is worth: with a lumped LGL
        # mass matrix the assembled rate M^-1 RHS IS what the integrator
        # advances, so the difference is the time-integration error and
        # nothing else — "it vanishes on an under-resolved solution exactly as
        # on a resolved one". The element form uses the ELEMENT RHS precisely
        # to avoid that. So nodal + "residual" is a BLIND sensor: nu ~ 0
        # everywhere, no shock capturing at all, and a Mach-7 bow shock has
        # nothing holding it. 231 steps is what that looks like.
        #
        # This does NOT condemn the nodal form. It is the right cure for what
        # the mu_dsgs field shows — the element kernel's staircase, which its
        # header records as "one wiggle per element in the smooth plateau" on
        # the Brio-Wu tube, measured — and it gives a C0 nu with no jump in
        # the diffusive flux at element interfaces. It just needs a sensor
        # that is not the assembled residual. The combination to try is
        # :ldsgs_nodal => true WITH :dsgs_sensor => "legacy", which is what
        # the MHD decks that exercise this path actually run: legacy makes
        # R ~ |dq/dt|, imperfect but not identically zero.
        #---------------------------------------------------------------------------
        :ldsgs_nodal          => false,
        :dsgs_sensor          => "residual",
        #
        # STARTUP HOLD OFF, as ffs_step has it. I turned it on for one run on
        # the argument that this case's initial field is smooth BY
        # CONSTRUCTION (the 5 mm blend of initialize.jl) and is therefore the
        # very condition the hold was written for. The run died SOONER — 398
        # steps against 887 — so the argument was wrong, and it is worth
        # writing down why.
        #
        # What the hold protects against is a sensor misreading a smooth
        # FIELD. What kills this case is a violent first few STEPS, and those
        # are violent no matter how smooth the field is: a 1568 m/s stream is
        # standing on a no-slip wall at t = 0 and a normal shock has to form
        # in front of it. That is ffs_step's own argument — "holding ν at zero
        # there integrates the most violent steps of the run with no
        # dissipation at all" — and it applies here for the same reason, which
        # I missed because I was looking at the initial condition instead of
        # at the first steps.
        :dsgs_hold_steps      => 0,
        :μ                    => [1.0, 1.0, 1.0, 1.0],
        #
        # :dsgs_Cmax => 0.03, not the ramp's 0.1, because μ_cap ∝ Δ and THIS
        # mesh is graded 34x (2.2 mm at the wall, 77 mm far field). The cap is
        # meant to be a ceiling reached at a shock, and per-cell it comes out
        #
        #     h = 2.2 mm  ->  ν_cap = 0.080 m²/s      (0.024 at Cmax = 0.03)
        #     h = 20  mm  ->  ν_cap = 0.72            (0.22)
        #     h = 77  mm  ->  ν_cap = 2.76            (0.83)
        #
        # against a molecular ν of 0.062 in the free stream. 0.03 was tried
        # and REVERTED: it went in the same run as the startup hold, both
        # changes cut dissipation, and the run died sooner. Back at the ramp's
        # 0.1. The coarse cells are then allowed 44x the molecular viscosity,
        # which is ugly but is not what is killing this case — the failure is
        # a tight cluster on the stagnation streamline, nowhere near the
        # coarse far field. If the far-field cap does become the problem, the
        # answer is less grading in the mesh (lc_far 0.06 -> 0.03 in
        # cylinder_M7.geo): μ_cap ∝ Δ cannot be undone from a deck.
        :dsgs_Cmax            => 0.1,
        :Pr                   => 0.1,
        :dsgs_norms           => "domain",
        #---------------------------------------------------------------------------
        # Mesh
        #
        # cylinder_M7.msh: unstructured all-quad, 14296 elements, minSICN
        # 0.609, no inverted cells, edge lengths 0.00223-0.0770 m. The size
        # field clusters 2.2 mm cells on the cylinder — about six LGL nodes
        # across the 2.8 mm boundary layer at :nop => 4 — and relaxes to
        # 60-77 mm in the far field. Near-square quads throughout: a modest
        # cluster, not a stretched y+ = 1 mesh.
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/CompEuler/shock_circle_M7/cylinder_M7.msh",
        # NOT a columnar partition. :lxy_partition bins cells into a uniform
        # nx x ny grid by centroid (_compute_xy_partition, mesh.jl), which is
        # what a 1D-implicit column solver — IMEX, HEVI — needs and nothing
        # else does. A uniform geometric bin only balances a uniform mesh, and
        # this one is graded 30x from the 2.2 mm wall cells to the 77 mm far
        # field: measured max/ideal load is 5.68x on 16 ranks, 6.66x on 32 and
        # 7.55x on 64, where one rank would own 1687 cells against an ideal
        # 223 and every other rank waits for it. (The same measurement on the
        # uniform ffs_step_M7 mesh gives 1.19-1.27x, and the unstructured
        # ffs_step_M7_round one 1.27-1.36x, which is why it never showed up
        # there.) The GLOBAL DEFAULT IS STILL true, because the same flag also
        # selects the mesh-READ strategy — see the note on :lxy_partition in
        # mod_inputs.jl — so a deck on a graded mesh has to say so itself, as
        # rampCaoEtAl2021 does.
        :lxy_partition        => false,
        #---------------------------------------------------------------------------
        # CURVE THE CYLINDER. This is not optional on a curved wall.
        #
        # gmsh writes a LINEAR grid, so the circle arrives as 64 straight
        # segments whose endpoints merely happen to lie on it. Filling those
        # elements with LGL nodes would put every high-order node on a CHORD,
        # and the wall the solver sees would stay a polygon however large
        # :nop is. On a NO-SLIP wall that is fatal twice: the polygon corners
        # shed spurious vorticity into the boundary layer, and the
        # wall-normal direction the temperature gradient — the heat flux — is
        # taken along is wrong by O(h) at every node.
        #
        # exact_geometry.jl snaps the high-order nodes of the "cylinder" edges
        # onto the true circle and blends the element interiors (Kopriva,
        # J. Sci. Comput. 26(3):301-327, 2006, §3: the linear-blending
        # transfinite map, which stays in P^N and so preserves the discrete
        # metric identities and free-stream preservation exactly).
        #
        # Explicit centre and radius rather than the :circle shorthand, which
        # fits them from the linear vertices and is refused on a refined grid
        # where new vertices sit at chord midpoints.
        #
        # Watch the first run for "# SNAP HIGH-ORDER NODES ONTO EXACT
        # GEOMETRY" and the wall distance it reports; shock_circle gets 3e-16.
        #---------------------------------------------------------------------------
        :exact_geometry       => Dict("cylinder" => (:circle, 1.0, 0.0, 0.2)),
        #---------------------------------------------------------------------------
        # Plotting
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        :output_dir           => "./output",
        :loutput_pert         => false,
        # Numerical schlieren from ρ. The field to look at for the bow shock,
        # its standoff (0.212 R = 42 mm at Mach 7) and the wake.
        :lschlieren           => true,
        :schlieren_k          => 20.0,
        #---------------------------------------------------------------------------
        # No refinement: the size field already puts the resolution where it
        # is needed. :init_refine_lvl => 1 would quarter every cell, and the
        # explicit :exact_geometry above is what keeps the snap working if you
        # do turn it on — but halve Δt with it.
        #---------------------------------------------------------------------------
        :linitial_refine      => false,
        :init_refine_lvl      => 1,
        :ladapt               => false,
    ) #Dict

    return inputs

end
