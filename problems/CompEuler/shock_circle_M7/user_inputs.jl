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
        # The VISCOUS limit is comfortable here and worth writing down, since
        # it is the one that bound ffs_step: the free-stream kinematic
        # viscosity is large at this density, ν = μ/ρ = 0.062 m²/s, but
        # νΔt/Δx² is only 0.006. The DynSGS coefficient, capped at
        # C_max = 0.1 below, adds about as much again.
        :Δt                   => 1.5e-8,
        :diagnostics_at_times => (0:1.0e-5:1.0e-3),
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
        :dsgs_sensor          => "residual",
        :dsgs_hold_steps      => 0,               # impulsive start: see ffs_step
        :μ                    => [1.0, 1.0, 1.0, 1.0],
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
