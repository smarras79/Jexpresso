function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # SoliWaveIslandDSGS: 2D non-linear shallow water equations, the
        # SoliWaveIsland case stabilized by the residual-based DynSGS
        # (DSGS_SW, kernel/physics/SGS.jl) instead of a constant viscosity.
        #
        # Section 5.5 of Marras, Kopera, Constantinescu, Suckale, Giraldo,
        # "A residual-based shock capturing scheme for the continuous/
        #  discontinuous spectral element solution of the 2D shallow water
        #  equations", Advances in Water Resources 114 (2018) 45-63.
        #
        # A right-propagating solitary wave runs up and runs down a circular
        # conical island that sits in the middle of a closed rectangular basin.
        # Solved with CG SEM (no Riemann flux).
        #
        # Conservative variables:
        #   q = [H, Hu, Hv]
        # where H is the local water depth, Hb(x,y) is the (time-independent)
        # bathymetry, and the still-water depth far from the island is h0.
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK54(),
        #:ode_solver           => CarpenterKennedy2N54(),
        :Δt                   => 0.01,
        :tinit                => 0.0,
        :tend                 => 12.0,
        :diagnostics_at_times => (0:1.0:25.0),
        :case                 => "soliwave_island",
        :lsource              => true,
        :SOL_VARS_TYPE        => TOTAL(),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # Stabilization: DynSGS for the shallow-water system (Marras et al.
        # 2018, in the Marras-Nazarov form of the Euler/MHD kernels): one
        # residual-based kinematic ν per element,
        #
        #     ν = min( C_max Δ (|v| + √(gH)),  C_R Δ² max_i ‖R_i‖/‖δq_i − ⟨δq_i⟩‖ ),
        #
        # from the BDF2 residual of the three equations normalized by the
        # spread of the departure δq = q − qe from the lake at rest (the cone
        # would otherwise set the scale of H), applied as ∇·(ν∇(H − He)) on
        # the continuity equation and in stress form on (Hu, Hv)
        # (user_primitives.jl). :μ are the per-slot multipliers. The sibling
        # SoliWaveIsland uses a constant μ = 0.05 on the three equations
        # instead; here ν vanishes where the solution is resolved and rises
        # to the first-order cap at the wave front and the wet/dry ring.
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :visc_model           => DSGS_SW(),
        :dsgs_sensor          => "residual",  # element-wise strong residual (DSGS.md §1.2)
        :μ                    => [1.0, 1.0, 1.0],
        :dsgs_CR              => 1.0,       # C_R,   residual viscosity
        :dsgs_Cmax            => 0.5,       # C_max, first-order cap
        # Background floor: 5 % of the first-order viscosity, ν ≥ 0.05·Δ·(|v|+√(gH))
        # ≈ 0.02 m²/s in the still basin (the sibling's constant is 0.05). The
        # residual cannot see the element-scale transverse mode that the CG
        # solution of the solitary wave develops when ν vanishes on it
        # (measured: Hv ripples along the whole wave front by t = 2 s with
        # C_min = 0, caught by the residual only once grown; that run still
        # completes). With 0.05 the mode never appears, the run to t = 25 s
        # is clean and the residual part still rises to 0.035 at the island
        # and 0.14 at the run-up.
        :dsgs_Cmin            => 0.05,
        :dsgs_norms           => "domain",  # residual normalized by the spread over the whole basin
        :ldsgs_nodal          => false,     # true: ν at every node (Dao & Nazarov form), :dsgs_Cl its C_l
        :dsgs_Cl              => 0.0,
        :dsgs_swe_g           => 9.81,      # = _G_SWE of user_flux.jl
        :dsgs_swe_hmin        => 1.0e-3,    # = _H_WET_SWE of user_flux.jl (wet/dry threshold)
        #---------------------------------------------------------------------------
        # CG filter: OFF. The filter acts on the full depth H (it only
        # subtracts qe from the momentum components), and at rest H is
        # cone-shaped with a kink at the wet/dry ring: re-projecting it
        # every step perturbs the lake-at-rest equilibrium that the
        # well-balanced flux/source split preserves. The constant
        # artificial viscosity above is enough to stabilise the fronts.
        #---------------------------------------------------------------------------
        :lfilter              => false,
        #---------------------------------------------------------------------------
        # Mesh
        # Generate with:
        #   gmsh -2 problems/ShallowWater/SoliWaveIslandDSGS/SoliWaveIsland.geo \
        #        -o problems/ShallowWater/SoliWaveIslandDSGS/SoliWaveIsland.msh
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/ShallowWater/SoliWaveIslandDSGS/SoliWaveIsland.msh",
        #---------------------------------------------------------------------------
        # Plotting / output: one PNG per variable at every diagnostic time
        # (H-it<n>.png, Hu-it<n>.png, Hv-it<n>.png) plus the DynSGS
        # coefficient actually applied, μ_dsgs_H-it<n>.png (the three slots
        # carry the same ν, so one panel says it all). Set :lplot_surf3d to
        # true for the Spline2D surface rendering instead of the nodal map,
        # or switch to "vtk" for ParaView output (fields H, Hu, Hv, mu_dsgs).
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :plot_dsgs            => true,
        :plot_dsgs_vars       => ["H"],
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        :output_dir           => "./output",
        :loutput_pert         => false,
        #---------------------------------------------------------------------------
        # init_refinement
        #---------------------------------------------------------------------------
        :linitial_refine     => false,
        :init_refine_lvl     => 1,
    ) #Dict

    return inputs

end
