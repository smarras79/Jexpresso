function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # SoliWaveIsland_amr: 2D non-linear shallow water equations WITH AMR.
        #
        # AMR-enabled analogue of problems/ShallowWater/SoliWaveIsland. The
        # physics (flux/source/bc/primitives) are identical; the only additions
        # are the adaptive-mesh-refinement switches below and the two flag
        # functions (user_get_adapt_flags! / user_get_preadapt_flags!) in
        # initialize.jl. See docs/amr_setup.md for the full AMR contract.
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
        #
        # AMR criterion: the mesh refines where the water-level (free-surface)
        # perturbation |η| = |H - He| exceeds a tolerance, i.e. it tracks the
        # solitary wave and its run-up/run-down on the island, and coarsens
        # back to the base mesh in the still water behind the front. See
        # user_get_adapt_flags! in initialize.jl.
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK54(),
        :Δt                   => 0.01,
        :tinit                => 0.0,
        :tend                 => 25.0,
        :diagnostics_at_times => (0:1.0:25.0),
        :case                 => "soliwave_island_amr",
        :lsource              => true,
        :SOL_VARS_TYPE        => TOTAL(),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # Physical parameters / artificial viscosity.
        # We use a small constant Laplacian viscosity on all three equations
        # to keep the CG solution stable around the wave front. This is a
        # placeholder for the Dyn-SGS scheme of Marras et al. (2018), which
        # would replace the constant μ with a residual-based coefficient.
        # Note: the continuity equation diffuses the depth perturbation
        # H - He (see user_primitives.jl), not the full cone-shaped depth,
        # so the lake at rest stays an exact equilibrium.
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :visc_model           => AV(),
        :μ                    => [0.05, 0.05, 0.05],
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
        # Mesh (coarse, level-0 base for AMR).
        # Generate with:
        #   gmsh -2 problems/ShallowWater/SoliWaveIsland_amr/SoliWaveIsland_amr.geo \
        #        -o meshes/gmsh_grids/SoliWaveIsland_amr.msh
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./meshes/gmsh_grids/SoliWaveIsland.msh",
        #---------------------------------------------------------------------------
        # Plotting / output. AMR writes a p4est-forest checkpoint alongside
        # each VTK diagnostic dump, so use VTK here (the tested AMR output /
        # restart path, same as CompEuler/theta_amr). Open simulation.pvd in
        # ParaView to see the adapting mesh.
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        :output_dir           => "./output",
        :loutput_pert         => false,
        #---------------------------------------------------------------------------
        # preadapt: criterion-driven static pre-refinement before t=0. Refines
        # the vertical band around the initial solitary-wave crest so the wave
        # is resolved from the very first timestep (before the first runtime AMR
        # pass at t = amr_freq*Δt). Uses geometry only (no solution field exists
        # yet) — see user_get_preadapt_flags! in initialize.jl and
        # docs/amr_setup.md section 3.
        #---------------------------------------------------------------------------
        :lpreadapt            => true,
        :preadapt_max_level   => 1,
        #---------------------------------------------------------------------------
        # AMR: refine/coarsen on the water-level perturbation every :amr_freq
        # base timesteps, up to :amr_max_level. See user_get_adapt_flags! in
        # initialize.jl. The runtime dt is auto-rescaled by the current max
        # refinement level (dt = Δt / 2^ad_lvl_max) to keep the CFL stable as
        # elements shrink — set :Δt for the level-0 resolution.
        #---------------------------------------------------------------------------
        :lamr                 => true,
        :amr_freq             => 20,     # AMR pass every 20*Δt = 0.2 s (wave moves ~0.35 m)
        :amr_max_level        => 2,      # local Δx down to 1/4 = 0.25 m on the wave
        #---------------------------------------------------------------------------
        # AMR restart (off by default). If enabled, resumes from the p4est
        # forest + VTK solution saved at :restart_vtk_iout (see
        # user_read_vtu_point_data! in user_primitives.jl and
        # docs/amr_setup.md section 4).
        #---------------------------------------------------------------------------
        :lrestart_amr         => false,
        # :restart_vtk_iout   => 10,
        #---------------------------------------------------------------------------
    ) #Dict

    return inputs

end
