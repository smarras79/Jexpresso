function user_inputs()

    day = 86400.0 #seconds in a day
    year = 360.0 * day #seconds in a year (model year of 360 days)

    inputs = Dict(
        #---------------------------------------------------------------------------
        # ReducedGravity_EkmanWind: the ReducedGravity deck (Choi et al. 2004
        # wind-driven double gyre, reduced gravity, CG SEM) with ONE change, in
        # user_source.jl: the wind stress is handed to the active layer over an
        # Ekman-depth scale H_E = 20 m (smoothstep in H) instead of over
        # 1 mm < H < 0.1 m. Everything else -- mesh, Δt, AV, filter, clipping
        # limiter, thin-film relaxation -- is identical to ReducedGravity, so a
        # side-by-side run isolates the effect of the wind forcing on the thin
        # (outcropping) layers. See user_source.jl and README.md.
        #
        # Conservative variables:
        #   q = [H, Hu, Hv]
        # where H is the local water depth.
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK54(), #CarpenterKennedy2N54(),
        :Δt                   => 600.0, 
        :tinit                => 0.0,
        :tend                 => 6.0 * year, #6 years in seconds
        :diagnostics_at_times => (0: 30*day : 6.0*year),
        # Debugging diagnostic: check the MPI-global minimum layer thickness
        # after every accepted timestep and print only when it is negative.
        # Increase the interval to throttle output after H becomes negative.
        :lmin_h_diagnostic    => true,
        :min_h_diagnostic_interval => 1,
        :case                 => "double_gyre_reduced_gravity_ekmanwind",
        :lsource              => true,
        :SOL_VARS_TYPE        => TOTAL(),
        # Enforce H >= swe_hmin after every SSPRK stage and again after each
        # accepted step. The latter is also the integrator-independent path
        # used if this case is changed to an algorithm without stage hooks.
        :swe_limiter          => "clip",
        :swe_limiter_schedule => "stage",
        :swe_hmin             => 1.0e-3,
        # Accumulate every stage clipping event and print one MPI-global
        # summary after each accepted timestep that required correction.
        :swe_limiter_diagnostic => true,
        :swe_limiter_diagnostic_interval => 1,
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
        :μ                    => [0.0, 300.0, 300.0],
        #---------------------------------------------------------------------------
        # CG filter: OFF. The filter acts on the full depth H (it only
        # subtracts qe from the momentum components), and at rest H is
        # cone-shaped with a kink at the wet/dry ring: re-projecting it
        # every step perturbs the lake-at-rest equilibrium that the
        # well-balanced flux/source split preserves. The constant
        # artificial viscosity above is enough to stabilise the fronts.
        #---------------------------------------------------------------------------
        :lfilter     => true,
        :filter_type => "erf",
        :mu_x        => 0.05,
        :mu_y        => 0.05,
        #---------------------------------------------------------------------------
        # Mesh
        # Generate with:
        #   gmsh -2 problems/ShallowWater/ReducedGravity/DoubleGyre.geo \
        #        -o problems/ShallowWater/ReducedGravity/DoubleGyre.msh
        #---------------------------------------------------------------------------

        #NOTE - need to create that geo file and mesh for the double gyre case, and check if it is appropriate for the time step I want to use
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/ShallowWater/ReducedGravity/DoubleGyre.msh",
        #---------------------------------------------------------------------------
        # Plotting / output: one PNG per variable at every diagnostic time
        # (H-it<n>.png, Hu-it<n>.png, Hv-it<n>.png). Set :lplot_surf3d to
        # true for the Spline2D surface rendering instead of the nodal map,
        # or switch back to "vtk" for ParaView output.
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => false,
        :lwrite_initial       => true,
        :lrestart_vtk         => false,
        :restart_vtk_input_dir => "./output/ShallowWater/ReducedGravity_EkmanWind/output/",
        :output_dir           => "./output",
        :loutput_pert         => true,
        #---------------------------------------------------------------------------
    ) #Dict

    return inputs

end
