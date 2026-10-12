function user_inputs()

    day  = 86400.0
    year = 360.0 * day

    inputs = Dict(
        #---------------------------------------------------------------------------
        # ReducedGravityDSGS: the wind-driven, reduced-gravity DoubleGyre case
        # from problems/ShallowWater/ReducedGravity, stabilized with the same
        # residual-based shallow-water DynSGS configuration used by
        # problems/ShallowWater/SoliWaveIslandDSGS.
        #
        # Conservative variables:
        #   q = [H, Hu, Hv]
        # where H is the active-layer thickness.
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK54(),
        :Δt                   => 600.0,
        :tinit                => 0.0,
        :tend                 => 6.0 * year,
        :diagnostics_at_times => (0:30*day:6.0*year),
        :lmin_h_diagnostic    => true,
        :min_h_diagnostic_interval => 1,
        :case                 => "double_gyre_reduced_gravity_dsgs",
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
        # Integration and quadrature properties: identical to ReducedGravity.
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # Stabilization: the DSGS_SW configuration from SoliWaveIslandDSGS.
        #
        # The model computes one residual-based kinematic viscosity per
        # element and applies it to [H, Hu, Hv]. The three entries in :μ are
        # dimensionless per-equation multipliers. The gravity-wave cap must
        # use the reduced gravity g'=0.015 m/s² from ReducedGravity/user_flux.jl,
        # rather than the 9.81 m/s² used by SoliWaveIslandDSGS.
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :visc_model           => DSGS_SW(),
        :dsgs_sensor          => "residual",
        :μ                    => [1.0, 1.0, 1.0],
        :dsgs_CR              => 1.0e-1,
        :dsgs_Cmax            => 0.5e-2,
        :dsgs_Cmin            => 0.05e-3,
        :dsgs_norms           => "domain",
        :ldsgs_nodal          => false,
        :dsgs_Cl              => 0.0,
        :dsgs_swe_g           => 0.015,
        :dsgs_swe_hmin        => 1.0e-2,
        #---------------------------------------------------------------------------
        # No additional modal filter: DynSGS supplies the stabilization.
        #---------------------------------------------------------------------------
        :lfilter              => false,
        #---------------------------------------------------------------------------
        # Reuse the DoubleGyre mesh from the ReducedGravity case.
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/ShallowWater/ReducedGravity/DoubleGyre.msh",
        #---------------------------------------------------------------------------
        # Output: retain the ReducedGravity settings and add the DynSGS field.
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :plot_dsgs            => true,
        :plot_dsgs_vars       => ["H"],
        :loverwrite_output    => false,
        :lwrite_initial       => true,
        :output_dir           => "./output",
        :loutput_pert         => true,
        :lrestart_vtk         => true,
        :restart_vtk_input_dir => "./output/ShallowWater/ReducedGravityDSGS/output-08Oct2026-123236",
        #---------------------------------------------------------------------------
    )

    return inputs
end
