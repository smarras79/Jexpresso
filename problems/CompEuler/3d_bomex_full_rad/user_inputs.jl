function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK54(),
        :Δt                   => 0.25,
        :tinit                => 0.0,
        :tend                 => 21600.0,  # 6-hour BOMEX simulation
        :diagnostics_at_times => (1, 10, 50, 100, (300:300:21600)...),
        :lsource              => true,
        :SOL_VARS_TYPE        => PERT(),
        #---------------------------------------------------------------------------
        # Restart options
        #---------------------------------------------------------------------------
        :lrestart             => false,
        :restart_time         => 0.0,
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # Physical parameters
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :visc_model           => SMAG(),
        :ivisc_equations      => [1, 2, 3, 4, 5, 6, 7],
        :μ                    => [0.0, 1.0, 1.0, 1.0, 3.0, 3.0, 3.0],
        :lmoist               => true,   # activates SAM microphysics + RT flux injection
        :lprecip              => false,  # BOMEX is non-precipitating shallow cumulus
        :energy_equation      => "energy",
        :lrichardson          => true,
        #---------------------------------------------------------------------------
        # Radiation coupling (disabled for dynamics-only verification)
        # Re-enable all entries below once dynamics are validated.
        #---------------------------------------------------------------------------
        #:RT_atmos_coupling    => true,
        #:RT_radiative_heating => true,
        #:radiation_time_step  => 300.0,
        #:extra_dimensions     => 2,
        #:adaptive_extra_meshes => false,
        #:RT_amr_threshold     => [0.99999],
        #:extra_dimensions_order   => 4,
        #:extra_dimensions_nelemx  => 4,
        #:extra_dimensions_nelemy  => 4,
        #:lcubed_sphere_angular_mesh => false,
        #:extra_dimensions_xmax => π,
        #:extra_dimensions_ymax => 2*π,
        #---------------------------------------------------------------------------
        # Mesh parameters
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./meshes/gmsh_grids/hexa_BOMEX_full_rad-8x8x6.msh",
        #---------------------------------------------------------------------------
        # Filter parameters
        #---------------------------------------------------------------------------
        #:lfilter             => true,
        #:mu_x                => 0.01,
        #:mu_y                => 0.01,
        #:filter_type         => "erf",
        #---------------------------------------------------------------------------
        # Output parameters
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :output_dir           => "./output/",
        :loutput_pert         => true,
        :lwrite_initial       => true,
        #---------------------------------------------------------------------------
        # AMR (off by default — enable for adaptive angular mesh)
        #---------------------------------------------------------------------------
        :ladapt               => false,
        :lamr                 => false,
        #---------------------------------------------------------------------------
    ) #Dict

    return inputs
end
