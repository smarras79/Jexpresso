function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK54(),
        :Δt                   => 0.05,
        :tinit                => 0.0,
        :tend                 => 5400,#21600.0,  # 6-hour BOMEX simulation
        :diagnostics_at_times => (0.25, 1, 10, 50, 100, (300:300:21600)...),
        :lsource              => true,
        :SOL_VARS_TYPE        => PERT(),
        #---------------------------------------------------------------------------
        # BOMEX forcing (Siebesma et al. 2003); each term can be switched off for tests
        #---------------------------------------------------------------------------
        :bomex_wind           => true,   # initial wind + Coriolis/geostrophic forcing
        :bomex_subsidence     => true,   # large-scale subsidence on hl, qt
        :bomex_radiation      => true,   # prescribed cooling (turn off when RT is coupled)
        :bomex_qt_ls          => true,   # large-scale drying
        :bomex_sponge         => true,   # top buffer layer, z >= 2250 m
        :bomex_noise          => true,   # initial θl/qt noise below 1600 m
        :bdy_fluxes           => true,   # surface fluxes (user_bc_neumann!)
        :lwall_model          => true,   # allocates the wall-flux arrays used by the surface fluxes
        #---------------------------------------------------------------------------
        # Restart options
        #---------------------------------------------------------------------------
        :lrestart             => false,
        :restart_time         => 300.0,
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
        :μ                    => [0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        :lmoist               => true,   # activates SAM microphysics + RT flux injection
        :lprecip              => false,  # BOMEX is non-precipitating shallow cumulus
        :energy_equation      => "energy",
        #:lrichardson          => true,
        #---------------------------------------------------------------------------
        # Radiation coupling (disabled for dynamics-only verification)
        # Re-enable all entries below once dynamics are validated.
        #---------------------------------------------------------------------------
        #:RT_atmos_coupling    => true,
        #:RT_radiative_heating => true,
        #:radiation_time_step  => 5.0,
        #:extra_dimensions     => 2,
        #:adaptive_extra_meshes => true,
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
        :gmsh_filename        => "./meshes/gmsh_grids/hexa_BOMEX_full_rad-16x16x19.msh",
        #---------------------------------------------------------------------------
        # Filter parameters
        #---------------------------------------------------------------------------
        # Needed: the non-dissipative CG scheme has a grid-scale instability in the
        # inversion that Smagorinsky does not control at small amplitude. Applied once
        # per time step in the orthogonal Legendre basis (see docs/modal_filter.md).
        :lfilter             => true,
        :mu_x                => 0.05,
        :mu_y                => 0.05,
        :mu_z                => 0.05,
        :filter_type         => "erf",
        :filter_basis        => "legendre",
        #---------------------------------------------------------------------------
        # Output parameters
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :output_dir           => "/scratch/tissaoui/",
        :loutput_pert         => true,
        :lwrite_initial       => true,
	:loverwrite_output   => true,
        #---------------------------------------------------------------------------
        # AMR (off by default — enable for adaptive angular mesh)
        #---------------------------------------------------------------------------
        :ladapt               => true,
        :lamr                 => false,
        #---------------------------------------------------------------------------
    ) #Dict

    return inputs
end
