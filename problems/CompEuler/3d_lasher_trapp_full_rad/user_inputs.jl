function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK54(),
        :Δt                   => 0.15,   # acoustic CFL ≈ 0.6 on 12x12x16 (min LGL spacing 87 m); scale with the mesh (0.5 was for 6x6x8)
        :tinit                => 0.0,
        :tend                 => 7200.0,   # 2-hour deep convection simulation
        :diagnostics_at_times => (1, 10, 15, 20, 30, 40, 50, 100, (300:300:7200)...),
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
        :μ                    => [0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],   # as BOMEX
        :lmoist               => true,   # activates SAM microphysics + RT flux injection
        :lprecip              => true,   # deep convection can produce ice and precipitation
        :energy_equation      => "energy",
        :lrichardson          => false,   # off for now
        :bdy_fluxes           => true,   # Gaussian surface heat flux (user_bc_neumann!)
        :lwall_model          => true,   # allocates the wall-flux arrays used by the surface flux
        :lt_noise_amp         => 0.5,    # m/s horizontal-wind perturbation below 1 km (paper: 0.5 coarse, 0.2 fine)
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
        # Sounding file (22 Jul 1995, 1545 UTC — Lasher-Trapp et al. 2001, Fig. 1)
        # Format: z(m)  P(Pa)  T(°C)  qv(g/kg)  u(m/s)  v(m/s)
        # u=v=0 per the paper's zero background wind state.
        #---------------------------------------------------------------------------
        :sounding_file        => "./data_files/LasherTrapp2001_CCOPE.dat",
        #---------------------------------------------------------------------------
        # Mesh parameters
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./meshes/gmsh_grids/hexa_LasherTrapp_full_rad-12x12x16.msh",
        #---------------------------------------------------------------------------
        # Filter parameters
        #---------------------------------------------------------------------------
        :lfilter             => true,
        :mu_x                => 0.05,
        :mu_y                => 0.05,
        :mu_z                => 0.05,
        :filter_type         => "erf",
        :filter_basis        => "legendre",   # applied once per time step (docs/modal_filter.md)
        #---------------------------------------------------------------------------
        # Output parameters
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :output_dir           => "./output/",
        :loutput_pert         => true,
        :lwrite_initial       => true,
        #---------------------------------------------------------------------------
        # AMR (off by default)
        #---------------------------------------------------------------------------
        :ladapt               => false,
        :lamr                 => false,
        #---------------------------------------------------------------------------
    ) #Dict

    return inputs
end
