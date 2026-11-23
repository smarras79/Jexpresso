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
        :visc_model           => DSGS(),  # residual-based; cleaner than SMAG + modal filter (no spurious cloud aloft)
        :dsgs_CR              => 0.1,     # residual coefficient; the default 1 damped the plume out (12x12x16)
        :μ                    => [0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],   # as BOMEX
        :lmoist               => true,   # activates SAM microphysics + RT flux injection
        :lprecip              => true,   # deep convection can produce ice and precipitation
        :energy_equation      => "energy",
        :lrichardson          => false,   # off for now
        :Pr_t                 => 0.44,    # turbulent Prandtl number (paper), heat (hl)
        :Sc_t                 => 0.44,    # same for total water and precipitation, so hl and qt mix alike
        :bdy_fluxes           => true,   # Gaussian surface heat flux (user_bc_neumann!)
        :lwall_model          => true,   # allocates the wall-flux arrays used by the surface flux
        :lt_noise_amp         => 0.2,    # m/s horizontal-wind perturbation below 1 km (paper: 0.5 coarse, 0.2 fine)
        :lt_noise_type        => "smooth",   # "smooth" (non-divergent Fourier modes, λ >= :lt_noise_lmin) or "pointwise"
        :lt_noise_lmin        => 1000.0, # m, shortest wavelength of the smooth noise
        :lt_heat_flux_max     => 300.0,  # W m-2, peak of the Gaussian surface heat flux (0 = off)
        :lt_heat_sigma        => 1000.0, # m, Gaussian half-width (paper gives none; cloud base ~2 km wide)
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
        # Sounding: SCMS95 ISS4 (Loughman Lake, FL), 22 Jul 1995 15:45 UTC (x7221545.is4), the one used by
        # Lasher-Trapp et al. (2001). Extracted pixel-by-pixel from the NCAR EOL archive skew-T
        # (data_files/scms95_iss4_tools/). Below 107 m the heated-sensor temperature is replaced by the mixed
        # layer. v3: the lowest 50 mb have a uniform mixing ratio of 17.1 g/kg (archive 50-mb mean), blended back
        # to the measured profile over the next 30 mb; moisture column written as specific humidity.
        # Archive analysis check (plain T): LCL 918 hPa (916), LFC 845 hPa (850), CIN -16 J/kg (-11).
        # (v1: mixing ratio, surface moisture cut, 0.4 g/kg too dry; v2: measured surface moisture with mixed
        # theta, statically unstable lowest 100 m.) z is hypsometric.
        # Format: z(m)  P(Pa)  T(°C)  qv(g/kg, specific)  u(m/s)  v(m/s)
        # u=v=0 per the paper's zero background wind state.
        #---------------------------------------------------------------------------
        :sounding_file        => "./data_files/LasherTrapp2001_SCMS_ISS4_22Jul1995_1545_v3.dat",
        #---------------------------------------------------------------------------
        # Mesh parameters
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./meshes/gmsh_grids/hexa_LasherTrapp_full_rad-12x12x16.msh",
        #---------------------------------------------------------------------------
        # Filter parameters
        #---------------------------------------------------------------------------
        :lfilter             => false,   # DSGS needs no filter; with SMAG() use true and mu 0.003 below
        :mu_x                => 0.003,   # applied every step: damping rate ~ mu/Δt (0.05 suppressed convection)
        :mu_y                => 0.003,
        :mu_z                => 0.003,
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
