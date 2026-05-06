function user_inputs()
  inputs = Dict(
      #---------------------------------------------------------------------------
      # User define your inputs below: the order doesn't matter
      #---------------------------------------------------------------------------
      :ode_solver           => CarpenterKennedy2N54(),
      :Δt                   => 0.2,
      :tinit                => 0.0,
      :tend                 => 36000.0,        # 1 hour — Jang & Hong Fig 1
      :diagnostics_at_times => (0:60:36000),
      #:case                 => "rtb",
      :lsource              => true,
      :lsponge              => true,
      :zsponge              => 15000.0,
      :SOL_VARS_TYPE        => PERT(),
      #---------------------------------------------------------------------------
      # Integration and quadrature properties
      #---------------------------------------------------------------------------
      :interpolation_nodes  => "lgl",
      :nop                  => 6,
      #---------------------------------------------------------------------------
      # Physical parameters/constants:
      #---------------------------------------------------------------------------
      :lvisc                => true,
      :ivisc_equations      => (2,3,4),
      :μ                    => [0.0, 50.0, 50.0, 0.0],
      #---------------------------------------------------------------------------
      # Mesh
      #---------------------------------------------------------------------------
      :lread_gmsh           => true,
      :gmsh_filename        => "./meshes/gmsh_grids/JH_periodic_6.msh",
      #---------------------------------------------------------------------------
      # Mountain parameters
      #---------------------------------------------------------------------------
      :lwarp                => true,
      :mount_type           => "agnesi",
      :a_mount              => 1000.0,        # half-width = 1km
      :h_mount              => 1000.0, #400.0,         # mountain height = 400m
      :c_mount              => 36000.0,       # mountain center at x=36km
      #---------------------------------------------------------------------------
      # Soundings and data files
      #---------------------------------------------------------------------------
      :sounding_file        => "./data_files/snwp_nonhydro.data",
      #---------------------------------------------------------------------------
      # Filter parameters
      #---------------------------------------------------------------------------
      :lfilter              => true,
      :mu_x                 => 0.25,
      :mu_y                 => 0.25,
      :filter_type          => "erf",
      #---------------------------------------------------------------------------
      # Plotting parameters
      #---------------------------------------------------------------------------
      :outformat            => "vtk",
     # :output_dir           => "./JaaaaH_10_0.25",
    # :output_dir           => "./JH_10_sponge=0.05",
  #  :output_dir           => "./JH_visc=10_400m",
  :output_dir           => "./JH_visc=50_1000m_nop6",
      :loverwrite_output    => true,
      :loutput_pert         => true,
      #---------------------------------------------------------------------------
  ) #Dict
# this is for ctop, and lateral coefficients = 0.1
  return inputs

end