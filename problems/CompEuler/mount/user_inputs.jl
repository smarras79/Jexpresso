function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(),
        :Δt                   => 0.3,
        :tinit                => 0.0,
        :tend                 => 21600.0,
        #:tinit                => 100.0,
        #:tend                 => 1000.0,
        #:lrestart             => true,
        #:restart_input_file_path => "./output/CompEuler/theta/output-19Nov2023-115126",
        :diagnostics_at_times => (0:60:43200),
        :case                 => "rtb",
        :lsource              => true, 
        :lmoist               => true,
        :lprecip              => true,
        :lsponge              => true,
        :zsponge              => 15000.0,
        :SOL_VARS_TYPE        => PERT(),
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc       => true, #false by default NOTICE: works only for Inexact
        #:visc_model  => SMAG()
        #:μ           => [0.0, 20.0, 20.0, 30.0, 30.0, 30.0], #horizontal viscosity constant for momentum
        :μ           => [0.0, 400.0, 400.0, 400.0, 400.0, 400.0], #horizontal viscosity constant for momentum
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_squall_line_coarse.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_squall_line_2D.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/squall2d.msh",
        :gmsh_filename       => "./meshes/gmsh_grids/m400_4.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_squall_line_LinEtAl.msh",
        #---------------------------------------------------------------------------
        # Mountain parameters
        #---------------------------------------------------------------------------
        :lwarp               => true,
        :mount_type          => "agnesi",
        :a_mount             => 10000.0,
        :h_mount  => 2500.0,   # currently 100.0 — needs to be 2.5 km
        :c_mount  => 40000.0,  # currently 0.0 — mountain center at 40 km   
        #---------------------------------------------------------------------------
        # Soundings and data files
        #---------------------------------------------------------------------------
      #  :sounding_file       => "./data_files/test_hill.data",
        #:sounding_file       => "./data_files/input_sounding",
        :sounding_file       => "./data_files/sounding-SAM-new.dat", 
        #---------------------------------------------------------------------------
        # Filter parameters
        #---------------------------------------------------------------------------
        :lfilter             => true,
        :mu_x                => 0.1,
        :mu_y                => 0.1,
        :filter_type         => "erf", #use "erf" for Boyd-Vandeven, "exp" for exponential filter, or "quad" for quadratic filter
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk", #"hdf5",
        #:output_dir          => "./output_prime",
        :output_dir          => "./mount_400_noqv",
        :loverwrite_output   => true,
        :loutput_pert        => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end