function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(), #ORK256(),#SSPRK33(), #SSPRK33(), #SSPRK54(),
        :Δt                   => 0.01,
        :tinit                => 0.0,
        :tend                 => 43200.0,#86400.0,
        :diagnostics_at_times => 0.0:600.0:43200.0,
        #:tinit                => 100.0,
        #:tend                 => 1000.0,
        #:lrestart             => true,
        #:restart_input_file_path => "./output/CompEuler/theta/output-19Nov2023-115126",
       # :diagnostics_at_times => (0.2, 100.0, 200.0, 400.0, 500.0, 600.0, 900.0, 1000.0, 1200.0, 1300, 1400, 1500, 1800, 2000, 2500, 3000, 3500, 4000, 4250, 4500, 4750, 5000, 7000),
        :case                 => "mountain",
        :lsource              => true, 
        :lmoist               => true,
        :lprecip              => true,
        :SOL_VARS_TYPE        => TOTAL(),
        #:backend              => MetalBackend(),
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => true, #false by default NOTICE: works only for Inexact
        :ivisc_equations      => [ 2, 3, 4, 5, 6],
        :μ                   => [0.0, 200.0, 200.0, 200.0, 200.0, 200.0],
        #:μ                   => [0.0, 150.0, 150.0, 150.0, 150.0, 150.0], #horizontal viscosity constant for momentum
        #:μ                   => [0.0, 200.0, 200.0, 300.0, 300.0, 300.0], #horizontal viscosity constant for momentum
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_squall_line_coarse.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_squall_line_2D.msh",
       # :gmsh_filename       => "./meshes/gmsh_grids/fine_hill_6.msh",
        :gmsh_filename       => "./meshes/gmsh_grids/medium_h_4.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/squall2d.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_squall_line_LinEtAl.msh",
        #---------------------------------------------------------------------------
        # Mountain parameters
        #---------------------------------------------------------------------------
        :lwarp               => true,
        :mount_type          => "gauss",
        :h_mount => 1000.0,  # Change from 100.0 to 1000.0 (1km mountain)
        :a_mount => 10000.0,  # Keep this (10km half-width)
        :c_mount             => 0.0,
        #---------------------------------------------------------------------------
        # Soundings and data files
        #---------------------------------------------------------------------------
        #:sounding_file       => "./data_files/input_sounding_hill_pressure.txt",
        :sounding_file       => "./data_files/test_sounding.data",
        #:sounding_file       => "./data_files/a.data",
        #:sounding_file       => "./data_files/sounding-SAM-new.dat", 
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
        :output_dir          => "./sm_hill_4/",
        :loverwrite_output   => true,
        :loutput_pert        => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
