function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(),
        :Δt                   => 1.25,
        :tinit                => 0.0,
        :tend                 => 10000.0,
        #:tinit                => 100.0,
        #:tend                 => 1000.0,
        #:lrestart             => true,
        #:restart_input_file_path => "./output/CompEuler/theta/output-19Nov2023-115126",
        :diagnostics_at_times => (0:10:10000),
        :case                 => "rtb",
        :lsource              => true, 
        :lmoist               => true,
        :lprecip              => true,
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
        :μ           => [0.0, 200.0, 200.0, 300.0, 300.0, 300.0], #horizontal viscosity constant for momentum
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_squall_line_coarse.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_squall_line_2D.msh",
        :gmsh_filename       => "./meshes/gmsh_grids/squall2d.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_squall_line_LinEtAl.msh",
        #---------------------------------------------------------------------------
        # Mountain parameters
        #---------------------------------------------------------------------------
        :lwarp               => false,
        :mount_type          => "agnesi",
        :a_mount             => 10000.0,
        :h_mount             => 100.0,
        :c_mount             => 0.0,
        :lsponge             => true,
        :zsponge             => 15000.0,
        #---------------------------------------------------------------------------
        # Soundings and data files
        #---------------------------------------------------------------------------
        :sounding_file       => "./data_files/test_sounding.data",
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
        :output_dir          => "./output",
        :loverwrite_output   => true,
        :loutput_pert        => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
