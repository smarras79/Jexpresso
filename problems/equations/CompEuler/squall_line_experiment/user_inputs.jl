function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK54(), #ORK256(),#SSPRK33(), #SSPRK33(), #SSPRK54(),
        :Δt                   => 0.2,
        :tinit                => 0.0,
        :tend                 => 7500.0,
        #:tinit                => 100.0,
        #:tend                 => 1000.0,
        #:lrestart             => true,
        #:restart_input_file_path => "./output/CompEuler/theta/output-19Nov2023-115126",
        :diagnostics_at_times => (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,1250, 1500,1750, 2000,2250, 2500,2750, 3000,3250, 3500,3750, 4000,4250, 5000, 6000, 6500, 7000, 7500),
        :case                 => "rtb",
        :lsource              => true, 
        :lmoist               => true,
        #:lprecip              => true,
        :SOL_VARS_TYPE        => PERT(),
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
        :ivisc_equations      => [1, 2, 3, 4, 5, 6, 7],
        :μ                   => [0.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0], #horizontal viscosity constant for momentum
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_squall_line.msh",
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_squall_line_coarse.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_squall_line_periodic_stretched.msh",
        #---------------------------------------------------------------------------
        # Mountain parameters
        #---------------------------------------------------------------------------
        #:lwarp               => true,
        :mount_type          => "agnesi",
        :a_mount             => 10000.0,
        :h_mount             => 100.0,
        :c_mount             => 0.0,
        #---------------------------------------------------------------------------
        # Soundings and data files
        #---------------------------------------------------------------------------
        :sounding_file       => "./data_files/test_sounding.data",
        #---------------------------------------------------------------------------
        # Filter parameters
        #---------------------------------------------------------------------------
        :lfilter             => false,
        :mu_x                => 0.05,
        :mu_y                => 0.05,
        :mu_z                => 0.05,
        :filter_type         => "erf", #use "erf" for Boyd-Vandeven, "exp" for exponential filter, or "quad" for quadratic filter
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk", #"hdf5",
        :output_dir          => "./output/",
        :loverwrite_output   => false,
        :loutput_pert        => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
