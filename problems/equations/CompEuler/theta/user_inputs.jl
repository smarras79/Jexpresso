function user_inputs()
    
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK33(), #ORK256(),#SSPRK33(), #SSPRK33(), #SSPRK54(),
        #:Δt                   => 0.02,
        :Δt                   => 0.2,
        :tinit                => 0.0,
        :tend                 => 1000.0,
        #:tinit                => 100.0,
        #:tend                 => 1000.0,
        #:lrestart             => true,
        :restart_input_file_path => "./output/CompEuler/theta/output-19Nov2023-115126",
        :diagnostics_at_times => (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000),
        :case                 => "rtb",
        :lsource              => true, 
        #:backend              => MetalBackend(),
        #:SOL_VARS_TYPE        => PERT(), #TOTAL() is default
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lvisc                => true, #false by default NOTICE: works only for Inexact
        :ivisc_equations      => (1, 2, 3, 4),
        :μ                   => (0.0, 20.0, 20.0, 60.0), #horizontal viscosity constant for momentum
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/square_UNSTR_20el.msh", #for nop=4
        #:gmsh_filename       => "./meshes/gmsh_grids/plate_hole.msh", #for nop=4
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_RTB20x20.msh", #for nop=4
        #---------------------------------------------------------------------------
        # Filter parameters
        #---------------------------------------------------------------------------
        #:lfilter             => true,
        #:mu_x                => 0.01,
        #:mu_y                => 0.01,
        #:filter_type         => "erf",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :loverwrite_output   => true,
        :lwrite_initial      => false,
        :output_dir          => "./output",
        #:output_dir          => "./test/CI-run",
        :loutput_pert        => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
