function user_inputs()
    
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(), #SSPRK54(), #ORK256(),#SSPRK33(), #SSPRK33(),
        :Δt                   => 0.5,
        :tinit                => 0.0,
        :tend                 => 1000.0,
        :diagnostics_at_times => (0:100:1000),
        :restart_time         => 500,
        :lrestart             => false,
        :restart_input_file_path => "/home/leon/njit/Jexpresso_gigales/Jexpresso/problems/equations/CompEuler/theta",
        :case                 => "rtb",
        :lsource              => true, 
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
        :visc_model           => AV(),
        :μ                   => [0.0, 125.0, 125.0, 125.0], #horizontal viscosity constant for momentum
        #:visc_model           => WALE(),
        #:visc_model           => VREM(),
        #:visc_model           => SMAG(),
        #:μ                   => [0.0, 1.0, 1.0, 1.0], #horizontal viscosity constant for momentum
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:lwarp => true,
        #:lstretch => true,
        :stretch_factor => 1.25,
        :mount_type => "agnesi",        
        :h_mount => 1000.0,
        :a_mount => 2000.0,
        :c_moint => 5000.0,
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_RTB20x20.msh", #for nop=4
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x10_periodic.msh", #for nop=4
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x10.msh", #for nop=4
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
        :lwrite_initial      => true,
        :output_dir          => "./output",
        #:output_dir          => "./test/CI-run",
        :loutput_pert        => true,  #this is only implemented for VTK for now
        #---------------------------------------------------------------------------
        # init_refinement
        #---------------------------------------------------------------------------
        :linitial_refine     => false,
        :init_refine_lvl     => 1,
        #---------------------------------------------------------------------------
        # AMR
        #---------------------------------------------------------------------------
        :ladapt              => false,
        #---------------------------------------------------------------------------
        # AMR parameters
        #---------------------------------------------------------------------------
        :amr_freq            => 200,
        :amr_max_level       => 2,
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
