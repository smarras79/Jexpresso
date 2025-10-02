function user_inputs()
    inputs = Dict(
    #---------------------------------------------------------------------------
    # User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------
    :ode_solver           => SSPRK54(), #ORK256(),#SSPRK33(), #SSPRK33(), #SSPRK54(),
        :Δt                   => 0.2,
        :tinit                => 0.0,
        :tend                 => 1000.0,
        #:tinit                => 100.0,
        #:tend                 => 1000.0,
        :lrestart             => true,
        :restart_input_file_path => "./output/CompEuler/3d/output/restart",
        :diagnostics_at_times => (0:100:1000),
        :restart_time         => 150.0,
        :lsource              => true,
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lwall_model          => true,
        :lvisc                => true, #false by default
        :visc_model           => AV(), #VREM(), #SMAG(),
        # smagorinsky, cs = 0.23, input cs^2 for momentum cs^2/Pr for other equations, where Pr = 1/3
        #:μ                    => [0.1587, 0.0529, 0.0529, 0.0529, 0.1587],
        :μ                    => [0.0, 60.0, 60.0, 60.0, 60.0],
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        #:lwarmup             => true,
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_2x1x1.msh",
        # :gmsh_filename       => "./meshes/gmsh_grids/2x2x2.msh",
        # :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x10x10.msh",
        # :gmsh_filename        => "./meshes/gmsh_grids/LESICP_stretched.msh",
        #:gmsh_filename_c      => "./meshes/gmsh_grids/LESICP_stretched.msh",
        :gmsh_filename        => "./meshes/gmsh_grids/hexa_TFI_10x1x10.msh",
        #:gmsh_filename_c      => "./meshes/gmsh_grids/hexa_TFI_10x1x10.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_RTB_periodic3D.msh",
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
        :output_dir          => "./output/",
        :loverwrite_output   => true,
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
        :amr                 => true,
        #---------------------------------------------------------------------------
        # AMR parameters
        #---------------------------------------------------------------------------
        :amr_freq            => 20,
        :amr_max_level       => 1,
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------
    
    return inputs    
end
