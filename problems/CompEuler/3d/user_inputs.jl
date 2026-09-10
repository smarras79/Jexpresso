function user_inputs()
    inputs = Dict(
    #---------------------------------------------------------------------------
    # User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(), #ORK256(),#SSPRK33(), #SSPRK33(), #SSPRK54(),
        :Δt                   => 0.65,
        :tinit                => 0.0,
        :tend                 => 1000.0,
        :diagnostics_at_times => (100:100:1000),
        # :diagnostics_at_times => (5, 100:100:1000...),        
        :lsource              => true,
        #---------------------------------------------------------------------------
        # restart options
        #---------------------------------------------------------------------------
        # set restart_time to enable write restart files every [restart_time] seconds 
        :restart_time         => 10000000.0, 
        # the default restart output dir is $(your_output_dir)/restart but you can always specify
        # :restart_output_file_path => "./output/CompEuler/3d/output/restart",
        :lrestart             => false,
        # the default restart input dir is $(your_output_dir)/restart but you can always specify
        # :restart_input_file_path => "./output/CompEuler/3d/output/restart",
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes =>"lgl",
        :nop                 => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        #:lwall_model          => true,
        :lvisc                => true, #false by default
        :visc_model           => VREM(),
        #:visc_model           => DSGS(),  # residual-based DynSGS, 3D θ form (DSGS.md §3.1).
                                           # :μ are multipliers of its ν, not a viscosity: with
                                           # [1,1,1,1,1] this case runs to t = 1000 s (ν ≈ 250 Pa·s
                                           # in the developed flow, 1.3e4 at start-up on the
                                           # bubble's cone edge). :μ[5] = 2 doubles the θ slot and
                                           # crosses the explicit parabolic limit at start-up
                                           # (ν·Δt/Δx_min² = 0.77; the run says so and dies at the
                                           # first step) — halve :Δt for that.
        #:visc_model           => AV(),
        #:visc_model           => SMAG(),
        # smagorinsky, cs = 0.23, input cs^2 for momentum cs^2/Pr for other equations, where Pr = 1/3
        #:μ                    => [0.1587, 0.0529, 0.0529, 0.0529, 0.1587],
        #:μ                    => [0.0, 125.0, 125.0, 125.0, 125.0],
        :μ                    => [0.0, 1.0, 1.0, 1.0, 2.0],
        :energy_equation     => "theta",
        #:lrichardson => true,
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
        #:lwarmup             => true,
        # Warping:
        :lwarp => false,
        :mount_type => "agnesi",
        :a_mount => 4000.0,
        :h_mount => 1000.0,
        :c_mount => 5000.0,

        # Stretching factors:
        :lstretch => false,
        :stretch_factor => 1.5,
        :stretch_type => "fixed_first_twoblocks_strong", #strong means that the top is constrained
        :first_zelement_size => 250.0,
        :zlevel_transition => 5000.0,
        
        # GMSH files:
        :lread_gmsh          => true, #If false, a 1D problem will be enforced
        #:gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_2x1x1.msh",
        # :gmsh_filename       => "./meshes/gmsh_grids/2x2x2.msh",
        # :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x10x10.msh",
        :gmsh_filename       => "./meshes/gmsh_grids/hexa_TFI_10x1x10.msh",
        #:gmsh_filename       => "./meshes/gmsh_grids/unstructured_xz_20x2x20.msh",
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
        :lamr                 => false,
        #---------------------------------------------------------------------------
        # AMR parameters
        #---------------------------------------------------------------------------
        :amr_freq            => 100,
        :amr_max_level       => 1,
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------
    
    return inputs    
end
