function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(), #ORK256(),#SSPRK33(), #SSPRK33(), #SSPRK54(),
        :Δt                   => 0.05,
        :tinit                => 0.0,
        :tend                 => 10800.0,
	:lrestart             => false,
	#:restart_output_file_path => "",
	:restart_time         => 500,
	#:diagnostics_at_times => (11500.0:10.0:15000.0),
	:diagnostics_at_times => (0.0:10.0:10800.0),
	#:diagnostics_at_times => (0:10:100..., 1250:500:5000..., 5000:100:8500...,  9000:10:10800.0...),
	:lsource              => true,
	:lsponge              => true,
	:zsponge              => 2500.0,
        :sounding_file        =>"./data_files/input_sounding_teamx_u00_flat_noheader.dat",
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  =>"lgl",
        :nop                  => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lwall_model          => true,
        :ifirst_wall_node_index=> 5, # This must be between 2 <= :first_wall_node_index <= nop+1
        :bdy_fluxes           => true,
        :lvisc                => true, #false by default
        :visc_model           => SMAG(),
        #:visc_model           => AV(),
        :μ                    => [0.0, 1.0, 1.0, 1.0, 1.0], #horizontal viscosity constant for momentum
        #:μ                    => [0.0, 10, 10, 10, 10], #horizontal viscosity constant for momentum	
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
	#:lwarmup          => true,
        :lread_gmsh       => true, #If false, a 1D problem will be enforced
	:gmsh_filename    => "./meshes/gmsh_grids/LESICP_16x16x36.msh",
	
        # Warping:
        :lwarp => false,
        :mount_type => "LESICP",
        :h_mount => 1000.0,
        :a_mount => 10240.0,
	:z_transition_start => -1000.0,
	:z_transition_end => 2200.0,

        # Stretching factors:
        :lstretch => false,
        :stretch_factor => 1.15,
        :stretch_type => "fixed_first_twoblocks_strong", #strong means that the top is constrained
        :first_zelement_size => 10.0,
        :zlevel_transition => 2000.0,
        
        #---------------------------------------------------------------------------
        # Filter parameters
        #---------------------------------------------------------------------------
        :lfilter             => true,
        :mu_x                => 0.1,
        :mu_y                => 0.1,
	:mu_z                => 0.1,
        :filter_type         => "erf",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :output_dir          => "/scratch/smarras/smarras/output/LESICP2_64x16x36_10kmX5kmX3dot5km-coarse",
        #:output_dir          => "./output",
        :loverwrite_output   => true,  #this is only implemented for VTK for now
        :lwrite_initial      => true,
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
