function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(), #ORK256(),#SSPRK33(), #SSPRK33(), #SSPRK54(),
        :Δt                   => 0.04,
        :tinit                => 0.0,
        :tend                 => 10800.0,
	:lrestart             => false,
	#:restart_output_file_path => "",
	:restart_time         => 250,
	:diagnostics_at_times => (0:5:10800.0),
        :lsource              => true,
	:lsponge              => true,
	:zsponge              => 2500.0,
        :sounding_file        =>"./data_files/input_sounding_teamx_u10_ridge1000_noheader.dat",
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
        #:μ                    => [0.0, 0.53, 0.53, 0.53, 1.6], #horizontal viscosity constant for momentum
        :μ                    => [0.0, 10, 10, 10, 15], #horizontal viscosity constant for momentum
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
	#:lwarmup          => true,
        :lread_gmsh       => true, #If false, a 1D problem will be enforced
        :gmsh_filename_c  => "./meshes/gmsh_grids/LESICP_64x32x36_10kmX5kmX3km.msh",
        #:gmsh_filename    => "./meshes/gmsh_grids/LESICP_32x16x18_10kmX5kmX3km.msh",
	#:gmsh_filename    => "./meshes/gmsh_grids/LESICP_64x64x36_10kmX10kmX3km.msh",
	#:gmsh_filename    => "./meshes/gmsh_grids/LESICP_64x32x36_10kmX5kmX3km.msh",
	:gmsh_filename    => "./meshes/gmsh_grids/LESICP_64x16x36_10kmX5kmX3dot5km.msh",

        # Warping:
        :lwarp => true,
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
        #:lfilter             => true,
        #:mu_x                => 0.01,
        #:mu_y                => 0.01,
        #:filter_type         => "erf",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :output_dir          => "/scratch/smarras/smarras/output/LESICP6_scaling-8nodes-64x16x36_10kmX10kmX3dot5km/",
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
