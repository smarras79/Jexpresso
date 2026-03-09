function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(), #ORK256(),#SSPRK33(), #SSPRK33(), #SSPRK54(),
        :Δt                   => 0.5,
        :tinit                => 0,
        :tend                 => 10000.0,
        # :tend                 => 43200.0,
	:lrestart             => false,
	#:restart_output_file_path => "",
	# :restart_time         => 500,
	# :statistics_time      => 100,
	# :diagnostics_at_times => (0:4:100),
	:diagnostics_at_times => (100:100:10000),
	# :diagnostics_at_times => (0:4:40..., 100:500:600..., 610:10:700...,  800:100:1000.0...),
        :lsource              => true,
        :lmoist               => true,
        :lprecip              => true,
        :SOL_VARS_TYPE        => TOTAL(),
        :LST                  => true,
	:lsponge              => true,
	:zsponge              => 19000.0,
        :sounding_file        =>"./data_files/GIGALES_GATE_IDEAL_sounding.dat",
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  =>"lgl",
        :nop                  => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lwall_model          => true,
        :ifirst_wall_node_index=> 2, # This must be between 2 <= :first_wall_node_index <= nop+1
        :bdy_fluxes           => true,
        :lvisc                => true, #false by default
        :visc_model           => SMAG(),
        :μ                    => [0.0, 1, 1, 1, 1, 1, 1], #horizontal viscosity constant for momentum
        # :visc_model           => AV(),
        # :μ           => [0.0, 100.0, 100.0, 100.0, 200.0, 200.0, 200.0], #horizontal viscosity constant for momentum
        :energy_equation      => "energy",
        # :lrichardson          => true,
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
	#:lwarmup          => true,
        :lread_gmsh       => true, #If false, a 1D problem will be enforced
        # :gmsh_filename_c    => "./meshes/gmsh_grids/LESICP_64x16x36_10kmX5kmX3dot5km.msh",
        #:gmsh_filename    => "./meshes/gmsh_grids/LESICP_32x16x18_10kmX5kmX3km.msh",
	#:gmsh_filename    => "./meshes/gmsh_grids/LESICP_64x32x36_10kmX5kmX3km.msh",
        :gmsh_filename    => "./meshes/gmsh_grids/hexa_TFI_giga_les_30kmx12kmx25km.msh",	
	# :gmsh_filename    => "./meshes/gmsh_grids/hexa_TFI_giga_les_128kmx128kmx25km_1600m.msh",
	# :gmsh_filename    => "./meshes/gmsh_grids/hexa_TFI_giga_les.msh",
	
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
        :mu_x                => 0.5,
        :mu_y                => 0.5,
	:mu_z                => 0.5,
        :filter_type         => "erf",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :output_dir          => "./output_gigales_energy_moist/",
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
        # :ladapt              => false,
        :lamr                 => true,
        #---------------------------------------------------------------------------
        # AMR parameters
        #---------------------------------------------------------------------------
        :amr_freq            => 100,
        :amr_max_level       => 1,
        :amr_start_time      => 3500.0
        #---------------------------------------------------------------------------
    ) #Dict
    #---------------------------------------------------------------------------
    # END User define your inputs below: the order doesn't matter
    #---------------------------------------------------------------------------

    return inputs
    
end
