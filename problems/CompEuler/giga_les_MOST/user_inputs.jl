function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(), #ORK256(),#SSPRK33(), #SSPRK33(), #SSPRK54(),
        :Δt                   => 0.4,
        :tinit                => 0,
        :tend                 => 10800.0,
	:lrestart             => false,
	#:restart_output_file_path => "",
	:restart_time         => 500,
	:diagnostics_at_times => (0:1:50..., 100:500:5000..., 5000:250:8500...,  9000:10:10800.0...),
        :lsource              => true,
        :lmoist               => true,
        :lprecip              => true,
        :SOL_VARS_TYPE        => PERT(),
        :LST                  => true,
        # :LST                  => false,
	:lsponge              => true,
	:zsponge              => 19000.0,
        :sounding_file        =>"./data_files/Cirrus.dat",
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  =>"lgl",
        :nop                  => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :lwall_model          => false,
        :ifirst_wall_node_index=> 5, # This must be between 2 <= :first_wall_node_index <= nop+1
        :bdy_fluxes           => false,
        :lvisc                => true, #false by default
        # :visc_model           => SMAG(),
        :ivisc_equations      => [1, 2, 3, 4, 5, 6, 7],
        :visc_model           => AV(),
        :μ           => [0.0, 200.0, 200.0, 200.0, 300.0, 300.0, 300.0], #horizontal viscosity constant for momentum
        # :μ                    => [0.0, 1, 1, 1, 2, 2, 2], #horizontal viscosity constant for momentum
        :energy_equation      => "total_energy",
        :lrichardson          => true,
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
	#:lwarmup          => true,
        :lread_gmsh       => true, #If false, a 1D problem will be enforced
        # :gmsh_filename_c    => "./meshes/gmsh_grids/LESICP_64x16x36_10kmX5kmX3dot5km.msh",
        #:gmsh_filename    => "./meshes/gmsh_grids/LESICP_32x16x18_10kmX5kmX3km.msh",
	#:gmsh_filename    => "./meshes/gmsh_grids/LESICP_64x32x36_10kmX5kmX3km.msh",
	:gmsh_filename    => "./meshes/gmsh_grids/hexa_TFI_giga_les.msh",
	
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
        :lfilter             => false,
        :mu_x                => 0.5,
        :mu_y                => 0.5,
	:mu_z                => 0.5,
        :filter_type         => "erf",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :output_dir          => "./output_gigales_test/",
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
        :amr                 => false,
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
