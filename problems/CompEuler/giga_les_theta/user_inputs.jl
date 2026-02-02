function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(), #ORK256(),#SSPRK33(), #SSPRK33(), #SSPRK54(),
        :Δt                   => 0.8,
        :tinit                => 0,
        # :tend                 => 100.0,
        :tend                 => 10800.0,
	:lrestart             => false,
	#:restart_output_file_path => "",
	:restart_time         => 500,
	# :diagnostics_at_times => (0:4:100),
	:diagnostics_at_times => (0:4:40..., 100:100:5000..., 5000:250:8500...,  9000:10:10800.0...),
	# :diagnostics_at_times => (0:4:40..., 100:500:600..., 610:10:700...,  800:100:1000.0...),
        :lsource              => true,
        :lmoist               => false,
        :lprecip              => false,
        :SOL_VARS_TYPE        => TOTAL(),
        :LST                  => false,
        # :LST                  => false,
	:lsponge              => true,
	:zsponge              => 17000.0,
        # :sounding_file        =>"./data_files/Cirrus_new.dat",
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
        :ifirst_wall_node_index=> 5, # This must be between 2 <= :first_wall_node_index <= nop+1
        :bdy_fluxes           => true,
        :lvisc                => true, #false by default
        :visc_model           => SMAG(),
        # :visc_model           => AV(),
        # :μ           => [0.0, 400.0, 400.0, 400.0, 800.0, 800.0, 800.0], #horizontal viscosity constant for momentum
        :μ                    => [0.0, 4, 4, 4, 4, 4, 4], #horizontal viscosity constant for momentum
        :energy_equation      => "theta",
        # :lrichardson          => true,
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
	#:lwarmup          => true,
        :lread_gmsh       => true, #If false, a 1D problem will be enforced
        # :gmsh_filename_c    => "./meshes/gmsh_grids/LESICP_64x16x36_10kmX5kmX3dot5km.msh",
        #:gmsh_filename    => "./meshes/gmsh_grids/LESICP_32x16x18_10kmX5kmX3km.msh",
	#:gmsh_filename    => "./meshes/gmsh_grids/LESICP_64x32x36_10kmX5kmX3km.msh",
	:gmsh_filename    => "./meshes/gmsh_grids/hexa_TFI_giga_les.msh",
	
  
        #---------------------------------------------------------------------------
        # Filter parameters -bigger weaker
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
        :output_dir          => "./output_gigales_energy_dry/",
        #:output_dir          => "./output",
        :loverwrite_output   => false,  #this is only implemented for VTK for now
        :lwrite_initial      => true,
        #---------------------------------------------------------------------------
        # init_refinement
        #---------------------------------------------------------------------------
        :linitial_refine     => true,
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
