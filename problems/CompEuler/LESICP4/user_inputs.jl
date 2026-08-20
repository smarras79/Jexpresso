function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(), #ORK256(),#SSPRK33(), #SSPRK33(), #SSPRK54(),
        :Δt                   => 0.04,
        :tinit                => 0,
        :tend                 => 10800.0,
	:lrestart             => false,
	#:restart_output_file_path => "",
	:restart_time         => 500,
	:diagnostics_at_times => (0:10:100..., 1250:250:5000..., 5000:100:8500...,  9000:10:10800.0...),
        :lsource              => true,
	:lsponge              => true,
	# TABLES specifies a 5000 m domain top with Rayleigh damping above 3000 m.
	# This mesh tops out at 3500 m and damps from 2500 m, so the sponge starts
	# 500 m lower than the protocol and occupies 29% of the column — gravity
	# waves radiated from the CBL top are absorbed closer in than in the other
	# models. Raising :zsponge alone would leave too thin a damping layer; this
	# needs a taller mesh (LESICP_*_10kmX10kmX5km.msh) and :zsponge => 3000.0.
	:zsponge              => 2500.0,
        # NOTE (ridge100 forced convection): this "_noheader" file is NOT a header-stripped copy of the
        # distributed TABLES sounding. All the u00 cases and u10_ridge1000 are
        # byte-identical to theirs; input_sounding_teamx_u10_ridge100.dat differs by up to
        # 2.2 K in theta and 4.8 m/s in u -- the protocol profile carries a
        # super-adiabatic surface layer, an Ekman spiral, and a ~200 m entrainment
        # zone, while this one is a uniform 300 K / u=10 idealisation down to z=0.
        # A correctly stripped copy is committed alongside; switch to it once the
        # intercomparison coordinators confirm which profile is the agreed IC:
        #:sounding_file        => "./data_files/input_sounding_teamx_u10_ridge100_tables_noheader.dat",
        :sounding_file        =>"./data_files/input_sounding_teamx_u10_ridge100_noheader.dat",
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  =>"lgl",
        :nop                  => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        :user_heatflux        => 0.12,
        :lwall_model          => true,
        :ifirst_wall_node_index=> 2, # This must be between 2 <= :first_wall_node_index <= nop+1
        :bdy_fluxes           => true,
        :lvisc                => true, #false by default
        :visc_model           => SMAG(),
        # Smagorinsky constant. Was PhysConst.C_s = 0.21; ABL LES runs 0.13-0.18
        # and nu_t goes as C_s^2, so 0.21 alone is ~1.7x Lilly.
        :C_s                  => 0.16,
        # Buoyancy correction on nu_t. Without it the full eddy diffusivity acts
        # across the capping inversion and smears it over a few hundred metres.
        :lrichardson          => true,
        # Near-wall limit l = min(C_s*Delta, kappa*z) on the mixing length.
        # false: :lwarp is on, so height above the domain floor is not the
        # distance to the wall over the ridge.
        :lwall_damping        => false,
        #:visc_model           => AV(),
        #:μ                    => [0.0, 0.53, 0.53, 0.53, 1.6], #horizontal viscosity constant for momentum
        # :μ is a 0/1 MASK under a dynamic SGS model, not a viscosity: it
        # multiplies the eddy viscosity the closure already computed. The old
        # values ([0.0, 5, 5, 5, 5]) were AV constants and inflated C_s by sqrt(μ).
        # Tune the closure through :C_s instead.
        :μ                    => [0.0, 1.0, 1.0, 1.0, 1.0],
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
	#:lwarmup          => true,
        :lread_gmsh       => true, #If false, a 1D problem will be enforced
        :gmsh_filename_c    => "./meshes/gmsh_grids/LESICP_64x16x36_10kmX5kmX3dot5km.msh",
        #:gmsh_filename    => "./meshes/gmsh_grids/LESICP_32x16x18_10kmX5kmX3km.msh",
	#:gmsh_filename    => "./meshes/gmsh_grids/LESICP_64x32x36_10kmX5kmX3km.msh",
	:gmsh_filename    => "./meshes/gmsh_grids/LESICP_64x64x36_10kmX10kmX3dot5km.msh",
	
        # Warping:
        :lwarp => true,
        :mount_type => "LESICP",
        :h_mount => 100.0,
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
        :output_dir          => "/scratch/smarras/smarras/output/LESICP4_scaling-8nodes-64x16x36_10kmX10kmX3dot5km/",
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
