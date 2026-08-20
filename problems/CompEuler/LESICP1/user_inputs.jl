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
	#:restart_time         => 10800,
	:diagnostics_at_times => (100, 9000:5:10800.0...),
        :lsource              => true,
        :sounding_file        => "./data_files/input_sounding_teamx_u00_flat_noheader.dat",
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
        :lwall_damping        => true,
        #:visc_model           => AV(),
        :ivisc_equations      => [1, 2, 3, 4, 5],
        # smagorinsky, cs = 0.23, input cs^2 for momentum cs^2/Pr for other equations, where Pr = 1/3
        #:μ                    => [0.0, 0.53, 0.53, 0.53, 1.6], #horizontal viscosity constant for momentum
        # :μ is a 0/1 MASK under a dynamic SGS model, not a viscosity: it
        # multiplies the eddy viscosity the closure already computed. The old
        # values ([0.0, 5, 5, 5, 5]) were AV constants and inflated C_s by sqrt(μ).
        # Tune the closure through :C_s instead.
        :μ                    => [0.0, 1.0, 1.0, 1.0, 1.0],
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
	# :lwarmup          => true,
        :lread_gmsh       => true, #If false, a 1D problem will be enforced
        :gmsh_filename_c  => "./meshes/gmsh_grids/scaling_32x32x32.msh",
        #:gmsh_filename    => "./meshes/gmsh_grids/LESICP_32x16x18_10kmX5kmX3km.msh",
	#:gmsh_filename    => "./meshes/gmsh_grids/LESICP_64x64x36_10kmX10kmX3km.msh",
	:gmsh_filename    => "./meshes/gmsh_grids/LESICP_64x32x36_10kmX5kmX3km.msh",
        #:gmsh_filename    => "./meshes/gmsh_grids/LESICP_80x40x10_10kmX1kmX3km.msh",

        # Warping:
        :lwarp => false,
        :mount_type => "agnesi",
        :a_mount => 4000.0,
        :h_mount => 1000.0,
        :c_mount => 5000.0,

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
        :output_dir          => "/scratch/smarras/smarras/output/LESICP1_scaling-8nodes-64x32x36_10kmX10kmX3km/",
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
