function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(), #ORK256(),#SSPRK33(), #SSPRK33(), #SSPRK54(),
        :Δt                   => 0.02,
        :tinit                => 0.0,
        :tend                 => 10000.0,
	:lrestart             => false,
	:restart_time         => 500,
	:diagnostics_at_times => (0:1:10..., 10:10:100..., 500:100:5000..., 5000:100:8500...,  9000:10:10800.0...),
	:lsource              => true,
	:lsponge              => true,
	:zsponge              => 2500.0,
        #:sounding_file        =>"./data_files/input_sounding_teamx_u00_flat_noheader.dat",
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  =>"lgl",
        :nop                  => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        #:user_heatflux          => 0.12,
        #:lwall_model            => true,
        :ifirst_wall_node_index => 5, # This must be between 2 <= :first_wall_node_index <= nop+1
        :bdy_fluxes             => false,
        :lvisc                  => true, #false by default
        #:visc_model             => SMAG(),
        :visc_model           => AV(),
        :μ                      => [0.0, 10, 10, 10, 10], #horizontal viscosity constant for momentum	
        #---------------------------------------------------------------------------
        #LES statistics
        #---------------------------------------------------------------------------
	:statistics_time      => (7200.0:10.0:10800.0),
        :lesprofile_vars      => ["u_mean", "v_mean", "w_mean", "t_mean", "p_mean"],
        :lesstress_vars       => ["upup_res", "upvp_res", "upwp_res", "vpvp_res", "vpwp_res", "wpwp_res",
                                   "tptp_res", "uptp_res", "vptp_res", "wptp_res",
                                   "upup_sfs", "upvp_sfs", "upwp_sfs", "vpvp_sfs", "vpwp_sfs", "wpwp_sfs",
                                   "tptp_sfs", "uptp_sfs", "vptp_sfs", "wptp_sfs",
                                   "uppp", "vppp", "wppp", "eps", "eps_t", "rho",
                                   "upupup", "upupvp", "upupwp",
                                   "vpvpup", "vpvpvp", "vpvpwp",
                                   "wpwpup", "wpwpvp", "wpwpwp",
                                   "upuptp", "vpvptp", "wpwptp"],
        :lesspectra_vars      => [],
        #---------------------------------------------------------------------------
        # Mesh paramters and files:
        #---------------------------------------------------------------------------
	#:lwarmup          => true,
        :lread_gmsh       => true, #If false, a 1D problem will be enforced
	#:gmsh_filename    => "./meshes/gmsh_grids/LESICP_16x16x36.msh",
        #:gmsh_filename    => "./meshes/gmsh_grids/LESICP_coarse_test.msh",
        :gmsh_filename    => "./meshes/gmsh_grids/Cuyo.msh",
        #:gmsh_filename    => "./meshes/gmsh_grids/Cuyo_MOST.msh",

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
        #:lfilter             => true,
        :mu_x                => 0.5,
        :mu_y                => 0.5,
	:mu_z                => 0.5,
        :filter_type         => "erf",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        #:output_dir          => "/scratch/smarras/smarras/output/LESICP2_64x16x36_10kmX5kmX3dot5km-coarse",
        :output_dir          => "./output",
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
