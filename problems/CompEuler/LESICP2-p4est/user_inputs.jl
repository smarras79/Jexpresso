function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(), #ORK256(),#SSPRK33(), #SSPRK33(), #SSPRK54(),
        :Δt                   => 0.02,
	# ------------------------------------------------------------------
	# HEVI (horizontally explicit / vertically implicit). Swap the two
	# lines above for the two below to take the vertical acoustic terms
	# implicitly; nothing else in this deck has to change. See
	# src/kernel/solvers/hevi/README.md.
	#
	# Worth about 2x on this mesh and no more: Δx=80 m against Δz=41.7 m
	# at nop=4 in both directions is an acoustic anisotropy of 1.92, and
	# HEVI removes exactly the vertical acoustic term. Run once with
	# :lcfl_report => true first -- if the report says SGS diffusion
	# rather than vertical sound is what binds Δt, HEVI will buy nothing
	# and implicit vertical diffusion is the fix instead.
	#
	#:ode_solver           => HEVI_ARK(:ARS343),
	#:Δt                   => 0.04,
	#:lcfl_report          => true,
        :tinit                => 0.0,
        :tend                 => 10800.0,
	:lrestart             => false,
	#:lrestart_vtk	      => true,
	#:restart_output_file_path => "",
	:restart_time         => 9000.0,
	#:diagnostics_at_times => (11500.0:10.0:15000.0),
	#:diagnostics_at_times => (0.0:50.0:10800.0),
	:diagnostics_at_times => (100:100:1000..., 1000:1000:9000.0...,10800),
	:lsource              => true,
	#:lsponge              => true,
	#:zsponge              => 2500.0, hard coded in user_source.jl
        :sounding_file        =>"./data_files/input_sounding_teamx_u10_flat_noheader.dat",
        #---------------------------------------------------------------------------
        #Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  =>"lgl",
        :nop                  => 4,      # Polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters/constants:
        #---------------------------------------------------------------------------
        # SURFACE LAYER -- what changed, and where the full story is. The MOST
        # and near-wall fixes apply to this deck; write-up in
        # LESICP2-64x64x60-imex/user_inputs.jl and docs/boundary_conditions.md
        # section 2.2.1.
        #   * :user_heatflux sets delta_hf = 1, so MOST's own w'theta' is
        #     discarded -- but it was still setting theta* -> L -> psi_m -> the
        #     DRAG from a free surface node nothing pins. CM_MOST! is now handed
        #     the imposed flux. THE THETA TERM IS UNCHANGED; ONLY THE DRAG MOVES.
        #   * obukhov_length gained its missing rho (|L| was ~20% small at sea
        #     level), and its near-neutral guard no longer reports a STABLE
        #     layer for an unstable flux.
        #   * zeta = z/L is bounded to [-5, 2] and the drag coefficient sees a
        #     0.1 m/s gustiness floor. Both are defaults; set :most_zeta_min /
        #     :most_zeta_max / :most_u_min here only to tune them. The negative
        #     bound is the CONVECTIVE side -- zeta's sign is the sign of L, not
        #     a height.
        #   * :lwall_damping now damps AT the wall. It used to return the
        #     UNDAMPED (C_s*Delta)^2 on every z = 0 node while damping their
        #     neighbours, so the eddy viscosity spiked on the node with the
        #     smallest h_z. mu_turb near the ground therefore changes.
        :user_heatflux        => 0.12,
	# MUST be true. With false the mesh is built through Gridap's
	# GmshDiscreteModel(parts, ...) branch instead of the rank-0 read +
	# _compute_xy_partition column split, and the solution injects energy
	# out of nothing: still air with every forcing term off reached
	# 196 m/s in 100 s, independent of mesh, dt, C_s and :lrichardson.
	# p4est path: the xy column partitioner is NOT used. mesh.jl only consults
	# :lxy_partition when :linitial_refine is false, so this must be false here
	# and the decomposition comes from p4est's space-filling curve instead.
	:lxy_partition          => false,
        :lwall_model          => true,
        :ifirst_wall_node_index=> 2, # This must be between 2 <= :first_wall_node_index <= nop+1
        :bdy_fluxes           => true,
        :lvisc                => true, #false by default
        :visc_model           => SMAG(),
        # Smagorinsky constant. ABL LES runs 0.13-0.18
        :C_s                  => 0.16,
        # Buoyancy correction on nu_t. Without it the full eddy diffusivity acts
        # across the capping inversion and smears it over a few hundred metres.
        :lrichardson          => true,
        # Near-wall limit l = min(C_s*Delta, kappa*z) on the mixing length.
        :lwall_damping        => true,
        #:visc_model           => AV(),
        #:μ                    => [0.0, 0.53, 0.53, 0.53, 1.6], #horizontal viscosity constant for momentum
        # :μ is a 0/1 MASK under a dynamic SGS model, not a viscosity: it
        # multiplies the eddy viscosity the closure already computed. The old
        # values ([0.0, 5, 5, 5, 5]) were AV constants and inflated C_s by sqrt(μ).
        # Tune the closure through :C_s instead.
        :μ                    => [0.0, 1.0, 1.0, 1.0, 2.0],
        #---------------------------------------------------------------------------
        #LES statistics
        #---------------------------------------------------------------------------
	:statistics_time      => (9000.0:10:10800.0),
	#:statistics_time      => (10.0:10.0:100),
        #:statistics_online_start    => 9000.0,
	#:statistics_online_interval => 0.2,
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
        :lread_gmsh       => true, #If false, a 1D problem will be enforce
	# COARSE mesh (16x16x15). :init_refine_lvl => 1 refines it to the 32x32x30
	# grid the run actually uses. Must be the COLUMN-MAJOR file produced by
	# tools/reorder_msh_columns.jl -- p4est partitions contiguous ranges of the
	# mesh's element order, so a layer-ordered file gives horizontal slabs and
	# leaves most ranks with no ground surface at all.
	:gmsh_filename    => "./problems/CompEuler/LESICP2-p4est/LESICP_coarse_16x16x15_cols.msh",
		
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
        :mu_x                => 0.25,
        :mu_y                => 0.25,
	:mu_z                => 0.25,
        :filter_type         => "erf",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
	:output_dir          => "/scratch/smarras/smarras/output_new/LESICP2_128x128x120_10240mX10240mX5000m/",
        #:output_dir          => "./output_new/coarse-LESICP2_16x4x120_10kmX10kmX5km/",
        :loverwrite_output   => true,  #this is only implemented for VTK for now
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
