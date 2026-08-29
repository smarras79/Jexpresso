function user_inputs()

    # Implicit vertical diffusion is opt-in per run: DBG_VDIFF=1. Read into a
    # local here rather than inline because two keys below have to agree about
    # it -- switching it on also switches the linearisation to :PS, without
    # which the implicit operator would carry no diffusion at all.
    _vdiff = parse(Bool, get(ENV, "DBG_VDIFF", "false"))

    #---------------------------------------------------------------------------
    # HEVI or fully explicit -- ONE switch.
    #
    #   lexplicit = false  ->  HEVI_ARK(:ARS232)        -> ./output_hevi/
    #   lexplicit = true   ->  CarpenterKennedy2N54()   -> ./output_explicit/
    #
    # One switch rather than two commented-out lines, because a wall-clock
    # comparison is only worth something if everything else is identical
    # between the two runs, and swapping comments by hand is how one of those
    # quietly ends up different.
    #---------------------------------------------------------------------------
    lexplicit = false

    # Step sizes -- each scheme at its own largest safe step.
    #
    # Timing them at a COMMON step would measure cost per step, which is not the
    # question: HEVI costs more per stage and pays for it by taking a bigger
    # step. Only wall-clock to a fixed physical time answers it.
    #
    # The explicit 0.01 is the value LESICP2-coarse ships. The HEVI value is set
    # from the joint IMEX limit the setup report measures on this mesh -- see
    # the DBG_DT override to sweep it, and the setup report's
    #   "Δt is NN% of that limit  (recommended: X.XXX s, i.e. 70%)"
    # line to see where a given choice sits.
    # MEASURED on this mesh by the setup report (1 rank):
    #
    #   node spacing      xi 138.1 m   eta 138.1 m   zeta 6.907 m   (20:1)
    #   acoustic limits   xi 0.3877 s                zeta 0.01994 s
    #   dominant term     vertical acoustics -- this is what HEVI is for
    #   joint IMEX        explicit half 5.1 1/s, implicit half 34.6 1/s
    #                     Delta t_max = 0.2410 s, recommended 0.1687 (70%)
    #
    # 0.16 is just under that recommendation, so there is margin for the rate
    # estimate, for the flow speeding up as the boundary layer develops, and for
    # the rank count. Against the explicit 0.01 the deck ships, that is 16x on
    # the step.
    Δt   = parse(Float64, get(ENV, "DBG_DT", lexplicit ? "0.01" : "0.16"))
    tend = parse(Float64, get(ENV, "DBG_TEND", "1000.0"))

    inputs = Dict(

        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => lexplicit ? CarpenterKennedy2N54() : HEVI_ARK(:ARS232),
        :Δt                   => Δt,
        :tinit                => 0.0,
        :tend                 => tend,
        #---------------------------------------------------------------------------
        # HEVI options (ignored when lexplicit = true)
        #
        # :hevi_verify        setup self-check AND the joint-stability guard that
        #                     refuses to start above the measured limit. Cheap;
        #                     leave it on. It is what stops a tableau/Δt pair
        #                     that grows a fraction of a percent per step --
        #                     slow enough to look like anything but a time
        #                     integrator problem.
        # :hevi_linearization :RS freezes the operator coefficients at the
        #                     reference state qe (default); :PS refreshes them
        #                     from the solution every :hevi_update_freq steps and
        #                     refactorises. On rtb_hevi :PS cost 12.6% per step
        #                     and bought nothing, because a small perturbation
        #                     never leaves qe behind. THIS case is a 1000 s
        #                     convective boundary layer with surface heating, so
        #                     it DOES leave qe behind -- :PS is the one worth
        #                     testing here. Switch with DBG_LIN=PS.
        # :hevi_wall_flux     zero the implicit vertical mass flux at floor/lid.
        # :lcfl_report        print the stability table at startup. Read the
        #                     "dominant term" line.
        #---------------------------------------------------------------------------
        # :implicit_vdiff  MAKE THE VERTICAL SGS DIFFUSION IMPLICIT AS WELL.
        #
        # HEVI and the 3D IMEX remove the acoustic limit. What is left is
        # whatever was second, and on a mesh refined in z that is the SGS
        # diffusion of the vertical derivative, nu/dz^2. It is invisible at
        # startup -- the CFL report is taken on a laminar sounding where nu_t
        # is ~0 -- and it arrives on its own schedule as the boundary layer
        # spins up, which is what a blow-up at a fixed MODEL time looks like.
        #
        # This puts d/dz(mu d/dz) on u, v, w and theta into the same column
        # operator that already carries the vertical acoustics, so it costs a
        # wider band and no new solve. The horizontal and cross terms stay
        # explicit; they are dz/dx of what is removed.
        #
        # It REQUIRES :PS under a dynamic closure -- the coefficient is read
        # from the closure, which returns its molecular floor at t = 0, so :RS
        # would freeze it there and the operator would carry no diffusion at
        # all. DBG_VDIFF=1 therefore switches the linearisation with it (and
        # DBG_LIN still overrides).
        #---------------------------------------------------------------------------
        :implicit_vdiff       => _vdiff,
        #---------------------------------------------------------------------------
        :hevi_verify          => parse(Bool, get(ENV, "DBG_VERIFY", "true")),
        :hevi_linearization   => Symbol(get(ENV, "DBG_LIN", _vdiff ? "PS" : "RS")),
        :hevi_update_freq     => parse(Int, get(ENV, "DBG_UPDFREQ", "5")),
        :hevi_wall_flux       => true,
        :lcfl_report          => true,
        # :lcfl_report_every  re-print that table every N steps. :lcfl_report shows it
        #                   ONCE, at startup, where nu_t is ~0 on a laminar sounding --
        #                   so its viscous row says diffusion never limits dt, and says
        #                   it truthfully, at t = 0. Watch that row cross the line as the
        #                   boundary layer spins up. Collective and cheap; 0 = off.
        :lcfl_report_every    => parse(Int, get(ENV, "DBG_CFL_EVERY", "0")),
	:lrestart             => false,
	#:lrestart_vtk	      => true,
	#:restart_output_file_path => "",
	:restart_time         => 9000.0,
	#:diagnostics_at_times => (11500.0:10.0:15000.0),
	#:diagnostics_at_times => (0.0:50.0:10800.0),
	:diagnostics_at_times => (0:100:tend),
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
	# Set EXPLICITLY here rather than left to the default, for two reasons.
	#
	# (1) HEVI wants it. The implicit operator is column-local, so a partition
	#     that keeps whole vertical columns on one rank makes the solve
	#     communication-free; the SFC partition splits columns and the solver
	#     then has to gather and scatter every stage. The setup report prints
	#     "split across ranks: N of M local column-instances" -- it should read
	#     0.0% here.
	#
	# (2) The original deck carried a comment saying this MUST be true because
	#     with false "the solution injects energy out of nothing: still air
	#     with every forcing term off reached 196 m/s in 100 s". That was the
	#     ghost-element double-counting in the GmshDiscreteModel(parts, ...)
	#     branch -- each rank kept Gridap's ghost cells and assemble_mpi! then
	#     summed them again, which put the lumped mass at 144% of the domain
	#     volume on a 4-rank theta run. Fixed since (the local model is
	#     restricted to its owned cells), so false is no longer poison. It is
	#     still the wrong choice for HEVI, for reason (1).
        :lxy_partition        => true,
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
        :lread_gmsh       => true, #If false, a 1D problem will be enforced
	#:gmsh_filename    => "./meshes/gmsh_grids/LESICP_16x16x36.msh",
        #:gmsh_filename    => "./meshes/gmsh_grids/LESICP_coarse_test.msh",
	#:gmsh_filename    => "./meshes/gmsh_grids/LESICP_128x128x125_10kmX10kmX5km.msh",
        :gmsh_filename    => "./problems/CompEuler/LESICP2-coarse-hevi/LESICP_8x2x60_6400mX1600mX5000m.msh",
        #:gmsh_filename    => "./meshes/gmsh_grids/LESICP_64x64x60_10kmX10kmX5km.msh",
		
        # Warping:
        #=:lwarp => true,
        :mount_type => "LESICP",
        :h_mount => 1000.0,
        :a_mount => 10240.0,
        :z_transition_start => 0.0,
        :z_transition_end => 3000.0,=#
        #=:lwarp => false,
        :mount_type => "LESICP",
        :h_mount => 1000.0,
        :a_mount => 10240.0,
	:z_transition_start => -1000.0,
	:z_transition_end => 2200.0,
        =#
        # Stretching factors:
        :lstretch => true,
        :stretch_factor => 1.15,
        :stretch_type => "two_block uniformish", #strong means that the top is constrained
        :first_zelement_size => 10.0,
        :zlevel_transition => 2000.0,
        
        #---------------------------------------------------------------------------
        # Filter parameters
        #---------------------------------------------------------------------------
        :lfilter             => true,
        :mu_x                => 0.25,
        :mu_y                => 0.25,
	:mu_z                => 0.25,
        :filter_type         => "erf",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        # Separate directories per scheme, so an explicit and a HEVI run of the
        # same case cannot overwrite each other's output and be compared by
        # accident. (iter_1 is the initial condition and is byte-identical
        # between the two by construction -- compare a LATER iterate.)
	:output_dir          => lexplicit ? "./output_explicit/" : "./output_hevi/",
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
        #:amr                 => true,
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
