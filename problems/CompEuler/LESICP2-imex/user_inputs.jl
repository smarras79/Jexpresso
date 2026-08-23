#==============================================================================
 LESICP2-imex -- LESICP2 with all three time integrators behind one switch.

 A copy of problems/CompEuler/LESICP2. Physics, closure, wall model, filter,
 stretching, sounding, statistics and MESH SETTINGS are byte-identical to that
 case; the only differences are which integrator runs, the step it takes, the
 output directory, and the per-scheme option blocks below.

   DBG_SCHEME=imex      IMEX_ARK(:ARS343)       all acoustics implicit, 3D
   DBG_SCHEME=hevi      HEVI_ARK(:ARS232)       vertical acoustics implicit
   DBG_SCHEME=explicit  CarpenterKennedy2N54()  nothing implicit   (default)

   for s in explicit hevi imex; do
     DBG_SCHEME=$s mpiexecjl -n 16 julia --project=. \
       -e 'using Jexpresso; Jexpresso.run_case("CompEuler","LESICP2-imex")'
   done

 READ THIS BEFORE THE FIRST IMPLICIT RUN
 ---------------------------------------
 THE STEP SIZES BELOW ARE PLACEHOLDERS. Both implicit arms default to the
 explicit deck's own Δt, which is safe but throws away the entire point of the
 scheme. The right value is a property of YOUR grid, and the setup report
 measures it on the first run and prints it:

   HEVI:    "joint IMEX stability: ... Δt_max here = X  (recommended: 0.7X)"
   IMEX3D:  "neutral up to Δt = X s; this run is at N% of it"

 Take the recommended figure, put it in DBG_DT (or edit the defaults below),
 and rerun. With :hevi_verify / :imex_verify on -- they are -- HEVI REFUSES to
 start above its measured limit and IMEX3D warns, so an over-large step fails
 loudly with the right number in the message rather than diverging quietly.

 This case's grid is not the one the numbers in LESICP2-coarse-imex were
 measured on, so nothing is assumed about it here.

 WHAT DECIDES WHETHER EITHER SPLIT IS WORTH IT ON THIS GRID
 ----------------------------------------------------------
 :lcfl_report => true is on below. Read its "dominant term" line first:

   vertical acoustics    -> HEVI is for exactly this
   horizontal acoustics  -> HEVI cannot help; IMEX3D removes it
   SGS diffusion         -> neither helps; the vertical diffusion wants to be
                            implicit, which is not implemented
   advection             -> already at the floor

 and its "Δt gain available, per scheme" table, which states what each split
 would buy on this grid before any per-step cost.

 HEVI's gain is capped by the acoustic ANISOTROPY (h_x/h_z) and nothing else.
 IMEX3D's step is capped by advection instead, but its stage solve is
 iterative and costs ~30 x (γΔt·c/h_x) Krylov iterations, so a coarse
 HORIZONTAL mesh is what makes it cheap. See
 src/kernel/solvers/hevi/README_IMEX3D.md.

 See docs/hevi/hevi.md for the HEVI split, and README_IMEX3D.md for the 3D one.
==============================================================================#
function user_inputs()
    #---------------------------------------------------------------------------
    # WHICH SCHEME -- one switch, three values. Defaults to :explicit so this
    # case behaves exactly like LESICP2 until you ask for something else.
    #---------------------------------------------------------------------------
    scheme = Symbol(get(ENV, "DBG_SCHEME", "explicit"))
    scheme in (:imex, :hevi, :explicit) ||
        error("LESICP2-imex: DBG_SCHEME must be imex, hevi or explicit; got $scheme")

    # PLACEHOLDER STEPS -- see the header. Both implicit arms start at the
    # explicit deck's Δt, which is safe and slow; raise them from the setup
    # report's measured limit.
    Δt = parse(Float64, get(ENV, "DBG_DT", "0.02"))

    inputs = Dict(
        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => scheme === :imex ? IMEX_ARK(:ARS343)      :
                                 scheme === :hevi ? HEVI_ARK(:ARS232)      :
                                                    CarpenterKennedy2N54(),
        #---------------------------------------------------------------------------
        # Time integration. Full option lists and defaults:
        #   docs/hevi/hevi.md                          (HEVI)
        #   src/kernel/solvers/hevi/README_IMEX3D.md   (IMEX3D)
        #---------------------------------------------------------------------------
        # The stability table, and the "dominant term" line that says whether
        # either split can help on this grid at all. One RHS evaluation.
        :lcfl_report          => true,
        # Per-step progress with a wall clock and an ETA to :tend, from rank 0,
        # nothing collective. Free, and the only thing that distinguishes "still
        # compiling" from "hung" on a first run.
        :lstep_heartbeat      => parse(Bool, get(ENV, "DBG_HEARTBEAT", "true")),
        #--- HEVI (ignored unless DBG_SCHEME=hevi) ---------------------------------
        # :hevi_verify is also the joint-stability guard: it REFUSES to start
        # above the limit it measures on this mesh, and names the recommended
        # step in the error. Leave it on.
        :hevi_verify          => parse(Bool, get(ENV, "DBG_VERIFY", "true")),
        :hevi_linearization   => Symbol(get(ENV, "DBG_LIN", "RS")),
        :hevi_update_freq     => parse(Int, get(ENV, "DBG_UPDFREQ", "5")),
        :hevi_wall_flux       => true,
        #--- IMEX3D (ignored unless DBG_SCHEME=imex) -------------------------------
        :imex_verify          => parse(Bool, get(ENV, "DBG_VERIFY", "true")),
        # THE ONE INPUT THAT CANNOT BE MEASURED AT SETUP: the largest flow speed
        # this run is expected to reach, m/s. IMEX3D's explicit half is
        # ADVECTION, which is ~0 at t = 0 and grows as the boundary layer
        # develops, so an estimate from the initial state would call any Δt safe.
        # Set it from the physics of YOUR case.
        :imex_umax            => parse(Float64, get(ENV, "DBG_UMAX", "15.0")),
        :imex_rtol            => parse(Float64, get(ENV, "DBG_RTOL", "1.0e-8")),
        # NOT the knob it looks like: this operator is skew, its spectrum is a
        # line segment, and restarting costs ~4% rather than the stagnation an
        # elliptic problem shows. Leave it at 20.
        :imex_restart         => parse(Int, get(ENV, "DBG_RESTART", "20")),
        :imex_maxiter         => parse(Int, get(ENV, "DBG_MAXITER", "200")),
        # :none is there to MEASURE what the column preconditioner buys. On a
        # 20:1 mesh it is 25x in iterations, so do not run production with it.
        :imex_precond         => Symbol(get(ENV, "DBG_PRECOND", "column")),
        :imex_lateral_walls   => Symbol(get(ENV, "DBG_LATWALL", "auto")),
        :imex_wall_flux       => true,
        :imex_linearization   => Symbol(get(ENV, "DBG_LIN", "RS")),
        :imex_update_freq     => parse(Int, get(ENV, "DBG_UPDFREQ", "5")),
        # The Krylov iteration count IS the cost of this scheme and it drifts as
        # the flow develops. Watch it on a first production run.
        :imex_monitor         => parse(Bool, get(ENV, "DBG_IMEXMON", "true")),
        :imex_monitor_every   => parse(Int, get(ENV, "DBG_IMEXMONEVERY", "200")),
        :Δt                   => Δt,
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
        :user_heatflux        => 0.12,
	# MUST be true. With false the mesh is built through Gridap's
	# GmshDiscreteModel(parts, ...) branch instead of the rank-0 read +
	# _compute_xy_partition column split, and the solution injects energy
	# out of nothing: still air with every forcing term off reached
	# 196 m/s in 100 s, independent of mesh, dt, C_s and :lrichardson.
	:lxy_partition          => true,
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
	#:gmsh_filename    => "./problems/CompEuler/LESICP2-coarse/LESICP_8x2x60_6400mX1600mX5000m.msh",
        #:gmsh_filename    => "/scratch/smarras/smarras/large_meshes/LESICP_128x128x120_10240mX10240mX5000m.msh
	:gmsh_filename    => "/scratch/smarras/smarras/large_meshes/LESICP_16x16x125_10240mX10240mX5000m.msh",
		
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
        :mu_x                => 0.25,
        :mu_y                => 0.25,
	:mu_z                => 0.25,
        :filter_type         => "erf",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
	# One tree per scheme, so an A/B/C leaves all three solutions on disk
	# rather than overwriting each other.
	:output_dir          => "/scratch/smarras/smarras/output_new/LESICP2_16x16x125_10240mX10240mX5000m_" * String(scheme) * "/",
	#:output_dir          => "/scratch/smarras/smarras/output_new/LESICP2_128x128x120_10240mX10240mX5000m/,"
        #:output_dir          => "./output_new/coarse-LESICP2_16x4x120_10kmX10kmX5km/",
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
