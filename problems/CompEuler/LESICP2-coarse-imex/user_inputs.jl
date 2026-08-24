#==============================================================================
 LESICP2-coarse-imex -- LESICP2-coarse, run with the fully 3D IMEX integrator.

 Originally a copy of problems/CompEuler/LESICP2-coarse-hevi with one more
 scheme wired in.

 THE GUARANTEE THIS DECK MAKES is that its own THREE ARMS are identical to each
 other: across DBG_SCHEME only :ode_solver, :Δt and :output_dir change, so any
 difference in the answer or the wall clock between them is the time integrator
 and nothing else. That is checked, not assumed.

 IT NO LONGER MATCHES THE STANDALONE -hevi DECK, and deliberately so -- the two
 have since been edited apart:

     key                     -coarse   -coarse-hevi   -coarse-imex (here)
     :lfilter                false     true           false
     :first_zelement_size    12.0      10.0           10.0

 So do NOT compare a run of this deck against a run of LESICP2-coarse-hevi:
 that comparison now carries a filter difference as well as an integrator one.
 Compare WITHIN this deck instead -- DBG_SCHEME=hevi gives the same HEVI
 configuration under the same physics as DBG_SCHEME=imex, which is what the
 standalone -hevi case was for before this one existed.

 THREE SCHEMES, ONE SWITCH
 -------------------------
   DBG_SCHEME=imex      IMEX_ARK(:ARS343)       all acoustics implicit, 3D
   DBG_SCHEME=hevi      HEVI_ARK(:ARS232)       vertical acoustics implicit
   DBG_SCHEME=explicit  CarpenterKennedy2N54()  nothing implicit

 so the A/B/C is three commands and nothing else changes between them:

   for s in imex hevi explicit; do
     DBG_SCHEME=$s mpiexecjl -n 16 julia --project=. \
       -e 'using Jexpresso; Jexpresso.run_case("CompEuler","LESICP2-coarse-imex")'
   done

 WHY THIS MESH IS INTERESTING FOR BOTH SPLITS
 --------------------------------------------
 The setup report measures, on this grid:

   node spacing    xi 138.1 m   eta 138.1 m   zeta 6.907 m      (20:1)
   dominant term   vertical acoustics

 20:1 anisotropy is a great HEVI case -- HEVI removes exactly the vertical
 acoustic term and its gain is capped by that ratio, which is why the sister
 case runs at 0.16 s against the explicit deck's 0.01 s, a 16x step.

 IMEX3D removes the horizontal acoustic terms as well, so its step is capped
 by ADVECTION instead. From the same measured rates, the wedge stability limit
 (ark_wedge_dt_max, ARS343) is 1.355 s, and this deck ships 0.9 s -- 90x the
 explicit step and 5.6x HEVI's.

 BUT A BIGGER STEP IS NOT AUTOMATICALLY A FASTER RUN, AND THIS CASE IS CLOSE
 ---------------------------------------------------------------------------
 The 3D stage solve is iterative, and its cost grows LINEARLY with Δt:
 measured on this grid's spacings,

   Δt     CFL_h = γΔt·c/h_x    GMRES iterations
   0.50        0.55                 12
   1.00        1.10                 24
   1.50        1.64                 38
   2.00        2.19                 54

 so at Δt = 0.9 expect roughly 22 iterations per stage solve, three solves per
 step. Writing a = A_vert/T_rhs (see below), cost per SIMULATED second in
 rhs!-equivalents comes out at

   explicit   500
   HEVI       18.8 + 31.8a
   IMEX3D      4.4 +  297a

 Both implicit schemes beat explicit by 20-25x here. Between them the answer
 flips at a ~ 0.055: below that IMEX3D wins, above it HEVI does. That is close
 enough that it is worth MEASURING rather than arguing about, which is what
 this deck is for.

   A_vert  = cost of one application of HEVI's vertical acoustic operator.
             The profiler's "f_imp" row under a HEVI run.
   T_rhs   = cost of one full nonlinear rhs! evaluation. The "rhs!" row.
   a       = A_vert / T_rhs, read straight off those two numbers.

 To measure it, run the hevi arm with the profiler on -- the environment
 variable has to reach every RANK, so set it inside the -e string rather than
 only in the launching shell:

   mpiexecjl -n 16 julia --project=. -e '
       ENV["JEXPRESSO_HEVI_PROFILE"]="1";
       using Jexpresso; Jexpresso.run_case("CompEuler","LESICP2-coarse-imex")'

 The breakdown prints from rank 0 every 50 steps (JEXPRESSO_HEVI_PROFILE_ALLRANKS=1
 for all ranks). Under the imex arm the same block reports the 3D operator and
 the Krylov solve instead, and :imex_monitor => true adds a running average of
 the iteration count.

 See src/kernel/solvers/hevi/README_IMEX3D.md for the scheme, and
 docs/hevi/hevi.md for the HEVI split it is built on.
==============================================================================#
function user_inputs()

    # Implicit vertical diffusion is opt-in per run: DBG_VDIFF=1. Read into a
    # local here rather than inline because two keys below have to agree about
    # it -- switching it on also switches the linearisation to :PS, without
    # which the implicit operator would carry no diffusion at all.
    _vdiff = parse(Bool, get(ENV, "DBG_VDIFF", "false"))

    #---------------------------------------------------------------------------
    # WHICH SCHEME -- ONE switch, three values.
    #
    #   :imex      IMEX_ARK(:ARS343)       -> ./output_imex/
    #   :hevi      HEVI_ARK(:ARS232)       -> ./output_hevi/
    #   :explicit  CarpenterKennedy2N54()  -> ./output_explicit/
    #
    # One switch rather than commented-out lines, because a wall-clock
    # comparison is only worth something if everything else is identical
    # between the runs, and swapping comments by hand is how one of those
    # quietly ends up different.
    #
    # Override from the shell with DBG_SCHEME=hevi (etc.) so the three arms can
    # be run from a loop without editing the file at all.
    #---------------------------------------------------------------------------
    scheme = Symbol(get(ENV, "DBG_SCHEME", "imex"))
    scheme in (:imex, :hevi, :explicit) ||
        error("LESICP2-coarse-imex: DBG_SCHEME must be imex, hevi or explicit; got $scheme")
    lexplicit = scheme === :explicit

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
    # Each scheme at its own largest safe step -- the only comparison that means
    # anything, since every implicit scheme here costs more per step and pays
    # for it by taking a bigger one.
    #
    #   explicit  0.01   what LESICP2-coarse ships
    #   hevi      0.16   just under the 0.1687 the setup report recommends (70%
    #                    of the measured joint limit 0.2410)
    #   imex      0.9    just under the 0.948 the wedge limit recommends (70% of
    #                    ark_wedge_dt_max = 1.355 at the rates this mesh
    #                    measures: explicit half 2.1 1/s, implicit half 38 1/s)
    #
    # The imex step is NOT set to the stability limit on purpose. Krylov
    # iterations grow linearly with Δt (see the header), so past the point where
    # the four rhs! evaluations stop dominating a step, a bigger Δt buys nothing
    # and eventually costs. Sweep it with DBG_DT and watch s/step in the
    # heartbeat against the iteration count in the profile block.
    Δt_default = scheme === :imex ? "0.9" : (scheme === :hevi ? "0.16" : "0.01")
    Δt   = parse(Float64, get(ENV, "DBG_DT", Δt_default))
    tend = parse(Float64, get(ENV, "DBG_TEND", "1000.0"))

    inputs = Dict(

        #---------------------------------------------------------------------------
        # User define your inputs below: the order doesn't matter
        #---------------------------------------------------------------------------
        :ode_solver           => scheme === :imex ? IMEX_ARK(:ARS343)      :
                                 scheme === :hevi ? HEVI_ARK(:ARS232)      :
                                                    CarpenterKennedy2N54(),
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
        #---------------------------------------------------------------------------
        # IMEX3D options (ignored unless scheme = :imex).
        # Full list and defaults: src/kernel/solvers/hevi/README_IMEX3D.md
        #---------------------------------------------------------------------------
        # :imex_verify      setup self-check: the 3D operator against the vertical
        #                   one, f_imp(qe) == 0 exactly, one full Krylov solve with
        #                   its applied residual and iteration count, and the
        #                   cross-rank single-valuedness of the answer. It also
        #                   prints the wedge stability limit for this tableau at
        #                   this Δt on this mesh. Leave it on.
        :imex_verify          => parse(Bool, get(ENV, "DBG_VERIFY", "true")),
        # :imex_umax        the largest flow speed this run is expected to reach,
        #                   m/s. THE ONE INPUT THAT CANNOT BE MEASURED AT SETUP.
        #                   HEVI's explicit half is dominated by horizontal SOUND,
        #                   known at t = 0; IMEX3D's is dominated by ADVECTION,
        #                   which grows as the boundary layer develops. Estimating
        #                   it from the initial sounding would say any Δt is safe.
        #                   15 m/s covers the u10 sounding plus convective gusts;
        #                   raise it if the run gets windier and the stability
        #                   warning starts firing.
        :imex_umax            => 15.0,
        # :imex_rtol        stage-solve tolerance. 1e-8 is comfortably below the
        #                   third-order step's own truncation error; loosening to
        #                   1e-6 saves an iteration or two, 1e-4 starts costing
        #                   order.
        :imex_rtol            => parse(Float64, get(ENV, "DBG_RTOL", "1.0e-8")),
        # :imex_warm_start  START EACH STAGE SOLVE FROM THE PREVIOUS ONE'S ANSWER
        #                   rather than from zero. On by default, and the single
        #                   largest lever on this scheme's cost.
        #
        #                   Consecutive ARK stages differ by O(dt*f) while the
        #                   right-hand side is the WHOLE deviation u - qe, so the
        #                   guess arrives several orders down. This operator is
        #                   skew and GMRES on it converges LINEARLY -- iterations
        #                   go as log(residual/tol) -- so those orders come
        #                   straight off the count. Measured end to end through
        #                   the ARK stepper: 5.0 iterations/solve against 21.7
        #                   cold, with the two states agreeing to 1.3e-10.
        #
        #                   It cannot change what the solve returns: the tolerance
        #                   is unchanged and is measured against the true residual.
        #                   DBG_WARM=0 turns it off to measure the difference.
        :imex_warm_start      => parse(Bool, get(ENV, "DBG_WARM", "true")),
        # :imex_restart     GMRES(m). Below CFL_h ~ 3 this is NOT the knob it looks
        #                   like: the operator is skew, its spectrum is a short line
        #                   segment, and restarting costs ~4% rather than the
        #                   stagnation an elliptic problem shows -- 71 iterations at
        #                   m=20 against 68 at m=240. ABOVE CFL_h ~ 3 that reverses
        #                   and it becomes one of the most important keys here:
        #                   measured on an isotropic grid at CFL_h = 7.4, m=20 and
        #                   m=40 do not converge at all and m=80 is needed.
        #                   CFL_h = gamma*dt*c/h_x is the HORIZONTAL acoustic Courant
        #                   number -- the column preconditioner removes the vertical
        #                   acoustics, so the horizontal is all the Krylov iteration
        #                   sees. The setup report prints it and advises an m.
        :imex_restart         => parse(Int, get(ENV, "DBG_RESTART", "20")),
        :imex_maxiter         => parse(Int, get(ENV, "DBG_MAXITER", "200")),
        # :imex_precond     :column is HEVI's own banded column solve, which is
        #                   what makes the 3D solve affordable. :none is there to
        #                   MEASURE that -- run both and compare the iteration
        #                   counts the setup report prints.
        :imex_precond         => Symbol(get(ENV, "DBG_PRECOND", "column")),
        # :imex_lateral_walls :auto reads the deck's own boundary faces, skipping
        #                   periodic ones -- so a laterally periodic case gets no
        #                   lateral walls automatically, which is right.
        :imex_lateral_walls   => Symbol(get(ENV, "DBG_LATWALL", "auto")),
        :imex_wall_flux       => true,
        # :imex_linearization same :RS / :PS switch as HEVI's, and the same
        #                   argument applies: this is a 1000 s convective boundary
        #                   layer with surface heating, so the solution DOES leave
        #                   qe behind and :PS is worth testing. DBG_LIN=PS.
        :imex_linearization   => Symbol(get(ENV, "DBG_LIN", _vdiff ? "PS" : "RS")),
        :imex_update_freq     => parse(Int, get(ENV, "DBG_UPDFREQ", "5")),
        # :imex_monitor     a running average of the Krylov iteration count every
        #                   N solves. That count IS the cost of this scheme, and
        #                   it drifts as the flow develops, so watch it on a first
        #                   production run.
        :imex_monitor         => parse(Bool, get(ENV, "DBG_IMEXMON", "true")),
        :imex_monitor_every   => parse(Int, get(ENV, "DBG_IMEXMONEVERY", "200")),
        :lcfl_report          => true,
        # Per-step progress with a wall clock and an ETA to :tend, printed from
        # :lcfl_report_every  re-print that table every N steps. :lcfl_report shows it
        #                   ONCE, at startup, where nu_t is ~0 on a laminar sounding --
        #                   so its viscous row says diffusion never limits dt, and says
        #                   it truthfully, at t = 0. Watch that row cross the line as the
        #                   boundary layer spins up. Collective and cheap; 0 = off.
        :lcfl_report_every    => parse(Int, get(ENV, "DBG_CFL_EVERY", "0")),
        # rank 0 only and doing nothing collective -- so it cannot hang a large
        # job. Effectively free (one time_ns() per step, a printf every 100),
        # and it is the only thing that separates "still compiling" from "hung"
        # on a first run. Also the timing instrument: the s/step it reports is
        # measured over the interval since the PREVIOUS heartbeat, so it is a
        # steady-state rate with the JIT excluded, which total wall time is not.
        :lstep_heartbeat      => parse(Bool, get(ENV, "DBG_HEARTBEAT", "true")),
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
        # Byte-identical to LESICP2-coarse-hevi's, deliberately: this case exists
        # to be compared against that one, and a comparison across two meshes is
        # not a comparison.
        :gmsh_filename    => "./problems/CompEuler/LESICP2-coarse-imex/LESICP_8x2x60_6400mX1600mX5000m.msh",
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
        :lfilter             => false,
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
	:output_dir          => scheme === :imex ? "./output_imex/" :
                                scheme === :hevi ? "./output_hevi/" : "./output_explicit/",
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
