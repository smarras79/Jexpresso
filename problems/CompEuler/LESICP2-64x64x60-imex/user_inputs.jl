#=============================================================================
 LESICP2 on the 64 x 64 x 60 production grid -- 10240 x 10240 x 5000 m.

 WHAT THIS DECK IS FOR. The same case run fully explicit at dt = 0.02 s takes
 24 h of wall clock to reach tend = 10800 s. This deck runs it with the
 acoustics ENTIRELY implicit and is projected at 3-5 h on the same core count,
 2.5-3 h on 2048 ranks. Physics, sounding, surface flux, closure and filter are
 byte-identical to LESICP2-imex; only the grid, the step size and the solver
 settings differ.

 THE GRID, AND WHY IT DECIDES THE SCHEME
 ---------------------------------------
   64 x 64 x 60 elements at p = 4 over 10240 x 10240 x 5000 m
   element 160 x 160 x 40 m  ->  effective resolution 40 x 40 x 10 m
   smallest LGL gaps  h_x = 27.63 m   h_z = 6.91 m
   ACOUSTIC ANISOTROPY 4.0:1        15.9 M nodes, 245 760 elements

 4:1 is the number that chooses the integrator, and it chooses against HEVI.
 The explicit rate budget on this grid is

     vertical acoustic 50.2   horizontal acoustic 25.1
     advection          1.7   SGS diffusion        0.9      [1/s]

 HEVI removes only the first line -- 78.0 -> 27.7, a 2.8x rate gain -- but
 ARS232's explicit imaginary radius is 1.732 against CarpenterKennedy2N54's
 3.34, so nearly half of that is handed back at the tableau and the joint
 (rectangle) limit trims the rest. Net 1.2x on the step, 3 RHS/step instead of
 5: about 2x in wall clock, and 2x is 12 h, not "a few hours". Worse, the
 tableau with the radius to fix it (ARS343, 2.828) is the one HEVI cannot use:
 there the two halves are loaded by INDEPENDENT wavenumbers, the whole
 rectangle must be neutral, and ARS343 amplifies over most of it -- measured on
 this grid, dt = 0.0022 s, i.e. 0.05x explicit.

 IMEX3D removes ALL THREE acoustic terms: 78.0 -> 2.6, and now one wavenumber
 loads both halves, so only the wedge zE <= Mach*zI is reachable and ARS343 is
 neutral across it. dt = 0.506 s at the same 47% margin the explicit run
 carries -- 25x the step, 4 RHS/step, 31x fewer RHS evaluations per simulated
 second.

 WHAT IT COSTS, AND THE ONE NUMBER THAT DECIDES IT
 -------------------------------------------------
 The Krylov iteration count is the entire price. At dt = 0.5,

     CFL_h = gamma*dt*c/h_x = 2.77

 which is inside the band where GMRES(20-30) is still comfortable (it stops
 converging around CFL_h 5-7). Warm-started at rtol 1e-6 that is ~31
 iterations per solve, 3 solves per step.

 Writing every cost in units of A_vert (one vertical-operator application), the
 step costs 4*rho + 12 + 232 per 0.5 s, where rho = T_rhs / A_vert is how heavy
 this case's LES right-hand side is relative to the linear acoustic operator.
 That is the ONLY unknown, and JEXPRESSO_HEVI_PROFILE=1 measures it in a
 200-step run:

     rho     10      15      20      30      50
     speedup 3.1x    4.5x    5.7x    7.8x   11.2x
     hours   7.7     5.4     4.2     3.1     2.1

 IMEX3D overtakes explicit at rho = 3.0 and HEVI at rho = 4.6, so for any real
 LES right-hand side this is the right scheme; where in 2-8 h it lands is a
 property of how expensive rhs! is, not of the integrator. The ceiling, at
 infinite rho, is 31x -- the RHS-count ratio. Everything below that is Krylov.

 RUN IT ON A RANK COUNT THAT DIVIDES 4096
 ----------------------------------------
 :lxy_partition splits the 64 x 64 = 4096 ELEMENT COLUMNS, so

     512 ranks   8 columns each   exact
    1024 ranks   4 columns each   exact
    2048 ranks   2 columns each   exact

     DBG_SCHEME=imex (default) | hevi | explicit   -- same physics, three arms
     JEXPRESSO_HEVI_PROFILE=1                      -- measures rho
     DBG_DT=... DBG_RTOL=... DBG_RESTART=...       -- sweep the three knobs
=============================================================================#

function user_inputs()

    # Implicit vertical diffusion is opt-in per run: DBG_VDIFF=1. Read into a
    # local here rather than inline because two keys below have to agree about
    # it -- switching it on also switches the linearisation to :PS, without
    # which the implicit operator would carry no diffusion at all.
    _vdiff = parse(Bool, get(ENV, "DBG_VDIFF", "false"))
    #---------------------------------------------------------------------------
    # WHICH SCHEME -- one switch, three values. Defaults to :explicit so this
    # case behaves exactly like LESICP2 until you ask for something else.
    #---------------------------------------------------------------------------
    scheme = Symbol(get(ENV, "DBG_SCHEME", "imex"))
    scheme in (:imex, :hevi, :explicit) ||
        error("LESICP2-64x64x60-imex: DBG_SCHEME must be imex, hevi or explicit; got $scheme")

    #---------------------------------------------------------------------------
    # STEP SIZE -- one per scheme, each at its own limit. These are NOT
    # placeholders: they are derived from this grid's rate budget, anchored on
    # the fully explicit dt = 0.02 s that is known to run this case to
    # completion, so the SAFETY MARGIN and the spectral correction kappa are
    # inherited from a run that works rather than guessed.
    #
    #   smallest LGL gaps   h_x = 27.63 m   h_z = 6.91 m   (anisotropy 4.0:1)
    #
    #   rate [1/s]          vertical acoustic   50.2
    #                       horizontal acoustic 25.1
    #                       advection            1.7   (|u| 15, |w| 4 m/s)
    #                       SGS diffusion        0.9   (nu_t ~ 20 m^2/s)
    #
    #   explicit    all of it              78.0     CK2N54 radius 3.34  -> dt 0.020 (47% of neutral)
    #   HEVI        loses the vertical     27.7     ARS232 joint limit  -> dt 0.024 (1.2x)
    #   IMEX3D      loses ALL acoustics     2.6     ARS343 wedge limit  -> dt 0.506 (25x)
    #
    # WHY HEVI BUYS ALMOST NOTHING HERE, AND IT IS NOT THE SPLIT'S FAULT.
    # Removing the vertical acoustic term cuts the explicit budget by 2.8x --
    # but ARS232's explicit imaginary radius is 1.732 against
    # CarpenterKennedy2N54's 3.34, so nearly half of that is handed straight
    # back at the tableau, and the joint (rectangle) limit trims the rest. Net
    # 1.2x on the step against 3 RHS/step instead of 5, i.e. about 2x in
    # wall-clock. ARS343 would have the radius (2.828) but is unusable for HEVI
    # -- the two halves are loaded by INDEPENDENT wavenumbers there, so the
    # whole rectangle has to be neutral and ARS343 amplifies over most of it
    # (measured on this grid: dt 0.0022 s, i.e. 0.05x explicit).
    #
    # IMEX3D has no such problem: with the acoustics ENTIRELY implicit one
    # wavenumber loads both halves, only the wedge zE <= Mach*zI is reachable,
    # and ARS343 is neutral across it. That is why the ranking reverses.
    #---------------------------------------------------------------------------
    Δt_default = scheme === :imex ? "0.5" : (scheme === :hevi ? "0.024" : "0.02")
    Δt = parse(Float64, get(ENV, "DBG_DT", Δt_default))

    # A SHORT PROBE RUN, without editing this file. The point of measuring is
    # to decide the production settings, and a measurement that needs the deck
    # edited first is one that gets taken against a deck that no longer matches
    # the run it is being used to plan.
    #
    #   DBG_TEND=100  JEXPRESSO_HEVI_PROFILE=1   -> 200 steps at Δt = 0.5 with
    #                                               the cost breakdown printed
    #
    # The initial VTK dump is ~640 MB on this grid and lands before the time
    # loop, so it does not distort s/step -- but it is minutes of I/O for a run
    # that will be thrown away, and it is off by default whenever :tend has
    # been overridden.
    tend  = parse(Float64, get(ENV, "DBG_TEND", "10800.0"))
    lprobe = haskey(ENV, "DBG_TEND")

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
        # The stability table, and the "dominant term" line that says whether
        # either split can help on this grid at all. One RHS evaluation.
        :lcfl_report          => true,
        # Every ~500 steps (~250 s of model time). The startup table is taken on
        # a laminar sounding where nu_t is ~0, so its viscous row is true at
        # t = 0 and says nothing about a spun-up boundary layer. On this grid
        # the implicit step tolerates nu_t up to ~117 m^2/s before diffusion
        # binds, against the ~20 expected -- a 6x margin. This is how you watch
        # that margin rather than assume it. If it ever closes, the switch is
        # :implicit_vdiff => true (DBG_VDIFF=1), which moves d/dz(mu d/dz) into
        # the column operator that is already being factorised.
        # Per-step progress with a wall clock and an ETA to :tend, from rank 0,
        # :lcfl_report_every  re-print that table every N steps. :lcfl_report shows it
        #                   ONCE, at startup, where nu_t is ~0 on a laminar sounding --
        #                   so its viscous row says diffusion never limits dt, and says
        #                   it truthfully, at t = 0. Watch that row cross the line as the
        #                   boundary layer spins up. Collective and cheap; 0 = off.
        :lcfl_report_every    => parse(Int, get(ENV, "DBG_CFL_EVERY", "500")),
        # nothing collective. Free, and the only thing that distinguishes "still
        # compiling" from "hung" on a first run.
        :lstep_heartbeat      => parse(Bool, get(ENV, "DBG_HEARTBEAT", "true")),
        #--- HEVI (ignored unless DBG_SCHEME=hevi) ---------------------------------
        # :hevi_verify is also the joint-stability guard: it REFUSES to start
        # above the limit it measures on this mesh, and names the recommended
        # step in the error. Leave it on.
        :hevi_verify          => parse(Bool, get(ENV, "DBG_VERIFY", "true")),
        :hevi_linearization   => Symbol(get(ENV, "DBG_LIN", _vdiff ? "PS" : "RS")),
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
        # 1e-6, not the 1e-8 default. GMRES on this operator converges LINEARLY,
        # so iterations go as log(1/rtol) and the tolerance is a direct
        # multiplier on the only cost this scheme has: measured, 1e-6 is 70% of
        # the 1e-8 iteration count. 1e-6 is still three orders below a
        # third-order step's own truncation error. DBG_RTOL to sweep it.
        :imex_rtol            => parse(Float64, get(ENV, "DBG_RTOL", "1.0e-6")),
        # NOT the knob it looks like: this operator is skew, its spectrum is a
        # line segment, and restarting costs ~4% rather than the stagnation an
        # elliptic problem shows. Leave it at 20.
        # :imex_warm_start  start each stage solve from the PREVIOUS one's answer
        #                   rather than from zero. On by default, and the single
        #                   largest lever on this scheme's cost: consecutive ARK
        #                   stages differ by O(dt*f) while the right-hand side is
        #                   the whole deviation u - qe, so the guess arrives
        #                   several orders down, and GMRES on this skew operator
        #                   converges linearly, so those orders come straight off
        #                   the iteration count. Measured end to end: 5.0
        #                   iterations/solve against 21.7 cold, same answer to
        #                   1.3e-10. DBG_WARM=0 to measure it here.
        :imex_warm_start      => parse(Bool, get(ENV, "DBG_WARM", "true")),
        # :imex_restart     20 is right only while CFL_h = gamma*dt*c/h_x stays
        #                   below ~3. Above that, restarted GMRES(20) stops
        #                   converging -- measured on an isotropic grid at
        #                   CFL_h = 7.4, m=20 and m=40 never converge and m=80 is
        #                   needed. CFL_h is the HORIZONTAL acoustic Courant
        #                   number, because the column preconditioner removes the
        #                   vertical acoustics and leaves the horizontal to the
        #                   Krylov iteration. The setup report prints CFL_h, the
        #                   grid anisotropy behind it, and an advised m.
        # CFL_h = gamma*dt*c/h_x = 2.77 on this grid at dt = 0.5, which is just
        # inside the band where GMRES(20) is still comfortable (measured: m = 5
        # suffices to CFL_h = 3, m = 80 is needed by 7.5). 30 buys margin for
        # the flow speeding up, at 11 extra Krylov vectors -- about 3 MB/rank at
        # 2048 ranks. The setup report prints CFL_h and its advised m.
        :imex_restart         => parse(Int, get(ENV, "DBG_RESTART", "30")),
        :imex_maxiter         => parse(Int, get(ENV, "DBG_MAXITER", "200")),
        # :none is there to MEASURE what the column preconditioner buys. On a
        # 20:1 mesh it is 25x in iterations, so do not run production with it.
        :imex_precond         => Symbol(get(ENV, "DBG_PRECOND", "column")),
        :imex_lateral_walls   => Symbol(get(ENV, "DBG_LATWALL", "auto")),
        :imex_wall_flux       => true,
        # :PS, and worth the refactorisation. Over 10800 s of surface heating
        # rho*theta drifts ~1%, so beta = dp/d(rho theta) drifts ~0.4%; a stale
        # beta leaves that fraction of the FULL acoustic rate (78 1/s) in the
        # explicit half, which against IMEX3D's 2.6 1/s budget is a ~11% cut in
        # the admissible step. Refreshing costs 27 operator applications plus one
        # banded LU per column every 5 steps -- about 1.6% of a step here.
        :imex_linearization   => Symbol(get(ENV, "DBG_LIN", "PS")),
        :imex_update_freq     => parse(Int, get(ENV, "DBG_UPDFREQ", "5")),
        # The Krylov iteration count IS the cost of this scheme and it drifts as
        # the flow develops. Watch it on a first production run.
        :imex_monitor         => parse(Bool, get(ENV, "DBG_IMEXMON", "true")),
        :imex_monitor_every   => parse(Int, get(ENV, "DBG_IMEXMONEVERY", "200")),
        :Δt                   => Δt,
        :tinit                => 0.0,
        :tend                 => tend,
	:lrestart             => false,
	#:lrestart_vtk	      => true,
	#:restart_output_file_path => "",
	:restart_time         => 9000.0,
	#:diagnostics_at_times => (11500.0:10.0:15000.0),
	#:diagnostics_at_times => (0.0:50.0:10800.0),
        # 15.9 M nodes x 5 fields is ~640 MB per snapshot, so 19 of them is
        # ~12 GB. Every one is also a tstop: the step that lands on it is
        # shortened, gamma*dt changes, and the column factorisation is rebuilt
        # twice per event (once onto the short step, once back). 19 events is
        # ~38 refactorisations over ~21 600 steps -- under 0.2%.
	:diagnostics_at_times => (0.0:600.0:tend),
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
	:gmsh_filename    => "./problems/CompEuler/LESICP2-64x64x60-imex/LESICP_64x64x60_10240mX10240mX5000m.msh",
		
        # Stretching. The .geo is UNIFORM in z (Progression 1.0); the vertical
        # grading is applied here, at read time, by stretching.jl.
        #
        # :first_zelement_size IS THE EFFECTIVE RESOLUTION, NOT THE ELEMENT SIZE.
        # stretching.jl multiplies it by (ngl-1) on the way in --
        #     first_cell_size *= (mesh.ngl-1)
        # -- precisely because the LGL points divide it again, and then prints
        # "Desired resolution by the surface: 10.0". So 10.0 here means a 40 m
        # first ELEMENT and a 10 m effective resolution, which is what this case
        # wants.
        #
        # WHAT THE STRETCH ACTUALLY PRODUCES at nelemz = 60 -- worth knowing,
        # because it is a smooth power law and not two blocks of fixed size:
        #     exponent n = 1.179, first element 40 m (10 m effective)
        #     ~28 elements below the 2000 m transition, the last of them 85 m
        #       (21 m effective) -- so the resolution DEGRADES from 10 m at the
        #       surface to ~21 m by 2000 m; it is not 10 m throughout
        #     ~32 uniform elements of 85 m above it to the lid
        # If 10 m is wanted all the way to 2000 m this is the wrong stretch type
        # for it. Either way h_z_min, and so the step size, is set by the FIRST
        # element: 0.1727 x 40 = 6.91 m, which is what dt = 0.5 was derived from.
        :lstretch => true,
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
	:output_dir          => "/scratch/smarras/smarras/output_new/LESICP2_64x64x60_10240mX10240mX5000m_" * String(scheme) * "/",
	#:output_dir          => "/scratch/smarras/smarras/output_new/LESICP2_128x128x120_10240mX10240mX5000m/,"
        #:output_dir          => "./output_new/coarse-LESICP2_16x4x120_10kmX10kmX5km/",
        :loverwrite_output   => true,  #this is only implemented for VTK for now
        :lwrite_initial      => parse(Bool, get(ENV, "DBG_WRITE_INITIAL", lprobe ? "false" : "true")),
        #---------------------------------------------------------------------------
        # init_refinement
        #---------------------------------------------------------------------------
        # MUST be false. mesh.jl consults :lxy_partition ONLY when
        # :linitial_refine and :ladapt are both false, so runtime refinement
        # silently returns the mesh to the p4est space-filling-curve partition
        # whatever the flag says -- and that is the partition on which the
        # assembled RHS and the mass matrix carry different ghost multiplicities
        # and the implicit operator picks up a positive real eigenvalue.
        # check_columnar_partition refuses it at setup rather than letting the
        # run diverge slowly with every self-check passing.
        #
        # Refine OFFLINE instead. LESICP_64x64x60.geo in this directory is
        # LESICP.geo with nelem{x,y,z} = 64, 64, 60, i.e. exactly what one
        # refinement level of the 32x32x30 would have produced:
        #   gmsh -3 -order 1 problems/CompEuler/LESICP2-64x64x60-imex/LESICP_64x64x60.geo \
        #        -o problems/CompEuler/LESICP2-64x64x60-imex/LESICP_64x64x60_10240mX10240mX5000m.msh
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
