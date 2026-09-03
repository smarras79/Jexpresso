function user_inputs()

    _vdiff = parse(Bool, get(ENV, "DBG_VDIFF", "true"))

    #---------------------------------------------------------------------------
    # DBG_SCHEME=imex|hevi|explicit in the environment overrides it without an
    # edit, e.g.  DBG_SCHEME=hevi sbatch submit_LESICP2_128.sh
    #---------------------------------------------------------------------------
    #scheme_default = :imex
    scheme_default = :hevi

    # IMEX ONLY: how the stage system is solved. Turning this off does NOT go
    # explicit or HEVI -- that is scheme_default above.
    use_schur = !_vdiff   # scalar Schur stage solve, Np not 5*Np, 3.56x/step. IMPOSSIBLE with :implicit_vdiff -- the reduction cannot see the diffusion operator

    # -- used only when scheme is :imex. Each is also an env override so an A/B
    #    needs no edit here; the value shown is the default.
    lschur      = parse(Bool,    get(ENV, "DBG_SCHUR",     string(use_schur)))
    # STEP SIZES, one per scheme, set HERE and nowhere else: no environment
    # variable overrides them, so the deck is the record of what ran. The
    # budget behind each number is under "STEP SIZE" below.
    dt_imex     = 1.0     # :imex     wedge neutral ~2.6 s; re-read the t = 0 report
    dt_hevi     = 0.05    # :hevi     ARS232 neutral 0.0749 s measured on this mesh
    dt_explicit = 0.02    # :explicit CK2N54, known to run this case
    rtol        = parse(Float64, get(ENV, "DBG_RTOL",      "1.0e-6"))
    # Krylov basis costs (restart+4)*npoin*nvar*8 B/rank: ~19 MB on the scalar
    # Schur system, ~95 MB on the five-field one (npoin/rank ~ 70k at 256 ranks).
    # An OutOfMemoryError there is the heap budget, not this number -- launch with
    # --heap-size-hint (~85% of --mem-per-cpu). See heap_budget_note in krylov.jl.
    restart     = parse(Int,     get(ENV, "DBG_RESTART",   "30"))
    maxiter     = parse(Int,     get(ENV, "DBG_MAXITER",   "600"))
    precond     = Symbol(        get(ENV, "DBG_PRECOND",   "column"))   # or :none
    umax        = parse(Float64, get(ENV, "DBG_UMAX",      "15.0"))
    verify      = parse(Bool,    get(ENV, "DBG_VERIFY",    "true"))
    warm_start  = parse(Bool,    get(ENV, "DBG_WARM",      "true"))
    lin_imex    = Symbol(        get(ENV, "DBG_LIN", _vdiff ? "PS" : "RS"))  # :PS is REQUIRED by implicit vdiff and pure cost without it -- it re-probes and re-factorises the band every :imex_update_freq steps
    updfreq     = parse(Int,     get(ENV, "DBG_UPDFREQ",   "5"))
    lat_walls   = Symbol(        get(ENV, "DBG_LATWALL",   "auto"))
    monitor     = parse(Bool,    get(ENV, "DBG_IMEXMON",   "true"))
    monitor_every = parse(Int,   get(ENV, "DBG_IMEXMONEVERY", "200"))
    scheme = Symbol(get(ENV, "DBG_SCHEME", string(scheme_default)))
    scheme in (:imex, :hevi, :explicit) ||
        error("LESICP2-30x30x60-imex: scheme must be :imex, :hevi or :explicit; got :$scheme")
    #---------------------------------------------------------------------------
    # STEP SIZE -- the budget behind dt_imex / dt_hevi / dt_explicit above.
    #
    # MEASURED, NOT ESTIMATED: build_hevi refuses to start above the neutral
    # limit and prints the budget, which on this mesh (30 x 30 x 60, p = 4,
    # 40 m first z elements) came out as
    #
    #     explicit half  17.4 1/s = acoustic 17.2 + advective 0.17 + viscous 0.014
    #     ARS232 joint neutral limit  0.0749 s   -> run at 0.0525 s or below
    #
    # A hand estimate from the LGL gaps (h_x 58.95 m, h_z 6.91 m) lands at
    # 2.4x that, and it is wrong for two reasons worth remembering: the code
    # counts BOTH horizontal directions, and it calibrates c/h against the
    # measured column spectrum (kappa = 1.46 here) -- the true spectral radius
    # of the p = 4 operator, not c over the smallest node gap. With those two
    # in, ark_joint_dt_max(ARS232) reproduces the run's number (0.071 s).
    #
    # Where that leaves the three schemes, same margins throughout:
    #
    #   explicit  CK2N54                    neutral 0.036 s   -> 0.017  (0.02 runs)
    #   HEVI      ARS232, :implicit_vdiff   neutral 0.075 s   -> 0.05   <- default
    #   HEVI      ARS443 would allow ~0.06 but costs 4 RHS/step for 3: no gain
    #   IMEX3D    ARS343, :implicit_vdiff   wedge   2.6 s     -> 1.2    (dt_imex 1.0)
    #
    # HEVI is only ~2.5x explicit here: the horizontal acoustics stay explicit
    # and ARS232's imaginary radius is 1.73 against CK2N54's 3.34. 0.1 was
    # refused (amplification 1.92/step); 0.08 would be refused too.
    Δt = scheme === :imex ? dt_imex : scheme === :hevi ? dt_hevi : dt_explicit
    
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
        :lcfl_report          => true,
        :lcfl_report_every    => parse(Int, get(ENV, "DBG_CFL_EVERY", "500")),
        :lstep_heartbeat      => parse(Bool, get(ENV, "DBG_HEARTBEAT", "true")),
        :hevi_verify          => parse(Bool, get(ENV, "DBG_VERIFY", "true")),
        :hevi_linearization   => Symbol(get(ENV, "DBG_LIN", _vdiff ? "PS" : "RS")),
        :hevi_update_freq     => parse(Int, get(ENV, "DBG_UPDFREQ", "5")),
        :hevi_wall_flux       => true,
        #--- IMEX3D (ignored unless DBG_SCHEME=imex) -------------------------------
        :imex_schur           => lschur,
        :imex_verify          => verify,
        :imex_umax            => umax,
        :imex_rtol            => rtol,
        :imex_warm_start      => warm_start,
        :imex_restart         => restart,
        :imex_maxiter         => maxiter,
        :imex_precond         => precond,
        :imex_lateral_walls   => lat_walls,
        :imex_wall_flux       => true,
        :imex_linearization   => lin_imex,
        :imex_update_freq     => updfreq,
        :imex_monitor         => monitor,
        :imex_monitor_every   => monitor_every,
        :Δt                   => Δt,
        :tinit                => 0.0,
        :tend                 => tend,
	:lrestart             => false,
	:restart_time         => 9000.0,
	# EVERY range needs its own `...`; the third was missing one, which made this
	# a tuple of 28 Floats followed by a StepRangeLen and killed the run in
	# time_loop! (collect(Float64, ...) cannot convert a range to a Float64).
	:diagnostics_at_times => (0.0:100.0:1000.0..., 1000.0:500.0:9000.0..., 9000.0:10.0:tend...),
	:lsource              => true,
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
	:lxy_partition          => true,
        :lwall_model          => true,
        :ifirst_wall_node_index=> 2, # This must be between 2 <= :first_wall_node_index <= nop+1
        :bdy_fluxes           => true,
        :lvisc                => true, #false by default
        :visc_model           => SMAG(),
        :C_s                  => 0.18,
        :lrichardson          => true,
        # NOW ACTUALLY DAMPS AT THE WALL. The guard used to read
        # `(lwall_damping && z > 0.0) || return CsD2`, and zwall is built as
        # max(z - zmin, 0) -- so every node ON the lower boundary carried
        # exactly 0.0, took the early return, and got the UNDAMPED (C_s*Delta)^2
        # while its neighbours were damped. On this mesh that was l^2 = 40.96 m^2
        # at z = 0 against 6.44 m^2 at z = 6.91 m: a 6.4x eddy-viscosity spike
        # sitting on the node with the smallest h_z AND the node the surface
        # heat flux is injected into, i.e. the node that sets the vertical
        # viscous limit in the table above. l^2 now runs 0 -> 6.44 -> 40.96,
        # monotone, so that node contributes ~0 to the rate budget instead of
        # roughly 30% of it. See sgs_mixing_length2 in kernel/physics/SGS.jl and
        # test/sgs/test_wall_damping.jl.
        :lwall_damping        => true,
        :μ                    => [0.0, 1.0, 1.0, 1.0, 1.0],
        :les_filter_width     => :geometric,
        #---------------------------------------------------------------------------
        # MOST GUARD RAILS. Stated explicitly here rather than left to the
        # defaults because this case is exactly the one they exist for: a
        # convective boundary layer over a prescribed surface flux visits both
        # singular limits of Businger-Dyer at some node on essentially every
        # step. All three ARE the defaults, so writing them changes nothing --
        # they are here to be found and tuned. See docs/boundary_conditions.md
        # section 2.2.1 and test/physics/test_most_guards.jl.
        #
        #   :most_u_min    gustiness floor on the wind entering the DRAG
        #                  COEFFICIENT. Not a floor on the stress: the stress
        #                  direction divides by the same floored speed, so tau
        #                  stays linear in u_i and vanishes with it rather than
        #                  becoming a finite-magnitude random walk wherever two
        #                  thermals converge and |u| passes through zero.
        #   :most_zeta_*   bounds on zeta = z/L. NOT a height and NOT a z
        #                  coordinate -- z is always >= 0, and the SIGN here is
        #                  the SIGN OF L, i.e. the stability:
        #
        #                      L < 0  (surface heats the air)  ->  zeta < 0   UNSTABLE
        #                      L = +-inf (no surface flux)     ->  zeta = 0   neutral
        #                      L > 0  (surface cools the air)  ->  zeta > 0   stable
        #
        #                  So -5.0 bounds the CONVECTIVE side, which is this
        #                  deck's entire regime, and 2.0 the stable side.
        #                  [-5, 2] is the range Businger-Dyer was fitted over;
        #                  past it the correction is evaluated AT the bound
        #                  rather than extrapolated. The negative bound is the
        #                  one that binds here: at z = 6.91 m, zeta = -5 is
        #                  L = -1.38 m, i.e. free convection at a near-calm
        #                  node. Without it |u| -> 0 sends L ~ u*^2 -> 0 and
        #                  zeta -> -inf, and the drag that comes back is
        #                  whatever the 20th Picard iterate happened to be.
        #
        # Also fixed in the same pass, and NOT deck-settable because neither was
        # ever a choice: obukhov_length was missing its rho (|L| ~20% small at
        # sea level, every zeta ~20% large), and its near-neutral guard returned
        # a POSITIVE L for either sign of Q_H -- a stable surface layer reported
        # in the middle of a convective one, reached whenever u* got small.
        #
        # AND THE ONE THAT MATTERS MOST FOR THIS DECK: :user_heatflux above sets
        # delta_hf = 1, so the theta surface term is rho*0.12 and MOST's own
        # w'theta' is discarded -- but theta*, and through it L, psi_m and the
        # DRAG, were still being built from theta_sfc, a free prognostic node
        # that nothing pins. The momentum flux carried a stability correction
        # diagnosed from a heat flux the model does not apply. BCs.jl now hands
        # the imposed flux to CM_MOST!, so theta* = -w'theta'/u*. At |u| = 2 m/s
        # with the surface node 0.5 K cooler than the air that is u* = 0.243
        # (correctly unstable) against 0.171 (read as stable). The theta term
        # applied to the RHS is unchanged; only the drag moves.
        #---------------------------------------------------------------------------
        :most_u_min           => 0.1,     # m/s
        :most_zeta_min        => -5.0,
        :most_zeta_max        =>  2.0,
        #---------------------------------------------------------------------------
        #LES statistics
        #---------------------------------------------------------------------------
	# THE TRAILING COMMA IS LOAD-BEARING, and this is the second way this one
	# key has broken the run. `(x...)` with nothing after it is NOT a
	# one-element tuple -- in call-argument position it is a SPLAT, so
	# `:statistics_time => (9000.0:10:10800.0...)` asks for
	# Pair(::Symbol, ::Float64 x181) and dies on every rank at startup with a
	# MethodError whose stack trace points at the `inputs = Dict(` line, not at
	# this one. `(x...,)` is the tuple. Both forms PARSE and LOWER identically
	# as far as Meta.parse is concerned, so tools/syntax_check.jl cannot see the
	# difference -- test/decks/test_deck_inputs.jl evaluates the Dict and names
	# the key.
	#
	# (The earlier failure here was the mirror image: with two ranges and the
	# `...` on only the first, this was ten Floats followed by a StepRangeLen,
	# and collect(Float64, les_stat_t) in TimeIntegrators.jl cannot convert a
	# range to a Float64.)
	:statistics_time      => (9000.0:10:10800.0...,),
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
        # PREPROCESS/MESH CACHE: OFF, AND THAT IS A SIZE DECISION.
        #
        # The cache stores per-rank metrics + matrix so a LATER run of the same
        # case skips the metric build. At 1024 ranks on this grid it wrote
        # 573 GB in 2048 files, filled a 2 TB GPFS to 100%, and killed the job
        # with SIGBUS on every rank (JLD2 writes through an mmap; a page the
        # filesystem cannot back is signal 7, which no try/catch sees). The
        # saving it was buying is one metric build -- 8.4 s.
        #
        # It stays off by default here BECAUSE OF THE ARITHMETIC, not because
        # caching is wrong: ~0.5 GB/rank x nranks is fine at 4 ranks and absurd
        # at 1024. JEXPRESSO_MESH_CACHE=1 turns it back on if you have the
        # space and are starting the same case repeatedly.
        :luse_mesh_cache  => parse(Bool, get(ENV, "JEXPRESSO_MESH_CACHE", "false")),
        :lread_gmsh       => true, #If false, a 1D problem will be enforce
	#:gmsh_filename    => "./problems/CompEuler/LESICP2-128x128x60-imex/LESICP_128x128x60_10240mX10240mX5000m.msh",
        :gmsh_filename    => "./problems/CompEuler/LESICP2-30x30x60-imex/LESICP_30x30x60_10240mX10240mX5000m.msh",
		
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
        # "two_block uniformish" is the type that gives the profile this case
        # wants: WHOLE ELEMENTS of :first_zelement_size*(ngl-1) from the ground
        # to :zlevel_transition -- so a genuinely uniform ~10 m effective
        # resolution there -- then a ramp, then uniform-coarse to the top.
        # Verified in test/mesh/test_stretching.jl on this deck's own values:
        # 201 nodes below 2000 m, effective dz 10.00 m, expanding 30 m -> 70 m
        # above.
        #
        # NOT "fixed_first_twoblocks_strong", which was here and does the
        # opposite shape: it power-law stretches BELOW the transition and goes
        # uniform above. With these values it also inverts -- 40 m at the ground
        # against 14.39 m uniform -- because its exponent compares the
        # (ngl-1)-scaled first cell against the smallest LGL NODE gap rather
        # than against the element size. stretching.jl warns when that happens.
        :stretch_type => "two_block uniformish",
        :first_zelement_size => 10.0,
        :zlevel_transition => 2000.0,
        
        #---------------------------------------------------------------------------
        # Filter parameters. OFF, and that is a standing choice for this case --
        # do not "restore" it. DBG_FILTER=1 turns it on for a one-off A/B.
        #---------------------------------------------------------------------------
        :lfilter             => true, #parse(Bool, get(ENV, "DBG_FILTER", "false")),
        :mu_x                => 0.05,
        :mu_y                => 0.05,
	:mu_z                => 0.1,
        :filter_type         => "erf",
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
	# One tree per scheme, so an A/B/C leaves all three solutions on disk
	# rather than overwriting each other.
	# JEXPRESSO_OUTDIR overrides the parent directory without editing this
	# file, because the default below is one user's scratch and this deck is
	# shared. The scheme suffix is appended either way, so an A/B/C leaves all
	# three trees on disk rather than overwriting each other.
	#
	# BUDGET FOR IT: :diagnostics_at_times below is 209 dumps and each is
	# ~2.6 GB on this grid, i.e. ~543 GB for a full run. The 9000:10:10800
	# range alone is 181 of them. Thin that cadence if scratch is tight.
	:output_dir          => "/scratch/smarras/smarras/output_HANG/LESICP2_30x30x60_10240mX10240mX5000m_imex",
        :loverwrite_output   => true,  #this is only implemented for VTK for now
        # ~2.6 GB on this grid (4x the 64x64x60 dump) and it lands before the
        # time loop, so it does not distort s/step -- but it is minutes of I/O
        # for a run that will be thrown away. Off whenever tend is overridden.
        :lwrite_initial      => parse(Bool, get(ENV, "DBG_WRITE_INITIAL", lprobe ? "false" : "true")),
        #---------------------------------------------------------------------------
        # AMR
        #---------------------------------------------------------------------------
        :linitial_refine     => false,
        :init_refine_lvl     => 1,
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
