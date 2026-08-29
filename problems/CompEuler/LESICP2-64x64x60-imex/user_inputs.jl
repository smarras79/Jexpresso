function user_inputs()

    # IMPLICIT VERTICAL DIFFUSION -- ON BY DEFAULT, AND WHY THAT CHANGED.
    #
    # This was opt-in (DBG_VDIFF=1) and the run died at t = 500 s -- REPRODUCED
    # at t = 550 s with the modal filter off and the Schur arm, so the filter
    # was never implicated and the explicit viscous rate is the whole story:
    # dt*2*(:mu[5]/Pr_t)*nu_t/h_z^2 <= 2 caps nu_t at 33 m^2/s here. A blow-up at
    # a fixed MODEL time, on a step size the acoustic analysis below says is
    # 47% of neutral, is the signature of a rate that is not in that analysis
    # and arrives on its own schedule -- and on this mesh there is exactly one:
    # the SGS diffusion of the vertical derivative, nu_eff/dz^2, which is
    # invisible at t = 0 (the CFL report is taken on a laminar sounding where
    # nu_t ~ 0) and grows with the boundary layer. See the rate table under
    # "STEP SIZE" below for the arithmetic; the short version is that at
    # dz_min = 6.91 m the vertical viscous rate reaches the whole explicit
    # budget by the time nu_t reaches ~30 m^2/s, which a convective boundary
    # layer does in the first few hundred seconds.
    #
    # :implicit_vdiff puts d/dz(mu d/dz) on u, v, w and theta into the same
    # column operator that already carries the vertical acoustics: a wider
    # band, no new solve. It is the only lever here that removes the term
    # rather than tiptoeing around it, which is why it is now the default.
    # DBG_VDIFF=0 restores the old behaviour for an A/B.
    #
    # Read into a local rather than inline because two keys below have to agree
    # about it -- switching it on also switches the linearisation to :PS,
    # without which the implicit operator would carry no diffusion at all.
    _vdiff = parse(Bool, get(ENV, "DBG_VDIFF", "true"))

    #---------------------------------------------------------------------------
    # TWO INDEPENDENT SWITCHES: the first picks the INTEGRATOR, the second only
    # picks how the IMEX stage system is solved. Turning the second off does
    # NOT go explicit.
    #---------------------------------------------------------------------------
    use_imex  = true      # true -> IMEX-ARK(ARS343), all acoustics implicit, dt 0.5;  false -> fully explicit, dt 0.02
    use_schur = !_vdiff   # IMEX only: scalar Schur stage solve, Np not 5*Np, 3.56x/step. IMPOSSIBLE with :implicit_vdiff -- the reduction cannot see the diffusion operator

    # -- used only when use_imex. Each is also an env override so an A/B
    #    needs no edit here; the value shown is the default.
    lschur      = parse(Bool,    get(ENV, "DBG_SCHUR",     string(use_schur)))
    dt_imex     = parse(Float64, get(ENV, "DBG_DT",        "0.5"))
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
    scheme = Symbol(get(ENV, "DBG_SCHEME", use_imex ? "imex" : "explicit"))
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
    #                       SGS diffusion        2.4   (nu_eff ~ 57 m^2/s, see below)
    #
    #   explicit    all of it              79.5     CK2N54 radius 3.34  -> dt 0.020 (47% of neutral)
    #   HEVI        loses the vertical     29.2     ARS232 joint limit  -> dt 0.024 (1.2x)
    #   IMEX3D      loses ALL acoustics     4.1     ARS343 wedge limit  -> dt 0.506 (25x)
    #
    # WHERE THE SGS ROW COMES FROM, AND WHY IT IS THE ONE THAT KILLED THE RUN.
    # nu_eff is the largest diffusivity ANY equation is handed, not nu_t. With
    # :mu = [0, 1, 1, 1, 2] and Pr_t = 0.7 the theta equation gets
    # 2/0.7 = 2.857 x nu_t, so nu_eff = 2.857 * nu_t, and the vertical rate is
    #
    #     2 nu_eff / h_z^2   with h_z = 6.91 m
    #
    #     nu_t [m^2/s]     5     10     20     40     80
    #     rate [1/s]      0.6    1.2    2.4    4.8    9.6
    #     dt * rate       0.3    0.6    1.2    2.4    4.8      at dt = 0.5
    #
    # ARS343's explicit tableau is neutral on the negative real axis out to
    # about 2, so this term alone takes the step unstable somewhere between
    # nu_t = 20 and 40 m^2/s -- which a convective boundary layer over a
    # 0.12 K m/s surface flux reaches within the first few hundred seconds.
    # That is the t = 500 s failure, and it is why :implicit_vdiff is now on by
    # default at the top of this function rather than opt-in: it is the only
    # lever that REMOVES this rate.
    #
    # WHAT DOES NOT FIX IT, RECORDED SO IT IS NOT TRIED AGAIN:
    #   * raising :C_s. nu_t goes as C_s^2, so 0.16 -> 0.20 multiplies every
    #     number in the table above by 1.56 and brings the failure FORWARD.
    #   * the filter. It damps the top Legendre mode; it does not change an
    #     eigenvalue of the semi-discrete operator. Worth having anyway (it is
    #     back on below), but it is not this.
    #   * a smaller dt. It works -- dt = 0.2 buys a factor of 2.5 -- and it
    #     costs the entire reason this deck exists.
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
    # dt_imex is set in the switch block above (and honours DBG_DT); the other
    # two schemes keep their own limits here.
    Δt = scheme === :imex ? dt_imex :
         parse(Float64, get(ENV, "DBG_DT", scheme === :hevi ? "0.024" : "0.02"))

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
	# BOTH ranges have to be splatted. With the `...` on only the first, this
	# is a tuple of ten Floats followed by a StepRangeLen, and
	# TimeIntegrators.jl does `collect(Float64, les_stat_t)` on it -- which
	# cannot convert a range to a Float64 and dies at startup with a
	# MethodError. Splatting both gives the 191 sample times intended.
	:statistics_time      => (9000.0:10:10800.0...),
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
	:output_dir          => "/scratch/smarras/smarras/output_new/filter_LESICP2_64x64x60_10240mX10240mX5000m_" * String(scheme) * "/",
        :loverwrite_output   => true,  #this is only implemented for VTK for now
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
