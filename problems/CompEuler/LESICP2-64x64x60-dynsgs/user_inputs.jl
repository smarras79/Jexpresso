#=============================================================================
 LESICP2 on the 64 x 64 x 60 production grid, with DynSGS instead of
 Smagorinsky -- 10240 x 10240 x 5000 m.

 This is LESICP2-64x64x60-imex with ONE physics change: the sub-grid closure.
 Grid, sounding, surface flux, wall model, filter, integrator and every solver
 setting are the same, deliberately, so the two decks are an A/B on the
 closure and nothing else. Read that deck's header for the grid and the IMEX
 argument; this one only explains what DynSGS is and what it does differently.

 WHAT DynSGS IS
 --------------
 Marras, Nazarov & Giraldo, JCP 301 (2015) 77-101. Instead of measuring the
 STRAIN RATE and assuming the sub-grid stress is proportional to it, DynSGS
 measures how badly the discrete solution fails to satisfy the PDE it is
 supposed to satisfy, and puts viscosity exactly there. Per element,

     R_i     = dq_i/dt + div F_i - s_i                      (the residual)
     mu_res  = C1 Delta^2 max_i ‖R_i‖_inf,e / ‖q'_i - <q'_i>‖_inf,Omega
     mu_max  = C2 Delta (‖v‖ + c)_inf,e
     nu      = max(0, min(mu_max, mu_res))                  [m^2/s]

 with dq_i/dt from a BDF2 stencil on the last three STEP states and div F - s
 read off the right-hand side the solver has just assembled. C1 = 1, C2 = 0.5;
 there is no tuned constant deciding WHERE the model is active, which is the
 whole point -- a Smagorinsky coefficient cannot tell a well-resolved shear
 layer from an under-resolved one and damps both.

 The coefficient is then handed to the SAME 3D stress tensor Smagorinsky uses
 (kernel/physics/SGS.jl, SGS_diffusion for NSD_3D):

     momentum   mu_t = rho * nu
     theta      kappa_t = mu_t / Pr_t          Pr_t = 0.7
     mass       0        (:mu[1] = 0 below)

 so everything downstream -- the CFL table's viscous row, :implicit_vdiff, the
 subfilter stresses in the LES statistics -- reports and uses the DynSGS
 coefficient without knowing it is not Smagorinsky.

 WHAT TO EXPECT, AND THE ONE THING TO WATCH
 ------------------------------------------
 DynSGS is a STABILISATION first and a physical closure second. In the
 surface layer of a wall-modelled PBL the sub-filter stress carries most of
 the momentum flux and has to be there whether or not the local solution is
 smooth -- that is what makes the log law -- and no residual sensor produces
 it on its own. So the first thing to look at in the output is not stability,
 it is the near-surface profiles:

   * mu_dsgs_rhou / mu_dsgs_rhotheta in the VTU are the per-element
     coefficients actually applied (piecewise constant per element).
   * the les_statistics *_sfs columns are the sub-filter stresses. If
     <u'w'>_sfs collapses in the first few hundred metres while <u'w'>_res
     does not rise to compensate, the surface layer is under-stressed and the
     wind profile will be too steep near the ground.

 If that happens, the lever is :dsgs_add_smagorinsky => true, which makes the
 eddy viscosity mu_smag + mu_dsgs rather than mu_dsgs alone: the closure keeps
 its wall behaviour and DynSGS adds dissipation only where the residual says
 the discretisation needs it. It is off here because whether the residual term
 alone is enough is precisely the question this deck exists to answer.

 STABILITY BUDGET, WHICH IS NOT THE SAME AS THE SMAGORINSKY DECK'S
 -----------------------------------------------------------------
 The reason LESICP2-64x64x60-imex died at t = 500 s is that its SGS diffusion
 is explicit and nu_eff grows with the boundary layer until 2 nu_eff/h_z^2
 crosses what ARS343's explicit tableau will carry. DynSGS does not exempt
 this deck from that: it produces an eddy viscosity like any other closure and
 that viscosity has the same parabolic limit. What changes is only WHERE the
 viscosity sits -- concentrated on the elements whose residual is large rather
 than spread over every element with a strain rate.

 So :implicit_vdiff is ON here too, for the same reason and by the same
 default. Watch the "viscous" term in the IMEX self-check at startup and in
 the periodic CFL report (:lcfl_report_every) rather than assuming it.

 DIFFERENCES FROM THE SMAGORINSKY DECK, IN FULL. There are four:
   :visc_model            SMAG() -> DSGS()
   :mu[5]                 2.0 -> 1.0     (see the note at :mu below)
   :dsgs_C1 / :dsgs_C2    new, the model's two coefficients
   :output_dir            ..._dynsgs_ so an A/B leaves both on disk

 RUN IT ON A RANK COUNT THAT DIVIDES 4096
 ----------------------------------------
 :lxy_partition splits the 64 x 64 = 4096 ELEMENT COLUMNS, so

     512 ranks   8 columns each   exact
    1024 ranks   4 columns each   exact
    2048 ranks   2 columns each   exact

     DBG_SCHEME=imex (default) | hevi | explicit   -- same physics, three arms
     DBG_FILTER=1                                  -- modal filter ON (off by default)
     DBG_DSGS_C1=0                                 -- DynSGS off (nu = 0), everything else identical
     DBG_DSGS_C1=... DBG_DSGS_C2=...               -- sweep the two coefficients
     DBG_DSGS_RES=strict                           -- true BDF2 residual instead of the tendency form
     DBG_DT=... DBG_RTOL=... DBG_RESTART=...       -- sweep the three knobs

 THE TWO IMEX ARMS. :imex_schur and :implicit_vdiff cannot both be on -- the
 Schur reduction inverts the momentum row pointwise and cannot see a vertical
 diffusion operator -- so `use_schur = !_vdiff` picks one:

     (default)     implicit vertical diffusion, FIVE-FIELD stage solve
                   stable as nu_t grows; ~3.56x slower per step

     DBG_VDIFF=0   explicit vertical diffusion, SCHUR stage solve
                   3.56x faster per step; the vertical viscous rate
                   2*mu[5]/Pr_t*nu_t/h_z^2 is back in the explicit budget and
                   ARS343 carries it only to nu_t ~ 40-80 m^2/s here (twice the
                   Smagorinsky deck's 20-40, because :mu[5] is 1.0 not 2.0)

 Both take ALL acoustics implicitly at dt = 0.5. Watch the CFL report's viscous
 row and JEXPRESSO_DSGS_MONITOR=1: the DBG_VDIFF=0 arm ends when nu_t reaches
 the number above, and whether DynSGS gets there is the question it answers.
=============================================================================#

function user_inputs()

    # IMPLICIT VERTICAL DIFFUSION -- ON BY DEFAULT, AND WHY THAT CHANGED.
    #
    # This was opt-in (DBG_VDIFF=1) and the run died at t = 500 s. A blow-up at
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
        error("LESICP2-64x64x60-dynsgs: DBG_SCHEME must be imex, hevi or explicit; got $scheme")

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
        #---------------------------------------------------------------------------
        # THE CLOSURE -- the only physics that differs from LESICP2-64x64x60-imex
        #---------------------------------------------------------------------------
        :visc_model           => DSGS(),
        # The two DynSGS coefficients, eq. (9). These ARE the paper's values;
        # they are written out rather than left to the defaults in
        # io/mod_inputs.jl so that an A/B that changes them leaves a trace in
        # the deck.
        #
        # C1 scales the residual viscosity. 1.0 is parameter-free in the sense
        # that matters -- it does not decide WHERE the model is active, only
        # how hard it acts once the residual has decided.
        # DBG_DSGS_C1=0 makes mu_res = 0 and therefore nu = 0: DynSGS becomes a
        # no-op with every other setting untouched, which is the one-knob test for
        # "is the coefficient what is killing this run?".
        :dsgs_C1              => parse(Float64, get(ENV, "DBG_DSGS_C1", "1.0")),
        # C2 scales the first-order-upwind cap mu_max = C2 Delta (‖v‖+c). At
        # Delta = 40 m with c ~ 340 against ‖v‖ ~ 10 that is 0.5*40*350 = 7000
        # m^2/s -- and THE CAP IS NOT A SAFETY NET AT THIS VALUE.
        #
        # WITH THE DIFFUSION EXPLICIT (DBG_VDIFF=0), what the step can carry is
        #     dt * 2 * (mu[5]/Pr_t) * nu_t / h_z^2  <=  ~2   (ARS343, real axis)
        # and at h_z = 6.91 m, mu[5] = 1.0, Pr_t = 0.7, dt = 0.5 that is
        #     nu_t <~ 65 m^2/s.
        # The cap sits 100x above it, so it can only bind long after the run is
        # already dead. If the monitor (JEXPRESSO_DSGS_MONITOR=1) shows nu
        # pinned at 7000, C2 ~ 0.004 brings the cap to ~56 m^2/s -- but that is
        # then a first-order-upwind viscosity wearing DynSGS's name, not a
        # residual model, and lowering :dsgs_C1 or :dsgs_residual => :strict is
        # the honest lever. On the implicit arm this rate is removed and the
        # ceiling is set by accuracy instead.
        :dsgs_C2              => parse(Float64, get(ENV, "DBG_DSGS_C2", "0.5")),
        # :tendency (default) rolls the history BEFORE the read, so R is ~1.5x
        # the physical TENDENCY, not a truncation error -- it is large wherever
        # the solution is changing fast, which in a spinning-up PBL is
        # everywhere. :strict rolls after, giving the literal BDF2 residual;
        # measured 10x smaller on sod1d. DBG_DSGS_RES=strict to try it.
        :dsgs_residual        => Symbol(get(ENV, "DBG_DSGS_RES", "tendency")),
        # Add the Smagorinsky viscosity to the residual one instead of
        # replacing it. OFF: see "WHAT TO EXPECT" in the header for the
        # near-surface diagnostic that decides whether it should be on.
        :dsgs_add_smagorinsky => false,
        # Read only when :dsgs_add_smagorinsky is on. Left at the Smagorinsky
        # deck's value so the two arms are comparable if it is switched on.
        :C_s                  => 0.16,
        # Likewise: the buoyancy stability function multiplies the SMAGORINSKY
        # part only. A residual is not a strain rate, and rescaling it by a
        # Richardson number would be asking the sensor a question it has
        # already answered. With :dsgs_add_smagorinsky => false these two keys
        # affect only the N2 and f_Ri diagnostic fields.
        :lrichardson          => true,
        :lwall_damping        => true,
        # :μ IS A PER-EQUATION MASK, NOT A VISCOSITY (same as under SMAG).
        #
        # NOTE slot 5 is 1.0 HERE and 2.0 in the Smagorinsky deck. That 2.0
        # gives the theta equation :μ[5]/Pr_t = 2.857 x the momentum
        # diffusivity -- an effective turbulent Prandtl number of 0.35 -- which
        # is a hand-tuning of the Smagorinsky arm with no counterpart here: the
        # DynSGS coefficient is derived from the residual of the theta equation
        # itself, so there is nothing for a factor of two to correct. It is
        # also the term that binds the explicit parabolic limit, so carrying it
        # over would import a stability problem for no modelling reason.
        :μ                    => [0.0, 1.0, 1.0, 1.0, 1.0],
        # Domain (rather than rank-local) norms for <q'> and ‖q' - <q'>‖. Two
        # extra Allreduce per RHS call, 8 per step under ARS343. Off by
        # default; switch it on when mu has to be identical across rank counts.
        # :ldsgs_global_norms   => true,
        #---------------------------------------------------------------------------
        #LES statistics
        #---------------------------------------------------------------------------
	# BOTH ranges have to be splatted. With the `...` on only the first, this
	# is a tuple of ten Floats followed by a StepRangeLen, and
	# TimeIntegrators.jl does `collect(Float64, les_stat_t)` on it -- which
	# cannot convert a range to a Float64 and dies at startup with a
	# MethodError. Splatting both gives the 191 sample times intended.
	:statistics_time      => (10.0:10:100.0..., 9000.0:10:10800.0...),
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
	# SHARED WITH LESICP2-64x64x60-imex ON PURPOSE. It is a 28 MB file and
	# the two decks are an A/B on the closure: a second copy is 28 MB of
	# disk whose only possible future is to drift out of step with the one
	# the other arm reads. LESICP_64x64x60.geo is copied here so the mesh
	# can be regenerated from this directory if that one ever goes away.
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
        :lfilter             => parse(Bool, get(ENV, "DBG_FILTER", "false")),
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
	# The ARM is in the path, not just the scheme: the two IMEX configurations
	# below are a deliberate A/B and writing both into one directory loses the
	# first one silently.
	#   _dynsgs_imex_schur   Schur stage solve, EXPLICIT vertical diffusion
	#   _dynsgs_imex_vdiff   five-field stage solve, IMPLICIT vertical diffusion
	:output_dir          => "/scratch/smarras/smarras/output_new/LESICP2_64x64x60_10240mX10240mX5000m_dynsgs_" *
	                        String(scheme) * (scheme === :imex ? (_vdiff ? "_vdiff" : "_schur") : "") * "/",
	#:output_dir          => "/scratch/smarras/smarras/output_new/LESICP2_128x128x120_10240mX10240mX5000m/,"
        #:output_dir          => "./output_new/coarse-LESICP2_16x4x120_10kmX10kmX5km/",
        :loverwrite_output   => true,  #this is only implemented for VTK for now
        :lwrite_initial      => parse(Bool, get(ENV, "DBG_WRITE_INITIAL", lprobe ? "false" : "true")),
        #---------------------------------------------------------------------------
        # AMR
        #---------------------------------------------------------------------------
        :linitial_refine     => false,
        :init_refine_lvl     => 1,
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
