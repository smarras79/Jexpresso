#---------------------------------------------------------------------------------
# FREE-STREAM PRESERVATION TEST — THE RUNG LADDER
#
# One property, one number, one flip at a time. The number is the RANGE of
# dp_rel in ParaView's Information tab at the last output frame. Every rung
# is ~500 steps, i.e. seconds of wall clock.
#
# THE RUNGS ARE CUMULATIVE: rung N means every block up to and including N
# is uncommented. Two of them live in user_bc.jl (the wall and the outflow)
# because user_bc_dirichlet! is not handed `inputs`; the rest are blocks at
# the bottom of this file, each one `inputs[:key] = value` so uncommenting
# can never collide with the base Dict.
#
#   RUNG  what it adds                        where            result
#   ----  ----------------------------------  ---------------  --------------
#     0   bare inviscid operator, no wall     (base deck)      PASS 5.8e-13
#     1   :lkep + ranocha flux differencing   base deck        PASS 2.5e-11
#     2   the no-slip isothermal wall         user_bc.jl       ?
#     3   Sutherland molecular viscosity      this file        ?
#     4   DynSGS, :dsgs_norms => "domain"     this file        ?
#     5   DynSGS, :dsgs_norms => "rank"       this file        ?
#    (+)  positivity, to get a coordinate     this file        optional
#
# RUNG 1 IS DONE AND IT PASSED: dp_rel came back at -1.9e-11 .. +2.5e-11,
# uniform, with no structure at the outflow, at the 15-degree kink or at the
# block junction. It is now part of the BASE DICT above (:lkep => true,
# :volume_flux => ranocha()), so every rung from here on carries it.
#
# 43x rung 0, and still pure round-off rather than a source. The budget says
# so -- 500 steps x 5 stages = 2500 RHS evaluations:
#
#   rung 0   5.8e-13 / 2500 = 2.3e-16 per evaluation =  1.0 eps
#   rung 1   2.5e-11 / 2500 = 1.0e-14 per evaluation = 45.0 eps
#
# 45 eps per evaluation is what flux differencing COSTS: instead of one
# pointwise flux per node it evaluates a two-point flux against every other
# node in the line, each with a logarithmic mean carrying a series expansion.
# A longer arithmetic chain accumulates more round-off; it does not
# manufacture a source. A broken FSP would show as STRUCTURE that GROWS, a
# pattern locked to the mesh and orders of magnitude larger. The gap still to
# be explained is 2.7e9x: the production first repair was at p = -52.1 Pa,
# i.e. dp_rel = -6.9e-2.
#
# SO THE WHOLE INVISCID PATH IS CLEARED -- volume operator, metrics, the
# "impose nothing" outflow, and flux differencing.
#
# RUNG 0 PASSED TOO: dp_rel came back at +/-5.8e-13 — machine
# zero — uniform over the whole domain including the outflow plane. So the
# metric terms are exact across the 15-degree kink and the two-block
# junction, the "impose nothing" outflow manufactures nothing, and the
# inviscid volume operator holds a constant state at nop = 4 for 500 steps.
#
# WHAT RUNG 0 CANNOT SEE, and it matters: a uniform field is invariant under
# ANY permutation of the nodes, so it cannot detect a halo-exchange or
# node-indexing bug. Rung 0 clears the operator and the metrics. It does not
# clear the parallel plumbing — for that, see THE MPI CHECK at the bottom.
#
# WHAT IS BEING HUNTED. In rampCaoEtAl2021_M7 the first negative pressure
# appeared 60 mm above the wall by step 200. The fastest signal in this flow
# is u + c = 1949 m/s, which needs 3.1e-5 s — step 31,000 — to cross 60 mm.
# That node was corrupted 150x earlier than anything could have reached it.
# So from rung 2 on, when a real wall makes dp_rel legitimately large, the
# question is NEVER how big it is. It is WHERE IT IS. A disturbance that
# appears at the top of the domain inside 500 steps did not travel there,
# and whatever put it there is the bug.
#
# READ THE LOCATION, NOT THE MAGNITUDE.
#---------------------------------------------------------------------------------
function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # RUNG 0 — the base deck. Everything that could confound the
        # measurement is off, and each "off" is load bearing:
        #
        #   :lvisc       => false   no molecular viscosity, no DynSGS. Either
        #                           would DAMP the error being measured.
        #   :lpositivity => false   MUST stay false through rung 5. The repair
        #                           exists to hide negative pressure; here
        #                           negative pressure is the signal.
        #   :lfilter     => false   same argument.
        #   :lsource     => false   nothing to manufacture a source.
        #   :lkep        => false   baseline is the production inviscid
        #                           operator; rung 1 turns it on.
        #
        # Mesh, :nop, state and :Δt are the production ones unchanged, so a
        # result here transfers directly to rampCaoEtAl2021_M7.
        #---------------------------------------------------------------------------
        :ode_solver           => CarpenterKennedy2N54(),
        :tinit                => 0.0,
        # 500 steps at the production Δt. The first positivity repair in the
        # production run happened within the first 1000 RHS calls, i.e. inside
        # step 200, so this is twice the window in which the real failure
        # first shows.
        :tend                 => 5.0e-7,
        :Δt                   => 1.0e-9,
        :diagnostics_at_times => (0:5.0e-8:5.0e-7),   # every 50 steps, 10 frames
        :lrestart             => false,
        :restart_time         => 0.0,
        :lsource              => false,
        :SOL_VARS_TYPE        => TOTAL(),
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        :energy_equation      => "energy",
        :lvisc                => false,
        :lfilter              => false,
        :lpositivity          => false,
        :volume_flux          => ranocha(),
        :lkep                 => true,        
        #---------------------------------------------------------------------------
        # The PRODUCTION mesh, read from the production directory. Swap for
        # ramp15.msh to ask the same question of the wall-stretched grid.
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/CompEuler/rampCaoEtAl2021/ramp15_uniform.msh",
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        :output_dir           => "./output",
        :loutput_pert         => false,
        :lschlieren           => false,
        #---------------------------------------------------------------------------
        :linitial_refine      => false,
        :ladapt               => false,
        :lxy_partition        => false,
    ) #Dict

    # RUNG 1 — ENTROPY-CONSERVATIVE FLUX DIFFERENCING
    #
    # rampCaoEtAl2021_M7 runs with this ON and rung 0 ran with it OFF, so it
    # is the largest untested difference between the test and production.
    #
    # WHY IT COULD BREAK FREE-STREAM PRESERVATION. Flux differencing does not
    # inherit FSP from the pointwise scheme. The two-point form is built on a
    # summation-by-parts operator and preserves a constant state only if the
    # discrete metric terms satisfy the metric identities in the SAME form the
    # flux differencing uses them. Rung 0 proved the metrics are exact for the
    # standard operator; it proves nothing about this one.
    #
    # EXPECTED: dp_rel stays at 1e-13. If it does not, flux differencing
    # manufactures a source on this mesh at Mach 7.7, and since the internal
    # energy is only 5.7% of the total energy here, the amplification into
    # the pressure is 17.5x. That would be the bug, and it would close the
    # campaign.
    #
    # user_flux.jl already carries user_fluxaux! and flux_turbo for all three
    # volume fluxes, so nothing else has to change. central_euler() is the
    # useful third data point: it is the plain central flux written in
    # flux-differencing FORM, so if ranocha() breaks FSP and central_euler()
    # does not, the entropy machinery is at fault rather than the form.
    # DONE -- these two are now set directly in the base Dict above, so this
    # block is spent. Kept only to name the alternatives: kennedy_gruber() is
    # kinetic-energy preserving only, and central_euler() is the plain central
    # flux written in flux-differencing FORM, so if ranocha() had broken FSP
    # and central_euler() had not, the entropy machinery would have been at
    # fault rather than the form.

    # RUNG 2 — THE NO-SLIP ISOTHERMAL WALL          (in user_bc.jl, not here)
    #
    # Set FSP_WALL = :noslip in user_bc.jl.
    #
    # From here on the uniform stream is NOT an exact solution — 1725 m/s
    # standing on a no-slip 293 K wall is a genuine discontinuity — so dp_rel
    # WILL be large near the wall and that is correct behaviour, not a
    # failure. THE MEASUREMENT CHANGES: stop reading the range, start reading
    # the map. In 500 steps the fastest signal travels 1949 m/s * 5e-7 s =
    # 0.97 mm. Everything above about 1 mm from the wall must still be at
    # 1e-13. If the top of the domain lights up, the disturbance did not
    # travel there and rung 2 has found the bug.
    #
    # Useful sub-test: the production case starts from a Pohlhausen boundary
    # layer precisely BECAUSE a uniform stream on a no-slip wall is illegal.
    # Rung 2 deliberately keeps the illegal start, because an illegal start
    # that stays LOCAL is still a passing result for the question being asked.


    # RUNG 3 — SUTHERLAND MOLECULAR VISCOSITY
    #
    # Real Navier-Stokes viscous terms, still with no artificial viscosity.
    # The specific thing worth watching is mu(T) ~ T^1.5: T here is the
    # specific internal energy, and if it ever goes negative that exponent is
    # a DomainError or a NaN, not a small error. Rung 3 is where that would
    # first appear.
    #
    # Requires RUNG 2 (a viscous run with no wall has nothing to be viscous
    # about).

    # inputs[:lvisc]            = true
    # inputs[:lsutherland]      = true
    # inputs[:sutherland_muref] = 1.716e-5
    # inputs[:sutherland_Tref]  = 273.15
    # inputs[:sutherland_S]     = 110.4
    # inputs[:Pr_lam]           = 0.71


    # RUNG 4 — DYNSGS WITH DOMAIN NORMS          (the production stabilisation)
    #
    # THIS IS THE RUNG THE CAUSALITY ARGUMENT POINTS AT. :dsgs_norms =>
    # "domain" normalises the residual by norms taken over the WHOLE domain,
    # which is an MPI.Allreduce. That is a LITERAL non-local channel: it is
    # the one mechanism in this deck that can carry a number from the
    # boundary layer to a node 60 mm away in zero time, which is exactly what
    # the production failure requires and what nothing physical can supply.
    #
    # RUN THIS ONLY WITH RUNG 2 ACTIVE. With the wall removed the field stays
    # identically equal to qe, so the departure norm in the denominator is
    # exactly 0 and you are measuring a 0/0, not the scheme.
    #
    # Settings are rampCaoEtAl2021's own, unchanged.

    # inputs[:lvisc]        = true              # required by :visc_model
    # inputs[:visc_model]   = DSGS()
    # inputs[:dsgs_sensor]  = "legacy"
    # inputs[:ldsgs_nodal]  = false
    # inputs[:dsgs_norms]   = "domain"
    # inputs[:μ]            = [1.0, 1.0, 1.0, 1.0]
    # inputs[:dsgs_Cmax]    = 0.1
    # inputs[:Pr]           = 0.1

    # RUNG 5 — THE SAME, WITH RANK NORMS                   (the discriminator)
    #
    # Uncomment on top of RUNG 4. "rank" reduces over each MPI rank instead of
    # the whole domain, so the non-local reach of the normalisation shrinks
    # from the domain to the partition.
    #
    # If rung 4 corrupts the far field and rung 5 corrupts it DIFFERENTLY —
    # different location, different time, or not at all — the normalisation
    # is carrying the corruption and :dsgs_norms is the bug.
    # If rungs 4 and 5 fail identically, the norm is innocent.
    #
    # NOTE that "rank" makes the answer partition-dependent BY DESIGN. That is
    # a defect in production and a feature here: it is the variable being
    # changed.
    # inputs[:dsgs_norms] = "rank"

    # OPTIONAL ADD-ON — POSITIVITY, FOR THE COORDINATE ONLY
    #
    # Turn this on ONLY as a second run of a rung that already failed, never
    # as part of the measurement: the repair alters the state and contaminates
    # dp_rel. What it buys is the single most useful diagnostic this campaign
    # has produced — the (x, y) of the FIRST negative pressure, printed with
    # the mass and energy it injected. That is the number that decoded the
    # production failure to the outflow plane in one line.
    #
    # Floors are 1e-6 of this free stream, as in rampCaoEtAl2021_M7.

    # inputs[:lpositivity]            = true
    # inputs[:positivity_rho_min]     = 2.1e-8
    # inputs[:positivity_p_min]       = 7.6e-4
    # inputs[:positivity_report]      = true
    # inputs[:positivity_report_every] = 200    # tighter than production: this
    #                                           # run is 2500 RHS calls total


    # THE MPI CHECK — no code change, and it is not optional
    #
    # Rung 0 passed, but a uniform field is invariant under any permutation of
    # the nodes, so it CANNOT detect a halo-exchange or node-indexing bug.
    # That check has to be done on a non-uniform field, which means on the
    # production case:
    #
    #   Run rampCaoEtAl2021_M7 on 1 (or 2) ranks and on 32, and compare
    #     - the printed CFL numbers at the same t
    #     - the positivity counts and its first (x, y)
    #     - the step at which it dies
    #
    # :dsgs_norms => "domain" exists precisely so the answer does NOT depend
    # on the rank count. If those numbers differ, the parallel path is
    # implicated and every Mach-7 result in this campaign is suspect.

    return inputs

end
