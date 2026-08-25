#===============================================================================
 CompEuler/rtb2d_schur -- 2D (slab) rising thermal bubble, 10 km x 10 km,
 20x1x80, nop 4. One element thick in y with a CYLINDRICAL bubble, because
 HEVI and IMEX3D are 3D-only by construction (README.md says why the guard is
 protective rather than lazy). The answer is exactly the 2D one.
 Same bubble, spacing, Delta t and Schur A/B as CompEuler/rtb3d_schur, on 1/16
 of the points -- a desktop run rather than a queue slot.

 EVERYTHING YOU NORMALLY CHANGE IS IN THE NUMBERED BLOCKS BELOW. The rationale
 -- why this mesh, what the A/B measures, what the numbers came out at, what to
 read in a profile -- is in README.md beside this file, not here.

 Run it with run_rtb2d.sh beside this file. NR defaults to 4 and NOT to your
 core count: the slab has 20 element columns and that is the rank ceiling.
===============================================================================#

function user_inputs()

    # Environment overrides, so the batch scripts can drive an A/B without
    # editing this file. The value written here is the default in every case.
    _b(k, d) = parse(Bool,    get(ENV, k, string(d)))
    _f(k, d) = parse(Float64, get(ENV, k, string(d)))
    _i(k, d) = parse(Int,     get(ENV, k, string(d)))

    #---------------------------------------------------------------------------
    # 1. TIME INTEGRATOR -- set exactly ONE of these true
    #---------------------------------------------------------------------------
    limex     = true      # IMEX-ARK  : ALL acoustics implicit (3D operator)
    lhevi     = false     # HEVI-ARK  : vertical acoustics implicit only
    lexplicit = false     # explicit RK, no implicit solve at all

    #---------------------------------------------------------------------------
    # 2. SCHUR REDUCTION -- ignored unless limex
    #---------------------------------------------------------------------------
    lschur = _b("DBG_SCHUR", true)   # solve the scalar Np system in P = β·θ
                                     # instead of the 5·Np one. 3.56x faster.

    #---------------------------------------------------------------------------
    # 3. RESOLUTION AND TIME
    #---------------------------------------------------------------------------
    mesh   = get(ENV, "DBG_MESH", "20x1x80")    # or "10x1x40" (smoke test)
    nop    = 4                                  # polynomial order
    Δt     = _f("DBG_DT", limex ? 0.6 : (lhevi ? 0.2 : 0.1))
    tend   = _f("DBG_TEND", 1000.0)             # 60.0 is enough for a timing run
    tinit  = 0.0

    #---------------------------------------------------------------------------
    # 4. IMPLICIT SOLVER -- ignored unless limex or lhevi
    #---------------------------------------------------------------------------
    precond       = :column     # :column | :none
    rtol          = _f("DBG_RTOL", 1.0e-8)
    restart       = _i("DBG_RESTART", 20)   # GMRES restart m
    maxiter       = 200
    warm_start    = true
    verify        = true        # setup self-check; off saves a startup solve
    umax          = 20.0        # largest expected |v|, for the stability report
    linearization = :RS         # :RS (frozen) | :PS (periodically refreshed)
    update_freq   = 5           # steps between refreshes, :PS only
    stab_guard    = :warn       # :warn | :error | :none
    implicit_vdiff = false      # vertical diffusion inside the implicit half
    lateral_walls  = :auto      # :auto | :bc | :box | :none
    monitor        = false      # per-solve Krylov iteration counts

    #---------------------------------------------------------------------------
    # 5. PHYSICS
    #---------------------------------------------------------------------------
    lvisc = true
    μ     = [0.0, 125.0, 125.0, 125.0, 125.0]   # per equation; ρ is inviscid

    #---------------------------------------------------------------------------
    # 6. OUTPUT
    #---------------------------------------------------------------------------
    # Plot dθ, NOT θ: θ is the total, 300 K of background with a 2 K bubble on
    # it, so it autoscales to the stratification and the bubble disappears.
    outdir    = lexplicit ? "./output_explicit/" :
                lhevi     ? "./output_hevi/"     :
                lschur    ? "./output_imex_schur/" : "./output_imex_full/"
    out_every = 100.0           # model seconds between snapshots
    heartbeat = true            # per-step wall time, JIT excluded

    #===========================================================================
     Nothing below here normally needs editing.
    ===========================================================================#

    count(x -> x, (limex, lhevi, lexplicit)) == 1 ||
        error("rtb2d_schur: set exactly ONE of limex / lhevi / lexplicit true.")

    # The switches in block 1 select the solver through :ode_solver, which is
    # what imex3d_enabled / hevi_enabled read. (:limex and :lhevi are accepted
    # as separate force-on keys; they build the solver but do not choose the
    # integrator, so setting one without a matching :ode_solver gets you the
    # setup cost of a scheme you are not running.)
    ode_solver = limex ? IMEX_ARK(:ARS343) :
                 lhevi ? HEVI_ARK(:ARS232) : CarpenterKennedy2N54()

    inputs = Dict(
        :ode_solver           => ode_solver,
        :Δt                   => Δt,
        :tinit                => tinit,
        :tend                 => tend,
        :diagnostics_at_times => (0:out_every:tend),
        :lstep_heartbeat      => heartbeat,
        :lsource              => true,
        :restart_time         => 1.0e7,
        :lrestart             => false,

        # -- IMEX (block 2 and 4) --
        :imex_schur           => lschur,
        # The fast scalar matvec for H, on. DBG_SCHUR_KERN=0 reaches the
        # reference form -- two full five-field applies, ~6x slower, kept only
        # as the independent statement of the same operator to debug against.
        # There is no reason to run it otherwise: they agree to 1.9e-16.
        :imex_schur_kernel    => _b("DBG_SCHUR_KERN", true),
        :implicit_vdiff       => implicit_vdiff,
        :imex_verify          => verify,
        :imex_warm_start      => warm_start,
        :imex_rtol            => rtol,
        :imex_restart         => restart,
        :imex_maxiter         => maxiter,
        :imex_precond         => precond,
        :imex_umax            => umax,
        :imex_lateral_walls   => lateral_walls,
        :imex_wall_flux       => true,
        :imex_linearization   => linearization,
        :imex_update_freq     => update_freq,
        :imex_stability_guard => stab_guard,
        :imex_monitor         => monitor,

        # -- HEVI (block 4, when lhevi) --
        :hevi_verify          => verify,
        :hevi_linearization   => linearization,
        :hevi_update_freq     => update_freq,
        :hevi_wall_flux       => true,

        :lcfl_report          => true,
        :lcfl_report_every    => 0,

        # -- discretisation --
        :interpolation_nodes  => "lgl",
        :nop                  => nop,

        # -- physics --
        :lvisc                => lvisc,
        :visc_model           => AV(),
        :μ                    => μ,
        :energy_equation      => "theta",

        # -- mesh. Rebuild with generate_mesh.jl; :lxy_partition makes the
        #    decomposition columnar, which IMEX3D requires.
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/CompEuler/rtb2d_schur/rtb2d_$(mesh).msh",
        :lxy_partition        => true,
        :lwarp                => false,
        :lstretch             => false,

        # -- output --
        :outformat            => "vtk",
        :output_dir           => outdir,
        :loverwrite_output    => true,
        :lwrite_initial       => true,

        # -- AMR: off. IMEX3D refuses it; the column topology, the factorised
        #    preconditioner and the gather/scatter plan are all built once.
        :linitial_refine      => false,
        :init_refine_lvl      => 1,
        :ladapt               => false,
    )
    return inputs
end
