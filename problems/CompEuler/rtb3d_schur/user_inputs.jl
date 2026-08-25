#===============================================================================
 CompEuler/rtb3d_schur -- 3D rising thermal bubble, 10 km cube, 20x20x80, nop 4.
 Built to measure the scalar Schur stage solve against the five-field one.

 EVERYTHING YOU NORMALLY CHANGE IS IN THE NUMBERED BLOCKS BELOW. The rationale
 -- why this mesh, what the A/B measures, what the numbers came out at, what to
 read in a profile -- is in README.md beside this file, not here.

 The 2D desktop analogue of this case is CompEuler/rtb2d_schur.
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
    lschur        = _b("DBG_SCHUR", false)      # scalar Np system in P = β·θ,
                                                # instead of the 5·Np one
    lschur_kernel = _b("DBG_SCHUR_KERN", true)  # fast matvec for H;
                                                # false = reference form (~6x slower)

    #---------------------------------------------------------------------------
    # 3. RESOLUTION AND TIME
    #---------------------------------------------------------------------------
    mesh   = get(ENV, "DBG_MESH", "20x20x80")   # or "10x10x40" (fits one rank)
    nop    = 4                                  # polynomial order
    Δt     = _f("DBG_DT", limex ? 0.6 : (lhevi ? 0.2 : 0.1))
    tend   = _f("DBG_TEND", 1000.0)             # 45.0 is enough for a timing run
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
        error("rtb3d_schur: set exactly ONE of limex / lhevi / lexplicit true.")

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
        :imex_schur_kernel    => lschur_kernel,
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
        :gmsh_filename        => "./problems/CompEuler/rtb3d_schur/rtb3d_$(mesh).msh",
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
