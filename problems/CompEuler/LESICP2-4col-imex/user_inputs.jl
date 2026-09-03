#=============================================================================
 LESICP2-4col-imex -- FOUR ELEMENT COLUMNS, ONE PER RANK. A LAPTOP TEST.

 WHAT IT IS FOR. This deck exists to exercise the column partition, the IMEX3D
 stage solve and the MOST/SGS surface machinery on four MPI ranks on a laptop,
 in seconds rather than hours. It is NOT a physics case: 2 x 2 elements over
 1280 x 1280 m is four columns of a large-eddy simulation, which resolves
 nothing. Do not read a profile off it.

 THE GRID, AND WHY IT IS EXACTLY THIS SHAPE
 ------------------------------------------
   2 x 2 x 60 elements at p = 4 over 1280 x 1280 x 5000 m
   element 640 x 640 m in the horizontal -> 160 m effective resolution
   smallest LGL gaps  h_x = h_y = 110.5 m   h_z = 6.907 m   (anisotropy 16:1)
   19 521 gridpoints total, 4 880 per rank on 4 ranks -- EXACTLY EQUAL

 :lxy_partition cuts the mesh into VERTICAL COLUMNS, so the usable parallelism
 is bounded by nelemx*nelemy and a rank count that does not divide the column
 grid leaves ranks owning zero elements (fatal). Four columns in a 2 x 2 grid
 divides exactly one way that uses four ranks:

     julia tools/pick_nranks.jl 2 2 60 4 4

     ranks  rank grid   block     cols/rank  pts/rank   halo/elem
     1      1 x 1       2 x 2     4          19521      2.0        <== RECOMMENDED
     2      1 x 2       2 x 1     2          9760       3.0
     4      2 x 2       1 x 1     1          4880       4.0

 pick_nranks recommends 1 rank because it optimises for EFFICIENCY -- it wants
 >= 50 000 points/rank and this whole mesh is 19 521, so every parallel row is
 marked "thin: comms-bound". That is the right answer for a production run and
 the wrong one here: this deck is a correctness test, and 4 ranks x 1 column
 x 4 880 points is a perfectly balanced decomposition. Run it on 4.

     mpiexec -n 4 julia --project=. src/Jexpresso.jl CompEuler LESICP2-4col-imex

 The vertical grid, the sounding, the surface flux, the closure and the filter
 are the production LESICP2-64x64x60-imex values, so the code paths taken are
 the same ones; only the horizontal extent, the step size and tend differ.

     DBG_TEND=...    model seconds (default 60 -- a smoke test, not a run)
     DBG_SCHEME=     imex (default) | hevi | explicit
     DBG_VDIFF=0     explicit vertical diffusion + the scalar Schur arm
     DBG_DT=...      step size, if you are sweeping it
=============================================================================#

function user_inputs()

    # Implicit vertical SGS diffusion, on by default, as in the production
    # deck: it is what removes nu_eff/dz^2 from the explicit budget. Read into
    # a local because :imex_linearization has to agree with it -- under a
    # dynamic closure the implicit operator carries no diffusion at all unless
    # the linearisation is :PS.
    _vdiff = parse(Bool, get(ENV, "DBG_VDIFF", "true"))

    use_imex  = true
    use_schur = !_vdiff   # IMPOSSIBLE with :implicit_vdiff -- the reduction cannot see the diffusion operator

    lschur      = parse(Bool,    get(ENV, "DBG_SCHUR",     string(use_schur)))
    rtol        = parse(Float64, get(ENV, "DBG_RTOL",      "1.0e-6"))
    restart     = parse(Int,     get(ENV, "DBG_RESTART",   "30"))
    maxiter     = parse(Int,     get(ENV, "DBG_MAXITER",   "600"))
    precond     = Symbol(        get(ENV, "DBG_PRECOND",   "column"))   # or :none
    umax        = parse(Float64, get(ENV, "DBG_UMAX",      "15.0"))
    verify      = parse(Bool,    get(ENV, "DBG_VERIFY",    "true"))
    warm_start  = parse(Bool,    get(ENV, "DBG_WARM",      "true"))
    lin_imex    = Symbol(        get(ENV, "DBG_LIN", _vdiff ? "PS" : "RS"))
    updfreq     = parse(Int,     get(ENV, "DBG_UPDFREQ",   "5"))
    lat_walls   = Symbol(        get(ENV, "DBG_LATWALL",   "auto"))
    monitor     = parse(Bool,    get(ENV, "DBG_IMEXMON",   "true"))
    monitor_every = parse(Int,   get(ENV, "DBG_IMEXMONEVERY", "20"))

    scheme = Symbol(get(ENV, "DBG_SCHEME", use_imex ? "imex" : "explicit"))
    scheme in (:imex, :hevi, :explicit) ||
        error("LESICP2-4col-imex: DBG_SCHEME must be imex, hevi or explicit; got $scheme")

    # STEP SIZE -- 0.5 s, THE PRODUCTION DECK'S VALUE, DELIBERATELY.
    #
    # This mesh could take more. The run's own CFL report at t = 0 gives
    #
    #     node spacing   h_x = h_y = 110.5 m   h_z = 6.907 m   anisotropy 16:1
    #     wedge neutral up to Δt = 0.8038 s
    #
    # against the 0.506 s LESICP2-64x64x60-imex derives, and the difference is
    # h_x, not h_z: the vertical grid is now identical (6.91 m both), but this
    # deck's elements are 640 m wide against production's 160 m. With ALL
    # acoustics implicit the step is set by the EXPLICIT half -- advection,
    # |u|/h_x -- so 4x the horizontal spacing is 4x less rate and a higher
    # limit. 0.5 s is 62% of it, comfortably inside.
    #
    # It is pinned to production's value anyway, because that is the whole
    # point of this deck: same tableau, same stage count, same Krylov work per
    # step as the case being smoke-tested is worth more than 10% of Δt. If you
    # raise it, raise it to measure something, and read the report.
    #
    # Do not copy LESICP2-coarse-imex's 0.9 s across: that is 112% of the
    # limit here.
    #
    # BOTH of those numbers move with :first_zelement_size, so re-read the
    # report if you touch the stretch -- see the note under :lstretch below.
    # With 40.0 there the grid inverts, h_z drops to 3.301 m, the neutral limit
    # falls to 0.4307 s, and at 0.9 s ark_wedge_amplification warns "amplifies
    # 4364% per step": the run reached t = 4.5 s and died in
    # perfectGasLaw_ρθtoP with a DomainError on a negative ρθ.
    # :imex_stability_guard => :error turns that warning into a refusal to
    # start, which is worth setting if you are sweeping Δt.
    Δt_default = scheme === :imex ? "0.5" : (scheme === :hevi ? "0.055" : "0.0033")
    Δt = parse(Float64, get(ENV, "DBG_DT", Δt_default))

    # VTK RESTART, opt-in. read_vtk_restart! (called from this case's own
    # initialize.jl) overwrites qn from a previous snapshot; qe stays as the
    # sounding built it. :restart_vtk_iout and :tinit are auto-detected from
    # the LAST entry in simulation.pvd, so this resumes from the newest dump
    # without an edit. SAME RANK COUNT AND SAME GRID -- each rank reads its own
    # piece iter_N/iter_N_<rank+1>.vtu and the npoin check is fatal otherwise.
    _restart_vtk = parse(Bool, get(ENV, "DBG_RESTART_VTK", "false"))
    restart_dir  = get(ENV, "DBG_RESTART_DIR",
                       "./output/4col-src/CompEuler/LESICP2-4col-imex/output")

    # 60 s, i.e. ~67 IMEX steps. This is a smoke test; DBG_TEND is honoured
    # (the production deck currently hardcodes :tend and silently ignores it).
    tend = parse(Float64, get(ENV, "DBG_TEND", "60.0"))

    inputs = Dict(
        #---------------------------------------------------------------------------
        # Time integration
        #---------------------------------------------------------------------------
        :ode_solver           => scheme === :imex ? IMEX_ARK(:ARS343)      :
                                 scheme === :hevi ? HEVI_ARK(:ARS232)      :
                                                    CarpenterKennedy2N54(),
        :implicit_vdiff       => _vdiff,
        :lcfl_report          => true,
        :lcfl_report_every    => parse(Int, get(ENV, "DBG_CFL_EVERY", "20")),
        :lstep_heartbeat      => parse(Bool, get(ENV, "DBG_HEARTBEAT", "true")),
        :hevi_verify          => verify,
        :hevi_linearization   => Symbol(get(ENV, "DBG_LIN", _vdiff ? "PS" : "RS")),
        :hevi_update_freq     => updfreq,
        :hevi_wall_flux       => true,
        #--- IMEX3D ----------------------------------------------------------------
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
        :lrestart_vtk          => _restart_vtk,
        :restart_vtk_input_dir => restart_dir,
        # A BARE RANGE, not a tuple. Every range inside a tuple needs its own
        # `...` and the tuple needs a trailing comma -- both have killed the
        # production deck at startup. One range needs neither.
        :diagnostics_at_times => (0.0:10.0:tend),
        :lsource              => true,
        :sounding_file        => "./data_files/input_sounding_teamx_u10_flat_noheader.dat",
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # Physical parameters -- the production LESICP2-64x64x60-imex values
        #---------------------------------------------------------------------------
        :user_heatflux        => 0.12,
        :lxy_partition        => true,   # REQUIRED: this is what makes it column-based
        :lwall_model          => true,
        :ifirst_wall_node_index => 2,    # must be 2 <= i <= nop+1
        :bdy_fluxes           => true,
        :lvisc                => true,
        :visc_model           => SMAG(),
        :C_s                  => 0.18,
        :lrichardson          => true,
        :lwall_damping        => true,
        :μ                    => [0.0, 1.0, 1.0, 1.0, 1.0],
        :les_filter_width     => :geometric,
        # MOST guard rails -- these ARE the defaults; written out because this
        # is the deck you reach for when one of them is suspected. See
        # docs/boundary_conditions.md 2.2.1 and test/physics/test_most_guards.jl.
        :most_u_min           => 0.1,     # m/s
        :most_zeta_min        => -5.0,    # sign is the sign of L: < 0 is UNSTABLE
        :most_zeta_max        =>  2.0,
        #---------------------------------------------------------------------------
        # LES statistics. A bare range again, and it starts after the first
        # third so a 60 s run still writes some.
        #---------------------------------------------------------------------------
        :statistics_time      => (20.0:10.0:tend),
        :lesprofile_vars      => ["u_mean", "v_mean", "w_mean", "t_mean", "p_mean"],
        # ALL 38, AND THE COUNT IS LOAD-BEARING. This list sizes the stress
        # cache, while user_les_profiles! in user_primitives.jl writes prof[1]
        # through prof[38] unconditionally. Shortening it to the 26 entries
        # before the triple moments is a BoundsError at the first statistics
        # time -- "attempt to access 26-element Vector{Float64} at index [27]",
        # thrown from a callback 20 s into the run, not at startup.
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
        # Mesh. Regenerate with:
        #   gmsh -3 problems/CompEuler/LESICP2-4col-imex/LESICP_2x2x60.geo \
        #        -o problems/CompEuler/LESICP2-4col-imex/LESICP_2x2x60_1280mX1280mX5000m.msh
        #---------------------------------------------------------------------------
        :lread_gmsh       => true,
        :gmsh_filename    => "./problems/CompEuler/LESICP2-4col-imex/LESICP_2x2x60_1280mX1280mX5000m.msh",
        # The .geo is UNIFORM in z; the grading is applied here at read time.
        # :first_zelement_size IS THE EFFECTIVE RESOLUTION -- stretching.jl
        # multiplies it by (ngl-1), so 40.0 means a 160 m first ELEMENT.
        # Identical to LESICP2-coarse-imex, which is where the step size above
        # comes from.
        # :first_zelement_size IS THE EFFECTIVE RESOLUTION, NOT THE ELEMENT
        # SIZE -- stretching.jl does `first_cell_size *= (mesh.ngl-1)`, so 10.0
        # means a 40 m first ELEMENT and 10 m effective resolution. This is the
        # production LESICP2-64x64x60-imex value and it must stay 10.0.
        #
        # 40.0 (which LESICP2-coarse-imex carries) INVERTS THE GRID on this
        # mesh, and the deck ran that way once. "two_block uniformish" compares
        # the scaled first cell against the mesh's own uniform element:
        #
        #   value  first element  squash_ratio  below 2000 m  top block   h_z
        #   10.0        40 m          0.48       50 elements    404 m    6.91 m
        #   40.0       160 m          1.92       12 elements     19 m    3.30 m
        #
        # squash_ratio > 1 means the "high resolution" block is COARSER than
        # the 83.3 m mesh it is built from, so the leftover elements pile up at
        # the LID: 19 m at 5000 m, 160 m at the surface, i.e. finest where the
        # ABL is not and coarsest where the surface flux goes in. It is not
        # fatal and nothing warns (warn_if_stretch_inverts guards the other
        # stretch type), it just shows up as a 33.5:1 anisotropy and a step
        # size a third of what the grid should allow.
        :lstretch => true,
        :stretch_factor => 1.15,
        :stretch_type => "two_block uniformish",
        :first_zelement_size => 10.0,
        :zlevel_transition => 2000.0,
        #---------------------------------------------------------------------------
        # Filter -- the production imex deck's settings
        #---------------------------------------------------------------------------
        :lfilter             => parse(Bool, get(ENV, "DBG_FILTER", "true")),
        :mu_x                => 0.05,
        :mu_y                => 0.05,
        :mu_z                => 0.1,
        :filter_type         => "erf",
        #---------------------------------------------------------------------------
        # Output -- LOCAL, not /scratch. One tree per scheme.
        #---------------------------------------------------------------------------
        :outformat           => "vtk",
        :output_dir          => "./output/LESICP2-4col-imex_" * String(scheme) * "/",
        :loverwrite_output   => true,
        :lwrite_initial      => parse(Bool, get(ENV, "DBG_WRITE_INITIAL", "true")),
        #---------------------------------------------------------------------------
        # AMR -- off
        #---------------------------------------------------------------------------
        :linitial_refine     => false,
        :init_refine_lvl     => 1,
        :ladapt              => false,
        :amr                 => false,
        :amr_freq            => 20,
        :amr_max_level       => 1,
    ) #Dict

    return inputs

end
