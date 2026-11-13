#---------------------------------------------------------------------------------
# fluxEmergenceLeake2010 — the 2.5D flux-emergence tests of
#
#   J. E. Leake, M. G. Linton & S. K. Antiochos, "Tests of dynamical flux
#   emergence as a mechanism for coronal mass ejection initiation",
#   ApJ 722:550 (2010), with the corrected resistivity of its Erratum,
#   ApJ 741:125 (2011),
#
# on the DynSGS configuration of problems/MHD/fluxEmergenceSon2025DSGS (the
# reference-weighted, conserved-form Marras-Nazarov operator, no limiter), plus
# the paper's explicit viscosity and anomalous resistivity. See README.md.
#
#   julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2010
#   mpiexec -n 64 julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2010
#
# The default is the erratum's Figure 3: the low-tension tube (tube 2) emerging
# into the medium quadrupole ("98 G"), to t = 2530 s, on LK10_ref2.msh. Other
# simulations of the paper are selected without editing this file:
#
#   JEXPRESSO_LK10_TUBE      1 | 2 | 3 (0: no tube, the background alone)        [2]
#   JEXPRESSO_LK10_FIELD     quadrupole | dipole | none                          [quadrupole]
#   JEXPRESSO_LK10_STRENGTH  strong | medium | weak: quadrupole q_d = 11.7, 5.85,
#                            3.9e7 G cm ("197, 98, 66 G"), dipole |B(0,0)| = 104,
#                            52, 26 G                                             [medium]
#   JEXPRESSO_LK10_REF       mesh 1 | 2 | 3 (LK10.geo; 3 must be generated)       [2]
#   JEXPRESSO_LK10_TEND      final time in SECONDS                               [2530]
#   JEXPRESSO_LK10_DOUT      snapshot interval in SECONDS (plus the figure times) [116 = 5 t0]
#   JEXPRESSO_LK10_DT        time step in t0 (default: c_h Δt/Δh = 0.25, below)
#   JEXPRESSO_LK10_MCZ       convection-zone polytropic index m (p ∝ T^(m+1))     [1.3]
#   JEXPRESSO_LK10_RELAX_YMIN energy relaxation only above this height, in L0    [-Inf]
#   JEXPRESSO_LK10_TUBE_EQ   thermal | relaxed: tube in thermal equilibrium (the paper),
#                            or on the relaxed background relation ρT₀(ρ) = p      [thermal]
#   JEXPRESSO_LK10_CMIN      DynSGS background floor C_min                        [0]
#   JEXPRESSO_LK10_FLOORS    1: positivity floors as the Runge-Kutta stage limiter [0]
#   JEXPRESSO_LK10_CURLB     0: DynSGS resistivity as the component Laplacian       [1: curl form]
#   JEXPRESSO_LK10_NU        dynamic | kinematic: the paper's dynamic ν below the
#                            photosphere (capped above), or kinematic everywhere   [dynamic]
#
# Times are in t0 = 23.19 s inside the code (user_flux.jl).
#---------------------------------------------------------------------------------
_lk10_env(k, d) = strip(get(ENV, "JEXPRESSO_LK10_" * k, d))
_lk10_num(k, d) = something(tryparse(Float64, _lk10_env(k, "")), d)

function _lk10_case()
    tube  = something(tryparse(Int, _lk10_env("TUBE", "")), 2)
    field = lowercase(_lk10_env("FIELD", "quadrupole"))
    str   = lowercase(_lk10_env("STRENGTH", "medium"))
    ref   = something(tryparse(Int, _lk10_env("REF", "")), 2)
    istr  = Dict("strong" => 1, "medium" => 2, "weak" => 3)
    haskey(istr, str) || error(" problems/MHD/fluxEmergenceLeake2010: JEXPRESSO_LK10_STRENGTH = $(str), expected strong, medium or weak.")
    ref in (1, 2, 3) || error(" problems/MHD/fluxEmergenceLeake2010: JEXPRESSO_LK10_REF = $(ref), expected 1, 2 or 3.")
    qd    = (11.7e7, 5.85e7, 3.9e7)[istr[str]]        # G cm (Sec. 2.2.3)
    Bdip  = (104.0, 52.0, 26.0)[istr[str]]            # G at (0, 0) (Fig. 3(b), Secs. 3.4-3.5)
    # Fastest wave of the initial state (code units; printed again by initialize.jl).
    # Field-free: the coronal sound speed is 15.9, but the emerged field reaches
    # several times that, hence the 60.
    cmax  = field == "quadrupole" ? (548.7, 274.7, 183.5)[istr[str]] :
            field == "dipole"     ? (258.2, 129.7, 66.2)[istr[str]]  : 60.0
    return (; tube, field, str, ref, qd, Bdip, cmax)
end

#---------------------------------------------------------------------------------
# Optional positivity floors (JEXPRESSO_LK10_FLOORS=1), as in
# problems/MHD/fluxEmergenceLeake2014: ρ ≥ f ρ_e(y) and p ≥ f p_e(y) of the
# reference atmosphere, p raised through E at fixed velocity and field. Not
# conservative where they fire; every hit is counted and reported by rank 0 each
# time the total passes a power of ten. Off by default: the parent case runs on
# DynSGS alone.
#---------------------------------------------------------------------------------
if !@isdefined(lk10_floor_frac)
    const lk10_floor_frac = Ref{Float64}(0.2)
    const lk10_floor_hits = Ref{Int}(0)
    const lk10_floor_next = Ref{Int}(1)
end

function lk10_positivity_limiter!(u, integrator, p, t)
    npoin = p.mesh.npoin
    qe    = p.qp.qe
    frac  = lk10_floor_frac[]
    γm1   = γ_mhd - 1.0
    nhit  = 0
    @inbounds for ip = 1:npoin
        iρ = ip; iρu = npoin + ip; iρv = 2npoin + ip; iE = 3npoin + ip; iρw = 4npoin + ip
        iBx = 5npoin + ip; iBy = 6npoin + ip; iBz = 7npoin + ip; iψ = 8npoin + ip
        ρmin = frac*qe[ip,1]
        pmin = frac*qe[ip,end]
        ρ = u[iρ]
        if ρ < ρmin
            u[iρ] = ρmin; ρ = ρmin; nhit += 1
        end
        ke = 0.5*(u[iρu]^2 + u[iρv]^2 + u[iρw]^2)/ρ
        me = 0.5*(u[iBx]^2 + u[iBy]^2 + u[iBz]^2) + 0.5*u[iψ]^2
        if γm1*(u[iE] - ke - me) < pmin
            u[iE] = pmin/γm1 + ke + me; nhit += 1
        end
    end
    comm = MPI.COMM_WORLD
    nall = MPI.Comm_size(comm) > 1 ? MPI.Allreduce(nhit, +, comm) : nhit
    lk10_floor_hits[] += nall
    if nall > 0 && lk10_floor_hits[] >= lk10_floor_next[]
        if MPI.Comm_rank(comm) == 0
            print(string(" # positivity floors: ", lk10_floor_hits[], " node hits so far on all ranks (",
                         nall, " in this stage) at t = ", round(t*LK10_T0, digits = 2), " s\n"))
            flush(stdout)
        end
        while lk10_floor_next[] <= lk10_floor_hits[]
            lk10_floor_next[] *= 10
        end
    end
    return nothing
end

# Explicit dissipation (Sec. 2.1), code units (user_flux.jl). The paper's
# ν = 2.6e4 g cm⁻¹ s⁻¹ is a DYNAMIC viscosity, 0.01326 in code units: below the
# photosphere it is used as such (kinematic ν/ρ_e(y), ρ_e = T^m of the
# polytrope), above it the kinematic value is held at its photospheric 0.01326
# (as a dynamic viscosity it would grow to 10⁶ in the corona; README).
# JEXPRESSO_LK10_NU=kinematic: 0.01326 everywhere. The anomalous resistivity
# is the :mhd_anomalous_eta pair below; no other resistivity (the erratum
# removed the local η̂ = 0.1 exp(−r̂²/0.66²) at the origin of the original runs).
lk10_nu_dynamic_cz() = !(lowercase(_lk10_env("NU", "dynamic")) in ("kinematic", "const", "constant"))
lk10_dissipation(x, y) = ((lk10_nu_dynamic_cz() && y < 0.0) ? LK10_NU/lk10_temperature(y)^lk10_mcz[] : LK10_NU, 0.0)

function user_inputs()

    c    = _lk10_case()
    mesh = ("LK10_ref1.msh", "LK10_ref2.msh", "LK10_ref3.msh")[c.ref]
    # c_h Δt/Δh = 0.25 on the smallest LGL gap of the mesh (0.1727 of the 2.5,
    # 1.25, 0.625 L0 core elements at :nop => 4); c_h is the fastest wave of the
    # initial state, at the top of the transition region (README, "Cost").
    Δh   = 0.17267*(2.5, 1.25, 0.625)[c.ref]
    Δt   = _lk10_num("DT", 0.25*Δh/c.cmax)
    tend = _lk10_num("TEND", 2530.0)/LK10_T0
    # Snapshots every 5 t0 (116 s) and at the times of the erratum's figures
    # (Figs. 1 and 3: 920, 1380, 1610, 1955, 2070, 2185, 2254, 2530 s)
    tfig = [920.0, 1380.0, 1610.0, 1955.0, 2070.0, 2185.0, 2254.0, 2530.0] ./ LK10_T0
    dout = _lk10_num("DOUT", 5.0*LK10_T0)/LK10_T0          # JEXPRESSO_LK10_DOUT: snapshot interval in s
    tout = sort(unique(round.(vcat(collect(0.0:dout:tend), filter(t -> t <= tend + 1e-9, tfig), [tend]), digits = 6)))
    floors = _lk10_env("FLOORS", "0") in ("1", "true", "yes", "on")

    inputs = Dict(
        #---------------------------------------------------------------------------
        # Time integration: the low-storage RK of the parent case, no limiter
        # (unless JEXPRESSO_LK10_FLOORS=1)
        #---------------------------------------------------------------------------
        :ode_solver           => floors ? CarpenterKennedy2N54(stage_limiter! = lk10_positivity_limiter!,
                                                               williamson_condition = false) :
                                          CarpenterKennedy2N54(),
        :Δt                   => Δt,
        :tinit                => 0.0,
        :tend                 => tend + Δt,   # one step past the last output (an output at tend itself is not written)
        :diagnostics_at_times => Tuple(tout),
        :conservation_every   => 100,
        :conservation_slots   => [1, 4],
        :restart_time         => 0.0,
        :lrestart             => false,
        :lsource              => true,       # gravity, energy relaxation, ψ damping, damping layers
        :lwell_balanced       => true,       # the reference (atmosphere + coronal field) is a discrete steady state
        :SOL_VARS_TYPE        => TOTAL(),
        :ode_adaptive_solver  => false,
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # The paper's problem (initialize.jl)
        #---------------------------------------------------------------------------
        :lk10_tube            => c.tube,
        :lk10_field           => c.field,
        :lk10_qd              => c.qd,
        :lk10_Bdip            => c.Bdip,
        :lk10_mcz             => _lk10_num("MCZ", 1.3),
        :lk10_relax           => true,       # Eqs. 5-6
        :lk10_tube_eq         => lowercase(_lk10_env("TUBE_EQ", "thermal")),   # tube in thermal equilibrium (the paper) or "relaxed"
        :lk10_relax_ymin      => _lk10_num("RELAX_YMIN", -Inf),   # relaxation only above this height (L0)
        :lk10_glm_alpha       => 1.0,
        :lk10_damp_xd         => 250.0,      # 37.5 Mm
        :lk10_damp_yd         => 500.0,      # 75 Mm
        :lk10_damp_rate       => 1.0,        # 1/t0 at the boundary
        :mhd_dissipation      => lk10_dissipation,
        :lk10_nu_mode         => lk10_nu_dynamic_cz() ? "dynamic below the photosphere, kinematic 0.01326 above" : "kinematic 0.01326 everywhere",
        :mhd_anomalous_eta    => (LK10_ETA0, LK10_JCRIT),
        #---------------------------------------------------------------------------
        # DynSGS: the parent case's reference-weighted conserved-form operator
        # (problems/MHD/fluxEmergenceSon2025DSGS), with γ = 5/3, fast-speed floors
        # and the thermal-energy residual of the low-β corona (as in
        # problems/MHD/fluxEmergenceLeake2014), no split energy flux, and no
        # background floor by default: C_min Δ c_f with c_f ~ 270 in the corona
        # would be a diffusivity of ~4 L0 v0 (README, "DynSGS").
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :μ                    => [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        :visc_model           => DSGS_MHD(),
        :dsgs_sensor          => "legacy",
        :dsgs_CR              => 1.0,
        :dsgs_Cmax            => 0.5,
        :dsgs_Cmin            => _lk10_num("CMIN", 0.0),
        :dsgs_gamma           => 5.0/3.0,
        :dsgs_Prt             => 1.0,
        :dsgs_nazarov_energy  => false,
        :dsgs_norms           => "element",
        :dsgs_local_rel       => 1.0,
        :dsgs_conserved       => true,
        :dsgs_ref_weight      => true,
        :dsgs_nodal_rho       => true,       # inactive with :dsgs_conserved
        :dsgs_fast_floors     => true,
        :dsgs_thermal_residual => true,
        # In-plane field diffused by the curl form −∇ₕ×(η j_z ẑ), which keeps the discrete ∇·B
        # (the component Laplacian does not: README, "Divergence of B"); JEXPRESSO_LK10_CURLB=0 for
        # the Laplacian of the parent case.
        :dsgs_curl_resistivity => !(_lk10_env("CURLB", "1") in ("0", "false", "no", "off")),
        :lrichardson          => false,
        :energy_equation      => "energy",
        :lkep                 => false,
        :entropy_variables    => false,
        #---------------------------------------------------------------------------
        # Mesh (LK10.geo): -45 ≤ x ≤ 45 Mm, -3 ≤ y ≤ 90 Mm, walls on all sides
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/MHD/fluxEmergenceLeake2010/" * mesh,
        :lfilter              => false,
        #---------------------------------------------------------------------------
        # Output: VTU (tools/plot_leake2010.py draws the paper's figures from it)
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :plot_dsgs            => true,
        :loverwrite_output    => false,
        :lwrite_initial       => true,
        :output_dir           => "./output/",
        :loutput_pert         => false,
        :linitial_refine      => false,
        :ladapt               => false,
    ) #Dict

    return inputs
end
