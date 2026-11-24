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
#   JEXPRESSO_LK10_REF       mesh 1 | 2 | 3 | 4 (LK10.geo; 3 must be generated; 4 is
#                            2's core spacing over the whole figure region)      [2]
#   JEXPRESSO_LK10_TEND      final time in SECONDS                               [2530]
#   JEXPRESSO_LK10_DOUT      snapshot interval in SECONDS (plus the figure times) [116 = 5 t0]
#   JEXPRESSO_LK10_DT        time step in t0 (default: c_h Δt/Δh = 0.25, below)
#   JEXPRESSO_LK10_MCZ       convection-zone polytropic index m (p ∝ T^(m+1))     [1.3]
#   JEXPRESSO_LK10_RELAX_YMIN energy relaxation only above this height, in L0    [-Inf]
#   JEXPRESSO_LK10_TUBE_EQ   thermal | relaxed: tube in thermal equilibrium (the paper),
#                            or on the relaxed background relation ρT₀(ρ) = p      [thermal]
#   JEXPRESSO_LK10_CMIN      DynSGS background floor C_min                        [0]
#   JEXPRESSO_LK10_FLOORS    0: no positivity floors (stage limiter)              [1]
#   JEXPRESSO_LK10_TFLOOR    temperature floor, in T_ph = 5700 K                 [0.1]
#   JEXPRESSO_LK10_ETA       induction | all: the DynSGS resistivity from the residuals
#                            of the induction equation alone, or from the max over
#                            all equations (one coefficient, the parent case)    [induction]
#   JEXPRESSO_LK10_RHOFAC    factor on the DynSGS mass diffusion (the :μ multiplier of the ρ
#                            slot), a test of the plasma it carries across the field  [1]
#   JEXPRESSO_LK10_RHOSENS   off | on | a0,a1,C: the mass slot of an element whose density has a
#                            node-to-node part (the highest Legendre mode of log10(ρ/ρ_e), from
#                            a0 to a1 dex) gets at least C Δ |(u, v)|, the first-order diffusivity
#                            of a contact carried by the flow; on = 0.03,0.15,1 (README, "The
#                            density sensor")                                   [off]
#   JEXPRESSO_LK10_HEAT      what the heat the DynSGS dissipates does: "shock": only the
#                            explicit ν, η and the DynSGS viscosity in compression heat;
#                            "viscous": all the viscous dissipation and the explicit η heat,
#                            the DynSGS resistivity not, as Lare's default tensor shock
#                            viscosity does (README, "The heat of the DynSGS"); or a fraction
#                            put into the plasma (1: all, energy conserving; 0: none) [shock]
#   JEXPRESSO_LK10_ETAFAC    factor on the DynSGS resistivity (the :μ multiplier of the B
#                            slots), a sensitivity test of the reconnection          [1]
#   JEXPRESSO_LK10_VISCFAC   factor on the DynSGS viscosity (the :μ multipliers of the ρu,
#                            ρv, ρw slots; after its caps, the case's ν untouched), a test
#                            of its drag on the plasma that drains from the dome     [1]
#   JEXPRESSO_LK10_VISC      momentum | velocity | nodal: the viscosity (DynSGS and ν)
#                            diffuses the momentum ρv/ρ_e, or acts on the velocity gradient,
#                            with the momentum of the mass the DynSGS diffuses; its stress
#                            weighs every node of an element with the lightest one's
#                            density (velocity) or each with its own, up to what the
#                            lightest takes at the stability limit (nodal) (README, "The
#                            viscosity on the velocity")                       [momentum]
#   JEXPRESSO_LK10_EINT      off | lowbeta | all: the energy equation as the total energy's
#                            own (off), or, where the internal energy is under 10⁻³…10⁻² of
#                            the total (lowbeta) or everywhere (all), as the changes of the
#                            kinetic and magnetic energies plus the internal energy's own
#                            equation, so that the pressure of a drained low-β region is not
#                            the remainder of two large energies (README, Positivity:
#                            "The internal energy where it is a small share")       [off]
#   JEXPRESSO_LK10_STABCAP   the DynSGS coefficient never above this viscous number, μ Δt/Δx²,
#                            which the explicit step still integrates (0: no cap)    [0.2]
#   JEXPRESSO_LK10_DSGS_STAGES  all | 1,5 | 1: the stages of the Runge-Kutta scheme at which
#                            the DynSGS terms are evaluated: all five; the first and the
#                            last, held in between (the step costs 0.52 of the five-stage one);
#                            or the first alone (0.47, and the default of
#                            JEXPRESSO_LK10_STABCAP becomes 0.08). The coefficient is then that
#                            of the first stage, with half the C_1 (README, "The DynSGS at
#                            two stages of the five")                                [all]
#   JEXPRESSO_LK10_C1        the DynSGS constant C_1 (the coefficient is C_1 Δ² × the largest
#                            normalized residual, up to its caps)                    [1]
#   JEXPRESSO_LK10_ETACAP    fast | flow: the DynSGS resistivity capped at C_max Δ (|v| + c_f)
#                            like the viscosity, or at C_max Δ |v|, the flow speed alone [fast]
#   JEXPRESSO_LK10_SENSOR    legacy | residual: the DynSGS coefficient from |∂ₜq| (the
#                            sensor the parent case was validated with), or from the
#                            element residual of the equations, once per step     [legacy]
#   JEXPRESSO_LK10_ENERGY    thermal | total: the DynSGS energy slot diffuses the
#                            internal energy, with the work and heat of the other
#                            slots (energy-consistent), or E as a whole          [thermal]
#   JEXPRESSO_LK10_CURLB     0: DynSGS resistivity as the component Laplacian       [1: curl form]
#   JEXPRESSO_LK10_NU        dynamic | kinematic: the paper's dynamic ν below the
#                            photosphere (capped above), or kinematic everywhere   [dynamic]
#   JEXPRESSO_LK10_PARTITION balanced | uniform: the ranks hold equal numbers of elements,
#                            or equal areas (the framework's default: on these graded
#                            meshes one of 128 ranks then holds 6.8 times the average
#                            and every step waits for it)                       [balanced]
#   JEXPRESSO_LK10_SUBNORMALS  flush | keep: subnormal numbers flushed to zero (they
#                            slow the run several times, initialize.jl), or IEEE  [flush]
#   JEXPRESSO_LK10_PROGRESS  minutes between the progress lines of the log (simulated time
#                            reached, hours so far and to go); 0: none (printed by the
#                            stage limiter, so not with JEXPRESSO_LK10_FLOORS=0)     [10]
#   JEXPRESSO_LK10_RESTART_DIR  continue an earlier run of this case from a snapshot:
#                            its output directory, or "latest" (the most recently
#                            written run of the same configuration); same switches,
#                            mesh and rank count (README, "Restart")             [none]
#   JEXPRESSO_LK10_RESTART_TIME  the snapshot nearest this time in SECONDS
#                            [the last complete one]
#
# Times are in t0 = 23.19 s inside the code (user_flux.jl).
#
# A JEXPRESSO_LK10_* variable that this deck does not read STOPS the run (lk10_check_env):
# a misspelled switch, or a switch newer than the checkout the job runs from, would
# otherwise be ignored without a word, and the run be the one it was meant to change.
#---------------------------------------------------------------------------------
if !@isdefined(_LK10_ASKED)
    const _LK10_ASKED = Set{String}(("SUBNORMALS",))     # every switch asked for (SUBNORMALS: in initialize.jl)
end
_lk10_env(k, d) = (push!(_LK10_ASKED, String(k)); strip(get(ENV, "JEXPRESSO_LK10_" * k, d)))
# switches that existed and were removed, with the reason the error gives
if !@isdefined(_LK10_WITHDRAWN)
    const _LK10_WITHDRAWN = Dict(
        "RHOFLOOR" => "the velocity-keeping density floor made every restart go non-finite and was withdrawn on 24 Nov 2026" *
                      " (README, Positivity: What the density floor adds); unset it, the floor is the one of before")
end
# the JEXPRESSO_LK10_* variables of the environment, "NAME=value, ..."; an error if one was not asked for
function lk10_check_env()
    set = sort!([String(k)[16:end] for k in keys(ENV) if startswith(String(k), "JEXPRESSO_LK10_")])
    for k in set
        haskey(_LK10_WITHDRAWN, k) && error(" problems/MHD/fluxEmergenceLeake2010: JEXPRESSO_LK10_", k, ": ", _LK10_WITHDRAWN[k], ".")
    end
    bad = [k for k in set if !(k in _LK10_ASKED)]
    isempty(bad) ||
        error(" problems/MHD/fluxEmergenceLeake2010: ", join("JEXPRESSO_LK10_" .* bad, ", "), length(bad) > 1 ? " are" : " is",
              " not read by this code: misspelled, or this checkout is older than the switch (git log -1; git pull where the job runs).",
              " It reads: ", join(sort!(collect(_LK10_ASKED)), ", "), ".")
    return join(("$(k)=$(ENV["JEXPRESSO_LK10_" * k])" for k in set), ", ")
end
_lk10_num(k, d) = something(tryparse(Float64, _lk10_env(k, "")), d)
# C_1 (:dsgs_C1, :lk10_C1); JEXPRESSO_LK10_CR, its name until Nov 2026, is still read, with a warning
function _lk10_c1()
    v = _lk10_env("C1", "")
    if isempty(v)
        old = _lk10_env("CR", "")        # asked for, so lk10_check_env accepts it
        isempty(old) || (@warn "JEXPRESSO_LK10_CR is now JEXPRESSO_LK10_C1"; v = old)
    end
    return something(tryparse(Float64, v), 1.0)
end
# one of the words `allowed` (the first is the default), in lower case; anything else is an error, not the default
function _lk10_word(k, allowed...)
    w = lowercase(_lk10_env(k, allowed[1]))
    w in allowed || error(" JEXPRESSO_LK10_$(k) must be one of $(join(allowed, ", ")) (got \"$(w)\").")
    return w
end
# JEXPRESSO_LK10_HEAT: "shock" (default), "viscous" or the fraction of the dissipated heat kept
# JEXPRESSO_LK10_DSGS_STAGES: nothing (all: every stage), (1, 5) or (1,)
function _lk10_stages()
    w = replace(lowercase(strip(_lk10_env("DSGS_STAGES", "all"))), " " => "")
    w == "all" && return nothing
    w in ("1,5", "15") && return (1, 5)
    w == "1" && return (1,)
    error(" problems/MHD/fluxEmergenceLeake2010: JEXPRESSO_LK10_DSGS_STAGES = \"$(w)\", expected all, 1,5 or 1.")
end

# JEXPRESSO_LK10_RHOSENS: nothing (off), or (a0, a1, C) (:dsgs_rho_sensor)
function _lk10_rhosens()
    v = lowercase(_lk10_env("RHOSENS", "off"))
    v == "off" && return nothing
    v == "on" && return (0.03, 0.15, 1.0)
    p = [tryparse(Float64, strip(x)) for x in split(v, ",")]
    (length(p) == 3 && all(x -> x !== nothing, p) && 0.0 <= p[1] < p[2] && p[3] >= 0.0) ||
        error(" JEXPRESSO_LK10_RHOSENS must be off, on or a0,a1,C with 0 <= a0 < a1 and C >= 0 (got \"$(v)\").")
    return (p[1], p[2], p[3])
end
function _lk10_heat()
    h = lowercase(_lk10_env("HEAT", "shock"))
    h in ("shock", "viscous") && return h
    v = tryparse(Float64, h)
    (v === nothing || !(0.0 <= v <= 1.0)) && error(" JEXPRESSO_LK10_HEAT must be shock, viscous or a number from 0 to 1 (got \"$(h)\").")
    return v
end

function _lk10_case()
    tube  = something(tryparse(Int, _lk10_env("TUBE", "")), 2)
    field = lowercase(_lk10_env("FIELD", "quadrupole"))
    str   = lowercase(_lk10_env("STRENGTH", "medium"))
    ref   = something(tryparse(Int, _lk10_env("REF", "")), 2)
    istr  = Dict("strong" => 1, "medium" => 2, "weak" => 3)
    haskey(istr, str) || error(" problems/MHD/fluxEmergenceLeake2010: JEXPRESSO_LK10_STRENGTH = $(str), expected strong, medium or weak.")
    ref in (1, 2, 3, 4) || error(" problems/MHD/fluxEmergenceLeake2010: JEXPRESSO_LK10_REF = $(ref), expected 1, 2, 3 or 4.")
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
# Positivity floors, as the Runge-Kutta stage limiter: on by default
# (JEXPRESSO_LK10_FLOORS=0 switches them off), as in the other Leake cases. Where the
# field is strong and β low, p is the small remainder E − ½ρv² − ½B² − ½ψ² and goes
# negative: a run of the default deck from ~200 s in the upper chromosphere under the
# quadrupole's strongest field (non-finite at 631 s), the field-free run at the top of
# the emerging tube by 1000 s (README, "Positivity"). Two floors, chosen not to act on
# what moves legitimately:
#   T = p/ρ ≥ T_min = 0.1 T_ph (570 K; JEXPRESSO_LK10_TFLOOR), below any temperature of
#     the model, through E at fixed velocity and field;
#   ρ ≥ f p_e(y)/T_cor, f = 0.2: f ρ_e in the corona, and in the chromosphere 150 times
#     below its density, so that the transition region can move (the other cases'
#     f ρ_e(y) and f p_e(y) would add mass and heat where it dips).
# Not conservative where they fire; every hit is counted and reported by rank 0 each
# time the total passes a power of ten.
#---------------------------------------------------------------------------------
if !@isdefined(lk10_floor_frac)
    const lk10_floor_frac = Ref{Float64}(0.2)
    const lk10_floor_hits = Ref{Int}(0)
    const lk10_floor_next = Ref{Int}(1)
end
if !@isdefined(lk10_floor_Tmin)
    const lk10_floor_Tmin = Ref{Float64}(0.1)
end
# Progress line (rank 0): the snapshots are 116 s of simulated time apart, an hour of wall
# time on the fine meshes, and nothing else is printed in between.
if !@isdefined(lk10_progress)
    const lk10_progress       = Ref{NTuple{3,Float64}}((0.0, 0.0, 0.0))   # wall time and t of the first call, wall time of the last line
    const lk10_progress_every = Ref{Float64}(600.0)                       # s of wall time between lines (0: none)
    const lk10_progress_tend  = Ref{Float64}(0.0)
end
function lk10_report_progress(t)
    now = time()
    w0, t0, wl = lk10_progress[]
    if w0 == 0.0
        lk10_progress[] = (now, t, now)
    elseif now - wl >= lk10_progress_every[] && t > t0
        left = (lk10_progress_tend[] - t)*(now - w0)/(t - t0)
        print(string(" # progress: t = ", round(t, digits = 2), " t0 = ", round(Int, t*LK10_T0), " s of ",
                     round(Int, lk10_progress_tend[]*LK10_T0), " s; ", round((now - w0)/3600, digits = 2),
                     " h so far, about ", round(max(left, 0.0)/3600, digits = 1), " h to go at this rate\n"))
        flush(stdout)
        lk10_progress[] = (w0, t0, now)
    end
    return nothing
end

function lk10_positivity_limiter!(u, integrator, p, t)
    npoin = p.mesh.npoin
    qe    = p.qp.qe
    frac  = lk10_floor_frac[]
    Tmin  = lk10_floor_Tmin[]
    γm1   = γ_mhd - 1.0
    nhit  = 0
    worst = (0.0, 0.0, 0.0, 0.0, 0.0)    # (depth below the floor, 1 = ρ / 2 = T, x, y, value/floor) on this rank
    X     = p.mesh.coords
    @inbounds for ip = 1:npoin
        iρ = ip; iρu = npoin + ip; iρv = 2npoin + ip; iE = 3npoin + ip; iρw = 4npoin + ip
        iBx = 5npoin + ip; iBy = 6npoin + ip; iBz = 7npoin + ip; iψ = 8npoin + ip
        ρmin = frac*qe[ip,end]/LK10_TCOR
        ρ = u[iρ]
        if ρ < ρmin
            d = 1.0 - ρ/ρmin
            d > worst[1] && (worst = (d, 1.0, X[1,ip], X[2,ip], ρ/ρmin))
            u[iρ] = ρmin; ρ = ρmin; nhit += 1
        end
        ke = 0.5*(u[iρu]^2 + u[iρv]^2 + u[iρw]^2)/ρ
        me = 0.5*(u[iBx]^2 + u[iBy]^2 + u[iBz]^2) + 0.5*u[iψ]^2
        pmin = Tmin*ρ
        pnow = γm1*(u[iE] - ke - me)
        if pnow < pmin
            d = 1.0 - pnow/pmin
            d > worst[1] && (worst = (d, 2.0, X[1,ip], X[2,ip], pnow/ρ))
            u[iE] = pmin/γm1 + ke + me; nhit += 1
        end
    end
    comm = MPI.COMM_WORLD
    nall = MPI.Comm_size(comm) > 1 ? MPI.Allreduce(nhit, +, comm) : nhit
    lk10_floor_hits[] += nall
    if nall > 0 && lk10_floor_hits[] >= lk10_floor_next[]
        # the deepest violation of this stage on all ranks: which floor, its value, where
        all = MPI.Comm_size(comm) > 1 ? MPI.Allgather(collect(worst), comm) : collect(worst)
        k   = argmax(all[1:5:end]); w = all[5k-4:5k]
        if MPI.Comm_rank(comm) == 0
            what = w[2] == 1.0 ? string("ρ = ", round(w[5], sigdigits = 3), " of its floor 0.2 p_e/T_cor") :
                                 string("T = ", round(w[5], sigdigits = 3), " T_ph (floor ", Tmin, ")")
            print(string(" # positivity floors: ", lk10_floor_hits[], " node hits so far on all ranks (",
                         nall, " in this stage) at t = ", round(t, digits = 3), " t0 = ",
                         round(t*LK10_T0, digits = 1), " s; deepest: ", what, " at (x, y) = (",
                         round(w[3]*0.15, digits = 2), ", ", round(w[4]*0.15, digits = 2), ") Mm\n"))
            flush(stdout)
        end
        while lk10_floor_next[] <= lk10_floor_hits[]
            lk10_floor_next[] *= 10
        end
    end
    (lk10_progress_every[] > 0 && MPI.Comm_rank(comm) == 0) && lk10_report_progress(t)
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
    mesh = ("LK10_ref1.msh", "LK10_ref2.msh", "LK10_ref3.msh", "LK10_ref4.msh")[c.ref]
    # c_h Δt/Δh = 0.25 on the smallest LGL gap of the mesh (0.1727 of the 2.5,
    # 1.25, 0.625 L0 core elements at :nop => 4); c_h is the fastest wave of the
    # initial state, at the top of the transition region (README, "Cost").
    Δh   = 0.17267*(2.5, 1.25, 0.625, 1.25)[c.ref]   # ref 4: ref 2's core elements, over the whole figure region
    Δt   = _lk10_num("DT", 0.25*Δh/c.cmax)
    tend = _lk10_num("TEND", 2530.0)/LK10_T0
    # Snapshots every 5 t0 (116 s) and at the times of the erratum's figures
    # (Figs. 1 and 3: 920, 1380, 1610, 1955, 2070, 2185, 2254, 2530 s)
    tfig = [920.0, 1380.0, 1610.0, 1955.0, 2070.0, 2185.0, 2254.0, 2530.0] ./ LK10_T0
    dout = _lk10_num("DOUT", 5.0*LK10_T0)/LK10_T0          # JEXPRESSO_LK10_DOUT: snapshot interval in s
    tout = sort(unique(round.(vcat(collect(0.0:dout:tend), filter(t -> t <= tend + 1e-9, tfig), [tend]), digits = 6)))
    floors = !(_lk10_env("FLOORS", "1") in ("0", "false", "no", "off"))
    # JEXPRESSO_LK10_DSGS_STAGES: the stages at which the DynSGS block is evaluated (rhs.jl, :dsgs_stages)
    stages = _lk10_stages()
    sensor = lowercase(_lk10_env("SENSOR", "legacy")) == "residual" ? "residual" : "legacy"
    lk10_floor_Tmin[] = _lk10_num("TFLOOR", 0.1)
    lk10_floor_hits[] = 0; lk10_floor_next[] = 1
    lk10_progress[] = (0.0, 0.0, 0.0); lk10_progress_tend[] = tend
    lk10_progress_every[] = 60.0*_lk10_num("PROGRESS", 10.0)

    inputs = Dict(
        #---------------------------------------------------------------------------
        # Time integration: the low-storage RK of the parent case, no limiter
        # but the positivity floors (above; JEXPRESSO_LK10_FLOORS=0: none)
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
        :lrestart             => false,      # the restart is the case's own (JEXPRESSO_LK10_RESTART_DIR, initialize.jl)
        :lk10_restart_dir     => String(_lk10_env("RESTART_DIR", "")),
        :lk10_floors          => floors,
        :lk10_floor_Tmin      => lk10_floor_Tmin[],   # T_ph units
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
        # (slots ρ, ρu, ρv, E, ρw, Bx, By, Bz, ψ; JEXPRESSO_LK10_VISCFAC scales the viscosity, _ETAFAC
        # the resistivity, _RHOFAC the mass diffusion)
        :μ                    => [_lk10_num("RHOFAC", 1.0), _lk10_num("VISCFAC", 1.0), _lk10_num("VISCFAC", 1.0), 1.0,
                                  _lk10_num("VISCFAC", 1.0), fill(_lk10_num("ETAFAC", 1.0), 3)..., 1.0],
        :visc_model           => DSGS_MHD(),
        # "legacy" (DSGS.md §1.2): the assembled RHS against a BDF2 of the RK stage state, which
        # reads -∂ₜq/2 at the first stage — R ≈ |∂ₜq|, so every moving structure is diffused,
        # the corona's waves included. "residual": the element residual, frozen once per step.
        :dsgs_sensor          => sensor,
        # The DynSGS coefficient is held at zero for the first steps of a run while the BDF2 history fills
        # (:dsgs_hold_steps, rhs.jl). At a restart that leaves a developed flow with no DynSGS for two steps
        # (a restart from 2185 s of a drained-dome run went non-finite at its fourth step, the coefficient
        # on its cap). The legacy sensor needs no hold there: on the history seeded with the restart state
        # its first reading is |∂ₜq| itself, the gradient it always reads. Restarts with it hold none.
        :dsgs_hold_steps      => (sensor == "legacy" && !isempty(strip(_lk10_env("RESTART_DIR", "")))) ? 0 : 2,
        # JEXPRESSO_LK10_DSGS_STAGES=1,5 or 1: the coefficient is that of the first stage, where the
        # legacy sensor reads twice its mean over the stages: half the C_1, same coefficient.
        :dsgs_stages          => stages,
        :lk10_C1              => _lk10_c1(),
        :dsgs_C1              => _lk10_c1()*((stages !== nothing && sensor == "legacy") ? 0.5 : 1.0),
        :dsgs_Cmax            => 0.5,
        :dsgs_Cmin            => _lk10_num("CMIN", 0.0),
        :dsgs_gamma           => 5.0/3.0,
        :dsgs_Prt             => 1.0,
        :dsgs_nazarov_energy  => false,
        # Energy-consistent DynSGS (SGS.jl, dsgs_thermal_energy; README, "Positivity"): the
        # E slot diffuses the internal energy and receives the work and the heat of the other
        # slots, so the DynSGS cannot drive p negative in the low-β chromosphere and corona.
        :dsgs_thermal_energy  => lowercase(_lk10_env("ENERGY", "thermal")) != "total",
        # ... of which only the DynSGS viscosity in compression and the explicit ν, η heat (shock),
        # or all the viscous dissipation, as Lare's default tensor shock viscosity heats (viscous)
        # (JEXPRESSO_LK10_HEAT; README, "The heat of the DynSGS")
        :dsgs_heating         => _lk10_heat(),
        # The mass diffusion where the density has a node-to-node part (SGS.jl, :dsgs_rho_sensor;
        # README, "The density sensor"): the contacts the residual sensor does not see
        :dsgs_rho_sensor      => _lk10_rhosens(),
        # Resistivity from the induction residuals alone (SGS.jl, linduction; README, "DynSGS in
        # the corona"): acoustic fronts in the low-β corona then do not make the field resistive.
        :dsgs_induction_sensor => _lk10_word("ETA", "induction", "all") == "induction",
        # The viscosity on the velocity gradient instead of the diffusion of ρv/ρ_e (SGS.jl,
        # _dsgs_velocity_viscosity_2d!; README, "The viscosity on the velocity"): opt-in while under test
        :dsgs_velocity_viscosity => _lk10_word("VISC", "momentum", "velocity", "nodal") != "momentum",
        # ... its stress weighing each node with its own density, up to what the element's lightest node
        # takes at the stability limit (VISC=nodal; SGS.jl, _dsgs_scap), or all with the lightest (velocity)
        :dsgs_velocity_weight => _lk10_word("VISC", "momentum", "velocity", "nodal") == "nodal" ? "node" : "smin",
        # The energy slot where the internal energy is a small share of the total: the chain rule
        # of the kinetic and magnetic energies plus ∂ₜe + ∇·(e u) + p∇·u = S_e (rhs.jl,
        # _mhd_internal_energy_2d!; README, Positivity: "The internal energy where it is a small share");
        # Lare's energy equation
        :mhd_internal_energy  => _lk10_word("EINT", "off", "lowbeta", "all"),
        # ... and capped at the first-order diffusivity of a field carried by the flow, C_max Δ |v|,
        # instead of the viscosity's C_max Δ (|v| + c_f) (SGS.jl, lflowcap): opt-in while under test
        :dsgs_eta_flow_cap    => lowercase(_lk10_env("ETACAP", "fast")) == "flow",
        # Never more diffusivity than the explicit step integrates (SGS.jl, :dsgs_stab_cap): the
        # first-order cap C_max Δ c_f follows the Alfvén speed, which doubles as the dome drains.
        # The scheme takes a viscous number of 0.42; with the DynSGS at stages 1 and 5, 0.6; at
        # stage 1 alone, 0.18 (rhs.jl, :dsgs_stages).
        :dsgs_stab_cap        => _lk10_num("STABCAP", stages == (1,) ? 0.08 : 0.2),
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
        # Rank partition (mesh.jl, _compute_xy_partition): the x-y bins with equal numbers of
        # elements. No ghost cells, so the solution is that of the uniform bins to round-off.
        :lxy_partition        => true,
        :lxy_partition_balanced => lowercase(_lk10_env("PARTITION", "balanced")) != "uniform",
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
    rt = _lk10_num("RESTART_TIME", NaN)              # seconds
    isnan(rt) || (inputs[:lk10_restart_time] = rt/LK10_T0)

    inputs[:lk10_env_set] = lk10_check_env()         # after every switch has been asked for
    return inputs
end
