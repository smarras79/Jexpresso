#---------------------------------------------------------------------------------
# sphere_time_loop.jl
#
# Time integration of the shallow water equations on the spherical shell.
#
# The scheme is whatever `inputs[:ode_solver]` names: this file builds an
# ODEProblem and hands it to OrdinaryDiffEq like every other Jexpresso case.
# The default deck asks for SSPRK33(), which is the Shu-Osher SSP-RK3 of
# Marras, Kopera & Giraldo (2015), QJRMS 141: 1727-1739, section 4.2:
#
#   q⁽¹⁾ = qⁿ + Δt R(qⁿ)
#   q⁽²⁾ = ¾ qⁿ + ¼ [ q⁽¹⁾ + Δt R(q⁽¹⁾) ]
#   qⁿ⁺¹ = ⅓ qⁿ + ⅔ [ q⁽²⁾ + Δt R(q⁽²⁾) ]
#
# THE LAGRANGE PROJECTION IS APPLIED AFTER EVERY STAGE, not only at the end of
# the step. Each stage is a full Euler update, so each one can push momentum off
# the shell; projecting only at the end would feed an off-manifold state into
# the next stage's flux evaluation. This is the discrete half of the constraint
# — the continuous half is the μx term in user_source.jl.
#
# HOW THE PER-STAGE HOOK IS OBTAINED
# ----------------------------------
# Every SSPRK integrator in OrdinaryDiffEq takes a `stage_limiter!` and a
# `step_limiter!`, and calls them at exactly the points this scheme needs. From
# OrdinaryDiffEqSSPRK's perform_step! for SSPRK33Cache:
#
#   u = uprev + dt*f(uprev)   ; stage_limiter!(u, integrator, p, t+dt)
#   u = (3uprev + u + dt*k)/4 ; stage_limiter!(u, integrator, p, t+dt/2)
#   u = (uprev + 2u + 2dt*k)/3; stage_limiter!(u, integrator, p, t+dt)
#                               step_limiter!(u, integrator, p, t+dt)
#
# So the mapping is exact and needs no hand-rolled loop:
#
#   stage_limiter! = the Lagrange projection P = I - xxᵀ/r²   (after each stage)
#   step_limiter!  = the modal filter, then one more projection (once per step)
#
# This is why `inputs[:ode_solver]` must name a limiter-capable integrator: the
# projection is part of the discretization, not an optional post-process, so an
# algorithm with nowhere to hang it is an error rather than a silent fallback.
#
# DIAGNOSTICS, printed every :ndiagnostics_prints steps:
#
#   mass     ∫φ dΩ = Σ M[ip] φ[ip]   — must be conserved to round-off; the
#                                      scheme is conservative and the closed
#                                      manifold has no boundary flux
#   energy   ∫(½φ|u|² + ½φ²) dΩ      — decays slowly with the filter on
#   drift    max|(φu)·x̂|             — what the projection removes each stage.
#                                      If it grows, the discretization is
#                                      leaving the manifold
#   dh       max|h - h(t=0)|         — for the balanced jet this IS the error
#
# S. Marras & contributors
#---------------------------------------------------------------------------------

export sphere_time_loop!
export sphere_cfl_dt
export sphere_diagnostics
export sphere_output_times
export sphere_dsgs_mu_bound
export sphere_dsgs_init_history!
export sphere_dsgs_roll_history!


#---------------------------------------------------------------------------------
# Δt from the CFL condition.
#
#   Δt = CFL min( Δmin / max(|u| + √φ) ,  Δmin² / ν )
#
# √φ = √(gh) is the gravity-wave speed, |u| the advective speed; Δmin is the
# smallest LGL node spacing, which for a spectral element clusters like 1/nop²
# towards the edges and is what actually limits the step.
#
# The SECOND branch is the diffusive limit of the artificial viscosity: it is a
# second derivative, so it is explicit-stable only for Δt ~ Δ²/ν and not Δ/u.
# For the Galewsky test at the shipped resolution the two are 1e2 s and 1e5 s,
# so the wave speed wins by three orders of magnitude and ν = 1e5 m²/s costs
# nothing — but ν is a free parameter, and raising it far enough silently turns
# a stable run into an exponentially growing one. Taking the min is what makes
# "the step is too big for this ν" impossible rather than mysterious.
#---------------------------------------------------------------------------------
function sphere_cfl_dt(q::AbstractMatrix{TF}, mesh::St_mesh,
                       metrics::St_sphere_metrics{TF}; cfl = 0.35, μmax = 0.0) where {TF}

    cmax = zero(TF)
    #
    # `npoin::Int`, not `npoin = Int(...)`. mesh.npoin is an ::Any field, so
    # `Int(mesh.npoin)` is a dynamic call and INFERS ::Any — which makes the
    # loop bound dynamic, widens cmax, and left this function returning ::Any.
    # Δt is what that ::Any became at the call site. The declaration forces the
    # conversion and hands the caller a Float64.
    #
    npoin::Int = Int(mesh.npoin)
    @inbounds for ip = 1:npoin
        φ = q[ip,1]
        φ > 0 || error(string(" # ERROR sphere_time_loop.jl: non-positive geopotential φ = ", φ,
                              " at node ", ip, ". The state has gone unphysical."))
        umag = sqrt(q[ip,2]^2 + q[ip,3]^2 + q[ip,4]^2)/φ
        cmax = max(cmax, umag + sqrt(φ))
    end

    Δt::TF = cfl*metrics.Δmin/cmax

    ν::TF = TF(μmax)
    if ν > zero(TF)
        Δt = min(Δt, cfl*metrics.Δmin*metrics.Δmin/ν)
    end

    return Δt
end


#---------------------------------------------------------------------------------
# sphere_output_times(inputs, tinit, tend) -> Vector{Float64}, ascending
#
# WHEN the run writes. Two schedules, in priority order, because Jexpresso has
# two ways of asking and the shell has to understand BOTH:
#
#   :ndiagnostics_outputs => n > 0    n equally spaced writes, the last at :tend.
#                                     What the decks in problems/ use.
#   :diagnostics_at_times => times    an explicit list. What the flat time loop
#                                     (TimeIntegrators.jl) has always used, and
#                                     what CI_MODE sets — see below.
#
# THIS FUNCTION IS WHY ShallowWater/SWsphere COULD NOT RUN UNDER CI. The shell
# loop used to read :ndiagnostics_outputs and nothing else, and the two keys are
# not independent: setting :diagnostics_at_times makes mod_inputs_user_inputs!
# zero :ndiagnostics_outputs (mod_inputs.jl, the else branch at ~line 560), on
# the reasoning that a deck naming explicit times does not also want a uniform
# cadence. CI_MODE sets :diagnostics_at_times => [tend] (run.jl) to force the
# single reference write. So under CI the shell saw
#
#     :ndiagnostics_outputs = 0   ⇒   nout = 0   ⇒   no write ever fires
#
# and, with CI_MODE's :lwrite_initial => false killing the initial write too,
# the case solved all the way to :tend and produced NOT ONE FILE. The failure
# surfaced as run_ci_case's "the solver returned without writing any HDF5
# output", which reads like a blow-up and is not one: the solve is fine, the
# schedule was empty.
#
# Times at or before :tinit are dropped — the initial state is written
# separately, by :lwrite_initial — and the list is sorted and de-duplicated so
# the monitor can walk it with a single index.
#---------------------------------------------------------------------------------
function sphere_output_times(inputs, tinit::Float64, tend::Float64)

    nout = Int(get(inputs, :ndiagnostics_outputs, 0))
    if nout > 0
        outdt = (tend - tinit)/nout
        return Float64[tinit + k*outdt for k = 1:nout]
    end

    raw = get(inputs, :diagnostics_at_times, nothing)
    raw === nothing && return Float64[]

    times = raw isa Number ? Float64[Float64(raw)] : collect(Float64, raw)
    # `> tinit`, strictly: mod_inputs_user_inputs! defaults
    # :diagnostics_at_times to :tend, so a deck that sets neither key still gets
    # exactly one write, at the end.
    return sort!(unique!(filter(x -> x > tinit, times)))
end


#---------------------------------------------------------------------------------
# DynSGS support: the ν bound, and the BDF2 history.
#---------------------------------------------------------------------------------
#
# sphere_dsgs_mu_bound(q, mesh, sp)
#
# The LARGEST ν the model can ever return on this state,
#
#   ν ≤ μ_max|e = C2 · (Δ_elem,e/(N+1)) · (‖u‖ + √φ)_{∞,e} ,
#
# which is what the diffusive branch of the CFL condition needs. DynSGS sizes ν
# per element and per stage, so there is no constant to hand sphere_cfl_dt the
# way a fixed-ν deck does; the model's OWN cap is the bound, and it is a genuine
# one — μ = max(0, min(μ_max, μ_res)) can never exceed it, whatever the residual
# does. Evaluated once, on the initial state, for the same reason Δt is: the
# step is fixed for the run.
#
# The bound is not tight (μ_res is normally an order of magnitude below the cap,
# see DSGS.md §7), so this over-estimates ν and therefore under-estimates the
# diffusive Δt. That is the safe direction, and on the shipped grid the
# advective limit still wins by a wide margin — the printout says by how much.
#
function sphere_dsgs_mu_bound(q::AbstractMatrix{TF}, mesh::St_mesh,
                              sp::St_sphere_params) where {TF}

    sp.ldsgs || return 0.0

    nelem::Int = Int(mesh.nelem)
    ngl::Int   = Int(mesh.ngl)
    conn       = mesh.connijk

    μmax = zero(TF)
    @inbounds for iel = 1:nelem
        Δ = sp.Δelem[iel]/ngl
        w = zero(TF)
        for j = 1:ngl, i = 1:ngl
            ip = conn[iel,i,j]
            φ  = q[ip,1]
            φ > 0 || continue
            umag = sqrt(q[ip,2]^2 + q[ip,3]^2 + q[ip,4]^2)/φ
            w = max(w, umag + sqrt(φ))
        end
        μmax = max(μmax, sp.C2*Δ*w)
    end

    return Float64(μmax)
end

#
# The BDF2 history the residual is built on.
#
# BOTH buffers start at the initial state, so the first residual is
# −M⁻¹·RHS(q⁰) — the actual ∂q/∂t of the initial condition — rather than the
# 3q⁰/(2Δt) an empty history would give. For the Galewsky jet that is nearly
# zero, because the balanced state is a steady solution.
#
# q.qn carries one column beyond the neqs solution slots (define_q allocates
# npoin × (neqs+1)), so these copies are column-wise rather than a whole-array
# `.=`, which would be a DimensionMismatch.
#
function sphere_dsgs_init_history!(sp::St_sphere_params, qn::AbstractMatrix,
                                   npoin::Int, neqs::Int)
    sp.ldsgs || return nothing
    @inbounds for ieq = 1:neqs, ip = 1:npoin
        sp.qnm1[ip,ieq] = qn[ip,ieq]
        sp.qnm2[ip,ieq] = qn[ip,ieq]
    end
    sp.thist[] = -1.0e30
    return nothing
end

#
# Roll it ONCE PER TIME STEP — not once per RK stage.
#
# The RHS is evaluated at every stage, and `t` sweeps t + c_i·Δt inside a step,
# so the gate below fires on the first stage of each step and is quiet for the
# rest. Rolling per stage instead would difference consecutive INTERMEDIATE
# states over a full Δt, which is not an approximation of ∂q/∂t at all — the
# defect the flat paths carried until the dsgs_qnm1/qnm2 pair was introduced
# for them (DSGS.md §4.4, §7.4). The 0.999 absorbs the round-off in c_i·Δt.
#
# After the roll, qnm2 holds qⁿ (the state this step starts from) and qnm1 holds
# qⁿ⁻¹, which is what compute_dsgs_viscosity! is handed as (q1, q2).
#
function sphere_dsgs_roll_history!(sp::St_sphere_params, u::AbstractMatrix,
                                   npoin::Int, neqs::Int, t::Float64)
    sp.ldsgs || return nothing
    Δt = sp.Δt[]
    Δt > 0 || return nothing
    t - sp.thist[] >= 0.999*Δt || return nothing

    @inbounds for ieq = 1:neqs, ip = 1:npoin
        sp.qnm1[ip,ieq] = sp.qnm2[ip,ieq]
        sp.qnm2[ip,ieq] = u[ip,ieq]
    end
    sp.thist[] = t
    return nothing
end


#---------------------------------------------------------------------------------
# Conserved quantities and the constraint drift.
#---------------------------------------------------------------------------------
function sphere_diagnostics(q::AbstractMatrix{TF}, mesh::St_mesh,
                            metrics::St_sphere_metrics{TF}) where {TF}

    mass  = zero(TF)
    ener  = zero(TF)
    npoin::Int = Int(mesh.npoin)     # ::Any field; see the note in sphere_cfl_dt
    M     = metrics.M
    @inbounds for ip = 1:npoin
        φ  = q[ip,1]
        m2 = q[ip,2]^2 + q[ip,3]^2 + q[ip,4]^2
        mass += M[ip]*φ
        ener += M[ip]*(0.5*m2/φ + 0.5*φ*φ)
    end
    # ::TF for the same reason: mesh.x/y/z are ::Any, so the call inside
    # sphere_normal_momentum is dynamic and its result infers ::Any.
    drift::TF = sphere_normal_momentum(q, mesh; ivar = 2)

    return mass, ener, drift
end


#---------------------------------------------------------------------------------
# What the RHS and the limiters get as OrdinaryDiffEq's `p`.
#
# CONCRETELY TYPED ON PURPOSE. `inputs` is a Dict{Symbol,Any}, so anything read
# out of it is inferred ::Any; every value the per-stage code touches is pulled
# out once here and stored in a parameterised field, which is what lets
# _sphere_ode_rhs! and the two limiters specialise instead of compiling against
# Any. This is the same function-barrier trick the shell kernels use, moved up
# to the integrator boundary.
#---------------------------------------------------------------------------------
struct St_sphere_ode_params{TMesh, TMetrics, TParams, TQe, TSVT}
    mesh::TMesh
    metrics::TMetrics
    sp::TParams
    qe::TQe
    SVT::TSVT
    lproject::Bool
    driftmax::Base.RefValue{Float64}
    npoin::Int
    neqs::Int
end


#
# R(q). OrdinaryDiffEq calls this once per stage.
#
function _sphere_ode_rhs!(du, u, p::St_sphere_ode_params, t)
    # DynSGS only: advance the BDF2 history on the FIRST stage of each step.
    # The gate lives here rather than in a limiter because the residual is read
    # inside sphere_rhs!, so the history has to be current before it runs.
    sphere_dsgs_roll_history!(p.sp, u, p.npoin, p.neqs, Float64(t))
    sphere_rhs!(du, u, p.qe, p.mesh, p.metrics, p.sp, p.SVT)
    return nothing
end


#
# STAGE LIMITER — the Lagrange projection, after every stage.
#
# Also accumulates max|(φu)·x̂| over the whole run, which is the "how far did the
# discretization try to leave the shell" diagnostic reported at the end.
#
function _sphere_stage_limiter!(u, integrator, p::St_sphere_ode_params, t)
    if p.lproject
        d = project_momentum_to_sphere!(u, p.mesh; ivar = 2)
        p.driftmax[] = max(p.driftmax[], d)
    end
    return nothing
end


#
# STEP LIMITER — the modal filter, then one more projection.
#
# The filter is an L² projection back onto the continuous space and so conserves
# mass exactly, but it does NOT preserve the tangency constraint, hence the
# second projection behind it.
#
function _sphere_step_limiter!(u, integrator, p::St_sphere_ode_params, t)
    sphere_filter!(u, p.mesh, p.metrics, p.sp)
    p.lproject && project_momentum_to_sphere!(u, p.mesh; ivar = 2)
    return nothing
end


#
# Rebuild whatever `inputs[:ode_solver]` names, with our two limiters installed.
#
# The deck stores an already-constructed algorithm (`SSPRK33()`), and the
# limiters are constructor arguments, so the only way to inject them is to
# reconstruct from the type. Every SSPRK integrator takes them positionally as
# `Alg(stage_limiter!, step_limiter!)`.
#
function _sphere_alg_with_limiters(alg)

    hasfield(typeof(alg), :stage_limiter!) && hasfield(typeof(alg), :step_limiter!) ||
        error(string(" # ERROR sphere_time_loop.jl: :ode_solver => ", nameof(typeof(alg)),
                     " has no stage_limiter!/step_limiter!, and the shell needs one.\n",
                     " #   The Lagrange projection P = I - xxᵀ/r² has to run after EVERY RK\n",
                     " #   stage or the next stage evaluates fluxes on an off-manifold state.\n",
                     " #   Use an SSPRK integrator (SSPRK33, SSPRK43, SSPRK54, ...), which is\n",
                     " #   also what the SSP property of the scheme requires."))

    Alg = Base.typename(typeof(alg)).wrapper
    return Alg(_sphere_stage_limiter!, _sphere_step_limiter!)
end


#---------------------------------------------------------------------------------
# The per-step monitor: diagnostics lines and VTK output.
#
# A callable struct rather than a closure so the callback is concretely typed
# and its buffers (ζ, ζ0) are allocated once. The GATING lives in the
# DiscreteCallback's condition (St_sphere_monitor_due below) so that the affect!
# body is entered only on the ~35 print steps and 24 output steps instead of all
# 6913 — the conditions themselves are the old hand-written loop's, unchanged.
#---------------------------------------------------------------------------------
mutable struct St_sphere_monitor{TQ, TMesh, TMetrics, TParams, TIn, TSVT}
    q::TQ
    mesh::TMesh
    metrics::TMetrics
    sp::TParams
    inputs::TIn
    SVT::TSVT
    OUTPUT_DIR::String
    ζ::Vector{Float64}
    ζ0::Vector{Float64}
    mass0::Float64
    ener0::Float64
    npoin::Int
    nsteps::Int
    nprint::Int
    touts::Vector{Float64}     # the output times, ascending (see sphere_output_times)
    iout::Int
    verbose::Bool
end

#
# Is this step a print step or an output step? Called once per step; keeping the
# test here rather than inside affect! is what stops the integrator from
# entering the callback body 6913 times to do nothing.
#
struct St_sphere_monitor_due{TMon}
    mon::TMon
end

function (due::St_sphere_monitor_due)(u, t, integrator)
    mon   = due.mon
    istep = integrator.stats.naccept
    return (mon.verbose && (istep % mon.nprint == 0 || istep >= mon.nsteps)) ||
           (mon.iout < length(mon.touts) &&
            t >= mon.touts[mon.iout+1] - 1.0e-9*integrator.dt)
end

function (mon::St_sphere_monitor)(integrator)

    u     = integrator.u
    t     = integrator.t
    istep = integrator.stats.naccept

    #--- diagnostics
    if mon.verbose && (istep % mon.nprint == 0 || istep >= mon.nsteps)
        mass, ener, drift = sphere_diagnostics(u, mon.mesh, mon.metrics)
        # max|ζ| is the instability indicator: for the Galewsky test the
        # height field barely moves while the perturbation grows, so |ζ|
        # rising by orders of magnitude is how the barotropic instability
        # announces itself long before it is visible in h.
        sphere_relative_vorticity!(mon.ζ, u, mon.mesh, mon.metrics, mon.sp)
        ζmax = maximum(abs, mon.ζ)
        dζ   = 0.0
        for ip = 1:mon.npoin
            dζ = max(dζ, abs(mon.ζ[ip] - mon.ζ0[ip]))
        end
        @printf(" #   step %6d  t = %10.1f s (%6.3f d)  δmass/mass = %9.2e  δE/E = %9.2e  |(φu)·x̂| = %9.2e  max|ζ| = %9.3e  max|ζ-ζ₀| = %9.3e\n",
                istep, t, t/86400, (mass-mon.mass0)/mon.mass0, (ener-mon.ener0)/mon.ener0,
                drift, ζmax, dζ)
        #
        # What the model is actually doing. Reported as the range over the
        # MOMENTUM slot, which is the one every configuration switches on, and
        # against the C2·Δ·(|u|+c) cap: μ saturating at the cap means DynSGS has
        # degraded to first-order-upwind dissipation everywhere and is no longer
        # discriminating, which is the failure mode to watch for. In the regime
        # the model is designed for it sits an order of magnitude below it.
        #
        if mon.sp.ldsgs
            νmom = @view mon.sp.μ_dsgs[:, min(2, size(mon.sp.μ_dsgs,2))]
            @printf(" #     DynSGS ν(momentum): mean = %9.3e  max = %9.3e m²/s  (cap = %9.3e)\n",
                    sum(νmom)/length(νmom), maximum(νmom),
                    sphere_dsgs_mu_bound(u, mon.mesh, mon.sp))
        end
        flush(stdout)
        isfinite(mass) || error(" # ERROR sphere_time_loop.jl: the solution has gone non-finite. Reduce :cfl, or switch on :lfilter and/or :lvisc (with :μ > 0).")
    end

    #--- output
    #
    # `while`, not `if`: two requested output times can fall inside one step
    # (:diagnostics_at_times is arbitrary, and the shell takes a fixed Δt with
    # no tstops), and skipping the second would silently drop a file the deck
    # asked for. Each pass writes its own numbered file, all with the state at
    # this step.
    #
    while mon.iout < length(mon.touts) &&
          t >= mon.touts[mon.iout+1] - 1.0e-9*integrator.dt
        mon.iout += 1
        # _sphere_write! reads q.qn (through user_uout!), and the integrator
        # works on its own copy of the state, so refresh it first.
        copyto!(mon.q.qn, u)
        sphere_relative_vorticity!(mon.ζ, u, mon.mesh, mon.metrics, mon.sp)
        _sphere_write!(mon.q, mon.mesh, mon.inputs, mon.OUTPUT_DIR, mon.iout, t, mon.SVT;
                       verbose = mon.verbose,
                       extra = ("vorticity" => mon.ζ,
                                sphere_dsgs_extra(mon.sp, mon.mesh, mon.q.qvars)...))
    end

    return nothing
end


#---------------------------------------------------------------------------------
# sphere_time_loop!(mesh, metrics, sp, q, inputs, OUTPUT_DIR)
#
# Advances q from :tinit to :tend and writes VTK output. Returns the final time.
#---------------------------------------------------------------------------------
function sphere_time_loop!(mesh::St_mesh,
                           metrics::St_sphere_metrics,
                           sp::St_sphere_params,
                           q,
                           inputs,
                           OUTPUT_DIR::String;
                           verbose = true)
    #
    # FUNCTION BARRIER. `inputs` is a Dict{Symbol,Any} (or a NamedTuple), so
    # inputs[:SOL_VARS_TYPE] is inferred as Any; passing it positionally into
    # _sphere_march! specializes the marching code on its concrete type instead
    # of compiling the whole loop against Any.
    #
    return _sphere_march!(mesh, metrics, sp, q, inputs, OUTPUT_DIR,
                          inputs[:SOL_VARS_TYPE]; verbose = verbose)
end


function _sphere_march!(mesh::St_mesh,
                        metrics::St_sphere_metrics,
                        sp::St_sphere_params,
                        q,
                        inputs,
                        OUTPUT_DIR::String,
                        SVT;
                        verbose = true)

    neqs::Int  = Int(sp.neqs)
    npoin::Int = Int(mesh.npoin)

    #
    # Read out of `inputs` ONCE, into DECLARED locals. `inputs` is a
    # Dict{Symbol,Any}, so `Float64(inputs[:tend])` is a dynamic call and infers
    # ::Any all by itself — the `::Float64` declaration is what actually pins it
    # down, not the conversion. Δt matters most: it multiplies the RHS in every
    # stage of every step, so an ::Any Δt boxes each of those products.
    #
    t::Float64    = Float64(get(inputs, :tinit, 0.0))
    tend::Float64 = Float64(inputs[:tend])

    lproject = get(inputs, :llagrange_projection, true) == true

    #--- time step
    #
    # :lcfl_dt (default true) selects the CFL step; set it false to use the
    # deck's :Δt. This is an EXPLICIT switch rather than "omit :Δt and you get
    # the CFL step", because mod_inputs_user_inputs! fills a missing :Δt with
    # 0.1 s (mod_inputs.jl) — so an absent :Δt is indistinguishable from a
    # deliberate 0.1, and taking it literally here silently turns a one-day run
    # into 864 000 steps.
    #
    cfl::Float64 = Float64(get(inputs, :cfl, 0.35))
    lcfl         = get(inputs, :lcfl_dt, true) == true

    # the largest artificial viscosity in play, for the diffusive branch of the
    # CFL condition (0 when :lvisc is off, and then that branch is skipped).
    #
    # Under DynSGS sp.μ is a dimensionless per-equation multiplier, not a
    # viscosity, so the bound has to come from the model's own wave-speed cap —
    # see sphere_dsgs_mu_bound.
    μmax::Float64 = sp.ldsgs ? maximum(sp.μ)*sphere_dsgs_mu_bound(q.qn, mesh, sp) :
                               Float64(maximum(sp.μ))

    local Δt::Float64
    if lcfl
        Δt = sphere_cfl_dt(q.qn, mesh, metrics; cfl = cfl, μmax = μmax)
    else
        haskey(inputs, :Δt) && Float64(inputs[:Δt]) > 0 ||
            error(" # ERROR sphere_time_loop.jl: :lcfl_dt => false requires a positive :Δt in user_inputs.jl.")
        Δt = Float64(inputs[:Δt])
    end

    nsteps::Int = max(1, ceil(Int, (tend - t)/Δt))
    Δt          = (tend - t)/nsteps     # land exactly on tend

    #
    # A guard, not a limit: an accidental Δt (see above) turns into hours of
    # silent grinding, which is far worse than an error that says what to do.
    #
    maxsteps = Int(get(inputs, :max_steps, 500_000))
    nsteps <= maxsteps ||
        error(string(" # ERROR sphere_time_loop.jl: ", nsteps, " steps requested (Δt = ", Δt,
                     " s to t = ", tend, " s), above :max_steps = ", maxsteps, ".\n",
                     " #   Set :lcfl_dt => true to take the CFL step (",
                     @sprintf("%.2f", sphere_cfl_dt(q.qn, mesh, metrics; cfl = cfl, μmax = μmax)),
                     " s here), or raise :max_steps if you really mean it."))

    #
    # DynSGS reads Δt (the BDF2 denominator) and the history buffers from sp
    # inside the RHS, so both have to be set before the first stage. Δt is
    # published only once the "land exactly on tend" adjustment above is done —
    # a residual built on the requested step rather than the taken one is wrong
    # by that rounding on every step of the run.
    #
    sp.Δt[] = Δt
    sphere_dsgs_init_history!(sp, q.qn, npoin, neqs)

    #--- the integrator, with the projection and the filter installed
    haskey(inputs, :ode_solver) ||
        error(" # ERROR sphere_time_loop.jl: :ode_solver is missing from user_inputs.jl.")
    alg = _sphere_alg_with_limiters(inputs[:ode_solver])

    #
    # Fixed-step, always. The SSPRK schemes carry no embedded error estimate, and
    # Δt here is the CFL step of an explicit scheme rather than an accuracy
    # choice, so there is nothing for a step controller to control.
    #
    get(inputs, :ode_adaptive_solver, false) == true &&
        @warn " # sphere_time_loop.jl: :ode_adaptive_solver is ignored on the shell; the CFL step is fixed."

    #--- output times
    touts  = sphere_output_times(inputs, t, tend)
    nprint = Int(get(inputs, :ndiagnostics_prints, max(1, nsteps ÷ 20)))

    if verbose
        println(" # ")
        println(" # TIME INTEGRATION (", nameof(typeof(alg)), ", projection every stage) ..........")
        @printf(" #   Δt = %.4f s ; %d steps to t = %.1f s (%.3f days) ; CFL = %.2f\n",
                Δt, nsteps, tend, tend/86400, cfl)
        println(" #   filter: ", sp.lfilter ? "ON" : "OFF")
        if sp.ldsgs
            @printf(" #   DynSGS: ON, C1 = %.3g, C2 = %.3g, per-equation multipliers = [%s]\n",
                    sp.C1, sp.C2, join((@sprintf("%.3g", m) for m in sp.μ), ", "))
            @printf(" #     Δ_elem = [%.3g, %.3g] m ; ν ≤ %.3g m²/s (the C2·Δ·(|u|+c) cap at t = 0)\n",
                    minimum(sp.Δelem), maximum(sp.Δelem), μmax)
        elseif sp.lvisc
            @printf(" #   artificial viscosity: ON, ν = [%s] m²/s per equation\n",
                    join((@sprintf("%.3g", ν) for ν in sp.μ), ", "))
        else
            println(" #   artificial viscosity: OFF")
        end
        sp.lfilter || sp.lvisc ||
            @warn " # sphere_time_loop.jl: BOTH :lfilter and :lvisc are off. The inviscid unfiltered shell solution grows grid-scale modes and is expected to blow up (Marras/Kopera/Giraldo 2015, section 4.2)."
        #
        # Say the schedule out loud. An empty one is a legitimate request
        # (:lgrid_only-style runs, a timing measurement) and an easy accident —
        # it is how ShallowWater/SWsphere silently produced nothing under CI —
        # so it is reported rather than left to be discovered by an empty
        # directory at the end of the run.
        #
        _fmt = get(inputs, :outformat, VTK())
        if isempty(touts)
            println(" #   output: none scheduled — no :ndiagnostics_outputs > 0 and no :diagnostics_at_times",
                    get(inputs, :lwrite_initial, true) == true ?
                        " (the initial state is still written)" : ", and :lwrite_initial is false")
        elseif !(_fmt isa VTK || _fmt isa HDF5)
            println(" #   output: ", length(touts), " time(s) scheduled but :outformat => ",
                    nameof(typeof(_fmt)), " writes nothing on the shell (VTK and HDF5 are implemented)")
        else
            @printf(" #   output: %d × %s, first at t = %.1f s, last at t = %.1f s\n",
                    length(touts), nameof(typeof(_fmt)), first(touts), last(touts))
        end
    end

    # relative vorticity: the field the Galewsky test is judged on. h barely
    # moves while the instability develops, so a height plot makes the roll-up
    # nearly invisible; ζ shows it.
    ζ  = zeros(Float64, npoin)
    # ζ at t=0, so the diagnostics can report the PERTURBATION vorticity.
    # max|ζ| alone is dominated by the jet's own shear (~1e-4 everywhere in the
    # band from the start) and barely moves; max|ζ-ζ₀| isolates what the
    # instability is actually doing and grows by orders of magnitude.
    ζ0 = zeros(Float64, npoin)

    mass0, ener0, _ = sphere_diagnostics(q.qn, mesh, metrics)

    # the initial condition
    sphere_relative_vorticity!(ζ, q.qn, mesh, metrics, sp)
    copyto!(ζ0, ζ)
    # μ_dsgs is all zeros at t = 0 (no RHS has been assembled yet), and it is
    # written anyway so the field exists on every file of the series and
    # ParaView can animate it without a gap at frame 0.
    _sphere_write!(q, mesh, inputs, OUTPUT_DIR, 0, t, SVT; verbose = verbose,
                   lwrite = get(inputs, :lwrite_initial, true) == true,
                   extra = ("vorticity" => ζ, sphere_dsgs_extra(sp, mesh, q.qvars)...))

    params = St_sphere_ode_params(mesh, metrics, sp, q.qe, SVT, lproject, Ref(0.0),
                                  npoin, neqs)

    monitor = St_sphere_monitor(q, mesh, metrics, sp, inputs, SVT, OUTPUT_DIR,
                                ζ, ζ0, mass0, ener0, npoin, nsteps, nprint,
                                touts, 0, verbose)

    # save_positions = (false,false): the monitor only reads the state, so there
    # is no need to snapshot it around the callback.
    cb = DiscreteCallback(St_sphere_monitor_due(monitor), monitor;
                          save_positions = (false, false))

    prob = ODEProblem{true, FullSpecialize}(_sphere_ode_rhs!,
                                            copy(q.qn), (t, tend), params)

    solution = solve(prob, alg;
                     dt             = Δt,
                     adaptive       = false,
                     callback       = cb,
                     save_everystep = false,
                     save_start     = false,
                     saveat         = Float64[])

    copyto!(q.qn, solution.u[end])
    tfinal = Float64(solution.t[end])

    #
    # DID IT ACTUALLY GET THERE?
    #
    # When the solution blows up, OrdinaryDiffEq does not throw: it stops early
    # and RETURNS, with a retcode saying why (:DtNaN, :Unstable, :MaxIters …).
    # Without this check the run then printed
    #
    #   # TIME INTEGRATION ......... DONE
    #   #   final: t = 173274.9 s ; δmass/mass = -6.643e+144 ; δE/E = NaN
    #
    # — the word DONE on a solution that had gone to NaN a day short of :tend.
    # The diagnostics callback only catches it if a print step happens to fall
    # after the blow-up and before the integrator gives up, which is a race:
    # here the last print was at 1.736 d, the run died at 2.005 d, and nothing
    # said so. An early stop is a FAILURE and has to be reported as one.
    #
    if !successful_retcode(solution)
        error(string(" # ERROR sphere_time_loop.jl: the integration stopped at t = ",
                     @sprintf("%.1f", tfinal), " s of the requested ", @sprintf("%.1f", tend),
                     " s (retcode :", solution.retcode, ").\n",
                     " #   The solution blew up. On the shell that almost always means too\n",
                     " #   little stabilisation rather than too large a step, so try, in order:\n",
                     " #     * :lfilter => true — the modal filter damps ALL the fields;\n",
                     " #     * :ivisc_equations => [1,2,3,4] — with the filter off, leaving the\n",
                     " #       continuity equation out gives φ no dissipation at all, and CG\n",
                     " #       has no upwinding to supply any;\n",
                     " #     * a larger :μ, or a smaller :cfl if the step really is the problem."))
    end

    if verbose
        mass, ener, drift = sphere_diagnostics(q.qn, mesh, metrics)
        println(" # TIME INTEGRATION ............................................ DONE")
        @printf(" #   final: t = %.1f s ; δmass/mass = %.3e ; δE/E = %.3e ; max drift removed = %.3e\n",
                tfinal, (mass-mass0)/mass0, (ener-ener0)/ener0, params.driftmax[])
    end

    return tfinal
end


#
# One output time. Numbered files, as the flat cases do.
#
# THE FORMAT IS THE DECK'S :outformat, not VTK unconditionally. VTK is what a
# human wants and what every deck in problems/ asks for, but it is not the only
# consumer: CI_MODE forces :outformat => "hdf5" because the reference comparison
# in test/ci_compare.jl reads HDF5, and a shell case that wrote .vtu anyway
# produced output the suite could not see. (That was the second half of why
# ShallowWater/SWsphere could not run under CI — see sphere_output_times for the
# first.) Anything else, including the default NONE(), writes nothing; a deck
# that names no format is asking for none.
#
function _sphere_write!(q, mesh::St_mesh, inputs, OUTPUT_DIR::String,
                        iout::Int, t::Float64, SVT; verbose = true, lwrite = true,
                        extra = nothing)

    lwrite || return nothing

    outformat = get(inputs, :outformat, VTK())

    if outformat isa HDF5
        _sphere_write_hdf5!(q, mesh, inputs, OUTPUT_DIR, iout, t; verbose = verbose)
        return nothing
    elseif !(outformat isa VTK)
        return nothing
    end

    # λ and φ let user_uout! write the velocity PROJECTED ONTO THE SHELL
    # (zonal, meridional, radial) instead of raw Cartesian components.
    for ip = 1:mesh.npoin
        user_uout!(ip, SVT, @view(q.qout[ip,:]), @view(q.qn[ip,:]), @view(q.qe[ip,:]);
                   lon = mesh.lon[ip], lat = mesh.lat[ip])
    end

    fname = iout == 0 ? "sphere_grid_ho" : @sprintf("sphere_%04d", iout)
    write_vtk_sphere_grid(mesh, fname, OUTPUT_DIR; q = q, extra = extra, verbose = false)

    verbose && @printf(" #   wrote %s.vtu at t = %.1f s\n", fname, t)
    return nothing
end


#
# The HDF5 branch: the CONSERVATIVE state, through the same write_hdf5 the flat
# cases use, so the files are the var_<ieq>_<rank>.h5 + t.h5 set that
# test/ci_compare.jl knows how to read and diff.
#
# Conservative (q.qn), not the shell-projected primitives of q.qout: this is a
# reference solution, so what belongs in it is what the scheme integrates, and
# it is also what the flat path writes (TimeIntegrators.jl hands write_output
# the solution vector with params.qp.qvars as the names).
#
# q.qn is npoin × (neqs+1) and write_hdf5 wants a flat npoin*nvar vector laid
# out variable by variable, so the extra trailing column is dropped here. The
# copy is deliberate: a @view of the first neqs columns is not contiguous, and
# `vec` of it would not be the layout write_hdf5 indexes into.
#
function _sphere_write_hdf5!(q, mesh::St_mesh, inputs, OUTPUT_DIR::String,
                             iout::Int, t::Float64; verbose = true)

    npoin = Int(mesh.npoin)
    neqs  = Int(q.neqs)

    title = @sprintf("Spherical shell solution at t=%.6f", t)
    write_hdf5(mesh.SD, mesh, vec(q.qn[1:npoin, 1:neqs]), q.qe, t, title,
               OUTPUT_DIR, inputs, q.qvars;
               iout = iout, nvar = neqs, case = string(get(inputs, :case, "")))

    verbose && @printf(" #   wrote var_1..%d_<rank>.h5 at t = %.1f s\n", neqs, t)
    return nothing
end
