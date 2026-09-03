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

    # The step has to be the SAME on every rank: cmax is a maximum over this
    # rank's nodes only, and ranks that disagree on Δt would march the same
    # equation at different speeds and exchange states from different times.
    comm = get_mpi_comm()
    MPI.Comm_size(comm) > 1 && (cmax = MPI.Allreduce(cmax, MPI.MAX, comm))

    Δt::TF = cfl*metrics.Δmin/cmax

    ν::TF = TF(μmax)
    if ν > zero(TF)
        Δt = min(Δt, cfl*metrics.Δmin*metrics.Δmin/ν)
    end

    return Δt
end


#---------------------------------------------------------------------------------
# Conserved quantities and the constraint drift.
#---------------------------------------------------------------------------------
function sphere_diagnostics(q::AbstractMatrix{TF}, mesh::St_mesh,
                            metrics::St_sphere_metrics{TF}) where {TF}

    comm   = get_mpi_comm()
    rank   = MPI.Comm_rank(comm)
    nparts = MPI.Comm_size(comm)

    # ∫φ and ∫E are sums over the OWNED nodes only. M carries the fully
    # assembled mass (build_sphere_metrics completes it across ranks), so a node
    # mirrored on three ranks would otherwise contribute three times and the
    # "conserved" mass would depend on the partition.
    mass, ener = _sphere_integrals(q, metrics.M, mesh.gip2owner, rank, Int(mesh.npoin))

    # ::TF for the same reason: mesh.x/y/z are ::Any, so the call inside
    # sphere_normal_momentum is dynamic and its result infers ::Any. A max needs
    # no ownership test — a mirrored node has the same value on both ranks.
    drift::TF = sphere_normal_momentum(q, mesh; ivar = 2)

    if nparts > 1
        mass  = MPI.Allreduce(mass,  MPI.SUM, comm)
        ener  = MPI.Allreduce(ener,  MPI.SUM, comm)
        drift = MPI.Allreduce(drift, MPI.MAX, comm)
    end

    return mass, ener, drift
end


#
# ∫φ dΩ and ∫(½φ|u|² + ½φ²) dΩ over this rank's own nodes. Typed barrier:
# mesh.gip2owner is an ::Any field and the loop would box on every node.
#
function _sphere_integrals(q::AbstractMatrix{TF}, M::AbstractVector{TF},
                           gip2owner, rank::Int, npoin::Int) where {TF}
    mass = zero(TF)
    ener = zero(TF)
    @inbounds for ip = 1:npoin
        gip2owner[ip] == rank || continue
        φ  = q[ip,1]
        m2 = q[ip,2]^2 + q[ip,3]^2 + q[ip,4]^2
        mass += M[ip]*φ
        ener += M[ip]*(0.5*m2/φ + 0.5*φ*φ)
    end
    return mass, ener
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
end


#
# R(q). OrdinaryDiffEq calls this once per stage.
#
function _sphere_ode_rhs!(du, u, p::St_sphere_ode_params, t)
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
# Then, on a forced run, the forcing for the NEXT step is drawn here, from the
# state that step will start from: the forcing is held fixed over a step
# (every stage sees the same u_F), which is what its amplitude assumes.
#
function _sphere_step_limiter!(u, integrator, p::St_sphere_ode_params, t)
    sphere_filter!(u, p.mesh, p.metrics, p.sp)
    p.lproject && project_momentum_to_sphere!(u, p.mesh; ivar = 2)
    p.sp.forcing === nothing ||
        sphere_forcing_update!(p.sp.forcing, u, integrator.dt, p.mesh, p.metrics, p.sp)
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
    nout::Int
    outdt::Float64
    tnext::Float64
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

#
# NOTE `mon.verbose` is deliberately NOT part of this test, although only rank 0
# ever prints. The body below reduces the diagnostics across ranks, so it is
# collective, and a condition that was true on rank 0 alone would leave that rank
# waiting inside MPI for ranks that never entered the callback. Whether a step is
# a print step is a property of the STEP, identical everywhere; who prints is
# decided inside.
#
function (due::St_sphere_monitor_due)(u, t, integrator)
    mon   = due.mon
    istep = integrator.stats.naccept
    return (istep % mon.nprint == 0 || istep >= mon.nsteps) ||
           (mon.iout < mon.nout && t >= mon.tnext - 1.0e-9*integrator.dt)
end

function (mon::St_sphere_monitor)(integrator)

    u     = integrator.u
    t     = integrator.t
    istep = integrator.stats.naccept

    #--- diagnostics. Collective on every rank; printed by whichever one was
    #    handed verbose = true.
    if istep % mon.nprint == 0 || istep >= mon.nsteps
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
        comm = get_mpi_comm()
        if MPI.Comm_size(comm) > 1
            ζmax = MPI.Allreduce(ζmax, MPI.MAX, comm)
            dζ   = MPI.Allreduce(dζ,   MPI.MAX, comm)
        end
        if mon.verbose
            @printf(" #   step %6d  t = %10.1f s (%6.3f d)  δmass/mass = %9.2e  δE/E = %9.2e  |(φu)·x̂| = %9.2e  max|ζ| = %9.3e  max|ζ-ζ₀| = %9.3e\n",
                    istep, t, t/86400, (mass-mon.mass0)/mon.mass0, (ener-mon.ener0)/mon.ener0,
                    drift, ζmax, dζ)
            flush(stdout)
        end
        # mass is a reduced quantity, so this fires on every rank at once.
        isfinite(mass) || error(" # ERROR sphere_time_loop.jl: the solution has gone non-finite. Reduce :cfl, or switch on :lfilter and/or :lvisc (with :μ > 0).")

        # On a forced run (Scott & Polvani 2007) the quantities the paper is
        # read in: the kinetic energy per unit mass, the velocity scale U and
        # the Rossby number Ro = U/(2aΩ) it implies, and how much energy the
        # forcing has actually been injecting relative to the ε₀ it was asked
        # for. Collective (reductions inside), so outside the verbose gate.
        if mon.sp.forcing !== nothing
            fd = sphere_forcing_diagnostics(u, mon.sp.forcing, mon.mesh, mon.metrics)
            if mon.verbose
                @printf(" #           forced: KE/mass = %9.3e m²/s²  U_rms = %8.3f m/s  Ro = %9.3e  ε_inj/ε₀ = %6.3f (mean)  %6.3f (now)  gain = %8.3f\n",
                        fd.KE, fd.U, fd.Ro, fd.Pmean_rel,
                        mon.sp.forcing.ε > 0 ? fd.Pinj/mon.sp.forcing.ε : NaN, fd.gain)
                flush(stdout)
            end
        end
    end

    #--- output
    if t >= mon.tnext - 1.0e-9*integrator.dt && mon.iout < mon.nout
        mon.iout += 1
        # _sphere_write! reads q.qn (through user_uout!), and the integrator
        # works on its own copy of the state, so refresh it first.
        copyto!(mon.q.qn, u)
        sphere_relative_vorticity!(mon.ζ, u, mon.mesh, mon.metrics, mon.sp)
        _sphere_write!(mon.q, mon.mesh, mon.inputs, mon.OUTPUT_DIR, mon.iout, t, mon.SVT;
                       verbose = mon.verbose,
                       extra   = _sphere_extra_fields(mon.ζ, mon.sp))
        _sphere_write_zonal_mean(u, mon.mesh, mon.metrics, mon.inputs, mon.OUTPUT_DIR,
                                 mon.iout, t, mon.sp; verbose = mon.verbose)
        mon.tnext += mon.outdt
    end

    return nothing
end


#
# The point fields written next to the state: the relative vorticity always,
# plus the vorticity forcing F = ∇ₛ²ψ_F on a forced run, so the scale and the
# strength of what drives the flow can be looked at on the same file.
#
_sphere_extra_fields(ζ, sp) = sp.forcing === nothing ?
    ("vorticity" => ζ,) :
    ("vorticity" => ζ, "vorticity_forcing" => sp.forcing.F)


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
    # CFL condition (0 when :lvisc is off, and then that branch is skipped)
    μmax::Float64 = Float64(maximum(sp.μ))

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
    nout   = Int(get(inputs, :ndiagnostics_outputs, 10))
    outdt  = nout > 0 ? (tend - t)/nout : Inf
    nprint = Int(get(inputs, :ndiagnostics_prints, max(1, nsteps ÷ 20)))

    if verbose
        println(" # ")
        println(" # TIME INTEGRATION (", nameof(typeof(alg)), ", projection every stage) ..........")
        @printf(" #   Δt = %.4f s ; %d steps to t = %.1f s (%.3f days) ; CFL = %.2f\n",
                Δt, nsteps, tend, tend/86400, cfl)
        println(" #   filter: ", sp.lfilter ? "ON" : "OFF")
        if sp.lvisc
            @printf(" #   artificial viscosity: ON, ν = [%s] m²/s per equation\n",
                    join((@sprintf("%.3g", ν) for ν in sp.μ), ", "))
        else
            println(" #   artificial viscosity: OFF")
        end
        sp.lfilter || sp.lvisc ||
            @warn " # sphere_time_loop.jl: BOTH :lfilter and :lvisc are off. The inviscid unfiltered shell solution grows grid-scale modes and is expected to blow up (Marras/Kopera/Giraldo 2015, section 4.2)."
        if sp.forcing !== nothing
            @printf(" #   forcing: ON (Scott & Polvani 2007), ε₀ = %.3e m²/s³, n ∈ [%d, %d], %s; drawn every step\n",
                    sp.forcing.ε, minimum(sp.forcing.n), maximum(sp.forcing.n),
                    sp.forcing.τ > 0 ? @sprintf("τ = %.3e s", sp.forcing.τ) : "white in time")
        end
    end

    # On a forced run the first step needs a forcing too: the step limiter draws
    # the forcing for the step AFTER the one that just completed, so the very
    # first draw happens here, from the initial state.
    sp.forcing === nothing ||
        sphere_forcing_update!(sp.forcing, q.qn, Δt, mesh, metrics, sp)

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
    _sphere_write!(q, mesh, inputs, OUTPUT_DIR, 0, t, SVT; verbose = verbose,
                   lwrite = get(inputs, :lwrite_initial, true) == true,
                   extra = _sphere_extra_fields(ζ, sp))
    get(inputs, :lwrite_initial, true) == true &&
        _sphere_write_zonal_mean(q.qn, mesh, metrics, inputs, OUTPUT_DIR, 0, t, sp;
                                 verbose = verbose)

    params = St_sphere_ode_params(mesh, metrics, sp, q.qe, SVT, lproject, Ref(0.0))

    monitor = St_sphere_monitor(q, mesh, metrics, sp, inputs, SVT, OUTPUT_DIR,
                                ζ, ζ0, mass0, ener0, npoin, nsteps, nprint,
                                nout, outdt, t + outdt, 0, verbose)

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

    # Collective, so outside the verbose gate: sphere_diagnostics reduces across
    # ranks and calling it on rank 0 alone would hang the others.
    mass, ener, drift = sphere_diagnostics(q.qn, mesh, metrics)
    driftmax = params.driftmax[]
    MPI.Comm_size(get_mpi_comm()) > 1 &&
        (driftmax = MPI.Allreduce(driftmax, MPI.MAX, get_mpi_comm()))

    if verbose
        println(" # TIME INTEGRATION ............................................ DONE")
        @printf(" #   final: t = %.1f s ; δmass/mass = %.3e ; δE/E = %.3e ; max drift removed = %.3e\n",
                tfinal, (mass-mass0)/mass0, (ener-ener0)/ener0, driftmax)
        # Where the results are. Repeated here so the location is on screen at
        # the end of the run instead of only next to the individual writes.
        println(" # Output written to: ", abspath(OUTPUT_DIR))
        ext = MPI.Comm_size(get_mpi_comm()) > 1 ? ".pvtu" : ".vtu"
        get(inputs, :lwrite_initial, true) == true &&
            println(" #   sphere_grid_ho", ext, "   grid + initial condition")
        @printf(" #   sphere_0001..%04d%s   %d output time%s\n",
                nout, ext, nout, nout == 1 ? "" : "s")
    end

    return tfinal
end


#
# Refresh the primitive output fields through the case's own user_uout! and
# write one VTK file. Numbered files, one per output time, as the flat cases do.
#
function _sphere_write!(q, mesh::St_mesh, inputs, OUTPUT_DIR::String,
                        iout::Int, t::Float64, SVT; verbose = true, lwrite = true,
                        extra = nothing)

    lwrite || return nothing

    # λ and φ let user_uout! write the velocity PROJECTED ONTO THE SHELL
    # (zonal, meridional, radial) instead of raw Cartesian components.
    for ip = 1:mesh.npoin
        user_uout!(ip, SVT, @view(q.qout[ip,:]), @view(q.qn[ip,:]), @view(q.qe[ip,:]);
                   lon = mesh.lon[ip], lat = mesh.lat[ip])
    end

    fname = iout == 0 ? "sphere_grid_ho" : @sprintf("sphere_%04d", iout)
    write_vtk_sphere_grid(mesh, fname, OUTPUT_DIR; q = q, extra = extra, verbose = false)

    # Print the FULL path, as the flat cases do in write_output (" # writing
    # <OUTPUT_DIR>/iter_N.pvtu ... DONE"). `abspath` because :output_dir is
    # typically relative ("./output"), so the bare string does not tell you
    # where the file actually landed.
    verbose && @printf(" #   writing %s%s at t = %.1f s ... DONE\n",
                       joinpath(abspath(OUTPUT_DIR), fname),
                       MPI.Comm_size(get_mpi_comm()) > 1 ? ".pvtu" : ".vtu", t)
    return nothing
end


#
# ZONAL MEANS, as a small text file next to each VTK output — one row per
# latitude band: lat [deg], ū (zonal, m/s), v̄ (meridional, m/s), h̄ = φ̄/g
# is NOT formed here (g belongs to the case), so φ̄ [m²/s²] is written instead,
# and the number of nodes in the band.
#
# Written when the deck sets :lzonal_mean => true, which the forced
# (Scott & Polvani) decks do: the paper's principal figures are exactly ū(φ,t),
# and ParaView will not produce a zonal average of an unstructured shell. The
# number of bands is :zonal_mean_nbins (default 90, i.e. 2°).
#
# Collective (sphere_zonal_mean reduces across ranks); rank 0 writes.
#
function _sphere_write_zonal_mean(u, mesh::St_mesh, metrics::St_sphere_metrics,
                                  inputs, OUTPUT_DIR::String, iout::Int, t::Float64, sp;
                                  verbose = true)

    get(inputs, :lzonal_mean, false) == true || return nothing

    nbins = Int(get(inputs, :zonal_mean_nbins, 90))
    lat, ū, v̄, φ̄, cnt = sphere_zonal_mean(u, mesh, metrics; nbins = nbins)

    MPI.Comm_rank(get_mpi_comm()) == 0 || return nothing

    isdir(OUTPUT_DIR) || mkpath(OUTPUT_DIR)
    fname = joinpath(OUTPUT_DIR, @sprintf("zonal_mean_%04d.dat", iout))
    open(fname, "w") do io
        @printf(io, "# zonal means at t = %.6e s (output %d)\n", t, iout)
        @printf(io, "# %10s %16s %16s %16s %8s\n", "lat[deg]", "u_zonal[m/s]", "v_merid[m/s]", "phi[m2/s2]", "nodes")
        for b = 1:nbins
            @printf(io, "  %10.4f %16.8e %16.8e %16.8e %8d\n", lat[b], ū[b], v̄[b], φ̄[b], cnt[b])
        end
    end
    verbose && @printf(" #   writing %s ... DONE\n", fname)
    return nothing
end
