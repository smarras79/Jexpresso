#---------------------------------------------------------------------------------
# sphere_time_loop.jl
#
# Time integration of the shallow water equations on the spherical shell.
#
# SSP-RK3 (Shu-Osher), the scheme of Marras, Kopera & Giraldo (2015), QJRMS
# 141: 1727-1739, section 4.2:
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
# WHY NOT OrdinaryDiffEq / Jexpresso's time_loop!
# -----------------------------------------------
# time_loop! is built around params_setup and a St_mesh; the shell path has
# neither (it bypasses sem_setup, for the reasons in sphere_mesh.jl). SSP-RK3 in
# thirty lines is also exactly what the paper specifies, and it gives the
# per-stage hook the projection needs, which a black-box integrator would only
# provide through a stage limiter. When the shell is folded into sem_setup this
# should become an ordinary Jexpresso integrator with a stage callback.
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
#   Δt = CFL Δmin / max(|u| + √φ)
#
# √φ = √(gh) is the gravity-wave speed, |u| the advective speed; Δmin is the
# smallest LGL node spacing, which for a spectral element clusters like 1/nop²
# towards the edges and is what actually limits the step.
#---------------------------------------------------------------------------------
function sphere_cfl_dt(q::AbstractMatrix{TF}, mesh::St_mesh_sphere{TF,TI},
                       metrics::St_sphere_metrics{TF}; cfl = 0.35) where {TF,TI}

    cmax = 0.0
    @inbounds for ip = 1:mesh.npoin
        φ = q[ip,1]
        φ > 0 || error(string(" # ERROR sphere_time_loop.jl: non-positive geopotential φ = ", φ,
                              " at node ", ip, ". The state has gone unphysical."))
        umag = sqrt(q[ip,2]^2 + q[ip,3]^2 + q[ip,4]^2)/φ
        cmax = max(cmax, umag + sqrt(φ))
    end
    return cfl*metrics.Δmin/cmax
end


#---------------------------------------------------------------------------------
# Conserved quantities and the constraint drift.
#---------------------------------------------------------------------------------
function sphere_diagnostics(q::AbstractMatrix{TF}, mesh::St_mesh_sphere{TF,TI},
                            metrics::St_sphere_metrics{TF}) where {TF,TI}

    mass = 0.0
    ener = 0.0
    @inbounds for ip = 1:mesh.npoin
        φ  = q[ip,1]
        m2 = q[ip,2]^2 + q[ip,3]^2 + q[ip,4]^2
        mass += metrics.M[ip]*φ
        ener += metrics.M[ip]*(0.5*m2/φ + 0.5*φ*φ)
    end
    drift = sphere_normal_momentum(q, mesh; ivar = 2)

    return mass, ener, drift
end


#---------------------------------------------------------------------------------
# sphere_time_loop!(mesh, metrics, sp, q, inputs, OUTPUT_DIR)
#
# Advances q from :tinit to :tend and writes VTK output. Returns the final time.
#---------------------------------------------------------------------------------
function sphere_time_loop!(mesh::St_mesh_sphere,
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
    # PERF NOTE, measured rather than assumed: sphere_rhs!, sphere_filter!,
    # project_momentum_to_sphere! and sphere_diagnostics are each allocation-free
    # once warm, and a 95-step march takes 0.91 s at npoin = 5402 (≈10 ms/step,
    # which is the three RHS evaluations). The loop as a whole still allocates
    # ~6 MB/step from somewhere not yet pinned down; GC accounts for ~6% of
    # runtime, so it is not currently worth chasing, but it is not zero either
    # and this barrier alone did NOT remove it.
    #
    return _sphere_march!(mesh, metrics, sp, q, inputs, OUTPUT_DIR,
                          inputs[:SOL_VARS_TYPE]; verbose = verbose)
end


function _sphere_march!(mesh::St_mesh_sphere,
                        metrics::St_sphere_metrics,
                        sp::St_sphere_params,
                        q,
                        inputs,
                        OUTPUT_DIR::String,
                        SVT;
                        verbose = true)

    neqs = sp.neqs
    npoin = Int(mesh.npoin)

    t    = Float64(get(inputs, :tinit, 0.0))
    tend = Float64(inputs[:tend])

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
    cfl  = Float64(get(inputs, :cfl, 0.35))
    lcfl = get(inputs, :lcfl_dt, true) == true

    if lcfl
        Δt = sphere_cfl_dt(q.qn, mesh, metrics; cfl = cfl)
    else
        haskey(inputs, :Δt) && Float64(inputs[:Δt]) > 0 ||
            error(" # ERROR sphere_time_loop.jl: :lcfl_dt => false requires a positive :Δt in user_inputs.jl.")
        Δt = Float64(inputs[:Δt])
    end

    nsteps = max(1, ceil(Int, (tend - t)/Δt))
    Δt     = (tend - t)/nsteps          # land exactly on tend

    #
    # A guard, not a limit: an accidental Δt (see above) turns into hours of
    # silent grinding, which is far worse than an error that says what to do.
    #
    maxsteps = Int(get(inputs, :max_steps, 500_000))
    nsteps <= maxsteps ||
        error(string(" # ERROR sphere_time_loop.jl: ", nsteps, " steps requested (Δt = ", Δt,
                     " s to t = ", tend, " s), above :max_steps = ", maxsteps, ".\n",
                     " #   Set :lcfl_dt => true to take the CFL step (",
                     @sprintf("%.2f", sphere_cfl_dt(q.qn, mesh, metrics; cfl = cfl)),
                     " s here), or raise :max_steps if you really mean it."))

    #--- output times
    nout   = Int(get(inputs, :ndiagnostics_outputs, 10))
    outdt  = nout > 0 ? (tend - t)/nout : Inf
    tnext  = t + outdt
    iout   = 0

    nprint = Int(get(inputs, :ndiagnostics_prints, max(1, nsteps ÷ 20)))

    if verbose
        println(" # ")
        println(" # TIME INTEGRATION (SSP-RK3, projection every stage) ..........")
        @printf(" #   Δt = %.4f s ; %d steps to t = %.1f s (%.3f days) ; CFL = %.2f\n",
                Δt, nsteps, tend, tend/86400, cfl)
        println(" #   filter: ", sp.lfilter ? "ON" : "OFF")
    end

    RHS = zeros(Float64, npoin, neqs)
    q1  = copy(q.qn)
    q2  = copy(q.qn)
    # relative vorticity: the field the Galewsky test is judged on. h barely
    # moves while the instability develops, so a height plot makes the roll-up
    # nearly invisible; ζ shows it.
    ζ   = zeros(Float64, npoin)
    # ζ at t=0, so the diagnostics can report the PERTURBATION vorticity.
    # max|ζ| alone is dominated by the jet's own shear (~1e-4 everywhere in the
    # band from the start) and barely moves; max|ζ-ζ₀| isolates what the
    # instability is actually doing and grows by orders of magnitude.
    ζ0  = zeros(Float64, npoin)
    qn  = q.qn
    qe  = q.qe

    mass0, ener0, _ = sphere_diagnostics(qn, mesh, metrics)
    h0max = maximum(@view qn[1:npoin, 1])

    # the initial condition
    sphere_relative_vorticity!(ζ, q.qn, mesh, metrics, sp)
    copyto!(ζ0, ζ)
    _sphere_write!(q, mesh, inputs, OUTPUT_DIR, iout, t, SVT; verbose = verbose,
                   lwrite = get(inputs, :lwrite_initial, true) == true,
                   extra = ("vorticity" => ζ,))

    driftmax = 0.0

    for istep = 1:nsteps

        #--- stage 1
        sphere_rhs!(RHS, qn, qe, mesh, metrics, sp, SVT)
        @inbounds for ieq = 1:neqs, ip = 1:npoin
            q1[ip,ieq] = qn[ip,ieq] + Δt*RHS[ip,ieq]
        end
        lproject && (driftmax = max(driftmax, project_momentum_to_sphere!(q1, mesh; ivar = 2)))

        #--- stage 2
        sphere_rhs!(RHS, q1, qe, mesh, metrics, sp, SVT)
        @inbounds for ieq = 1:neqs, ip = 1:npoin
            q2[ip,ieq] = 0.75*qn[ip,ieq] + 0.25*(q1[ip,ieq] + Δt*RHS[ip,ieq])
        end
        lproject && (driftmax = max(driftmax, project_momentum_to_sphere!(q2, mesh; ivar = 2)))

        #--- stage 3
        sphere_rhs!(RHS, q2, qe, mesh, metrics, sp, SVT)
        @inbounds for ieq = 1:neqs, ip = 1:npoin
            qn[ip,ieq] = (qn[ip,ieq] + 2.0*(q2[ip,ieq] + Δt*RHS[ip,ieq]))/3.0
        end
        lproject && (driftmax = max(driftmax, project_momentum_to_sphere!(qn, mesh; ivar = 2)))

        #--- stabilization
        sphere_filter!(qn, mesh, metrics, sp)
        lproject && project_momentum_to_sphere!(qn, mesh; ivar = 2)

        t += Δt

        #--- diagnostics
        if verbose && (istep % nprint == 0 || istep == nsteps)
            mass, ener, drift = sphere_diagnostics(qn, mesh, metrics)
            # max|ζ| is the instability indicator: for the Galewsky test the
            # height field barely moves while the perturbation grows, so |ζ|
            # rising by orders of magnitude is how the barotropic instability
            # announces itself long before it is visible in h.
            sphere_relative_vorticity!(ζ, qn, mesh, metrics, sp)
            ζmax = maximum(abs, ζ)
            dζ   = 0.0
            for ip = 1:npoin
                dζ = max(dζ, abs(ζ[ip] - ζ0[ip]))
            end
            @printf(" #   step %6d  t = %10.1f s (%6.3f d)  δmass/mass = %9.2e  δE/E = %9.2e  |(φu)·x̂| = %9.2e  max|ζ| = %9.3e  max|ζ-ζ₀| = %9.3e\n",
                    istep, t, t/86400, (mass-mass0)/mass0, (ener-ener0)/ener0, drift, ζmax, dζ)
            flush(stdout)
            isfinite(mass) || error(" # ERROR sphere_time_loop.jl: the solution has gone non-finite. Reduce :cfl, or switch :lfilter on.")
        end

        #--- output
        if t >= tnext - 1.0e-9*Δt && iout < nout
            iout += 1
            sphere_relative_vorticity!(ζ, qn, mesh, metrics, sp)
            _sphere_write!(q, mesh, inputs, OUTPUT_DIR, iout, t, SVT; verbose = verbose,
                           extra = ("vorticity" => ζ,))
            tnext += outdt
        end
    end

    if verbose
        mass, ener, drift = sphere_diagnostics(qn, mesh, metrics)
        println(" # TIME INTEGRATION ............................................ DONE")
        @printf(" #   final: t = %.1f s ; δmass/mass = %.3e ; δE/E = %.3e ; max drift removed = %.3e\n",
                t, (mass-mass0)/mass0, (ener-ener0)/ener0, driftmax)
    end

    return t
end


#
# Refresh the primitive output fields through the case's own user_uout! and
# write one VTK file. Numbered files, one per output time, as the flat cases do.
#
function _sphere_write!(q, mesh::St_mesh_sphere, inputs, OUTPUT_DIR::String,
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

    verbose && @printf(" #   wrote %s.vtu at t = %.1f s\n", fname, t)
    return nothing
end
