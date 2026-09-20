#---------------------------------------------------------------------------------
# pod.jl — POD of a running Jexpresso case: what to sample, when to sample it,
# how to reduce it across ranks, and what to leave on disk.
#
# The decomposition itself is in pod_core.jl, which knows nothing about meshes,
# cases or MPI. This file is the part that does.
#
# FOR ANY PROBLEM, FROM ONE LINE IN THE DECK. `:lpod => true` is the whole of
# the minimum: the fields default to every variable of the solution vector,
# whatever the case called them, the sampling window to the whole run, the
# sampling rate to :ndiagnostics_outputs, and the outputs to all of them. The
# case does not have to know POD exists — no user_*.jl is involved and no case
# supplies anything. Everything else is a :pod_* override, listed in
# pod_settings below and in docs/POD.md.
#
# WHAT GEOMETRY IT IS ON is decided here, not by the deck, because the deck
# already says: the number of space dimensions and :lspherical_shell. It selects
# how the modes are WRITTEN, which is the only place the difference shows —
#
#   1-D       a CSV of (x, mean, mode_001, …), sorted by x, plus a line plot.
#             A .vtu of a line grid is a worse format than a table for data
#             whose natural plot is a curve.
#   2-D       an unstructured .vtu of the quads the mesh is made of, plus a
#             contour map on an (x,y) raster.
#   3-D       an unstructured .vtu of hexahedra. No PNG: a 3-D mode needs a
#             slice or an isosurface, which is what ParaView is for, and any
#             projection this code picked would be the wrong one.
#   MANIFOLD  the .vtu of the shell itself (so the modes land on the geometry
#             they were computed on), plus an equirectangular map.
#
# The spectrum, the temporal coefficients, the CSVs and the .jld2 basis are the
# same in every dimension: they are properties of the decomposition, not of the
# grid.
#
# WHY IN THE RUN AND NOT AFTERWARDS FROM THE OUTPUT FILES. Two reasons, and both
# are about getting the right answer rather than about convenience. First, the
# inner product: POD is only optimal under the L² product (Eq. (4) of
# pod_core.jl), which needs the mass matrix — that is, the metrics of the run
# that produced the data. It is in no output file. Second, the sampling: POD
# wants snapshots uniformly spaced in time and dense enough to resolve the
# dynamics, while the output cadence is chosen to keep a movie small.
# :pod_nsnapshots therefore carries its own timer and the snapshots never touch
# the disk — unless :pod_write_snapshots asks for them, which is how a POD is
# recomputed offline over a different window without re-running the case.
#
# WHAT IT COSTS. One array of npoin × ncomp × nsnap Float64 per requested field,
# held for the run. :pod_max_memory_gb refuses, with the arithmetic in the
# message, anything that would not fit. The decomposition itself is O(N K²) and
# runs once, at the end.
#
# UNDER MPI the field stays distributed and nothing is gathered. The correlation
# matrix is a sum over nodes, so each rank forms its own K×K partial and one
# Allreduce completes it; every rank then solves the same small eigenproblem and
# builds its own slice of the modes. Two details make the result independent of
# the rank count: a node shared by several ranks is counted by its OWNER only
# (mesh.gip2owner, exactly as sphere_diagnostics does for the conserved
# integrals), and the sign convention of pod_core.jl is resolved globally.
#
# S. Marras & contributors
#---------------------------------------------------------------------------------

export St_pod_settings, St_pod_recorder, St_pod_target
export pod_settings, pod_recorder, pod_due, pod_record!, pod_record_flat!,
       pod_reset!, pod_finalize!
export pod_weights, pod_targets, pod_target, pod_save, pod_load


#---------------------------------------------------------------------------------
# What the deck asked for.
#---------------------------------------------------------------------------------
struct St_pod_settings
    lpod        ::Bool
    fields      ::Vector{Any}
    nsnap       ::Int        # sampling INTERVALS; there is one more snapshot than this
    tstart      ::Float64
    tend        ::Float64
    nmodes      ::Int
    nmodes_plot ::Int
    lmean       ::Bool
    method      ::Symbol
    nlon        ::Int        # equirectangular raster (manifold)
    nlat        ::Int
    nx          ::Int        # (x,y) raster (2-D)
    ny          ::Int
    lvtk        ::Bool
    lpng        ::Bool
    ldata       ::Bool
    lsnapshots  ::Bool
    tscale      ::Float64
    tlabel      ::String
    cmap        ::Symbol
    maxmem      ::Float64
end


"""
    pod_settings(inputs, tinit, tend) -> St_pod_settings

Read the POD block of `user_inputs.jl`. **`:lpod => true` is the only required
entry**; every other key below is an override of a default that works for any
case.

| key | default | meaning |
|:--|:--|:--|
| `:lpod`             | `false`                | run the decomposition at all |
| `:pod_fields`       | `[:all]`               | which fields — a variable name, `:all`, `:state`, or one of the shell's derived fields (see `pod_targets`) |
| `:pod_nsnapshots`   | the output cadence     | sampling intervals over the POD window (`:ndiagnostics_outputs`, else the number of `:diagnostics_at_times` intervals, else 50) |
| `:pod_tstart`       | `:tinit`               | start of the window — skip a transient with it |
| `:pod_tend`         | `:tend`                | end of the window |
| `:pod_nmodes`       | `0` (all)              | modes kept and saved |
| `:pod_nmodes_plot`  | `6`                    | modes drawn in the mode figure |
| `:pod_subtract_mean`| `true`                 | decompose the fluctuation, not the field |
| `:pod_method`       | `:auto`                | `:svd` serial, `:snapshot` under MPI |
| `:pod_nlon/:pod_nlat`| `720 / 360`           | equirectangular raster (manifold cases) |
| `:pod_nx/:pod_ny`   | `400 / 400`            | (x,y) raster (2-D cases) |
| `:pod_write_vtk/png/data` | `true`           | which outputs to produce |
| `:pod_write_snapshots` | `false`             | also dump the raw snapshots, for an offline re-run |
| `:pod_time_scale`   | `1.0`                  | multiplies t on the coefficient plots |
| `:pod_time_label`   | `"t"`                  | its axis label |
| `:pod_cmap`         | `:balance`             | diverging colour map for the modes |
| `:pod_max_memory_gb`| `4.0`                  | refuse a snapshot set larger than this |
"""
#
# How many snapshots, when the deck does not say.
#
# :ndiagnostics_outputs is the natural answer, but mod_inputs_user_inputs! sets
# it to ZERO for every case that asks for :diagnostics_at_times instead (the two
# are alternatives, and most flat cases use the second) — and "sample zero
# times" is not what such a deck means by leaving the POD block at one line. Fall
# back to the number of output times it did ask for, and then to a round number,
# so that `:lpod => true` on its own always means something sensible.
#
function _pod_default_nsnap(inputs)
    n = Int(get(inputs, :ndiagnostics_outputs, 0))
    n >= 1 && return n
    dat = get(inputs, :diagnostics_at_times, nothing)
    if dat !== nothing
        m = length(collect(dat)) - 1
        m >= 1 && return m
    end
    return 50
end


function pod_settings(inputs, tinit::Real, tend::Real)

    lpod = get(inputs, :lpod, false) == true

    rawf   = get(inputs, :pod_fields, [:all])
    fields = Any[f for f in (rawf isa AbstractVector ? rawf : [rawf])]

    nsnap  = Int(get(inputs, :pod_nsnapshots, _pod_default_nsnap(inputs)))
    tstart = Float64(get(inputs, :pod_tstart, tinit))
    tstop  = Float64(get(inputs, :pod_tend,   tend))

    if lpod && !(tstop > tstart)
        error(string(" # ERROR pod.jl: the POD window is empty — :pod_tstart = ", tstart,
                     " s is not before the end of the window, ", tstop, " s.\n",
                     " #   :pod_tstart has to fall inside the run, i.e. before :tend."))
    end
    if lpod && nsnap < 1
        error(" # ERROR pod.jl: :pod_nsnapshots must be at least 1 (it counts sampling INTERVALS, so 1 gives two snapshots).")
    end

    return St_pod_settings(lpod, fields, nsnap, tstart, tstop,
                           Int(get(inputs, :pod_nmodes, 0)),
                           Int(get(inputs, :pod_nmodes_plot, 6)),
                           get(inputs, :pod_subtract_mean, true) == true,
                           Symbol(get(inputs, :pod_method, :auto)),
                           Int(get(inputs, :pod_nlon, 720)),
                           Int(get(inputs, :pod_nlat, 360)),
                           Int(get(inputs, :pod_nx, 400)),
                           Int(get(inputs, :pod_ny, 400)),
                           get(inputs, :pod_write_vtk,  true) == true,
                           get(inputs, :pod_write_png,  true) == true,
                           get(inputs, :pod_write_data, true) == true,
                           get(inputs, :pod_write_snapshots, false) == true,
                           Float64(get(inputs, :pod_time_scale, 1.0)),
                           String(get(inputs, :pod_time_label, "t")),
                           Symbol(get(inputs, :pod_cmap, :balance)),
                           Float64(get(inputs, :pod_max_memory_gb, 4.0)))
end


#---------------------------------------------------------------------------------
# WHAT can be decomposed.
#
# A target is a name, a component count, and where to read those components
# from. Vector-valued targets are decomposed JOINTLY (see pod_from_snapshots),
# which is why :state is one target with neqs components and not neqs targets
# with one.
#---------------------------------------------------------------------------------
struct St_pod_target
    sym   ::Symbol            # :var (generic) or a derived shell field
    name  ::String
    comps ::Vector{String}
    ncomp ::Int
    icols ::Vector{Int}       # which columns of the source array (:var only)
    src   ::Symbol            # :q solution variables | :qout output variables
end

Base.show(io::IO, tg::St_pod_target) =
    print(io, "St_pod_target(", tg.name, tg.ncomp == 1 ? "" : string(" ×", tg.ncomp),
          ", ", tg.sym === :var ? string(tg.src, tg.icols) : string(":", tg.sym), ")")


"""
    pod_targets(fields, qvars, qoutvars, neqs, lshell) -> Vector{St_pod_target}

Turn `:pod_fields` into something the recorder can read. An entry may be

  * **a variable name** — `"rho"`, `:theta`, `"q"` — matched first against the
    case's solution variables (`qvars`, what the equations are written in) and
    then against its output variables (`qoutvars`, what `user_uout!` derives for
    a human). Case-insensitive.
  * **an integer** — that column of the solution vector, for a case whose
    variables are unnamed.
  * `:all` — one scalar target per SOLUTION variable. This is the default, and
    it is the one that needs the user to know nothing.
  * `:allout` — one scalar target per OUTPUT variable instead.
  * `:state` — every solution variable, stacked into ONE vector-valued target.
    This is the basis a Galerkin ROM of the system is projected onto: the modes
    are states of the whole system rather than unrelated scalars, and they share
    one set of temporal coefficients.

On a spherical shell (`lshell`) four derived fields are available besides, for
the shallow water system `q = [φ, φu, φv, φw]` with Cartesian momentum:

| symbol | components | what it is |
|:--|:--|:--|
| `:vorticity` | `ζ`        | relative vorticity `n̂·(∇ₛ×u)`, supplied by the time loop |
| `:h`         | `h`        | fluid depth `φ/g` |
| `:u`, `:v`   | `u_λ`, `u_φ` | zonal and meridional velocity, in the local tangent basis |
| `:velocity`  | `u_λ, u_φ` | the horizontal velocity as ONE vector target |

The velocity components are the tangent-basis ones and not the Cartesian ones,
because the Cartesian components of a flow on a sphere are a property of the
frame and not of the flow: their modes would show a rigid zonal jet as a dipole
straddling the prime meridian.
"""
function pod_targets(fields, qvars, qoutvars, neqs::Int, lshell::Bool)

    qn  = String[String(v) for v in qvars   if v !== nothing]
    on  = String[String(v) for v in qoutvars if v !== nothing]
    out = St_pod_target[]

    _scalar(name, icol, src) = St_pod_target(:var, name, [name], 1, [icol], src)

    for f in fields
        if f isa Integer
            1 <= f <= neqs ||
                error(string(" # ERROR pod.jl: :pod_fields asks for column ", f,
                             " of a solution vector with ", neqs, " variables."))
            push!(out, _scalar(f <= length(qn) ? qn[f] : string("q", f), Int(f), :q))
            continue
        end

        sym = f isa Symbol ? f : Symbol(String(f))

        if sym === :all
            isempty(qn) && error(" # ERROR pod.jl: :pod_fields => :all, but the case named no solution variables.")
            for (i, v) in enumerate(qn); push!(out, _scalar(v, i, :q)); end
            continue
        elseif sym === :allout
            isempty(on) && error(" # ERROR pod.jl: :pod_fields => :allout, but the case named no output variables.")
            for (i, v) in enumerate(on); push!(out, _scalar(v, i, :qout)); end
            continue
        elseif sym === :state
            isempty(qn) && error(" # ERROR pod.jl: :pod_fields => :state, but the case named no solution variables.")
            push!(out, St_pod_target(:var, "state", qn, length(qn), collect(1:length(qn)), :q))
            continue
        end

        if lshell
            sym === :vorticity && (push!(out, St_pod_target(sym, "vorticity", ["vorticity"], 1, Int[], :q)); continue)
            sym === :h         && (push!(out, St_pod_target(sym, "h",         ["h"],         1, Int[], :q)); continue)
            sym === :u         && (push!(out, St_pod_target(sym, "u",         ["u"],         1, Int[], :q)); continue)
            sym === :v         && (push!(out, St_pod_target(sym, "v",         ["v"],         1, Int[], :q)); continue)
            sym === :velocity  && (push!(out, St_pod_target(sym, "velocity",  ["u", "v"],    2, Int[], :q)); continue)
        end

        # a plain variable name, solution variables first
        name = String(sym)
        i    = findfirst(v -> lowercase(v) == lowercase(name), qn)
        if i !== nothing
            push!(out, _scalar(qn[i], i, :q)); continue
        end
        i = findfirst(v -> lowercase(v) == lowercase(name), on)
        if i !== nothing
            push!(out, _scalar(on[i], i, :qout)); continue
        end

        error(string(" # ERROR pod.jl: :pod_fields names \"", name,
                     "\", which is not a field of this case.\n",
                     " #   solution variables: ", isempty(qn) ? "(none named)" : join(qn, ", "), "\n",
                     " #   output variables  : ", isempty(on) ? "(none named)" : join(on, ", "), "\n",
                     " #   always available  : :all, :allout, :state",
                     lshell ? ", and on the shell :vorticity, :h, :u, :v, :velocity" : ""))
    end

    isempty(out) && error(" # ERROR pod.jl: :lpod => true but :pod_fields selected nothing.")

    # Two targets writing to files of the same name would overwrite each other.
    seen = Dict{String,Int}()
    for (k, tg) in enumerate(out)
        n = get(seen, tg.name, 0)
        seen[tg.name] = n + 1
        n == 0 && continue
        out[k] = St_pod_target(tg.sym, string(tg.name, "_", n+1), tg.comps, tg.ncomp, tg.icols, tg.src)
    end

    return out
end

"""
    pod_target(sym) -> St_pod_target

One shell field by symbol, for a caller that is not going through a deck.
"""
pod_target(sym::Symbol) = pod_targets([sym], String[], String[], 0, true)[1]

pod_needs_vorticity(tg::St_pod_target) = tg.sym === :vorticity
pod_needs_qout(tg::St_pod_target)      = tg.src === :qout


#
# Read one snapshot of one target.
#
#   uaux  npoin × neqs, the solution variables
#   qout  npoin × noutvar, the case's user_uout! variables (or `nothing`)
#   ζ     the relative vorticity, supplied by the shell time loop
#
function pod_extract!(dest::AbstractMatrix, tg::St_pod_target,
                      uaux::AbstractMatrix, qout, ζ, mesh, g::Float64)

    npoin = size(dest, 1)

    if tg.sym === :var
        src = tg.src === :qout ? qout : uaux
        src === nothing &&
            error(string(" # ERROR pod.jl: the target \"", tg.name,
                         "\" reads the output variables, but none were computed for this snapshot."))
        @inbounds for c = 1:tg.ncomp
            ic = tg.icols[c]
            for ip = 1:npoin
                dest[ip,c] = src[ip,ic]
            end
        end

    elseif tg.sym === :vorticity
        ζ === nothing && error(" # ERROR pod.jl: the vorticity target was asked for a snapshot without a vorticity field.")
        @inbounds for ip = 1:npoin
            dest[ip,1] = ζ[ip]
        end

    elseif tg.sym === :h
        @inbounds for ip = 1:npoin
            dest[ip,1] = uaux[ip,1]/g
        end

    else   # :u, :v, :velocity — the tangent-basis components
        @inbounds for ip = 1:npoin
            φ  = uaux[ip,1]
            ux, uy, uz = uaux[ip,2]/φ, uaux[ip,3]/φ, uaux[ip,4]/φ
            sλ, cλ = sincos(mesh.lon[ip])
            sφ, cφ = sincos(mesh.lat[ip])
            uλ = -sλ*ux + cλ*uy                        # u·e_λ  eastward
            uφ = -sφ*cλ*ux - sφ*sλ*uy + cφ*uz          # u·e_φ  northward
            if tg.ncomp == 2
                dest[ip,1] = uλ
                dest[ip,2] = uφ
            else
                dest[ip,1] = tg.sym === :v ? uφ : uλ
            end
        end
    end
    return dest
end


#---------------------------------------------------------------------------------
# The recorder: the snapshot buffers and the sampling clock.
#---------------------------------------------------------------------------------
mutable struct St_pod_recorder{TF}
    set      ::St_pod_settings
    targets  ::Vector{St_pod_target}
    X        ::Vector{Array{TF,3}}     # one npoin × ncomp × nsnapmax per target
    t        ::Vector{TF}
    nsnap    ::Int
    nsnapmax ::Int
    dt       ::TF
    tnext    ::TF
    npoin    ::Int
    neqs     ::Int
    nsd      ::Int                     # 1, 2 or 3 space dimensions
    lshell   ::Bool                    # …on a manifold?
    g        ::TF
    lvort    ::Bool
    lqout    ::Bool
    uaux     ::Array{TF,2}             # scratch: the flat state, unpacked
    qout     ::Array{TF,2}             # scratch: the output variables
end


"""
    pod_recorder(inputs, mesh, tinit, tend; kwargs...) -> St_pod_recorder | nothing

Allocate the snapshot buffers, or return `nothing` when `:lpod` is off — which
is what every call site tests, so that a case that does not ask for POD carries
no cost and no branch beyond a `=== nothing`.

  * `qvars`, `qoutvars`, `neqs` — the case's own variable names, so that
    `:pod_fields` can be written in them (and so that the default, every
    solution variable, needs no deck entry at all).
  * `nsd`, `lshell` — the geometry, which decides only how the modes are
    written out.
  * `dt_step` — the time step, used to warn when the deck asks to be sampled
    faster than the case is integrated.
"""
function pod_recorder(inputs, mesh, tinit::Real, tend::Real;
                      verbose::Bool = true, dt_step::Real = 0.0,
                      qvars = String[], qoutvars = String[], neqs::Int = 0,
                      nsd::Int = 2, lshell::Bool = false)

    set = pod_settings(inputs, tinit, tend)
    set.lpod || return nothing

    #
    # ADAPTIVITY AND POD DO NOT COMPOSE, and the failure is silent rather than
    # loud: a refined mesh has a different npoin, so snapshot k and snapshot
    # k+1 are vectors in different spaces. There is no correlation matrix to
    # form. (A ROM on an adapting grid needs the snapshots interpolated onto a
    # common reference mesh first — a different piece of machinery.)
    #
    get(inputs, :lamr, false) == true &&
        error(" # ERROR pod.jl: :lpod and :lamr => true cannot both be on. Adaptive refinement changes npoin between snapshots, so the snapshots do not live in the same space and there is no decomposition to compute.")

    targets  = pod_targets(set.fields, qvars, qoutvars, neqs, lshell)
    npoin    = Int(mesh.npoin)
    nsnapmax = set.nsnap + 1                       # the window's ends are both sampled
    ncomptot = sum(tg.ncomp for tg in targets)

    bytes = 8.0*npoin*ncomptot*nsnapmax
    bytes <= set.maxmem*2^30 ||
        error(string(" # ERROR pod.jl: the snapshot set would need ",
                     @sprintf("%.2f", bytes/2^30), " GB (", npoin, " nodes × ", ncomptot,
                     " components × ", nsnapmax, " snapshots × 8 B), above :pod_max_memory_gb = ",
                     set.maxmem, ".\n",
                     " #   Lower :pod_nsnapshots, ask for fewer :pod_fields, or raise the limit."))

    X  = [zeros(Float64, npoin, tg.ncomp, nsnapmax) for tg in targets]
    dt = (set.tend - set.tstart)/set.nsnap
    g  = Float64(get(inputs, :galewsky_gravity, get(inputs, :sphere_gravity, 9.80616)))

    lqout   = any(pod_needs_qout, targets)
    noutvar = max(length(qoutvars), 1)

    #
    # Nothing can be sampled faster than it is computed. A deck that asks for
    # more snapshots than there are time steps in the window gets the steps, and
    # is told so rather than being left to wonder why :pod_nsnapshots was not
    # honoured (pod_record! does the catching up).
    #
    if dt_step > 0 && dt < dt_step
        @warn string("POD: :pod_nsnapshots => ", set.nsnap, " asks for a snapshot every ",
                     @sprintf("%.4g", dt), " s, which is shorter than the time step (",
                     @sprintf("%.4g", dt_step), " s). At most ",
                     floor(Int, (set.tend - set.tstart)/dt_step) + 1,
                     " snapshots can be taken, one per step.")
    end

    if verbose
        println(" # ")
        println(" # POD — Proper Orthogonal Decomposition ........................")
        @printf(" #   fields          : %s\n", join((tg.ncomp == 1 ? tg.name :
                                                     string(tg.name, " (", join(tg.comps, ","), ")")
                                                     for tg in targets), ", "))
        @printf(" #   snapshots       : %d, every %.6g s, t = %.6g … %.6g s\n",
                nsnapmax, dt, set.tstart, set.tend)
        @printf(" #   geometry        : %s\n",
                lshell ? "spherical shell (2-D manifold)" : string(nsd, "-D"))
        @printf(" #   inner product   : the SEM mass matrix (discrete L²)\n")
        @printf(" #   mean            : %s\n", set.lmean ? "subtracted (POD of the fluctuation)" : "kept (POD of the field)")
        @printf(" #   memory          : %.1f MB held for the run\n", bytes/2^20)
    end

    return St_pod_recorder{Float64}(set, targets, X, zeros(Float64, nsnapmax), 0, nsnapmax,
                                    dt, Float64(set.tstart), npoin, neqs, nsd, lshell, g,
                                    any(pod_needs_vorticity, targets), lqout,
                                    zeros(Float64, npoin, max(neqs, 1)),
                                    zeros(Float64, lqout ? npoin : 0, lqout ? noutvar : 0))
end


"""
    pod_due(rec, t, dt) -> Bool

Is `t` a snapshot time? The `1e-9·dt` slack is the one the output callbacks use:
the sampling times are exact multiples of the POD interval, the integrator lands
on multiples of Δt, and without the slack a snapshot is missed whenever the two
differ in the last bit.
"""
pod_due(rec::Nothing, t, dt) = false
pod_due(rec::St_pod_recorder, t, dt) =
    rec.nsnap < rec.nsnapmax && t >= rec.tnext - 1.0e-9*abs(dt)


"""
    pod_reset!(rec)

Throw away everything recorded and rewind the clock to the start of the window.
Called after the integrator warm-up in `time_loop!`: the warm-up runs a
throw-away step with the REAL callback set, so without this the first snapshot
would be the warm-up's and the sampling clock would start one interval in.
"""
pod_reset!(rec::Nothing) = nothing
function pod_reset!(rec::St_pod_recorder)
    rec.nsnap = 0
    rec.tnext = rec.set.tstart
    return rec
end


"""
    pod_record!(rec, uaux, t, mesh; ζ = nothing, qout = nothing)

Take one snapshot of every target from the solution variables `uaux`
(`npoin × neqs`).
"""
pod_record!(rec::Nothing, args...; kwargs...) = nothing

function pod_record!(rec::St_pod_recorder, uaux::AbstractMatrix, t::Real, mesh;
                     ζ = nothing, qout = nothing)

    rec.nsnap < rec.nsnapmax || return nothing

    k = rec.nsnap + 1
    for (m, tg) in enumerate(rec.targets)
        pod_extract!(view(rec.X[m], :, :, k), tg, uaux, qout, ζ, mesh, rec.g)
    end
    rec.t[k] = Float64(t)
    rec.nsnap = k

    #
    # The next sampling time, STRICTLY after this one. The loop matters when the
    # sampling interval is finer than the time step, or when a step overshot the
    # interval: without it the recorder would fall a snapshot further behind at
    # every step and never catch up, silently stretching the POD window.
    #
    rec.tnext += rec.dt
    while rec.tnext <= t
        rec.tnext += rec.dt
    end
    return nothing
end


"""
    pod_record_flat!(rec, u, t, params)

The same, from the integrator's own state vector — which in the generic time
loop is FLAT (`npoin·neqs`, the layout SciML marches) rather than
`npoin × neqs`. Unpacks it into the recorder's own scratch, exactly as
`write_output` does, and derives the case's output variables when a target asks
for one. Its own buffers, and not `params.uaux`, because that one holds whatever
the last RHS evaluation left in it — which is a stage of the step, not the step.
"""
function pod_record_flat!(rec::St_pod_recorder, u, t::Real, params)

    rec.nsnap < rec.nsnapmax || return nothing

    npoin = rec.npoin
    nvar  = rec.neqs

    q = u
    if params.inputs[:backend] != CPU()
        q = KernelAbstractions.allocate(CPU(), TFloat, Int64(npoin*nvar))
        KernelAbstractions.copyto!(CPU(), q, u)
    end
    u2uaux!(rec.uaux, q, nvar, npoin)

    qout = nothing
    if rec.lqout
        call_user_uout(rec.qout, rec.uaux, params.qp.qe, params.mp,
                       params.inputs[:SOL_VARS_TYPE], npoin, nvar, size(rec.qout, 2))
        qout = rec.qout
    end

    return pod_record!(rec, rec.uaux, t, params.mesh; qout = qout)
end

pod_record_flat!(rec::Nothing, u, t, params) = nothing


#---------------------------------------------------------------------------------
# The inner product, and the two parallel hooks.
#---------------------------------------------------------------------------------
"""
    pod_weights(M, mesh) -> w

The diagonal mass matrix, with the nodes this rank does not OWN zeroed out. A
node on a partition seam lives on every rank that touches it and carries the
same value on each; counting it once per rank would weight the seams by the
partition, and the decomposition would then depend on the rank count. This is
the same ownership test `_sphere_integrals` applies to the conserved integrals.

A case built with `:lexact_integration` and no lumping has a FULL mass matrix;
its diagonal is taken, which is the lumped approximation of (4). Everything the
decomposition claims still holds, in that inner product.
"""
function pod_weights(M, mesh)
    npoin = Int(mesh.npoin)
    w     = Vector{Float64}(undef, npoin)

    if M isa AbstractMatrix
        @inbounds for ip = 1:npoin
            w[ip] = Float64(M[ip,ip])
        end
    else
        length(M) >= npoin ||
            error(string(" # ERROR pod.jl: the mass matrix has ", length(M),
                         " entries for ", npoin, " nodes."))
        @inbounds for ip = 1:npoin
            w[ip] = Float64(M[ip])
        end
    end

    comm = get_mpi_comm()
    if MPI.Comm_size(comm) > 1
        rank  = MPI.Comm_rank(comm)
        owner = mesh.gip2owner
        length(owner) >= npoin ||
            error(" # ERROR pod.jl: mesh.gip2owner is missing, so owned nodes cannot be told from mirrored ones.")
        @inbounds for ip = 1:npoin
            owner[ip] == rank || (w[ip] = 0.0)
        end
    end
    return w
end

_pod_reducer(comm) =
    MPI.Comm_size(comm) > 1 ? (A -> (MPI.Allreduce!(A, MPI.SUM, comm); A)) : _pod_noreduce!

function _pod_signreduce(comm)
    MPI.Comm_size(comm) > 1 || return _pod_localsign
    return function (maxabs::Float64, val::Float64)
        gmax = MPI.Allreduce(maxabs, MPI.MAX, comm)
        # Ranks that do not hold the global extremum abstain. A node mirrored on
        # several ranks makes several of them agree, which only changes the
        # magnitude of the sum, never its sign.
        c = (maxabs == gmax && maxabs > 0) ? (val < 0 ? -1.0 : 1.0) : 0.0
        s = MPI.Allreduce(c, MPI.SUM, comm)
        return s < 0 ? -1.0 : 1.0
    end
end


#---------------------------------------------------------------------------------
# The decomposition, and everything it leaves on disk.
#---------------------------------------------------------------------------------
"""
    pod_finalize!(rec, mesh, M, OUTPUT_DIR; verbose = true) -> Vector{St_pod}

Decompose every recorded field and write the results. Called once, at the end of
the run; returns the decompositions so that a caller (a ROM driver, a test) can
go on using them in memory. `M` is the case's mass matrix — `metrics.M` on the
shell, `params.M` elsewhere.
"""
pod_finalize!(rec::Nothing, mesh, M, OUTPUT_DIR; verbose = true) = St_pod{Float64}[]

function pod_finalize!(rec::St_pod_recorder, mesh, M, OUTPUT_DIR::String;
                       verbose::Bool = true)

    comm    = get_mpi_comm()
    nparts  = MPI.Comm_size(comm)
    lserial = nparts == 1
    rank0   = MPI.Comm_rank(comm) == 0

    if rec.nsnap < 2
        verbose && @warn string("POD: only ", rec.nsnap, " snapshot(s) were taken, so there is nothing ",
                                "to decompose. The run probably ended before :pod_tstart = ",
                                rec.set.tstart, " s.")
        return St_pod{Float64}[]
    end

    w    = pod_weights(M, mesh)
    red! = _pod_reducer(comm)
    sgn  = _pod_signreduce(comm)

    if verbose
        println(" # ")
        println(" # POD ..........................................................")
        @printf(" #   %d snapshots, t = %.6g … %.6g s\n",
                rec.nsnap, rec.t[1], rec.t[rec.nsnap])
    end

    out = St_pod{Float64}[]

    for (m, tg) in enumerate(rec.targets)

        # Only the snapshots actually taken: the buffer is sized for the whole
        # window, and a run stopped early leaves zeros behind it, which are not
        # data and would be decomposed as if they were.
        X = @view rec.X[m][:, :, 1:rec.nsnap]

        #
        # A field that never moved has no fluctuation, hence no correlation
        # matrix with a positive eigenvalue and no modes — POD of it is not a
        # small answer, it is undefined. That happens for real: a case whose
        # forcing is switched off integrates a state at rest, and a deck that
        # asks for the POD of every variable then asks for the decomposition of
        # a constant. Say so and move on to the next field rather than killing a
        # run that was otherwise fine.
        #
        if _pod_is_constant(X, comm)
            verbose && @warn string("POD: ", tg.name, " does not vary over the snapshots ",
                                    "(every one equals the mean to round-off), so there is ",
                                    "nothing to decompose. Skipping it.")
            continue
        end

        P = pod_from_snapshots(X, w, @view rec.t[1:rec.nsnap];
                               name       = tg.name,
                               comps      = tg.comps,
                               nmodes     = rec.set.nmodes,
                               lmean      = rec.set.lmean,
                               method     = rec.set.method,
                               lparallel  = !lserial,
                               reducer!   = red!,
                               signreduce = sgn)
        push!(out, P)

        verbose && _pod_report(P)

        rec.set.lvtk  && pod_write_modes(P, rec, mesh, OUTPUT_DIR; verbose = verbose)
        rec.set.ldata && rank0 && pod_write_data(P, OUTPUT_DIR; verbose = verbose)
        if rec.set.lsnapshots && rank0
            pod_write_snapshots(P, X, w, OUTPUT_DIR; verbose = verbose)
        end

        if rec.set.lpng
            if lserial
                pod_plot(P, mesh, rec, OUTPUT_DIR; verbose = verbose)
            elseif verbose
                println(" #   PNG: skipped — the raster needs the whole domain on one rank.")
                println(" #        The .vtu carries the same modes; re-run on one rank for the maps.")
            end
        end
    end

    if verbose
        println(" # POD ..................................................... DONE")
        println(" #   Output written to: ", abspath(OUTPUT_DIR))
    end

    return out
end


#
# Does this field move at all? Measured as the largest deviation from the FIRST
# snapshot, relative to the largest value the field takes — a scale-free test,
# so that it means the same for a geopotential of 1e5 and a vorticity of 1e-5.
#
function _pod_is_constant(X::AbstractArray{<:Real,3}, comm; rtol = 1.0e-14)
    npoin, ncomp, nsnap = size(X)
    dev = 0.0
    mag = 0.0
    @inbounds for k = 1:nsnap, c = 1:ncomp, ip = 1:npoin
        x   = Float64(X[ip,c,k])
        mag = max(mag, abs(x))
        dev = max(dev, abs(x - Float64(X[ip,c,1])))
    end
    if MPI.Comm_size(comm) > 1
        dev = MPI.Allreduce(dev, MPI.MAX, comm)
        mag = MPI.Allreduce(mag, MPI.MAX, comm)
    end
    return dev <= rtol*max(mag, 1.0e-300)
end


#
# The spectrum, as a table. This is the thing to read first: how many modes the
# flow actually has, and therefore how big a ROM has to be.
#
function _pod_report(P::St_pod)
    r = length(P.λ)
    println(" # ")
    @printf(" #   %s: %d modes from %d snapshots (%s), Σλ = %.6e\n",
            P.name, r, P.nsnap, P.method, P.total)
    @printf(" #     orthonormality  max|ΦᵀMΦ - I| = %.2e\n", P.orthoerr)
    println(" #     mode        λ            E [%]      ΣE [%]")
    for i = 1:min(r, 10)
        @printf(" #     %4d   %12.5e   %8.3f   %8.3f\n",
                i, P.λ[i], 100*P.energy[i], 100*P.cumenergy[i])
    end
    r > 10 && @printf(" #     …  (%d more)\n", r - 10)
    for frac in (0.90, 0.99)
        n = pod_rank_for_energy(P, frac)
        @printf(" #     %2.0f%% of the energy is in the first %d mode%s\n",
                100*frac, n, n == 1 ? "" : "s")
    end
    return nothing
end


#---------------------------------------------------------------------------------
# THE MODES, WRITTEN OUT — the one place the geometry matters.
#---------------------------------------------------------------------------------
"""
    pod_write_modes(P, rec, mesh, OUTPUT_DIR; verbose = true)

The mean and the modes as a field on the grid they were computed on:

  * 1-D — `pod_<field>_modes.csv`, a table of `x, mean, mode_001, …` sorted by
    `x`. A .vtu of a line grid carries the same numbers in a format nothing
    reads as easily.
  * 2-D — `pod_<field>.vtu`, the quads of the mesh, one point-data array per
    mode.
  * 3-D — `pod_<field>.vtu`, the hexahedra of the mesh, likewise.
  * shell — `pod_<field>.vtu` through `write_vtk_sphere_grid`, so the modes land
    on the same spherical geometry the solution does and can be overlaid on it.
"""
function pod_write_modes(P::St_pod, rec::St_pod_recorder, mesh, OUTPUT_DIR::String;
                         verbose::Bool = true)

    isdir(OUTPUT_DIR) || mkpath(OUTPUT_DIR)

    if rec.nsd == 1
        return _pod_write_modes_1d(P, mesh, OUTPUT_DIR; verbose = verbose)
    end

    fields = Pair{String,Vector{Float64}}[]
    for c = 1:P.ncomp
        cn = P.ncomp == 1 ? "" : string("_", P.comps[c])
        P.lmean && push!(fields, string("pod_mean", cn) => collect(view(P.q̄, :, c)))
        for i = 1:length(P.λ)
            push!(fields, @sprintf("pod_mode_%03d%s", i, cn) => collect(view(P.Φ, :, c, i)))
        end
    end

    fname = string("pod_", P.name)
    if rec.lshell
        write_vtk_sphere_grid(mesh, fname, OUTPUT_DIR; extra = fields, verbose = false)
    else
        _pod_write_vtk(mesh, rec.nsd, fname, OUTPUT_DIR, fields)
    end

    verbose && @printf(" #     %s%s   mean + %d modes\n",
                       joinpath(abspath(OUTPUT_DIR), fname),
                       MPI.Comm_size(get_mpi_comm()) > 1 ? ".pvtu" : ".vtu", length(P.λ))
    return nothing
end


#
# 1-D: a table, not a mesh file. The nodes are written in x order (the SEM
# numbering is elementwise and would draw a zig-zag), and duplicated element
# interfaces are kept — they carry the same value, and dropping them would
# misrepresent the discretization.
#
function _pod_write_modes_1d(P::St_pod, mesh, OUTPUT_DIR::String; verbose::Bool = true)

    npoin = P.npoin
    x     = [Float64(mesh.x[ip]) for ip = 1:npoin]
    perm  = sortperm(x)

    fout = joinpath(OUTPUT_DIR, string("pod_", P.name, "_modes.csv"))
    open(fout, "w") do io
        println(io, "# POD modes of ", P.name, " — ", length(P.λ), " modes from ", P.nsnap, " snapshots")
        println(io, "# the modes are orthonormal under the mass matrix: sum_ip M_ip phi_i phi_j = delta_ij")
        print(io, "x")
        P.lmean && print(io, ",mean")
        for i = 1:length(P.λ); @printf(io, ",mode_%03d", i); end
        println(io)
        for ip in perm
            @printf(io, "%.16e", x[ip])
            P.lmean && @printf(io, ",%.16e", P.q̄[ip,1])
            for i = 1:length(P.λ); @printf(io, ",%.16e", P.Φ[ip,1,i]); end
            println(io)
        end
    end
    verbose && @printf(" #     %s   mean + %d modes\n", abspath(fout), length(P.λ))
    return nothing
end


#
# 2-D and 3-D: the unstructured grid of the sub-elements, with one point-data
# array per mode. Built here rather than through write_output's writers because
# those take a solution vector and its variable names; this takes an arbitrary
# list of named nodal fields, and a POD basis is exactly that.
#
function _pod_write_vtk(mesh, nsd::Int, fname::String, OUTPUT_DIR::String,
                        fields::Vector{Pair{String,Vector{Float64}}})

    ngl   = Int(mesh.ngl)
    nelem = Int(mesh.nelem)
    npoin = Int(mesh.npoin)

    cells = if nsd == 2
        cs = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef, nelem*(ngl-1)^2)
        isel = 1
        for iel = 1:nelem, i = 1:ngl-1, j = 1:ngl-1
            cs[isel] = MeshCell(VTKCellTypes.VTK_QUAD,
                                Int64[mesh.connijk[iel,i,  j  ], mesh.connijk[iel,i+1,j  ],
                                      mesh.connijk[iel,i+1,j+1], mesh.connijk[iel,i,  j+1]])
            isel += 1
        end
        cs
    else
        cs = Vector{MeshCell{VTKCellType, Vector{Int64}}}(undef, nelem*(ngl-1)^3)
        isel = 1
        for iel = 1:nelem, i = 1:ngl-1, j = 1:ngl-1, k = 1:ngl-1
            cs[isel] = MeshCell(VTKCellTypes.VTK_HEXAHEDRON,
                                Int64[mesh.connijk[iel,i,  j,  k  ], mesh.connijk[iel,i+1,j,  k  ],
                                      mesh.connijk[iel,i+1,j+1,k  ], mesh.connijk[iel,i,  j+1,k  ],
                                      mesh.connijk[iel,i,  j,  k+1], mesh.connijk[iel,i+1,j,  k+1],
                                      mesh.connijk[iel,i+1,j+1,k+1], mesh.connijk[iel,i,  j+1,k+1]])
            isel += 1
        end
        cs
    end

    xs = Float64[mesh.x[ip] for ip = 1:npoin]
    ys = Float64[mesh.y[ip] for ip = 1:npoin]
    zs = nsd == 3 ? Float64[mesh.z[ip] for ip = 1:npoin] : zeros(Float64, npoin)

    fout_name = string(OUTPUT_DIR, "/", fname)
    comm      = get_mpi_comm()
    nparts    = MPI.Comm_size(comm)
    part      = MPI.Comm_rank(comm) + 1        # WriteVTK parts are 1-indexed

    vtkf = nparts > 1 ?
        pvtk_grid(fout_name, xs, ys, zs, cells, compress = false;
                  part = part, nparts = nparts, ismain = (part == 1)) :
        vtk_grid(fout_name, xs, ys, zs, cells, compress = false)

    for (name, vals) in fields
        vtkf[name, VTKPointData()] = vals
    end
    vtk_save(vtkf)
    return nothing
end


"""
    pod_write_data(P, OUTPUT_DIR; verbose = true)

The numbers behind the pictures: the spectrum and the truncation error as CSV,
the temporal coefficients as CSV, and the basis itself as JLD2.

The JLD2 file is the hand-off to a reduced-order model. It is written as plain
arrays rather than as a serialised `St_pod` on purpose — the same reason
sem_setup.jl gives for its metric cache: a stored struct stops loading the day
its definition changes, and a basis is worth more than the struct that held it.
Read it back with [`pod_load`](@ref).
"""
function pod_write_data(P::St_pod, OUTPUT_DIR::String; verbose::Bool = true)

    isdir(OUTPUT_DIR) || mkpath(OUTPUT_DIR)
    r = length(P.λ)
    e = pod_truncation_error(P)

    fspec = joinpath(OUTPUT_DIR, string("pod_", P.name, "_spectrum.csv"))
    open(fspec, "w") do io
        println(io, "# POD spectrum of ", P.name, " — ", P.nsnap, " snapshots, method = ", P.method)
        println(io, "# total energy <|q'|^2> = ", P.total, " ; max|Phi^T M Phi - I| = ", P.orthoerr)
        println(io, "# truncation_error[i] is the relative L2 error of the reconstruction with i modes")
        println(io, "mode,lambda,energy_fraction,cumulative_energy,truncation_error")
        for i = 1:r
            @printf(io, "%d,%.16e,%.16e,%.16e,%.16e\n",
                    i, P.λ[i], P.energy[i], P.cumenergy[i], e[i+1])
        end
    end

    fcoef = joinpath(OUTPUT_DIR, string("pod_", P.name, "_coefficients.csv"))
    open(fcoef, "w") do io
        println(io, "# POD temporal coefficients a_i(t) of ", P.name)
        print(io, "t")
        for i = 1:r; print(io, ",a", i); end
        println(io)
        for k = 1:P.nsnap
            @printf(io, "%.16e", P.t[k])
            for i = 1:r; @printf(io, ",%.16e", P.a[k,i]); end
            println(io)
        end
    end

    pod_save(P, joinpath(OUTPUT_DIR, string("pod_", P.name, ".jld2")))

    verbose && @printf(" #     %s{_spectrum.csv,_coefficients.csv,.jld2}\n",
                       joinpath(abspath(OUTPUT_DIR), string("pod_", P.name)))
    return nothing
end


"""
    pod_write_snapshots(P, X, w, OUTPUT_DIR; verbose = true)

The raw snapshot set and the quadrature weights, under `:pod_write_snapshots`.
This is what lets a decomposition be REDONE offline — over a shorter window,
about a different mean, with a different rank — without re-running the case,
which is the expensive half of building a reduced-order model:

    d = JLD2.load("output/pod_<field>_snapshots.jld2")
    P = pod_from_snapshots(d["snapshots"], d["weights"], d["t"]; nmodes = 8)
"""
function pod_write_snapshots(P::St_pod, X, w, OUTPUT_DIR::String; verbose::Bool = true)
    fout = joinpath(OUTPUT_DIR, string("pod_", P.name, "_snapshots.jld2"))
    JLD2.jldsave(fout;
                 name = P.name, comps = P.comps,
                 snapshots = Array{Float64}(X), weights = Vector{Float64}(w),
                 t = P.t, format = "jexpresso-pod-snapshots-1")
    verbose && @printf(" #     %s   %d snapshots, raw\n", abspath(fout), P.nsnap)
    return nothing
end


"""
    pod_save(P, path) ; pod_load(path) -> St_pod

Round-trip a basis through JLD2. `pod_load` is the entry point of a reduced-order
model built on a decomposition computed by an earlier run: it returns the same
`St_pod` the solver had, so `pod_project`/`pod_reconstruct` work unchanged.
"""
function pod_save(P::St_pod, path::String)
    dir = dirname(path)
    isempty(dir) || isdir(dir) || mkpath(dir)
    JLD2.jldsave(path;
                 name = P.name, comps = P.comps,
                 npoin = P.npoin, ncomp = P.ncomp, nsnap = P.nsnap,
                 t = P.t, mean = P.q̄, modes = P.Φ, lambda = P.λ, coefficients = P.a,
                 energy = P.energy, cumenergy = P.cumenergy, total = P.total,
                 lmean = P.lmean, method = String(P.method), orthoerr = P.orthoerr,
                 format = "jexpresso-pod-1")
    return path
end

function pod_load(path::String)
    isfile(path) || error(string(" # ERROR pod.jl: no POD basis at ", path, "."))
    d = JLD2.load(path)
    get(d, "format", "") == "jexpresso-pod-1" ||
        error(string(" # ERROR pod.jl: ", path, " is not a Jexpresso POD basis (format = ",
                     get(d, "format", "absent"), ")."))
    return St_pod{Float64}(d["name"], d["comps"], d["npoin"], d["ncomp"], d["nsnap"],
                           d["t"], d["mean"], d["modes"], d["lambda"], d["coefficients"],
                           d["energy"], d["cumenergy"], d["total"], d["lmean"],
                           Symbol(d["method"]), d["orthoerr"])
end
