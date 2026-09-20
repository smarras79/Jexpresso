#---------------------------------------------------------------------------------
# pod.jl — POD of a running Jexpresso simulation: what to sample, when to sample
# it, how to reduce it across ranks, and what to leave on disk.
#
# The decomposition itself is in pod_core.jl, which knows nothing about meshes,
# cases or MPI. This file is the part that does: it collects snapshots from the
# integrator as it runs, hands them to pod_core.jl under the SEM mass matrix as
# inner product, and writes the result out in the four forms it is wanted in —
#
#   pod_<field>.vtu                  the modes ON THE SPHERE, for ParaView
#   pod_<field>_modes.png            the modes on an equirectangular map
#   pod_<field>_spectrum.png         the energy spectrum and its cumulative sum
#   pod_<field>_coefficients.png     a_i(t), and the (a_1,a_2) phase portrait
#   pod_<field>_spectrum.csv         λ_i, E_i, ΣE_i, and the truncation error
#   pod_<field>_coefficients.csv     a_i(t_k)
#   pod_<field>.jld2                 the basis itself — what a ROM reads back
#
# WHY IN THE RUN AND NOT AFTERWARDS FROM THE VTK FILES. Two reasons, and both
# are about getting the right answer rather than about convenience. First, the
# inner product: POD is only optimal under the L² product (Eq. (4) of
# pod_core.jl), which needs the mass matrix — that is, the surface Jacobian and
# the LGL weights of the run that produced the data. It is on `metrics` here and
# in no output file. Second, the sampling: POD wants snapshots uniformly spaced
# in time and dense enough to resolve the dynamics, while the VTK cadence is
# chosen to keep a movie small. :pod_nsnapshots therefore carries its own timer,
# independent of :ndiagnostics_outputs, and the snapshots never touch the disk.
#
# WHAT IT COSTS. One array of npoin × ncomp × nsnap Float64 per requested field,
# held for the run: 49 snapshots of the vorticity on the shipped cubed sphere
# (npoin = 15 002) is 5.9 MB, and the default 101 snapshots of vorticity and h
# together is 24 MB. :pod_max_memory_gb refuses, with the arithmetic in the
# message, anything that would not fit. The decomposition itself is O(N K²) and
# runs once, at the end: sub-second at these sizes.
#
# UNDER MPI the field stays distributed and nothing is gathered. The correlation
# matrix is a sum over nodes, so each rank forms its own K×K partial and one
# Allreduce completes it; every rank then solves the same small eigenproblem and
# builds its own slice of the modes. Two details make the result independent of
# the rank count: a node shared by several ranks is counted by its OWNER only
# (mesh.gip2owner, exactly as sphere_diagnostics does for the conserved
# integrals), or the seam nodes would be weighted by how many ranks happen to
# touch them; and the sign convention of pod_core.jl is resolved globally.
#
# S. Marras & contributors
#---------------------------------------------------------------------------------

export St_pod_settings, St_pod_recorder
export pod_settings, pod_recorder, pod_due, pod_record!, pod_finalize!
export pod_weights, pod_save, pod_load


#---------------------------------------------------------------------------------
# What the deck asked for.
#---------------------------------------------------------------------------------
struct St_pod_settings
    lpod        ::Bool
    fields      ::Vector{Symbol}
    nsnap       ::Int        # sampling INTERVALS; there is one more snapshot than this
    tstart      ::Float64
    tend        ::Float64
    nmodes      ::Int
    nmodes_plot ::Int
    lmean       ::Bool
    method      ::Symbol
    nlon        ::Int
    nlat        ::Int
    lvtk        ::Bool
    lpng        ::Bool
    ldata       ::Bool
    tscale      ::Float64
    tlabel      ::String
    cmap        ::Symbol
    maxmem      ::Float64
end


"""
    pod_settings(inputs, tinit, tend) -> St_pod_settings

Read the POD block of `user_inputs.jl`. Every key is optional; the defaults are
the ones the SWsphere deck documents.

| key | default | meaning |
|:--|:--|:--|
| `:lpod`             | `false`                | run the decomposition at all |
| `:pod_fields`       | `[:vorticity, :h]`     | which fields (see `pod_target`) |
| `:pod_nsnapshots`   | `:ndiagnostics_outputs`| sampling intervals over the POD window |
| `:pod_tstart`       | `:tinit`               | start of the window — skip a transient with it |
| `:pod_nmodes`       | `0` (all)              | modes kept and saved |
| `:pod_nmodes_plot`  | `6`                    | modes drawn in the mode figure |
| `:pod_subtract_mean`| `true`                 | decompose the fluctuation, not the field |
| `:pod_method`       | `:auto`                | `:svd` serial, `:snapshot` under MPI |
| `:pod_nlon/:pod_nlat`| `720 / 360`           | equirectangular raster size |
| `:pod_write_vtk/png/data` | `true`           | which outputs to produce |
| `:pod_time_scale`   | `1/86400`              | multiplies t on the coefficient plots |
| `:pod_time_label`   | `"t [days]"`           | its axis label |
| `:pod_cmap`         | `:balance`             | diverging colour map for the modes |
| `:pod_max_memory_gb`| `4.0`                  | refuse a snapshot set larger than this |
"""
function pod_settings(inputs, tinit::Real, tend::Real)

    lpod = get(inputs, :lpod, false) == true

    rawf   = get(inputs, :pod_fields, [:vorticity, :h])
    fields = Symbol[Symbol(f) for f in (rawf isa AbstractVector ? rawf : [rawf])]

    nsnap  = Int(get(inputs, :pod_nsnapshots, get(inputs, :ndiagnostics_outputs, 48)))
    tstart = Float64(get(inputs, :pod_tstart, tinit))
    tstop  = Float64(get(inputs, :pod_tend,   tend))

    if lpod && !(tstop > tstart)
        error(string(" # ERROR pod.jl: the POD window is empty — :pod_tstart = ", tstart,
                     " s is not before the end of the window, ", tstop, " s.\n",
                     " #   :pod_tstart has to fall inside the run, i.e. before :tend."))
    end
    if lpod && Int(get(inputs, :pod_nsnapshots, get(inputs, :ndiagnostics_outputs, 48))) < 1
        error(" # ERROR pod.jl: :pod_nsnapshots must be at least 1 (it counts sampling INTERVALS, so 1 gives two snapshots).")
    end

    return St_pod_settings(lpod, fields, nsnap, tstart, tstop,
                           Int(get(inputs, :pod_nmodes, 0)),
                           Int(get(inputs, :pod_nmodes_plot, 6)),
                           get(inputs, :pod_subtract_mean, true) == true,
                           Symbol(get(inputs, :pod_method, :auto)),
                           Int(get(inputs, :pod_nlon, 720)),
                           Int(get(inputs, :pod_nlat, 360)),
                           get(inputs, :pod_write_vtk,  true) == true,
                           get(inputs, :pod_write_png,  true) == true,
                           get(inputs, :pod_write_data, true) == true,
                           Float64(get(inputs, :pod_time_scale, 1.0/86400.0)),
                           String(get(inputs, :pod_time_label, "t [days]")),
                           Symbol(get(inputs, :pod_cmap, :balance)),
                           Float64(get(inputs, :pod_max_memory_gb, 4.0)))
end


#---------------------------------------------------------------------------------
# WHAT can be decomposed.
#
# A target is a name, a component count, and a rule for reading those components
# off the conservative state. Vector-valued targets are decomposed jointly (see
# pod_from_snapshots), which is why :velocity is one target with two components
# and not two targets with one.
#---------------------------------------------------------------------------------
struct St_pod_target
    sym   ::Symbol
    name  ::String
    comps ::Vector{String}
    ncomp ::Int
end

"""
    pod_target(sym) -> St_pod_target

The fields `:pod_fields` may name, for the shallow water system on the shell
(`q = [φ, φu, φv, φw]`, Cartesian momentum):

| symbol       | components                | what it is |
|:-------------|:--------------------------|:--|
| `:vorticity` | `ζ`                       | relative vorticity `n̂·(∇ₛ×u)` — the field the Galewsky test is judged on, and the one whose modes are worth looking at |
| `:h`         | `h`                       | fluid depth `φ/g` |
| `:phi`       | `φ`                       | geopotential, i.e. `h` unscaled |
| `:u`         | `u_λ`                     | zonal velocity |
| `:v`         | `u_φ`                     | meridional velocity |
| `:velocity`  | `u_λ, u_φ`                | the horizontal velocity as ONE vector target |
| `:state`     | `φ, φu, φv, φw`           | the conservative state itself — the basis a Galerkin ROM of this system is projected onto |
"""
function pod_target(sym::Symbol)
    sym === :vorticity && return St_pod_target(sym, "vorticity", ["vorticity"], 1)
    sym === :h         && return St_pod_target(sym, "h",         ["h"],         1)
    sym === :phi       && return St_pod_target(sym, "phi",       ["phi"],       1)
    sym === :u         && return St_pod_target(sym, "u",         ["u"],         1)
    sym === :v         && return St_pod_target(sym, "v",         ["v"],         1)
    sym === :velocity  && return St_pod_target(sym, "velocity",  ["u", "v"],    2)
    sym === :state     && return St_pod_target(sym, "state",     ["phi", "phiu", "phiv", "phiw"], 4)
    error(string(" # ERROR pod.jl: :pod_fields names :", sym, ", which is not a POD target.\n",
                 " #   Available: :vorticity, :h, :phi, :u, :v, :velocity, :state."))
end

pod_needs_vorticity(tg::St_pod_target) = tg.sym === :vorticity


#
# Read one snapshot of one target out of the conservative state.
#
# The velocity components are the ZONAL and MERIDIONAL ones, projected onto the
# local tangent basis exactly as user_primitives.jl does for the VTK output:
# the raw Cartesian components of a flow on a sphere are a property of the frame
# and not of the flow, and their POD modes would be dominated by that frame — a
# rigid zonal jet would come out as a dipole straddling the prime meridian.
#
function pod_extract!(dest::AbstractMatrix, tg::St_pod_target,
                      u::AbstractMatrix, ζ, mesh, g::Float64)

    npoin = size(dest, 1)

    if tg.sym === :vorticity
        ζ === nothing && error(" # ERROR pod.jl: the vorticity target was asked for a snapshot without a vorticity field.")
        @inbounds for ip = 1:npoin
            dest[ip,1] = ζ[ip]
        end

    elseif tg.sym === :phi
        @inbounds for ip = 1:npoin
            dest[ip,1] = u[ip,1]
        end

    elseif tg.sym === :h
        @inbounds for ip = 1:npoin
            dest[ip,1] = u[ip,1]/g
        end

    elseif tg.sym === :state
        @inbounds for ip = 1:npoin, c = 1:4
            dest[ip,c] = u[ip,c]
        end

    else   # :u, :v, :velocity — the tangent-basis components
        iu = tg.sym === :v ? 2 : 1                     # which component goes to column 1
        @inbounds for ip = 1:npoin
            φ  = u[ip,1]
            ux, uy, uz = u[ip,2]/φ, u[ip,3]/φ, u[ip,4]/φ
            sλ, cλ = sincos(mesh.lon[ip])
            sφ, cφ = sincos(mesh.lat[ip])
            uλ = -sλ*ux + cλ*uy                        # u·e_λ  eastward
            uφ = -sφ*cλ*ux - sφ*sλ*uy + cφ*uz          # u·e_φ  northward
            if tg.ncomp == 2
                dest[ip,1] = uλ
                dest[ip,2] = uφ
            else
                dest[ip,1] = iu == 1 ? uλ : uφ
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
    g        ::TF
    lvort    ::Bool
end


"""
    pod_recorder(inputs, mesh, tinit, tend; verbose = true) -> St_pod_recorder | nothing

Allocate the snapshot buffers, or return `nothing` when `:lpod` is off — which
is what every call site tests, so that a case that does not ask for POD carries
no cost and no branch beyond a `=== nothing`.
"""
function pod_recorder(inputs, mesh, tinit::Real, tend::Real;
                      verbose::Bool = true, dt_step::Real = 0.0)

    set = pod_settings(inputs, tinit, tend)
    set.lpod || return nothing

    targets  = St_pod_target[pod_target(f) for f in set.fields]
    isempty(targets) && error(" # ERROR pod.jl: :lpod => true but :pod_fields is empty.")

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
    g  = Float64(get(inputs, :galewsky_gravity, 9.80616))

    #
    # Nothing can be sampled faster than it is computed. A deck that asks for
    # more snapshots than there are time steps in the window gets the steps, and
    # is told so rather than being left to wonder why :pod_nsnapshots was not
    # honoured (pod_record! does the catching up).
    #
    if dt_step > 0 && dt < dt_step
        @warn string("POD: :pod_nsnapshots => ", set.nsnap, " asks for a snapshot every ",
                     @sprintf("%.3f", dt), " s, which is shorter than the time step (",
                     @sprintf("%.3f", dt_step), " s). At most ",
                     floor(Int, (set.tend - set.tstart)/dt_step) + 1,
                     " snapshots can be taken, one per step.")
    end

    if verbose
        println(" # ")
        println(" # POD — Proper Orthogonal Decomposition ........................")
        @printf(" #   fields          : %s\n", join((tg.name for tg in targets), ", "))
        @printf(" #   snapshots       : %d, every %.1f s (%.3f d), t = %.1f … %.1f s\n",
                nsnapmax, dt, dt/86400, set.tstart, set.tend)
        @printf(" #   inner product   : the SEM mass matrix (discrete L² on the shell)\n")
        @printf(" #   mean            : %s\n", set.lmean ? "subtracted (POD of the fluctuation)" : "kept (POD of the field)")
        @printf(" #   memory          : %.1f MB held for the run\n", bytes/2^20)
    end

    return St_pod_recorder{Float64}(set, targets, X, zeros(Float64, nsnapmax), 0, nsnapmax,
                                    dt, Float64(set.tstart), npoin, g,
                                    any(pod_needs_vorticity, targets))
end


"""
    pod_due(rec, t, dt) -> Bool

Is `t` a snapshot time? The `1e-9·dt` slack is the one the output monitor uses:
the sampling times are exact multiples of the POD interval, the integrator lands
on multiples of Δt, and without the slack a snapshot is missed whenever the two
differ in the last bit.
"""
pod_due(rec::Nothing, t, dt) = false
pod_due(rec::St_pod_recorder, t, dt) =
    rec.nsnap < rec.nsnapmax && t >= rec.tnext - 1.0e-9*abs(dt)


"""
    pod_record!(rec, u, ζ, t, mesh)

Take one snapshot of every target from the state `u`. `ζ` may be `nothing` when
no target needs the vorticity.
"""
pod_record!(rec::Nothing, u, ζ, t, mesh) = nothing

function pod_record!(rec::St_pod_recorder, u::AbstractMatrix, ζ, t::Real, mesh)

    rec.nsnap < rec.nsnapmax || return nothing

    k = rec.nsnap + 1
    for (m, tg) in enumerate(rec.targets)
        pod_extract!(view(rec.X[m], :, :, k), tg, u, ζ, mesh, rec.g)
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


#---------------------------------------------------------------------------------
# The inner product, and the two parallel hooks.
#---------------------------------------------------------------------------------
"""
    pod_weights(mesh, metrics) -> w

The diagonal mass matrix, with the nodes this rank does not OWN zeroed out. A
node on a partition seam lives on every rank that touches it and carries the
same value on each; counting it once per rank would weight the seams by the
partition, and the decomposition would then depend on the rank count. This is
the same ownership test `_sphere_integrals` applies to the conserved integrals.
"""
function pod_weights(mesh, metrics)
    npoin = Int(mesh.npoin)
    w     = Vector{Float64}(undef, npoin)
    @inbounds for ip = 1:npoin
        w[ip] = Float64(metrics.M[ip])
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
    pod_finalize!(rec, mesh, metrics, OUTPUT_DIR; verbose = true) -> Vector{St_pod}

Decompose every recorded field and write the results. Called once, at the end of
the run; returns the decompositions so that a caller (a ROM driver, a test) can
go on using them in memory.
"""
pod_finalize!(rec::Nothing, mesh, metrics, OUTPUT_DIR; verbose = true) = St_pod{Float64}[]

function pod_finalize!(rec::St_pod_recorder, mesh, metrics, OUTPUT_DIR::String;
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

    w    = pod_weights(mesh, metrics)
    red! = _pod_reducer(comm)
    sgn  = _pod_signreduce(comm)

    if verbose
        println(" # ")
        println(" # POD ..........................................................")
        @printf(" #   %d snapshots, t = %.1f … %.1f s (%.3f … %.3f d)\n",
                rec.nsnap, rec.t[1], rec.t[rec.nsnap], rec.t[1]/86400, rec.t[rec.nsnap]/86400)
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
        # asks for the POD of h then asks for the decomposition of a constant.
        # Say so and move on to the next field rather than killing a run that
        # was otherwise fine.
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

        rec.set.lvtk  && pod_write_vtk(P, mesh, OUTPUT_DIR; verbose = verbose)
        rec.set.ldata && rank0 && pod_write_data(P, OUTPUT_DIR; verbose = verbose)

        if rec.set.lpng
            if lserial
                pod_plot(P, mesh, rec.set, OUTPUT_DIR; verbose = verbose)
            elseif verbose
                println(" #   PNG: skipped — the equirectangular raster needs the whole sphere on one rank.")
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


"""
    pod_write_vtk(P, mesh, OUTPUT_DIR; verbose = true)

The mean and the modes as point data on the spherical grid, in ONE file per
field (`pod_<field>.vtu`, or `.pvtu` under MPI). Written through the same
`write_vtk_sphere_grid` the solution uses, so the modes land on exactly the
geometry they were computed on and can be overlaid on it in ParaView.
"""
function pod_write_vtk(P::St_pod, mesh, OUTPUT_DIR::String; verbose::Bool = true)

    extras = Pair{String,Vector{Float64}}[]
    for c = 1:P.ncomp
        cn = P.ncomp == 1 ? "" : string("_", P.comps[c])
        P.lmean && push!(extras, string("pod_mean", cn) => collect(view(P.q̄, :, c)))
        for i = 1:length(P.λ)
            push!(extras, @sprintf("pod_mode_%03d%s", i, cn) => collect(view(P.Φ, :, c, i)))
        end
    end

    fname = string("pod_", P.name)
    write_vtk_sphere_grid(mesh, fname, OUTPUT_DIR; extra = extras, verbose = false)

    verbose && @printf(" #     %s%s   mean + %d modes\n",
                       joinpath(abspath(OUTPUT_DIR), fname),
                       MPI.Comm_size(get_mpi_comm()) > 1 ? ".pvtu" : ".vtu", length(P.λ))
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

    fbas = joinpath(OUTPUT_DIR, string("pod_", P.name, ".jld2"))
    JLD2.jldsave(fbas;
                 name = P.name, comps = P.comps,
                 npoin = P.npoin, ncomp = P.ncomp, nsnap = P.nsnap,
                 t = P.t, mean = P.q̄, modes = P.Φ, lambda = P.λ, coefficients = P.a,
                 energy = P.energy, cumenergy = P.cumenergy, total = P.total,
                 lmean = P.lmean, method = String(P.method), orthoerr = P.orthoerr,
                 format = "jexpresso-pod-1")

    verbose && @printf(" #     %s{_spectrum.csv,_coefficients.csv,.jld2}\n",
                       joinpath(abspath(OUTPUT_DIR), string("pod_", P.name)))
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
