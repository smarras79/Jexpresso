#=============================================================================
 columns.jl -- vertical column topology and the MPI redistribution that makes
 a column-local implicit solve possible under an arbitrary partition.

 THE PROBLEM
 -----------
 A HEVI split is only cheap because the implicit operator couples degrees of
 freedom exclusively within a vertical column: the Jacobian is block-banded
 per column and a direct banded LU costs nothing next to one explicit RHS.
 That argument holds only if each column lives entirely on one rank.

 Jexpresso's large LES runs do NOT satisfy that. They go through p4est
 precisely because the alternative -- reading the full 1.97M-element mesh on
 rank 0 and partitioning it in x-y with :lxy_partition -- is what ran the
 nodes out of memory. p4est orders leaves by (tree, Morton-within-tree) and
 partitions contiguous ranges of that order, so a rank's subdomain is a run
 of whole trees plus a Morton sub-cube at each end. Column locality therefore
 holds only at *tree-column* granularity, and only when the rank count
 divides the number of tree columns (see tools/pick_nranks.jl --columns).

 Rather than make HEVI conditional on a partition the production run cannot
 use, this module makes it work under any partition:

   * columns that happen to be rank-local are solved in place (fast path,
     no communication at all);
   * columns split across ranks are gathered onto an owner, solved, and
     scattered back.

 The owner is the rank holding the most levels of the column, so under a
 sensible partition the volume that actually moves is the two partial columns
 at each end of a rank's SFC range -- kilobytes, not megabytes. The report
 printed at setup states what fraction is split, so a partition that is
 accidentally terrible is visible rather than merely slow.

 HOW A COLUMN IS IDENTIFIED
 --------------------------
 By a *global structured index*, not by a floating-point hash, because two
 ranks holding the same node compute its coordinates from different elements
 and can disagree in the last bits.

   1. Every rank contributes its distinct x values; they are Allgathered and
      clustered with a tolerance into one global catalogue. Same for y. The
      clustering is deterministic given the same global multiset, so every
      rank builds a bit-identical catalogue. Column id = (ix, iy).
   2. Column bottom/top heights are reduced globally, and each node's
      normalised height (z - z_bot)/(z_top - z_bot) is clustered the same
      way into a level catalogue. Normalising rather than using z directly
      is what lets a terrain-following (warped) mesh work: the physical z of
      level k differs from column to column, the normalised height does not.

 Periodic images (x = 0 and x = xmax carry the same gip but different
 coordinates) come out as two distinct columns holding identical data. They
 are solved redundantly and give identical answers -- wasteful for the ~0.4%
 of columns on a periodic face, and not worth special-casing.
=============================================================================#

#-----------------------------------------------------------------------------
# Tracing
#
# "It hangs on the first step" is not a diagnosis: on a problem this size the
# first solve is also where Julia compiles perform_step!, the whole implicit
# chain and the integrator specialised on the real CallbackSet, which can take
# minutes and looks identical from outside. These probes make the difference
# visible in one run.
#
#     JEXPRESSO_HEVI_TRACE=1 julia --project=. ...
#
# Off by default and gated on a single Bool load, so the production path pays
# nothing. Only the first few steps are traced -- after that the interesting
# question is answered and the output would just be noise.
#-----------------------------------------------------------------------------
const HEVI_TRACE      = Ref(false)
const HEVI_TRACE_STEPS = Ref(3)
const HEVI_STEP_COUNT = Ref(0)
const _HEVI_T0        = Ref(0.0)

"""
    hevi_trace_init!() -> Bool

Read `JEXPRESSO_HEVI_TRACE` and start the clock. Called once from `build_hevi`.
"""
function hevi_trace_init!()
    on(k) = lowercase(strip(get(ENV, k, "0"))) in ("1", "true", "yes", "on")
    HEVI_TRACE[] = on("JEXPRESSO_HEVI_TRACE")
    HEVI_PROF[]  = HEVI_TRACE[] || on("JEXPRESSO_HEVI_PROFILE")
    HEVI_PROF_EVERY[] = parse(Int, get(ENV, "JEXPRESSO_HEVI_PROFILE_EVERY", "50"))
    HEVI_PROF_SKIP[]  = parse(Int, get(ENV, "JEXPRESSO_HEVI_PROFILE_SKIP", "10"))
    HEVI_MONITOR[]    = on("JEXPRESSO_HEVI_MONITOR")
    HEVI_STEP_COUNT[] = 0
    hevi_profile_reset!()
    (HEVI_TRACE[] || HEVI_PROF[]) && (_HEVI_T0[] = time())
    return HEVI_TRACE[]
end

@inline function hevi_trace(msg...)
    HEVI_TRACE[] || return nothing
    Base.print("   [hevi +", lpad(string(round(time() - _HEVI_T0[]; digits = 2)), 8), " s] ",
               msg..., "\n")
    Base.flush(stdout)
    return nothing
end

#-----------------------------------------------------------------------------
# Cost breakdown
#
# The question a timing comparison actually raises is not "is HEVI slower" but
# "WHERE does a HEVI step spend its time". A step is 4 rhs! evaluations, 4
# applications of the implicit operator and 3 banded column solves; if those
# do not add up to the measured step time, the cost is somewhere else
# entirely, and no amount of tuning the column solver will touch it.
#
# So this accumulates each phase separately and, crucially, reports the GAP
# between their sum and the wall time of the step. A large gap is itself the
# finding.
#
#     JEXPRESSO_HEVI_PROFILE=1   (implied by JEXPRESSO_HEVI_TRACE=1)
#-----------------------------------------------------------------------------
mutable struct HEVIProfile
    t_fimp::Float64;  n_fimp::Int
    t_solve::Float64; n_solve::Int
    t_rhs::Float64;   n_rhs::Int
    t_refac::Float64; n_refac::Int
    t_step::Float64;  n_step::Int
end
const HEVI_PROF       = Ref(false)
const HEVI_PROF_EVERY = Ref(50)
# Steps skipped before accumulating. The first few steps are where Julia
# compiles the RHS chain, so counting them inflates every row -- and rhs! most
# of all, which is exactly the row the other numbers are judged against.
const HEVI_PROF_SKIP  = Ref(10)
const HEVI_PROFILE    = HEVIProfile(0.0, 0, 0.0, 0, 0.0, 0, 0.0, 0, 0.0, 0)

@inline hevi_prof_on() = HEVI_PROF[] && HEVI_STEP_COUNT[] > HEVI_PROF_SKIP[]
@inline hevi_tic() = hevi_prof_on() ? time_ns() : UInt64(0)

function hevi_profile_reset!()
    P = HEVI_PROFILE
    P.t_fimp = 0.0;  P.n_fimp = 0
    P.t_solve = 0.0; P.n_solve = 0
    P.t_rhs = 0.0;   P.n_rhs = 0
    P.t_refac = 0.0; P.n_refac = 0
    P.t_step = 0.0;  P.n_step = 0
    return nothing
end

"""
    hevi_profile_report!()

Print the per-step cost breakdown and reset the counters.

The last two lines are the point. "accounted" is everything HEVI knows it
spent; "measured" is the wall time of the same steps. When they disagree, the
time is going somewhere outside the integrator's own work -- callbacks, saving,
allocation, GC -- and that is where to look next.
"""
function hevi_profile_report!()
    P = HEVI_PROFILE
    P.n_step == 0 && return nothing
    ns = P.n_step
    acc = (P.t_fimp + P.t_solve + P.t_rhs + P.t_refac) / ns
    tot = P.t_step / ns
    row(name, t, n) = Base.print("   │ ", rpad(name, 16),
        lpad(string(round(n / ns; digits = 2)), 6), " /step ",
        lpad(string(round(n > 0 ? t / n : 0.0; digits = 5)), 10), " s each ",
        lpad(string(round(t / ns; digits = 4)), 9), " s/step ",
        lpad(string(round(100 * t / ns / max(tot, eps()); digits = 1)), 6), " %\n")
    Base.print("   ┌── HEVI cost breakdown over ", ns, " steps (first ",
               HEVI_PROF_SKIP[], " skipped: JIT)\n")
    row("rhs!",        P.t_rhs,   P.n_rhs)
    row("f_imp",       P.t_fimp,  P.n_fimp)
    row("column solve", P.t_solve, P.n_solve)
    row("refactorise", P.t_refac, P.n_refac)
    Base.print("   │ ", rpad("accounted", 16), " " ^ 30,
               lpad(string(round(acc; digits = 4)), 9), " s/step\n")
    Base.print("   │ ", rpad("measured step", 16), " " ^ 30,
               lpad(string(round(tot; digits = 4)), 9), " s/step",
               "   gap ", round(tot - acc; digits = 4), " s\n")
    Base.print("   └", "─" ^ 60, "\n")
    Base.flush(stdout)
    hevi_profile_reset!()
    return nothing
end

#-----------------------------------------------------------------------------
# Instability monitor
#
# "It blew up at t = 15 s" is not localising. WHERE it starts is: a mode that
# grows from the floor or the lid points at the boundary treatment, one that
# grows in the interior points at the step size, and one confined to a single
# column points at the redistribution. This reports the largest deviation from
# the reference state per implicit variable, and the height it sits at, so the
# growth has a location before it becomes a DomainError.
#
#     :hevi_monitor => true      (or JEXPRESSO_HEVI_MONITOR=1)
#-----------------------------------------------------------------------------
const HEVI_MONITOR       = Ref(false)
const HEVI_MONITOR_EVERY = Ref(50)

# Fraction of the neutral joint-stability limit a run should sit at. The limit
# itself is where |R| first exceeds 1, and AT it the scheme is exactly neutral:
# the explicit rate that goes into it is a node-spacing estimate calibrated
# against a single measured spectrum, it assumes the flow speeds at t = 0, and
# the assembly at rank-shared nodes is mildly partition-dependent -- so a few
# percent either way decides between neutral and growing. rtb_hevi shipped at
# 97% of the limit and was stable on 1 and 3 ranks and divergent on 10.
const HEVI_DT_SAFETY = Ref(0.70)

# Scratch for the single-valuedness check, built on first use because it needs
# a live operator and is only ever wanted when the monitor is on.
const _HEVI_MULT = Ref{Union{Nothing, Matrix{Float64}}}(nothing)
const _HEVI_SCR  = Ref{Union{Nothing, Matrix{Float64}}}(nothing)

"""
    hevi_single_valued_error(u, p) -> (relerr, ip_worst)

How far the implicit variables are from being single-valued at nodes that more
than one rank holds.

THE POINT OF THIS CHECK
-----------------------
Every other HEVI self-check is self-consistent: `hevi_verify_band` probes the
operator and compares against a matrix built from the same probes, and
`hevi_verify_solve` inverts the matrix it was just given. Both pass whatever
the gather election does, because both feed it a field that is single-valued
by construction. Neither can see the one failure that only exists in parallel:
the implicit stage writing *different* values into the two copies of a node
that sits on a rank interface. A continuous-Galerkin state that has gone
double-valued is not a small error -- the next explicit RHS differentiates
across the jump, and the jump grows.

The measurement needs no serial reference. `assemble_mpi!` leaves every holder
of a node with the SUM over all holders, so a field that is single-valued with
value v comes back as m*v where m is the multiplicity -- and m is measured
once by assembling a field of ones. Anything else is a genuine disagreement
between ranks, in the field's own units.

Returns the max relative disagreement over the implicit variables (already
reduced across ranks) and the local point where it occurred.
"""
function hevi_single_valued_error(u, p)

    op = p.hevi.op
    op.dss_cache === nothing && return (0.0, 0)

    npoin = Int(p.mesh.npoin)
    nimp  = op.nimp

    if _HEVI_MULT[] === nothing || size(_HEVI_MULT[]) != (npoin, nimp)
        m = ones(Float64, npoin, nimp)
        assemble_mpi!(m, op.dss_cache)
        _HEVI_MULT[] = m
        _HEVI_SCR[]  = zeros(Float64, npoin, nimp)
    end
    mult = _HEVI_MULT[]::Matrix{Float64}
    w    = _HEVI_SCR[]::Matrix{Float64}

    @inbounds for (q, ieq) in enumerate(op.vars)
        off = (ieq - 1) * npoin
        for ip = 1:npoin
            w[ip, q] = u[off + ip]
        end
    end
    assemble_mpi!(w, op.dss_cache)

    err = 0.0; ipw = 0; scl = 0.0
    @inbounds for (q, ieq) in enumerate(op.vars)
        off = (ieq - 1) * npoin
        for ip = 1:npoin
            v = u[off + ip]
            scl = max(scl, abs(v))
            d = abs(w[ip, q] - mult[ip, q] * v)
            d > err && (err = d; ipw = ip)
        end
    end

    comm = p.hevi.cc.comm
    if MPI.Comm_size(comm) > 1
        err = MPI.Allreduce(err, MPI.MAX, comm)
        scl = MPI.Allreduce(scl, MPI.MAX, comm)
    end
    return (scl > 0 ? err / scl : err), ipw
end

"""
    hevi_monitor_step!(u, p)

Report max|u - qe| per implicit variable and the height where it occurs.

Printed every `HEVI_MONITOR_EVERY` steps, and immediately if any deviation is
non-finite or has grown past a plainly unphysical size -- which is the point:
the report lands BEFORE the equation of state throws, when the location still
means something.
"""
function hevi_monitor_step!(u, p)
    hevi  = p.hevi
    op    = hevi.op
    mesh  = p.mesh
    npoin = Int(mesh.npoin)
    qe    = p.qp.qe
    z     = @view mesh.coords[3, :]

    bad = false
    msg = IOBuffer()
    print(msg, "   [hevi monitor] ")
    @inbounds for (q, ieq) in enumerate(op.vars)
        off  = (ieq - 1) * npoin
        dmax = 0.0; ipmax = 1
        for ip = 1:npoin
            d = abs(u[off + ip] - qe[ip, ieq])
            (isfinite(d) && d <= dmax) && continue
            isfinite(d) || (bad = true)
            d > dmax && (dmax = d; ipmax = ip)
        end
        qref = 0.0
        for ip = 1:npoin
            qref = max(qref, abs(qe[ip, ieq]))
        end
        rel = qref > 0 ? dmax / qref : dmax
        rel > 0.5 && (bad = true)
        Base.print(msg, "eq", ieq, ": |δ|=", round(dmax; sigdigits = 3),
                   " (", round(100 * rel; sigdigits = 3), "% of qe) at z=",
                   round(z[ipmax]; digits = 1), "m   ")
    end
    # Parallel-only failure mode, invisible to every other self-check: the
    # implicit stage leaving the two copies of a rank-interface node holding
    # different values. Reported the moment it exceeds round-off, because by
    # the time it is visible in |u - qe| it has already been amplified.
    sv, ipsv = hevi_single_valued_error(u, p)
    if sv > 1.0e-12
        bad = true
        Base.print(msg, "|| CG state NOT single-valued across ranks: ",
                   round(sv; sigdigits = 3), " (rel) at z=",
                   round(z[min(max(ipsv, 1), npoin)]; digits = 1), "m ")
    end

    (bad || HEVI_STEP_COUNT[] % HEVI_MONITOR_EVERY[] == 0) &&
        (Base.println(String(take!(msg))); Base.flush(stdout))
    return bad
end

"""
    ColumnTopology

The rank-local view of the global column structure.

  `nlev`      levels per column, uniform and global
  `ncol`      number of columns this rank touches at all
  `gcol`      `gcol[icol]` -> global column id
  `node`      `node[lev, icol]` -> local point index, or 0 where this rank
              does not hold that level
  `nheld`     `nheld[icol]` -> number of levels held locally
  `ncolg`     total number of global columns
  `nsplit`    number of local columns this rank holds only part of
"""
struct ColumnTopology
    nlev::Int
    ncol::Int
    ncolg::Int
    gcol::Vector{Int}
    node::Matrix{Int}
    nheld::Vector{Int}
    nsplit::Int
end

#-----------------------------------------------------------------------------
# Global coordinate catalogues
#-----------------------------------------------------------------------------

"""
    _local_distinct(v, tol) -> Vector{Float64}

Sorted distinct values of `v`, merging anything closer together than `tol`.
"""
function _local_distinct(v::AbstractVector, tol::Real)
    s = sort!(collect(Float64, v))
    isempty(s) && return Float64[]
    out = Float64[s[1]]
    @inbounds for x in s
        if x - out[end] > tol
            push!(out, x)
        end
    end
    return out
end

"""
    global_catalogue(localvals, comm, tol; label, maxsize) -> Vector{Float64}

Build the globally-agreed sorted catalogue of distinct values.

Every rank contributes its own distinct values, they are Allgathered, and the
union is clustered once more with the same tolerance. Because the clustering
runs on the same global multiset with the same rule on every rank, the
resulting catalogue is bit-identical everywhere -- which is the whole point:
a column index derived from it is a genuine global identifier, not a
per-rank guess.

`maxsize` guards against being handed an unstructured mesh, where the number
of distinct horizontal positions is O(npoin) and this whole scheme is
inapplicable. Hitting it is a hard error with an explanation, not a silent
degradation.
"""
function global_catalogue(localvals::AbstractVector, comm, tol::Real;
                          label::AbstractString = "coordinate",
                          maxsize::Int = 100_000)

    loc = _local_distinct(localvals, tol)
    n   = length(loc)

    if MPI.Comm_size(comm) == 1
        cat = loc
    else
        counts = MPI.Allgather(Int32(n), comm)
        total  = sum(Int, counts)
        total > maxsize * MPI.Comm_size(comm) &&
            error("HEVI: rank-local distinct $label count is unmanageably large ($total). ",
                  "The column solver needs a vertically extruded mesh.")
        allv = Vector{Float64}(undef, total)
        MPI.Allgatherv!(loc, MPI.VBuffer(allv, counts), comm)
        cat = _local_distinct(allv, tol)
    end

    length(cat) > maxsize &&
        error("HEVI: found $(length(cat)) distinct $label values, more than the $maxsize limit. ",
              "This mesh does not look like a vertical extrusion of a horizontal grid; ",
              "the column-based implicit solver cannot be used on it.")
    return cat
end

"""
    catalogue_index(cat, x) -> Int

Index of the catalogue entry nearest to `x`. Binary search; `cat` is sorted.
"""
function catalogue_index(cat::Vector{Float64}, x::Real)
    i = searchsortedfirst(cat, x)
    i <= 1            && return 1
    i > length(cat)   && return length(cat)
    return (x - cat[i-1]) <= (cat[i] - x) ? i - 1 : i
end

#-----------------------------------------------------------------------------
# Topology construction
#-----------------------------------------------------------------------------

"""
    build_column_topology(mesh, comm; rtol = 1.0e-6) -> ColumnTopology

Group every local node into a (global column, level) slot.

`rtol` is relative to the domain extent and only has to separate genuinely
distinct grid positions from round-off between two ranks' evaluations of the
same node -- the gap between real neighbours is ~1e-3 of the extent, the
disagreement between ranks is ~1e-15, so anything in between works.
"""
function build_column_topology(mesh, comm; rtol::Real = 1.0e-6)

    npoin  = Int(mesh.npoin)
    coords = mesh.coords
    x = @view coords[1, 1:npoin]
    y = @view coords[2, 1:npoin]
    z = @view coords[3, 1:npoin]

    Lx = max(Float64(mesh.xmax - mesh.xmin), eps())
    Ly = max(Float64(mesh.ymax - mesh.ymin), eps())
    tolx = rtol * Lx
    toly = rtol * Ly

    xcat = global_catalogue(x, comm, tolx; label = "x")
    ycat = global_catalogue(y, comm, toly; label = "y")
    nx, ny = length(xcat), length(ycat)
    ncolg  = nx * ny

    # Global column id of every local node.
    gcol_of = Vector{Int}(undef, npoin)
    @inbounds for ip = 1:npoin
        ix = catalogue_index(xcat, x[ip])
        iy = catalogue_index(ycat, y[ip])
        gcol_of[ip] = (iy - 1) * nx + ix
    end

    # Column bottom and top, globally. A rank that holds only the upper part
    # of a column does not know its own column's ground level; the reduction
    # is what makes the normalised height below agree across the split.
    BIG  = 1.0e30
    zbot = fill( BIG, ncolg)
    ztop = fill(-BIG, ncolg)
    @inbounds for ip = 1:npoin
        c = gcol_of[ip]
        zz = z[ip]
        zz < zbot[c] && (zbot[c] = zz)
        zz > ztop[c] && (ztop[c] = zz)
    end
    if MPI.Comm_size(comm) > 1
        MPI.Allreduce!(zbot, MPI.MIN, comm)
        MPI.Allreduce!(ztop, MPI.MAX, comm)
    end

    # Normalised height -> level catalogue. On a flat mesh this is just z
    # rescaled; on a terrain-following mesh it is the computational level,
    # which is what actually lines up between columns.
    ζ = Vector{Float64}(undef, npoin)
    @inbounds for ip = 1:npoin
        c = gcol_of[ip]
        span = ztop[c] - zbot[c]
        ζ[ip] = span > 0 ? (z[ip] - zbot[c]) / span : 0.0
    end
    ζcat = global_catalogue(ζ, comm, rtol; label = "normalised height")
    nlev = length(ζcat)

    # Local columns, in ascending global id so the ordering is reproducible.
    present = sort!(unique(gcol_of))
    ncol    = length(present)
    slot    = Dict{Int,Int}(g => i for (i, g) in enumerate(present))

    node = zeros(Int, nlev, ncol)
    @inbounds for ip = 1:npoin
        ic = slot[gcol_of[ip]]
        il = catalogue_index(ζcat, ζ[ip])
        if node[il, ic] != 0 && node[il, ic] != ip
            error("HEVI: two distinct local nodes ($(node[il,ic]) and $ip) landed in the same ",
                  "(column $(present[ic]), level $il) slot. The mesh is not a clean vertical ",
                  "extrusion, or rtol=$rtol is too coarse for its horizontal spacing.")
        end
        node[il, ic] = ip
    end

    nheld  = [count(!iszero, @view node[:, ic]) for ic = 1:ncol]
    nsplit = count(<(nlev), nheld)

    return ColumnTopology(nlev, ncol, ncolg, present, node, nheld, nsplit)
end

#-----------------------------------------------------------------------------
# Ownership and the gather/scatter plan
#-----------------------------------------------------------------------------

"""
    assign_column_owners(topo, comm) -> (owner::Vector{Int}, own::Vector{Int})

Decide which rank solves each column, and list the ones this rank owns.

The owner is the rank holding the most levels; ties go to the lowest rank, so
the assignment is deterministic and independent of message arrival order.
Both criteria ride in a single global MAX-reduction over one integer per
column: `nheld * nranks + (nranks - 1 - rank)` is monotone in `nheld` first
and in `-rank` second.

`owner` is indexed by global column id and is meaningful for every column,
not just local ones. `own` lists the global ids this rank owns.
"""
function assign_column_owners(topo::ColumnTopology, comm)

    rank   = MPI.Comm_rank(comm)
    nranks = MPI.Comm_size(comm)

    score = zeros(Int, topo.ncolg)
    @inbounds for ic = 1:topo.ncol
        score[topo.gcol[ic]] = topo.nheld[ic] * nranks + (nranks - 1 - rank)
    end
    nranks > 1 && MPI.Allreduce!(score, MPI.MAX, comm)

    owner = Vector{Int}(undef, topo.ncolg)
    @inbounds for g = 1:topo.ncolg
        owner[g] = score[g] == 0 ? -1 : (nranks - 1 - (score[g] % nranks))
    end

    own = [g for g = 1:topo.ncolg if owner[g] == rank]
    return owner, own
end

"""
    ColumnComm

Everything needed to move the implicit variables between their home ranks and
the ranks that own the columns, and back.

The owner-side working array `X` is `(nimp*nlev, nown)` with the row index
running variables-fastest, `(il-1)*nimp + m`. That is exactly the ordering
that makes the vertical coupling banded with half-width `nimp*ngl - 1`, and
it makes `view(X, :, ic)` a contiguous vector -- the right-hand side LAPACK's
banded solver wants, with no copy.

Four index lists drive the two directions:

  gather  `g_send_ip[k]`   local points whose values are sent to `g_ranks[k]`
          `g_recv_slot[k]` owner-side slots filled from `g_ranks_in[k]`
  scatter `s_send_slot[k]` owner-side slots sent back to `s_ranks[k]`
          `s_recv_ip[k]`   local points written from `s_ranks_in[k]`

plus `self_*` for the part that never leaves the rank. The gather side is
*elected*: exactly one holder (the lowest-ranked one) contributes each node,
so the assembled column is bit-reproducible rather than depending on which
message landed last.
"""
struct ColumnComm
    comm::MPI.Comm
    rank::Int
    nimp::Int
    nlev::Int
    nown::Int
    X::Matrix{Float64}

    g_ranks::Vector{Int};    g_send_ip::Vector{Vector{Int}};    g_send_buf::Vector{Vector{Float64}}
    g_ranks_in::Vector{Int}; g_recv_slot::Vector{Vector{Int}};  g_recv_buf::Vector{Vector{Float64}}

    s_ranks::Vector{Int};    s_send_slot::Vector{Vector{Int}};  s_send_buf::Vector{Vector{Float64}}
    s_ranks_in::Vector{Int}; s_recv_ip::Vector{Vector{Int}};    s_recv_buf::Vector{Vector{Float64}}

    self_gather_ip::Vector{Int};  self_gather_slot::Vector{Int}
    self_scatter_ip::Vector{Int}; self_scatter_slot::Vector{Int}

    # column-local metadata, owner side
    own_gcol::Vector{Int}
    nsplit_global::Int
    nmoved_global::Int
end

# Slot index used by the gather/scatter plans. `Xf`, a reshape of X to
# (nimp, nlev*nown), addresses one node as Xf[m, slot]; the same memory is
# X[(il-1)*nimp + m, ic].
@inline _xslot(nlev, il, ic) = (ic - 1) * nlev + il

"""
    build_column_comm(topo, owner, own, comm, nimp) -> ColumnComm

Build the gather/scatter plan. One setup-time exchange of index lists; after
that every solve is two rounds of point-to-point messages with preallocated
buffers and no allocation.
"""
function build_column_comm(topo::ColumnTopology, owner::Vector{Int}, own::Vector{Int},
                           comm, nimp::Int)

    rank   = MPI.Comm_rank(comm)
    nranks = MPI.Comm_size(comm)
    nlev   = topo.nlev
    nown   = length(own)
    own_slot = Dict{Int,Int}(g => i for (i, g) in enumerate(own))

    #--- What this rank holds, sorted by destination owner --------------------
    # Report format is a flat [gcol, lev, gcol, lev, ...] Int vector.
    report   = [Int[] for _ = 1:nranks]
    report_ip = [Int[] for _ = 1:nranks]
    self_ip   = Int[];  self_slot = Int[]

    @inbounds for ic = 1:topo.ncol
        g = topo.gcol[ic]
        o = owner[g]
        for il = 1:nlev
            ip = topo.node[il, ic]
            ip == 0 && continue
            if o == rank
                push!(self_ip, ip)
                push!(self_slot, _xslot(nlev, il, own_slot[g]))
            else
                push!(report[o+1], g, il)
                push!(report_ip[o+1], ip)
            end
        end
    end

    #--- Exchange the reports -------------------------------------------------
    send_sizes = Int32[length(report[r]) for r = 1:nranks]
    recv_sizes = MPI.Alltoall(MPI.UBuffer(send_sizes, 1), comm)
    recv = [Vector{Int}(undef, Int(recv_sizes[r])) for r = 1:nranks]

    reqs = MPI.Request[]
    for r = 0:nranks-1
        send_sizes[r+1] > 0 && push!(reqs, MPI.Isend(report[r+1], r, 11, comm))
        recv_sizes[r+1] > 0 && push!(reqs, MPI.Irecv!(recv[r+1], r, 11, comm))
    end
    !isempty(reqs) && MPI.Waitall(reqs)

    #--- Election: lowest-ranked holder contributes each (column, level) ------
    # Candidates are every reporting rank plus this rank's own copies.
    elected = Dict{Int,Int}()          # X-slot -> winning rank
    for r = 0:nranks-1
        v = recv[r+1]
        @inbounds for e = 1:2:length(v)
            g, il = v[e], v[e+1]
            s = _xslot(nlev, il, own_slot[g])
            w = get(elected, s, typemax(Int))
            r < w && (elected[s] = r)
        end
    end
    @inbounds for (k, s) in enumerate(self_slot)
        w = get(elected, s, typemax(Int))
        rank < w && (elected[s] = rank)
    end

    #--- Owner-side lists -----------------------------------------------------
    g_ranks_in = Int[];  g_recv_slot = Vector{Int}[]
    s_ranks    = Int[];  s_send_slot = Vector{Int}[]
    for r = 0:nranks-1
        v = recv[r+1]
        isempty(v) && continue
        take = Int[]; give = Int[]
        @inbounds for e = 1:2:length(v)
            g, il = v[e], v[e+1]
            s = _xslot(nlev, il, own_slot[g])
            push!(give, s)                    # everything reported goes back
            elected[s] == r && push!(take, s) # only the elected one comes in
        end
        if !isempty(take)
            push!(g_ranks_in, r); push!(g_recv_slot, take)
        end
        push!(s_ranks, r); push!(s_send_slot, give)
    end

    # The owner's own copies: it contributes only the elected ones, but it
    # writes back every one it holds.
    self_gather_ip = Int[];  self_gather_slot = Int[]
    @inbounds for (k, s) in enumerate(self_slot)
        if elected[s] == rank
            push!(self_gather_ip, self_ip[k]); push!(self_gather_slot, s)
        end
    end

    #--- Home-side lists ------------------------------------------------------
    # Symmetric to the owner side: this rank must be told, per destination,
    # which of its reported entries it actually has to send. The owner knows
    # the election result; rather than a second exchange, each rank recomputes
    # it -- it cannot, because election needs the global holder set. So the
    # owner replies with a per-entry take/skip mask.
    mask_send = [Int8[] for _ = 1:nranks]
    for (k, r) in enumerate(s_ranks)
        v = recv[r+1]
        m = Vector{Int8}(undef, length(v) >> 1)
        @inbounds for e = 1:2:length(v)
            g, il = v[e], v[e+1]
            s = _xslot(nlev, il, own_slot[g])
            m[(e+1) >> 1] = elected[s] == r ? Int8(1) : Int8(0)
        end
        mask_send[r+1] = m
    end
    mask_recv = [Vector{Int8}(undef, length(report_ip[r])) for r = 1:nranks]
    reqs2 = MPI.Request[]
    for r = 0:nranks-1
        !isempty(mask_send[r+1]) && push!(reqs2, MPI.Isend(mask_send[r+1], r, 12, comm))
        !isempty(mask_recv[r+1]) && push!(reqs2, MPI.Irecv!(mask_recv[r+1], r, 12, comm))
    end
    !isempty(reqs2) && MPI.Waitall(reqs2)

    g_ranks = Int[]; g_send_ip = Vector{Int}[]
    s_ranks_in = Int[]; s_recv_ip = Vector{Int}[]
    for r = 0:nranks-1
        ips = report_ip[r+1]
        isempty(ips) && continue
        push!(s_ranks_in, r); push!(s_recv_ip, copy(ips))   # all of them come back
        m = mask_recv[r+1]
        sel = [ips[k] for k = 1:length(ips) if m[k] == 1]
        if !isempty(sel)
            push!(g_ranks, r); push!(g_send_ip, sel)
        end
    end

    # Every (owned column, level) slot must have exactly one elected
    # contributor, or the gather would leave stale values in X and the
    # implicit solve would quietly use them. On a conforming vertical
    # extrusion this is guaranteed; the check is here because the failure is
    # otherwise silent and would look like a bad implicit operator.
    covered = length(self_gather_ip) + sum(length, g_recv_slot; init = 0)
    covered == nown * nlev ||
        error("HEVI: the gather plan covers $covered of the $(nown*nlev) (column, level) ",
              "slots this rank owns. Some level of some column is held by no rank at all, ",
              "which means the mesh is not a conforming vertical extrusion of its ",
              "horizontal grid.")

    mkbuf(lists) = [Vector{Float64}(undef, nimp * length(l)) for l in lists]

    nmoved = sum(length, g_send_ip; init = 0)
    nsplit_g = nranks > 1 ? MPI.Allreduce(topo.nsplit, MPI.SUM, comm) : topo.nsplit
    nmoved_g = nranks > 1 ? MPI.Allreduce(nmoved,      MPI.SUM, comm) : nmoved

    return ColumnComm(comm, rank, nimp, nlev, nown,
                      zeros(Float64, nimp * nlev, nown),
                      g_ranks, g_send_ip, mkbuf(g_send_ip),
                      g_ranks_in, g_recv_slot, mkbuf(g_recv_slot),
                      s_ranks, s_send_slot, mkbuf(s_send_slot),
                      s_ranks_in, s_recv_ip, mkbuf(s_recv_ip),
                      self_gather_ip, self_gather_slot,
                      self_ip, self_slot,
                      copy(own), nsplit_g, nmoved_g)
end

#-----------------------------------------------------------------------------
# Runtime data movement
#-----------------------------------------------------------------------------

"""
    column_gather!(cc, field, vars)

Fill `cc.X` from `field[ip, vars[m]]`, moving whatever is not rank-local.

`field` is any `npoin x ncols` array (u reshaped, an RHS, a state deviation);
`vars` selects which of its columns are the implicit variables.
"""
function column_gather!(cc::ColumnComm, field::AbstractMatrix, vars)

    nimp = cc.nimp
    X    = cc.X
    Xf   = reshape(X, nimp, cc.nlev * cc.nown)

    @inbounds for (k, ips) in enumerate(cc.g_send_ip)
        buf = cc.g_send_buf[k]
        for (e, ip) in enumerate(ips), m = 1:nimp
            buf[(e-1)*nimp + m] = field[ip, vars[m]]
        end
    end

    reqs = MPI.Request[]
    for (k, r) in enumerate(cc.g_ranks)
        push!(reqs, MPI.Isend(cc.g_send_buf[k], r, 21, cc.comm))
    end
    for (k, r) in enumerate(cc.g_ranks_in)
        push!(reqs, MPI.Irecv!(cc.g_recv_buf[k], r, 21, cc.comm))
    end

    # Local part overlaps with the messages in flight.
    @inbounds for (e, ip) in enumerate(cc.self_gather_ip)
        s = cc.self_gather_slot[e]
        for m = 1:nimp
            Xf[m, s] = field[ip, vars[m]]
        end
    end

    !isempty(reqs) && MPI.Waitall(reqs)

    @inbounds for (k, slots) in enumerate(cc.g_recv_slot)
        buf = cc.g_recv_buf[k]
        for (e, s) in enumerate(slots), m = 1:nimp
            Xf[m, s] = buf[(e-1)*nimp + m]
        end
    end
    return nothing
end

"""
    column_scatter!(cc, field, vars)

Write `cc.X` back into `field[ip, vars[m]]` on every rank that holds the node,
including the duplicated copies at rank interfaces -- which is what keeps the
continuous-Galerkin solution single-valued after the implicit stage.
"""
function column_scatter!(cc::ColumnComm, field::AbstractMatrix, vars)

    nimp = cc.nimp
    Xf   = reshape(cc.X, nimp, cc.nlev * cc.nown)

    @inbounds for (k, slots) in enumerate(cc.s_send_slot)
        buf = cc.s_send_buf[k]
        for (e, s) in enumerate(slots), m = 1:nimp
            buf[(e-1)*nimp + m] = Xf[m, s]
        end
    end

    reqs = MPI.Request[]
    for (k, r) in enumerate(cc.s_ranks)
        push!(reqs, MPI.Isend(cc.s_send_buf[k], r, 22, cc.comm))
    end
    for (k, r) in enumerate(cc.s_ranks_in)
        push!(reqs, MPI.Irecv!(cc.s_recv_buf[k], r, 22, cc.comm))
    end

    @inbounds for (e, ip) in enumerate(cc.self_scatter_ip)
        s = cc.self_scatter_slot[e]
        for m = 1:nimp
            field[ip, vars[m]] = Xf[m, s]
        end
    end

    !isempty(reqs) && MPI.Waitall(reqs)

    @inbounds for (k, ips) in enumerate(cc.s_recv_ip)
        buf = cc.s_recv_buf[k]
        for (e, ip) in enumerate(ips), m = 1:nimp
            field[ip, vars[m]] = buf[(e-1)*nimp + m]
        end
    end
    return nothing
end
