#!/usr/bin/env julia
#
# pick_nranks.jl -- choose the MPI rank count for a Jexpresso run that uses
#                   :lxy_partition => true (column decomposition).
#
#   julia tools/pick_nranks.jl <nelemx> <nelemy> <nelemz> <nop> <max_cores> [Lx] [Ly]
#
#   e.g.  julia tools/pick_nranks.jl 128 128 120 4 2048
#         julia tools/pick_nranks.jl  16  16 125 4  256
#
# WHY THIS EXISTS
# ---------------
# With :lxy_partition the mesh is cut into VERTICAL COLUMNS: z is never
# partitioned, so the usable parallelism is bounded by nelemx*nelemy, and the
# split _compute_xy_partition picks is aspect-ratio aware. That combination
# makes the safe rank counts non-obvious: asking for more ranks than the chosen
# nx*ny grid can fill leaves ranks owning ZERO elements, which is fatal
# downstream (it shows up as "Min elements: 0" in the Load Balance Analysis
# block, then a crash after the high-order node phase).
#
# Rather than approximate that, this script REPLICATES _compute_xy_partition
# (src/kernel/mesh/mesh.jl) exactly -- including the centroid-span denominator
# and the clamp -- and reports what each candidate rank count would actually
# produce.
#
# THE CLOSED FORM, for the common square case (nelemx == nelemy == N):
#
#     use nparts = k^2  where k divides N
#
#   and pick the largest such k satisfying BOTH
#       k^2 <= max_cores                       (your allocation)
#       total_gridpoints / k^2 >= pts_target   (enough work per rank)
#
# Everything else here is that rule, generalised and checked against the real
# partitioner.

# ---- exact replication of _compute_xy_partition (mesh.jl:1720) -------------
function partition_counts(nelemx, nelemy, nparts, Lx, Ly)
    cx = [(i + 0.5) * Lx / nelemx for i in 0:nelemx-1]
    cy = [(j + 0.5) * Ly / nelemy for j in 0:nelemy-1]
    lx = maximum(cx) - minimum(cx) + 1e-10
    ly = maximum(cy) - minimum(cy) + 1e-10

    divisors  = [d for d in 1:nparts if nparts % d == 0]
    target_nx = sqrt(nparts * lx / ly)
    nx = divisors[argmin(abs.(divisors .- target_nx))]
    ny = nparts ÷ nx

    xmin, ymin = minimum(cx), minimum(cy)
    counts = zeros(Int, nparts)
    for x in cx, y in cy
        xi = clamp(floor(Int, (x - xmin) / lx * nx), 0, nx - 1)
        yi = clamp(floor(Int, (y - ymin) / ly * ny), 0, ny - 1)
        counts[xi * ny + yi + 1] += 1
    end
    return nx, ny, counts
end

# ──────────────────────────────────────────────────────────────────────────
# p4est / space-filling-curve mode.
#
#   julia tools/pick_nranks.jl --p4est <nelem> <nop> <max_cores>
#
# Completely different constraint structure from the column mode below, and
# the difference is worth stating because it inverts the usual advice.
#
# The column partitioner bins cells by centroid into an nx x ny grid, so the
# rank count must factor to match the mesh or ranks come out EMPTY and the run
# dies. p4est has no such requirement: it cuts the space-filling curve into
# nparts contiguous chunks of floor(N/P) or ceil(N/P) leaves. Any count works,
# and the element imbalance is at most ONE cell -- negligible except when the
# chunks themselves get tiny. Divisibility buys essentially nothing here; do
# not chase it.
#
# What does bind:
#
#   work per rank      A rank holding too few elements spends its time in halo
#                      exchange rather than in the RHS. This is the real cap on
#                      useful parallelism.
#
#   surface-to-volume  A Morton/Hilbert chunk of e elements is compact, roughly
#                      a cube, so its halo scales as e^(2/3) and halo per
#                      element as ~6/e^(1/3). Halving e per rank multiplies
#                      communication per unit work by 2^(1/3).
#
# Gridpoints are estimated as nelem * nop^3: in a structured mesh each element
# contributes nop^3 unique nodes, the rest being shared with its neighbours.
# The estimate is exact only in the limit of many elements per direction, and
# it always UNDERSTATES, because the outermost layer of nodes is shared with
# nothing. Measured against exact counts at nop=4:
#
#     128x128x120   125.8M vs 126.6M    -0.6%
#      32x32x30      1.97M vs  2.01M    -2.4%
#      16x16x15       246k vs   258k    -4.6%
#       8x2x60         61k vs    72k   -14.2%
#
# Thin or small meshes are where it drifts, since the boundary is a larger
# fraction of them. Erring low means the recommendation is conservative -- it
# credits each rank with slightly less work than it really has, and so asks for
# slightly fewer ranks. If your mesh is small or has a short direction, pass
# the exact gridpoint count via the column mode instead, which computes it.
function p4est_mode(nelem, nop, maxc)
    ptsmin = parse(Int, get(ENV, "JEXPRESSO_PTS_PER_RANK", "50000"))
    dof    = nelem * nop^3
    println("partitioning  : p4est, space-filling curve (any rank count is valid)")
    println("elements      : $(nelem)   nop=$(nop)  ->  ~$(dof) gridpoints")
    println("budget        : $(maxc) cores;  target >= $(ptsmin) gridpoints/rank\n")

    # Largest rank count that still keeps each rank above the work floor...
    rcap = max(1, min(maxc, dof ÷ max(ptsmin, 1)))
    # ...then round DOWN to a multiple of 8. The unrounded cap is usually an
    # awkward number (39, 157, 613): SLURM allocates whole cores, so 39 ranks
    # leaves most of a 64-core node idle and bills you for it either way. A
    # multiple of 8 packs into any common node width, and giving up a few ranks
    # of the cap costs nothing -- the work floor is a soft threshold, not a
    # cliff. p4est is indifferent to the exact number, so this is free.
    rec = rcap >= 8 ? (rcap ÷ 8) * 8 : rcap

    cands = Int[]
    r = 1
    while r <= maxc
        push!(cands, r); r *= 2
    end
    rec  in cands || push!(cands, rec)
    rcap in cands || push!(cands, rcap)
    maxc in cands || push!(cands, maxc)
    sort!(unique!(cands))

    println(rpad("ranks",8), rpad("elem/rank",12), rpad("dof/rank",12),
            rpad("imbalance",11), rpad("halo/elem",11), "note")
    for n in cands
        n > maxc && continue
        lo, hi = fld(nelem, n), cld(nelem, n)
        lo == 0 && continue
        imb  = hi / (nelem / n)
        halo = 6 / cbrt(nelem / n)
        note = n == rec ? "  <== RECOMMENDED" :
               (dof / n < ptsmin ? "  (thin: comms-bound)" : "")
        println(rpad(n,8), rpad("$(lo)-$(hi)",12), rpad(round(Int, dof/n),12),
                rpad(round(imb, digits=4),11), rpad(round(halo, digits=2),11), note)
    end

    lo, hi = fld(nelem, rec), cld(nelem, rec)
    println("\nRECOMMENDED: $(rec) ranks")
    println("   $(lo)-$(hi) elements/rank, ~$(round(Int, dof/rec)) gridpoints/rank,")
    println("   halo/elem ~$(round(6/cbrt(nelem/rec), digits=2)), imbalance $(round(hi/(nelem/rec), digits=4))x")
    rec == rcap || println("   ($(rcap) is the raw cap at this work floor; rounded down to pack into nodes)")
    nodes = cld(rec, 64)
    println("\n  #SBATCH --nodes=$(nodes)")
    println("  #SBATCH --ntasks-per-node=$(cld(rec, nodes))")
    println("\nAny other count up to $(maxc) also runs -- p4est does not care about")
    println("divisibility. Going above $(rcap) keeps working but each rank drops")
    println("below $(ptsmin) gridpoints, so more of its time goes to halo exchange")
    println("than to the RHS. Raise the floor with JEXPRESSO_PTS_PER_RANK to be")
    println("more conservative, or lower it to trade efficiency for wall-clock.")
    return nothing
end


#-----------------------------------------------------------------------------
# --columns : rank counts that keep vertical columns on one rank
#
# HEVI's implicit stage is a column solve, and it is cheapest when a column
# lives entirely on one rank. Under p4est that is not automatic. Leaves are
# ordered by (tree, Morton-within-tree) and the partition cuts contiguous
# ranges of that order, so:
#
#   * a partial tree is a Morton sub-cube -- it spans only part of the tree's
#     vertical extent, so column locality does NOT survive a cut inside a tree;
#   * whole trees do keep their columns, so if the rank count divides the
#     number of TREE COLUMNS (nelemx * nelemy of the COARSE mesh) and each rank
#     gets a whole number of them, every column is rank-local.
#
# The refinement level does not change this: refining multiplies the leaves per
# tree, not the tree layout. What it does change is how many elements a rank
# ends up with, which is the other half of the trade -- one tree column of a
# 16x16x15 coarse mesh refined three times is 7680 elements, which is a lot to
# put on one rank.
#
# This is guidance, not a guarantee: the coarse mesh must also be ordered
# column-major (tools/reorder_msh_columns.jl) for tree columns to be
# contiguous in the first place. The ground truth is the split percentage the
# HEVI setup prints at the top of a run -- if that says 0%, the partition is
# column-local whatever this tool predicted.
#-----------------------------------------------------------------------------
function columns_mode(cx, cy, cz, lvl, maxc)
    ncol   = cx * cy
    per    = 2^lvl
    trees  = cx * cy * cz
    leaves = trees * per^3
    println("coarse mesh   : $(cx) x $(cy) x $(cz) trees  ($(trees) trees)")
    println("refine level  : $(lvl)  ->  $(cx*per) x $(cy*per) x $(cz*per) elements ",
            "($(leaves) leaves)")
    println("tree columns  : $(ncol)")
    println("budget        : $(maxc) ranks\n")
    println("  ranks   tree cols/rank   elements/rank   note")
    any = false
    for n = 1:maxc
        ncol % n == 0 || continue
        any = true
        cpr = ncol ÷ n
        epr = leaves ÷ n
        note = epr > 4000 ? "heavy: consider a finer coarse mesh" :
               epr < 200  ? "light: halo may dominate the RHS"    : ""
        println(lpad(n, 7), lpad(cpr, 16), lpad(epr, 16), "   ", note)
    end
    any || println("  none -- $(ncol) tree columns has no divisor <= $(maxc)")
    println()
    println("Any rank count NOT in this list still runs: columns that straddle a rank")
    println("boundary are gathered onto an owner for the solve and scattered back. The")
    println("cost is proportional to how many straddle, which the run reports at setup.")
    return nothing
end

function main()
    if "--columns" in ARGS
        a = filter(x -> x != "--columns", ARGS)
        if length(a) < 5
            println("usage: julia tools/pick_nranks.jl --columns <coarse_nx> <coarse_ny> ",
                    "<coarse_nz> <refine_lvl> <max_ranks>")
            return
        end
        return columns_mode(parse.(Int, a[1:5])...)
    end
    if "--p4est" in ARGS
        a = filter(x -> x != "--p4est", ARGS)
        if length(a) < 3
            println("usage: julia tools/pick_nranks.jl --p4est <nelem> <nop> <max_cores>")
            return
        end
        return p4est_mode(parse.(Int, a[1:3])...)
    end
    if length(ARGS) < 5
        println("usage: julia tools/pick_nranks.jl <nelemx> <nelemy> <nelemz> <nop> <max_cores> [Lx] [Ly]")
        return
    end
    nelemx, nelemy, nelemz, nop, maxc = parse.(Int, ARGS[1:5])
    Lx = length(ARGS) >= 6 ? parse(Float64, ARGS[6]) : Float64(nelemx)
    Ly = length(ARGS) >= 7 ? parse(Float64, ARGS[7]) : Float64(nelemy)

    # Target gridpoints per rank. Below ~20k a spectral-element CPU rank spends
    # more time in halo exchange than in the RHS; 50k keeps it comfortably
    # compute-bound. Override with JEXPRESSO_PTS_PER_RANK.
    pts_target = parse(Int, get(ENV, "JEXPRESSO_PTS_PER_RANK", "50000"))

    totpts = (nelemx*nop + 1) * (nelemy*nop + 1) * (nelemz*nop + 1)
    println("mesh      : $(nelemx) x $(nelemy) x $(nelemz) elements, nop=$(nop)")
    println("domain    : Lx/Ly = $(round(Lx/Ly, digits=3))")
    println("gridpoints: $(totpts)")
    println("budget    : $(maxc) cores;  target >= $(pts_target) points/rank\n")

    rows = NamedTuple[]
    for n in 1:maxc
        nx, ny, c = partition_counts(nelemx, nelemy, n, Lx, Ly)
        mn, mx = minimum(c), maximum(c)
        mn == 0 && continue                       # empty rank => crash
        push!(rows, (n=n, nx=nx, ny=ny, mn=mn, mx=mx,
                     imb = mx/(sum(c)/n),
                     pts = totpts/n,
                     halo = 2*(nx/nelemx + ny/nelemy)))
    end

    # Valid = no empty ranks AND every rank owning the SAME number of columns.
    # Anything else wastes cores waiting on the heaviest rank at every halo
    # exchange, and the near-miss counts (nparts prime, or not factorable to
    # match the element grid) are never worth the extra cores.
    ok = filter(r -> r.mn == r.mx, rows)
    best = filter(r -> r.pts >= pts_target, ok)
    rec  = isempty(best) ? (isempty(ok) ? nothing : ok[1]) : best[end]

    println("COLUMN HEADINGS")
    println("  rank grid   : how the RANKS are arranged (nx x ny), NOT your element grid")
    println("  block       : element columns each rank owns, as (nelemx/nx) x (nelemy/ny)")
    println("  cols/rank   : total columns per rank = the two block numbers multiplied")
    println("  halo/elem   : halo faces per owned element -- pure communication overhead\n")
    println(rpad("ranks",7), rpad("rank grid",12), rpad("block",10), rpad("cols/rank",11),
            rpad("pts/rank",11), rpad("halo/elem",11), "note")
    for r in ok
        note = rec !== nothing && r.n == rec.n ? "  <== RECOMMENDED" :
               (r.pts < pts_target ? "  (thin: comms-bound)" : "")
        println(rpad(r.n,7), rpad("$(r.nx) x $(r.ny)",12),
                rpad("$(nelemx÷r.nx) x $(nelemy÷r.ny)",10),
                rpad(r.mn == r.mx ? "$(r.mn)" : "$(r.mn)-$(r.mx)", 11),
                rpad(round(Int, r.pts), 11),
                rpad(round(r.halo, digits=2), 11), note)
    end

    if rec === nothing
        println("\nNo valid rank count found at or below $(maxc).")
    else
        println("\nRECOMMENDED: $(rec.n) ranks")
        println("   ranks arranged $(rec.nx) x $(rec.ny) over your $(nelemx) x $(nelemy) element columns,")
        println("   so each rank owns a $(nelemx÷rec.nx) x $(nelemy÷rec.ny) block of columns " *
                "($(rec.mn) columns, full height),")
        println("   $(round(Int, rec.pts)) gridpoints/rank, halo/elem $(round(rec.halo, digits=2)).")
        if rec.nx != rec.ny
            println("\n   NOTE: $(rec.nx) x $(rec.ny) is not square because $(rec.n) is not a perfect")
            println("   square, so the blocks are elongated and the halo is larger than it")
            println("   needs to be. If you would rather have square blocks, use the largest")
            println("   k*k below with k dividing $(nelemx).")
        end
        nodes = cld(rec.n, 64)
        println("\n  #SBATCH --nodes=$(nodes)")
        println("  #SBATCH --ntasks-per-node=$(cld(rec.n, nodes))")
        println("\n  (64 ranks/node keeps per-node memory low during the mesh")
        println("   broadcast, which is the peak of the whole run.)")
    end
    println("\nRank counts NOT listed above would leave some ranks with zero")
    println("elements, or badly imbalanced. Do not use them.")
end

main()
