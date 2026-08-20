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

function main()
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
