#!/usr/bin/env julia
#
# plan_run.jl -- work out the whole configuration for a Jexpresso LES run:
#                coarse mesh, decomposition path, rank count, SBATCH geometry.
#
#   julia tools/plan_run.jl <NX> <NY> <NZ> <NOP> <LVL> <MAXCORES> <LX> <LY> [--machine]
#
#     NX,NY,NZ   the FINE grid the simulation should run on
#     LVL        p4est initial-refinement levels (0 = read the fine mesh directly)
#     MAXCORES   upper bound on ranks
#
# Emits a human summary, or KEY=VALUE lines for a shell to eval with --machine.
#
# The two paths have DIFFERENT rank rules, which is the entire reason this
# exists -- getting it wrong costs a queue cycle and a crash, and the two
# failure modes look nothing alike:
#
#   direct (:lxy_partition => true, :linitial_refine => false)
#       Cells are binned by centroid into an nx x ny grid chosen as the divisor
#       of nparts nearest sqrt(nparts*lx/ly). A rank whose bin catches no cells
#       owns nothing and the run dies. Valid rank counts are found by
#       replicating that binning, not by dividing.
#
#   p4est (:linitial_refine => true)
#       Each rank gets a CONTIGUOUS range of the space-filling curve over the
#       coarse trees. With a column-major coarse mesh a rank's range is a whole
#       number of columns exactly when nparts divides the coarse COLUMN count
#       (cnx*cny) -- the vertical extent cancels. Anything else splits a column
#       and only its bottom piece owns ground.

# ---- exact replication of _compute_xy_partition (src/kernel/mesh/mesh.jl) ----
function direct_counts(nelemx, nelemy, nparts, Lx, Ly)
    cx = [(i + 0.5) * Lx / nelemx for i in 0:nelemx-1]
    cy = [(j + 0.5) * Ly / nelemy for j in 0:nelemy-1]
    lx = maximum(cx) - minimum(cx) + 1e-10
    ly = maximum(cy) - minimum(cy) + 1e-10
    divisors  = [d for d in 1:nparts if nparts % d == 0]
    nx = divisors[argmin(abs.(divisors .- sqrt(nparts * lx / ly)))]
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
    if length(ARGS) < 8
        println("usage: julia tools/plan_run.jl <NX> <NY> <NZ> <NOP> <LVL> <MAXCORES> <LX> <LY> [--machine]")
        exit(1)
    end
    NX, NY, NZ, NOP, LVL, MAXC = parse.(Int, ARGS[1:6])
    LX, LY = parse.(Float64, ARGS[7:8])
    machine = "--machine" in ARGS
    ptsmin  = parse(Int, get(ENV, "JEXPRESSO_PTS_PER_RANK", "50000"))

    totpts = (NX*NOP + 1) * (NY*NOP + 1) * (NZ*NOP + 1)
    f = 2^LVL
    if LVL > 0 && (NX % f != 0 || NY % f != 0 || NZ % f != 0)
        println(stderr, "ERROR: fine grid $(NX)x$(NY)x$(NZ) is not divisible by 2^$(LVL) = $f.")
        println(stderr, "       Pick a level whose factor divides all three, or change the grid.")
        exit(1)
    end
    cnx, cny, cnz = NX ÷ f, NY ÷ f, NZ ÷ f

    # candidate rank counts, per path
    cands = NamedTuple[]
    if LVL == 0
        for n in 1:MAXC
            nx, ny, c = direct_counts(NX, NY, n, LX, LY)
            (minimum(c) == maximum(c) && minimum(c) > 0) || continue
            push!(cands, (n=n, nx=nx, ny=ny, ground=(NX*NY) ÷ n,
                          halo=2*(nx/NX + ny/NY)))
        end
    else
        cols = cnx * cny                      # coarse columns
        for n in 1:MAXC
            cols % n == 0 || continue         # must be whole columns
            colsper = cols ÷ n
            push!(cands, (n=n, nx=0, ny=0, ground=colsper * 4^LVL,
                          halo=2*sqrt(1/colsper)))   # ~ square block estimate
        end
    end
    isempty(cands) && (println(stderr, "ERROR: no valid rank count <= $MAXC."); exit(1))

    fat = filter(c -> totpts / c.n >= ptsmin, cands)
    rec = isempty(fat) ? cands[1] : fat[end]

    ranks = rec.n
    tpn   = min(ranks, 64)                    # <=64 ranks/node keeps setup memory sane
    nodes = cld(ranks, tpn)
    tpn   = cld(ranks, nodes)

    if machine
        # NB: DECOMP_PATH, never PATH -- the caller eval's these lines and
        # "PATH=..." would clobber the shell's command search path.
        println("DECOMP_PATH=", LVL == 0 ? "direct" : "p4est")
        println("FINE_NX=$NX"); println("FINE_NY=$NY"); println("FINE_NZ=$NZ")
        println("COARSE_NX=$cnx"); println("COARSE_NY=$cny"); println("COARSE_NZ=$cnz")
        println("REFINE_LVL=$LVL")
        println("RANKS=$ranks"); println("NODES=$nodes"); println("NTASKS_PER_NODE=$tpn")
        println("GROUND_PER_RANK=$(rec.ground)")
        println("PTS_PER_RANK=$(round(Int, totpts / ranks))")
        println("COLUMNS=$(LVL == 0 ? NX*NY : cnx*cny)")
        println("VALID_RANKS=\"", join((c.n for c in cands), " "), "\"")
    else
        println("fine grid      : $(NX) x $(NY) x $(NZ)  (nop=$NOP, $(totpts) gridpoints)")
        println("path           : ", LVL == 0 ? "direct read + lxy_partition" :
                                     "p4est, $(LVL) refinement level(s)")
        LVL > 0 && println("coarse mesh    : $(cnx) x $(cny) x $(cnz)  ($(cnx*cny*cnz) elements)")
        println("recommended    : $(ranks) ranks  ->  --nodes=$(nodes) --ntasks-per-node=$(tpn)")
        println("  gridpoints/rank : $(round(Int, totpts / ranks))")
        println("  ground cells/rank: $(rec.ground)   (must be > 0 on every rank)")
        println("valid rank counts: ", join((c.n for c in cands), ", "))
    end
end

main()
