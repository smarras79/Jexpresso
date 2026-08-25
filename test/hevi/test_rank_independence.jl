#=============================================================================
 test_rank_independence.jl -- does the answer depend on how many ranks
 computed it?

 WHY THIS IS SEPARATE FROM EVERY OTHER CHECK IN THIS DIRECTORY
 -------------------------------------------------------------
 The existing self-checks are all SELF-consistency at one rank count:
 the assembled band matches the operator it was probed from, the factorisation
 inverts the matrix it was handed, the gather/scatter round-trips, f_imp
 matches the vertical part of rhs!. Every one of them can pass on 1 rank and
 on 8 ranks while the two produce DIFFERENT ANSWERS, because each is checked
 against something that was itself computed on the same partition.

 Rank-independence is the property that actually matters to a user -- "is the
 science the same on 1024 cores as on 32" -- and nothing here was testing it.
 This does: it computes a field, gathers it keyed by GLOBAL node id, and
 writes it to a file. Run at one rank count, then another, and diff.

     julia --project=<env> test/hevi/test_rank_independence.jl          # writes
     mpiexecjl -n 3 julia --project=<env> test/hevi/test_rank_independence.jl
     # -> compares against the reference the 1-rank run wrote

 The 2x2x5 mock mesh on 3 ranks cuts through columns, which is the case the
 redistribution exists for and the one a single-rank run cannot reach.
=============================================================================#

using Test, MPI, LinearAlgebra, Printf, OrdinaryDiffEq
using Logging: with_logger, NullLogger

MPI.Initialized() || MPI.Init()
const COMM = MPI.COMM_WORLD
const RANK = MPI.Comm_rank(COMM)
const NR   = MPI.Comm_size(COMM)

const SRC = joinpath(@__DIR__, "..", "..", "src", "kernel", "solvers", "hevi")
include(joinpath(@__DIR__, "mock_sem.jl"))
include(joinpath(SRC, "columns.jl"))
include(joinpath(SRC, "vdiffusion.jl"))
include(joinpath(SRC, "operator.jl"))
include(joinpath(SRC, "factorize.jl"))
include(joinpath(SRC, "krylov.jl"))
include(joinpath(SRC, "acoustic.jl"))
include(joinpath(SRC, "ark.jl"))
include(joinpath(SRC, "hevi.jl"))
include(joinpath(SRC, "schur.jl"))
include(joinpath(SRC, "schur_precond.jl"))
include(joinpath(SRC, "schur_stage.jl"))
include(joinpath(SRC, "imex3d.jl"))

say(a...) = RANK == 0 && (println(a...); flush(stdout))

const NELX, NELY, NELZ, P = 2, 2, 5, 2
params = build_mock_params(nelx = NELX, nely = NELY, nelz = NELZ, p = P,
                           Lz = 1000.0, comm = COMM)
const NPOIN = Int(params.mesh.npoin)
const GIP   = params.mesh.ip2gip
const GN    = MPI.Allreduce(maximum(GIP), MPI.MAX, COMM)

# Gather a per-node field to rank 0, keyed by GLOBAL id, so the result is
# independent of the local numbering each rank happens to use.
function gather_global(X::AbstractMatrix)
    nq = size(X, 2)
    # -Inf, NOT zero. The reduction is a MAX over ranks and a rank contributes
    # this array for EVERY global id, including the ones it does not hold. Seeded
    # with zeros, a node whose true value is negative reduces to
    # max(negative, 0, 0) = 0 -- so on more than one rank the gathered field is
    # silently clamped at zero wherever it is negative, and these fields are
    # negative about half the time.
    #
    # That produced the perfect tell: every quantity differed from its 1-rank
    # reference by a relative 1.000e+00 exactly, which is what
    # max|0 - ref| / max|ref| is, not what a partition-dependent solver looks
    # like. A 1-rank run never sees it, because there is no second rank to
    # contribute the zero.
    G  = fill(-Inf, GN, nq)
    @inbounds for (ip, g) in enumerate(GIP), q = 1:nq
        G[g, q] = X[ip, q]          # duplicated nodes agree, or the test fails
    end
    MPI.Allreduce!(G, MPI.MAX, COMM)
    return G
end

# A deterministic field that is a function of the GLOBAL id alone, so every
# rank count starts from bit-identical data.
function global_field(nq)
    X = zeros(Float64, NPOIN, nq)
    @inbounds for ip = 1:NPOIN, q = 1:nq
        X[ip, q] = sinpi(1.0e-2 * (Float64(GIP[ip]) + 17.0q)) + 0.25cospi(3.0e-2 * GIP[ip])
    end
    return X
end

const REF = joinpath(@__DIR__, "rank_independence_reference.txt")

say("\n=== rank-independence: $NR rank(s), $(NELX)x$(NELY)x$(NELZ) elements, p=$P ===")

results = Dict{String,Matrix{Float64}}()

topo = build_column_topology(params.mesh, COMM)
vars = hevi_choose_vars(params.metrics, COMM)
op   = build_hevi_operator(params, topo, vars; lwall_flux = true)
owner, own = assign_column_owners(topo, COMM)
cc   = build_column_comm(topo, owner, own, COMM, length(vars))

@testset "rank independence ($NR ranks)" begin

    # 1. THE OPERATOR. f_imp is a pure function of the field; if the halo
    #    exchange or the DSS weighting is partition-dependent, this is where it
    #    shows, before any solve can average the error away.
    V = global_field(op.nimp)
    W = similar(V); fill!(W, 0.0)
    hevi_apply_A!(W, V, params, op)
    results["apply_A"] = gather_global(W)

    # 2. THE COLUMN SOLVE. Adds the gather/scatter and the banded factors, i.e.
    #    everything that moves data between ranks in a HEVI stage.
    #
    # THIS TESTSET HAD NEVER RUN. As first written it called
    #
    #     hevi_column_solve!(X, params, op, cc, fac)
    #
    # against a signature that has been `(dst, src, params, gdt::Real)` since
    # before this file existed -- five arguments where there are four, and the
    # operator, the column plan and the factors reached through `params.hevi`
    # rather than passed. It threw MethodError on every invocation, so the one
    # test in the repository whose whole purpose is "does the answer depend on
    # the rank count" has never answered it. Julia resolves that at CALL time,
    # which is why nothing static caught it.
    #
    # It also wants the FLAT layout and subtracts qe itself, so it is handed
    # `qe + deviation` rather than the deviation.
    gdt  = 0.05
    fac  = build_column_factorization(params, op, cc, topo, gdt)
    hev  = HEVICache(topo, op, cc, fac, hevi_fimp!, hevi_column_solve!, :ARS343,
                     gdt, :RS, 5)
    ph   = merge(params, (hevi = hev,))
    neqs = Int(params.neqs)
    qe   = params.qp.qe
    Bm   = global_field(neqs)
    srcv = zeros(Float64, neqs * NPOIN)
    @inbounds for ieq = 1:neqs, ip = 1:NPOIN
        srcv[(ieq-1)*NPOIN + ip] = qe[ip, ieq] + Bm[ip, ieq]
    end
    dstv = similar(srcv)
    hevi_column_solve!(dstv, srcv, ph, gdt)
    # back to the operator's packed deviation, which is what gather_global and
    # the single-valued check both expect.
    X = zeros(Float64, NPOIN, op.nimp)
    @inbounds for (q, ieq) in enumerate(op.vars), ip = 1:NPOIN
        X[ip, q] = dstv[(ieq-1)*NPOIN + ip] - qe[ip, ieq]
    end
    results["column_solve"] = gather_global(X)

    # 3. Duplicated nodes must already agree ACROSS ranks before the gather.
    #    gather_global takes a MAX, so a disagreement would be silently hidden
    #    in 1 and 2 above; this measures it directly.
    sv = imex_single_valued_error(X, op, COMM)
    @test sv < 1.0e-12
    say(@sprintf("  single-valued across ranks: %.2e", sv))

    # 4. THE SCALAR SCHUR STAGE SOLVE, which is the reason this file was
    #    revisited. It is the path with the most ways to be partition-dependent
    #    and the least prior coverage of it:
    #
    #      * schur_setup_rhs! caches q_b and m_b PER NODE, and schur_recover!
    #        reads them back AFTER a distributed Krylov solve -- so anything
    #        that makes those caches disagree at a rank-shared node produces a
    #        wrong state that no residual check would see, because the residual
    #        is measured on the scalar system the cache is not part of;
    #      * its preconditioner gathers each column onto one owner, and on this
    #        mesh at 3 ranks columns ARE split, which a 1-rank run cannot reach.
    #
    #    test_schur_stage.jl compares Schur against the five-field solve at 1
    #    and at 3 ranks, which is a strong check and still not this one: a fault
    #    that moved BOTH answers together would pass it. This compares the Schur
    #    answer against the SAME Schur answer computed on a different partition.
    opf   = build_hevi_fast_operator(params, topo; lwall_flux = true,
                                     theta_advective = true)
    inner = build_distributed_inner(opf.dss_cache, NPOIN, COMM)
    sinp  = Dict{Symbol,Any}(:imex_rtol => 1.0e-12, :imex_restart => 40,
                             :imex_maxiter => 400, :imex_precond => :column)
    sch   = build_imex3d_schur(params, topo, COMM, sinp, opf, gdt, inner)
    sws   = GMRESWorkspace(NPOIN, opf.nimp, inner; m = 40, maxiter = 400,
                           rtol = 1.0e-12, atol = 1.0e-30)
    sic   = IMEX3DCache(topo, opf, nothing, sch, sws,
                        zeros(Float64, NPOIN, opf.nimp),
                        zeros(Float64, NPOIN, opf.nimp),
                        imex3d_fimp!, imex3d_solve_schur!, :ARS343, gdt,
                        :RS, 1, false, COMM)
    ps    = merge(params, (imex = sic,))
    sdst  = similar(srcv)
    imex3d_solve_schur!(sdst, srcv, ps, gdt)
    XS = zeros(Float64, NPOIN, opf.nimp)
    @inbounds for (q, ieq) in enumerate(opf.vars), ip = 1:NPOIN
        XS[ip, q] = sdst[(ieq-1)*NPOIN + ip] - qe[ip, ieq]
    end
    results["schur_solve"] = gather_global(XS)
    svs = imex_single_valued_error(XS, opf, COMM)
    @test svs < 1.0e-12
    say(@sprintf("  schur solve single-valued across ranks: %.2e", svs))

    if RANK == 0
        if NR == 1
            open(REF, "w") do io
                for (k, G) in sort(collect(results); by = first)
                    println(io, "# ", k, " ", size(G, 1), " ", size(G, 2))
                    for i = 1:size(G, 1), q = 1:size(G, 2)
                        @printf(io, "%d %d %.17e\n", i, q, G[i, q])
                    end
                end
            end
            println("  wrote 1-rank reference -> ", REF)
            println("  now run: mpiexecjl -n 3 julia --project=<env> ", basename(@__FILE__))
        elseif isfile(REF)
            ref = Dict{String,Matrix{Float64}}(); key = ""; rows = 0; cols = 0
            for l in eachline(REF)
                if startswith(l, "#")
                    p = split(l); key = p[2]; rows = parse(Int, p[3]); cols = parse(Int, p[4])
                    ref[key] = zeros(Float64, rows, cols)
                else
                    p = split(l)
                    ref[key][parse(Int, p[1]), parse(Int, p[2])] = parse(Float64, p[3])
                end
            end
            for (k, G) in sort(collect(results); by = first)
                haskey(ref, k) || continue
                d = maximum(abs, G .- ref[k])
                s = max(maximum(abs, ref[k]), eps())
                @printf("  %-14s max|%d rank - 1 rank| = %.3e  (relative %.3e)\n",
                        k, NR, d, d/s)
                # Round-off only. Assembly order differs between partitions, so
                # bit-identity is not the bar; anything above this is a real
                # partition dependence.
                @test d/s < 1.0e-12
            end
        else
            println("  no reference file -- run at 1 rank first")
        end
    end
end

MPI.Barrier(COMM)
say("=== rank-independence checks done on $NR rank(s) ===\n")
