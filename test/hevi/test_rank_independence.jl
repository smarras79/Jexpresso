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
    G  = zeros(Float64, GN, nq)
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
    gdt = 0.05
    fac = build_column_factorization(params, op, cc, topo, gdt)
    B = global_field(op.nimp)
    X = copy(B)
    hevi_column_solve!(X, params, op, cc, fac)
    results["column_solve"] = gather_global(X)

    # 3. Duplicated nodes must already agree ACROSS ranks before the gather.
    #    gather_global takes a MAX, so a disagreement would be silently hidden
    #    in 1 and 2 above; this measures it directly.
    sv = imex_single_valued_error(X, op, COMM)
    @test sv < 1.0e-12
    say(@sprintf("  single-valued across ranks: %.2e", sv))

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
