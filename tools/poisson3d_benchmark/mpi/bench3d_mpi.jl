#=============================================================================
 tools/poisson3d_benchmark/mpi/bench3d_mpi.jl — run configurations of the
 MPI 3D periodic Poisson solver comparison (poisson3d_mpi.jl) on all ranks
 and write results.csv (rank 0).

 Same protocol as ../bench3d.jl: per configuration a discarded warm-up run
 on a small grid (same solver and N), the garbage collector, then the
 recorded run. Timings are those of the slowest rank (every phase starts and
 ends at a barrier).

 USAGE
   mpiexec -n 128 julia --project=tools/poisson3d_benchmark/mpi \
       tools/poisson3d_benchmark/mpi/bench3d_mpi.jl \
       --solver boomeramg --nop 4 --nes 4,8,16,32 --outdir ppb3d_mpi/parts/boomeramg --resume

 OPTIONS
   --solver S / --solvers S1,S2     mumps boomeramg jacobi
   --ne n / --nes n1,n2             elements per direction (unknowns n = (ne*N)^3)
   --nop N / --nops N1,N2           SEM order
   --r 0.5                          Fourier decay rate of the exact solution
   --rtol 1e-12                     CG tolerance (relative, preconditioned residual)
   --ordering metis|auto            MUMPS fill-reducing ordering (sequential METIS
                                    nested dissection, as the serial benchmark)
   --theta 0.5                      BoomerAMG strong threshold
   --outdir DIR   --resume          (--resume: skip configurations with an "ok" row)

 COLUMNS: those of ../bench3d.jl (rhs is part of assembly here; maxrss_gb is
 the sum over ranks of each rank's peak memory, the running maximum of the
 session), plus nranks, the process grid, MUMPS's own memory figure (INFOG(22),
 all ranks) and the largest per-rank peak memory.
=============================================================================#
using MPI
MPI.Init()
using HYPRE
HYPRE.Init()
using Printf
include(joinpath(@__DIR__, "poisson3d_mpi.jl"))
using .P3DMPI

const COLS = (:solver, :d, :r, :ne, :nop, :Ng, :n, :solved, :linf, :l2rel, :assembly, :rhs, :setup,
              :solve, :total, :iters, :nnz, :factor_nnz, :skeleton_nnz, :ordering,
              :julia_threads, :blas_threads, :nranks, :grid, :mumps_mem_gb, :maxrss_rank_gb,
              :maxrss_gb, :status)

function parse_args(args)
    o = Dict{Symbol, Any}(:r => P3DMPI.R_DEFAULT, :ordering => :metis, :rtol => 1e-12, :theta => 0.5,
                          :outdir => "ppb3d_mpi", :resume => false)
    ints(s) = parse.(Int, split(s, ','))
    i = 1
    nxt() = (i += 1; args[i])
    while i <= length(args)
        a = args[i]
        if a in ("--solver", "--solvers"); o[:solvers] = Symbol.(split(nxt(), ','))
        elseif a in ("--ne", "--nes");     o[:nes] = ints(nxt())
        elseif a in ("--nop", "--nops");   o[:nops] = ints(nxt())
        elseif a == "--r";                 o[:r] = parse(Float64, nxt())
        elseif a == "--ordering";          o[:ordering] = Symbol(nxt())
        elseif a == "--rtol";              o[:rtol] = parse(Float64, nxt())
        elseif a == "--theta";             o[:theta] = parse(Float64, nxt())
        elseif a == "--outdir";            o[:outdir] = nxt()
        elseif a == "--resume";            o[:resume] = true
        else error("unknown option $a (see the header of bench3d_mpi.jl)")
        end
        i += 1
    end
    haskey(o, :solvers) && haskey(o, :nes) && haskey(o, :nops) || error("need --solver(s), --ne(s), --nop(s)")
    return o
end

function write_csv(path, rows)
    tmp = path * ".tmp"
    open(tmp, "w") do io
        println(io, join(COLS, ","))
        for r in rows
            println(io, join((replace(string(get(r, c, "")), "," => ";") for c in COLS), ","))
        end
    end
    mv(tmp, path; force = true)
end

function main(args)
    comm = MPI.COMM_WORLD; me = MPI.Comm_rank(comm); np = MPI.Comm_size(comm)
    o = parse_args(args)
    csv = joinpath(o[:outdir], "results.csv")
    rows = Dict{Symbol, Any}[]
    if me == 0
        mkpath(o[:outdir])
        if o[:resume] && isfile(csv)                  # keep what finished; failed rows are rerun
            lines = readlines(csv); hdr = Symbol.(split(lines[1], ","))
            for l in lines[2:end]
                r = Dict{Symbol, Any}(zip(hdr, split(l, ",")))
                get(r, :status, "") == "ok" && push!(rows, r)
            end
        end
    end
    done = MPI.bcast(Set((string(r[:solver]), parse(Int, string(r[:ne])), parse(Int, string(r[:nop]))) for r in rows), 0, comm)
    pmax = maximum(MPI.Dims_create(np, [0, 0, 0]))   # every direction needs >= pmax nodes
    kw = (; r = o[:r], rtol = o[:rtol], ordering = o[:ordering], theta = o[:theta])
    for ne in o[:nes], N in o[:nops], s in o[:solvers]
        (string(s), ne, N) in done && continue
        n = (ne * N)^3
        me == 0 && @printf("P3DMPI  run  %-9s ne=%d N=%d n=%d on %d ranks ...\n", s, ne, N, n, np)
        row = Dict{Symbol, Any}(c => "" for c in COLS)
        try
            ne * N >= pmax || error("$np ranks need at least $pmax nodes per direction; ne*N = $(ne * N)")
            nw = max(2, cld(pmax, N))                 # warm-up grid, discarded
            run_config(s, min(ne, nw), N; kw...)
            GC.gc(); GC.gc(); MPI.Barrier(comm)
            res = run_config(s, ne, N; kw...)
            for (k, v) in pairs(res); row[k] = v; end
            row[:status] = "ok"
        catch e
            e isa InterruptException && rethrow()
            msg = sprint(showerror, e)
            me == 0 && @warn "P3DMPI: $s ne=$ne N=$N failed" exception = (e, catch_backtrace())
            merge!(row, Dict(:solver => s, :d => 3, :r => o[:r], :ne => ne, :nop => N, :Ng => ne * N, :n => n,
                             :nranks => np, :status => "error: " * first(replace(msg, '\n' => ' ', ',' => ';'), 200)))
        end
        row[:julia_threads] = 1; row[:blas_threads] = 1
        if me == 0
            push!(rows, row)
            write_csv(csv, rows)
            @printf("P3DMPI  %-9s ne=%d N=%d n=%d  L∞=%s  assembly=%s s  setup=%s s  solve=%s s  its=%s  mem=%s GB  %s\n",
                    s, ne, N, n, string(row[:linf]), string(row[:assembly]), string(row[:setup]), string(row[:solve]),
                    string(row[:iters]), string(row[:maxrss_gb]), row[:status])
        end
    end
    me == 0 && println("wrote ", csv)
end

main(copy(ARGS))
