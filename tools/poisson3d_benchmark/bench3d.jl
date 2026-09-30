#=============================================================================
 tools/poisson3d_benchmark/bench3d.jl — run configurations of the 3D (or 2D)
 periodic Poisson solver comparison (poisson3d.jl) and write results.csv.

 One Julia session per call. Every configuration (solver, ne, N) follows the
 protocol of the 2D benchmark: a first run that is discarded (compilation,
 first-call effects), the garbage collector, then the recorded SECOND run.
 With --warmup small the discarded run uses 2^d elements (same solver and N)
 so a large configuration is not solved twice; with --warmup same it is the
 identical configuration.

 USAGE
   one configuration (what each SLURM job runs):
     julia --project=. -t 8 tools/poisson3d_benchmark/bench3d.jl \
           --solver sem_amg --ne 32 --nop 4 --outdir ppb3d/parts/sem_amg_ne32_N4
   a local sweep (one session, rows appended, resumable):
     julia --project=. -t 4 tools/poisson3d_benchmark/bench3d.jl \
           --solvers sem,sem_amg,sc_amg,fft --nes 4,8,16 --nops 2,4 --outdir ppb3d_local --resume

 OPTIONS
   --solver S / --solvers S1,S2     sem sem_amg sem_jacobi sc_direct sc_amg ps fft
   --ne n / --nes n1,n2             elements per direction
   --nop N / --nops N1,N2           SEM order (Fourier grids: ne*N points per direction)
   --d 3                            dimension (2 or 3)
   --r 0.5                          Fourier decay rate of the exact solution
   --ordering metis|amd             fill-reducing ordering of the Cholesky factorisations
   --rtol 1e-12                     CG tolerance (preconditioned residual)
   --warmup small|same
   --blas-threads n                 BLAS threads (CHOLMOD supernodes, dense kernels);
                                    Julia threads (element loops) are set with `julia -t`
   --outdir DIR   --resume

 COLUMNS of results.csv: the row of P3D.run_config (errors; assembly, rhs,
 setup, solve and total seconds; CG iterations; nnz of K, of the Cholesky
 factor, of the skeleton matrix) plus julia_threads, blas_threads, the peak
 resident memory of the process (maxrss_gb: per configuration when each runs
 in its own process, as under SLURM; the running maximum in a local sweep),
 and status ("ok", or the error of a configuration that failed).
=============================================================================#
using Jexpresso, Printf, LinearAlgebra
include(joinpath(@__DIR__, "poisson3d.jl"))
using .P3D

const COLS = (:solver, :d, :r, :ne, :nop, :Ng, :n, :solved, :linf, :l2rel, :assembly, :rhs, :setup,
              :solve, :total, :iters, :nnz, :factor_nnz, :skeleton_nnz, :ordering,
              :julia_threads, :blas_threads, :maxrss_gb, :status)

function parse_args(args)
    o = Dict{Symbol, Any}(:d => 3, :r => nothing, :ordering => :metis, :rtol => 1e-12, :warmup => :small,
                          :outdir => "ppb3d_results", :resume => false)
    ints(s) = parse.(Int, split(s, ','))
    i = 1
    nxt() = (i += 1; args[i])
    while i <= length(args)
        a = args[i]
        if a in ("--solver", "--solvers"); o[:solvers] = Symbol.(split(nxt(), ','))
        elseif a in ("--ne", "--nes");     o[:nes] = ints(nxt())
        elseif a in ("--nop", "--nops");   o[:nops] = ints(nxt())
        elseif a == "--d";                 o[:d] = parse(Int, nxt())
        elseif a == "--r";                 o[:r] = parse(Float64, nxt())
        elseif a == "--ordering";          o[:ordering] = Symbol(nxt())
        elseif a == "--rtol";              o[:rtol] = parse(Float64, nxt())
        elseif a == "--warmup";            o[:warmup] = Symbol(nxt())
        elseif a == "--blas-threads";      BLAS.set_num_threads(parse(Int, nxt()))
        elseif a == "--outdir";            o[:outdir] = nxt()
        elseif a == "--resume";            o[:resume] = true
        else error("unknown option $a (see the header of bench3d.jl)")
        end
        i += 1
    end
    haskey(o, :solvers) && haskey(o, :nes) && haskey(o, :nops) || error("need --solver(s), --ne(s), --nop(s)")
    o[:r] === nothing && (o[:r] = P3D.R_DEFAULT[o[:d]])
    return o
end

function write_csv(path, rows)
    tmp = path * ".tmp"
    open(tmp, "w") do io
        println(io, join(COLS, ","))
        for r in rows
            println(io, join((replace(string(r[c]), "," => ";") for c in COLS), ","))
        end
    end
    mv(tmp, path; force = true)
end

function main(args)
    o = parse_args(args)
    mkpath(o[:outdir]); csv = joinpath(o[:outdir], "results.csv")
    rows = Dict{Symbol, Any}[]
    if o[:resume] && isfile(csv)                      # keep what is there
        lines = readlines(csv); hdr = Symbol.(split(lines[1], ","))
        for l in lines[2:end]
            push!(rows, Dict{Symbol, Any}(zip(hdr, split(l, ","))))
        end
    end
    done = Set((string(r[:solver]), parse(Int, string(r[:ne])), parse(Int, string(r[:nop]))) for r in rows)
    d = o[:d]
    kw = (; r = o[:r], rtol = o[:rtol], ordering = o[:ordering])
    for ne in o[:nes], N in o[:nops], s in o[:solvers]
        (string(s), ne, N) in done && continue
        n = (ne * N)^d
        @printf("P3D  run  %-10s d=%d ne=%d N=%d n=%d (threads: julia %d, blas %d) ...\n",
                s, d, ne, N, n, Threads.nthreads(), BLAS.get_num_threads())
        row = Dict{Symbol, Any}(c => "" for c in COLS)
        try
            run_config(s, d, o[:warmup] === :same ? ne : min(ne, 2), N; kw...)   # warm-up, discarded
            GC.gc(); GC.gc()
            res = run_config(s, d, ne, N; kw...)                                    # recorded
            for (k, v) in pairs(res); row[k] = v; end
            row[:status] = "ok"
        catch e
            e isa InterruptException && rethrow()
            msg = sprint(showerror, e)
            @warn "P3D: $s ne=$ne N=$N failed" exception = (e, catch_backtrace())
            merge!(row, Dict(:solver => s, :d => d, :r => o[:r], :ne => ne, :nop => N, :Ng => ne * N, :n => n,
                             :status => "error: " * first(replace(msg, '\n' => ' ', ',' => ';'), 200)))
        end
        row[:julia_threads] = Threads.nthreads(); row[:blas_threads] = BLAS.get_num_threads()
        row[:maxrss_gb] = round(Sys.maxrss() / 2^30, digits = 3)
        push!(rows, row)
        write_csv(csv, rows)
        @printf("P3D  %-10s ne=%d N=%d n=%d  L∞=%s  setup=%s s  solve=%s s  total=%s s  its=%s  maxrss=%.2f GB  %s\n",
                s, ne, N, n, string(row[:linf]), string(row[:setup]), string(row[:solve]), string(row[:total]),
                string(row[:iters]), row[:maxrss_gb], row[:status])
    end
    println("wrote ", csv)
end

main(copy(ARGS))
