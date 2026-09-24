#=============================================================================
 tools/periodic_poisson_benchmark/pipeline.jl

 Benchmark of every solver of the doubly periodic Poisson problem
 Elliptic/poisson_periodic_sem (-∇²u = f on [0,2π]², exact u known):

     :sem   direct SEM solve on 16×16 spectral elements of order N (:nop)
     :ps    pseudo-spectral (Fourier collocation, Kopriva) on a 16N × 16N grid
     :fft   FFT (FFTW) on a 16N × 16N grid

 i.e. at the SAME number of unknowns, (16N)², for N = 2..8.

 TIMING PROTOCOL. Every configuration is run TWICE in the same Julia session
 with run_case(...; inputs = overrides) and only the SECOND run is recorded,
 so compilation and first-call effects never enter a number. Before the
 recorded run the garbage collector is run, so the first run's garbage is
 not collected inside the second run's timers. Each phase is a single
 wall-clock measurement (time_ns) of that second run — no repetition, no
 BenchmarkTools minimum (:lbenchmark_solve => false). The mesh / SEM
 preprocess caches are switched off (:luse_mesh_cache => false), so the
 second run BUILDS its infrastructure instead of loading it from disk, and
 output files are switched off (:outformat => "none").

 WHAT IS MEASURED (per-phase timers recorded by the driver and the solvers
 in Jexpresso.JX_TIMINGS):
   solve     the solve step alone: triangular solves (SEM), four dense
             N×N products (pseudo-spectral), rfft/scale/brfft (FFT)
   setup     the solver's own infrastructure: periodic reduction + sparse Cholesky
             factorisation (SEM), 1-D eigen-decompositions (pseudo-spectral),
             FFTW plan (FFT)
   rhs       sampling / assembling the right-hand side
   sem_setup mesh read + SEM infrastructure (basis, metrics, mass and
             Laplacian assembly) — needed by the SEM only
   total     time-to-solution including all the infrastructure the METHOD
             needs:  SEM  = sem_setup + rhs + setup + solve
                     ps / fft = rhs + setup + solve
             (the Jexpresso driver still runs sem_setup before dispatching to
             a spectral solver; the spectral methods do not use it, so it is
             not charged to them — it is in `wall`)
   wall      the whole second run_case call, for reference

 USAGE — from the REPL (the intended way):
     julia --project=.
     julia> using Jexpresso
     julia> include("tools/periodic_poisson_benchmark/pipeline.jl")
     julia> rows = run_periodic_poisson_benchmark()            # N = 2..8
 or as a script (also one Julia session, same protocol):
     julia --project=. tools/periodic_poisson_benchmark/pipeline.jl

 OUTPUT (in this directory): results.csv, results.md (the table), and the
 figures in assets/ drawn by plot.py (error vs DOFs, vs order, vs solve time,
 vs total time-to-solution), when python3 is available.
=============================================================================#
using Jexpresso, Printf

const PPB_DIR    = @__DIR__
const PPB_EQS    = "Elliptic"
const PPB_CASE   = "poisson_periodic_sem"
const PPB_LABELS = Dict(:sem => "SEM", :ps => "pseudo-spectral", :fft => "FFT")

function _ppb_overrides(solver::Symbol, nop::Int, nel::Int)
    return Dict{Symbol, Any}(
        :nop              => nop,
        :fft_N            => nel * nop,          # the pseudo-spectral grid defaults to it
        :lfft             => solver === :fft,
        :lpseudospectral  => solver === :ps,
        :luse_mesh_cache  => false,              # build, never load, the infrastructure
        :lbenchmark_solve => false,              # single-shot solve timer
        :outformat        => "none",
    )
end

function _ppb_run(ov)
    wall = @elapsed Jexpresso.run_case(PPB_EQS, PPB_CASE; inputs = ov)
    return (wall = wall, t = copy(Jexpresso.JX_TIMINGS), err = Jexpresso.JX_LAST_SOLVE_ERR[])
end

"""
    run_periodic_poisson_benchmark(; nops = 2:8, nel = 16, solvers = (:sem, :ps, :fft),
                                   outdir = <this directory>, plot = true) -> rows

Run the benchmark (see the file header), write results.csv / results.md into
`outdir`, draw the figures, and return the rows.
"""
function run_periodic_poisson_benchmark(; nops = 2:8, nel::Int = 16,
                                        solvers = (:sem, :ps, :fft),
                                        outdir::AbstractString = PPB_DIR,
                                        plot::Bool = true)
    rows = NamedTuple[]
    for nop in nops, solver in solvers
        ov = _ppb_overrides(solver, nop, nel)
        _ppb_run(ov)                             # 1st run: warm-up, discarded
        GC.gc()
        r = _ppb_run(ov)                         # 2nd run: recorded
        g(k) = get(r.t, k, 0.0)
        method_total = g(:rhs) + g(:setup) + g(:solve) + (solver === :sem ? g(:sem_setup) : 0.0)
        row = (solver = solver, nop = nop, N = nel * nop, dofs = (nel * nop)^2,
               linf = r.err.linf, l2rel = r.err.l2rel,
               solve = g(:solve), setup = g(:setup), rhs = g(:rhs),
               sem_setup = solver === :sem ? g(:sem_setup) : 0.0,
               reduce = g(:reduce), factorize = g(:factorize),
               total = method_total, wall = r.wall)
        push!(rows, row)
        @printf("PPB  %-15s N=%d  dofs=%6d  L∞=%.2e  solve=%.3e s  total=%.3e s  wall=%.3e s\n",
                PPB_LABELS[solver], nop, row.dofs, row.linf, row.solve, row.total, row.wall)
    end
    mkpath(outdir)
    _ppb_write_csv(joinpath(outdir, "results.csv"), rows)
    _ppb_write_md(joinpath(outdir, "results.md"), rows, nel)
    if plot
        try
            run(`python3 $(joinpath(PPB_DIR, "plot.py")) $(joinpath(outdir, "results.csv"))`)
        catch e
            @warn "plot.py failed; the CSV and the table are written" exception = e
        end
    end
    return rows
end

const _PPB_COLS = (:solver, :nop, :N, :dofs, :linf, :l2rel, :solve, :setup, :rhs,
                   :sem_setup, :reduce, :factorize, :total, :wall)

function _ppb_write_csv(path, rows)
    open(path, "w") do io
        println(io, join(_PPB_COLS, ","))
        for r in rows
            println(io, join((getproperty(r, c) for c in _PPB_COLS), ","))
        end
    end
    println("wrote ", path)
end

_ppb_t(x) = x == 0 ? "—" : x < 1e-3 ? @sprintf("%.3g µs", 1e6x) : x < 1 ? @sprintf("%.3g ms", 1e3x) : @sprintf("%.3g s", x)
_ppb_e(x) = @sprintf("%.1e", x)
_ppb_n(x) = replace(string(x), r"(\d)(?=(\d{3})+$)" => s"\1 ")

function _ppb_write_md(path, rows, nel)
    open(path, "w") do io
        println(io, "| method | SEM order N | unknowns (grid) | ‖e‖∞ | relative ‖e‖₂ | solve | setup | RHS | SEM infrastructure | time-to-solution | run_case wall-clock |")
        println(io, "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for r in sort(rows; by = r -> (findfirst(==(r.solver), (:sem, :ps, :fft)), r.nop))
            println(io, "| ", PPB_LABELS[r.solver], " | ", r.solver === :sem ? string(r.nop) : "—", " | ", _ppb_n(r.dofs),
                    " (", r.N, "²) | ", _ppb_e(r.linf), " | ", _ppb_e(r.l2rel), " | ",
                    _ppb_t(r.solve), " | ", _ppb_t(r.setup), " | ", _ppb_t(r.rhs), " | ",
                    _ppb_t(r.sem_setup), " | **", _ppb_t(r.total), "** | ", _ppb_t(r.wall), " |")
        end
    end
    println("wrote ", path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_periodic_poisson_benchmark()
end
