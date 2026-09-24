#=============================================================================
 tools/sem_pconvergence/sweep.jl — SEM vs pseudo-spectral vs FFT on the
 doubly periodic Poisson problem Elliptic/poisson_periodic_sem

 For polynomial orders N = 2..8 on the deck's 16×16-element mesh, runs
   • the direct SEM solve                      ((16N)² periodic unknowns)
   • the pseudo-spectral (Fourier collocation) solve  on a 16N × 16N grid
   • the FFT solve                                    on a 16N × 16N grid
 i.e. the three at the SAME number of unknowns, and writes their L∞ /
 relative L² errors against the exact solution and their solve times:

     julia --project=. tools/sem_pconvergence/sweep.jl [out.csv]

 then draw the README figures with
     python3 tools/sem_pconvergence/plot.py

 Times: `solve_s` is each driver's SOLVER TIMING (BenchmarkTools minimum):
   SEM              sparse LU factorisation + solve of the periodic system
   pseudo-spectral  the solve only (four dense N×N products)
   FFT              the solve only (rfft, scale, brfft)
 `setup_s` is the one-time setup the two spectral solves exclude (the 1-D
 eigen-decompositions, the FFTW plan), measured here separately.
=============================================================================#
using Jexpresso

const OUT  = isempty(ARGS) ? joinpath(@__DIR__, "periodic_poisson_solvers.csv") : ARGS[1]
const DECK = joinpath(dirname(dirname(@__DIR__)), "problems", "Elliptic", "poisson_periodic_sem", "user_inputs.jl")
const NOPS = 2:8
const NEL  = 16                 # elements per direction in the deck's mesh

const NOP_LINE = ":nop                  => 4,"
const LFFT     = ":lfft                 => false,"
const LPS      = ":lpseudospectral      => false,"
const FFTN     = ":fft_N                => 64,"     # the spectral grid size (the pseudo-spectral one defaults to it)

function deck(orig; nop = 4, solver = :sem, N = 64)
    # Edit the deck's own lines: a key repeated in the Dict literal would be
    # overridden by the deck's later entry, silently.
    s = replace(orig, NOP_LINE => ":nop                  => $nop,", FFTN => ":fft_N                => $N,")
    solver === :fft && (s = replace(s, LFFT => ":lfft                 => true,"))
    solver === :ps  && (s = replace(s, LPS  => ":lpseudospectral      => true,"))
    return s
end

best(f; reps = 3) = minimum(@elapsed(f()) for _ in 1:reps)

orig = read(DECK, String)
for line in (NOP_LINE, LFFT, LPS, FFTN)
    occursin(line, orig) || error("sweep.jl: could not find `$line` in $DECK")
end
rows = ["solver,nop,unknowns,linf,l2rel,solve_s,setup_s"]
try
    for n in NOPS
        N = NEL * n
        for solver in (:sem, :ps, :fft)
            write(DECK, deck(orig; nop = n, solver = solver, N = N))
            Jexpresso.run_case("Elliptic", "poisson_periodic_sem")   # (overwrites ARGS)
            e = Jexpresso.JX_LAST_SOLVE_ERR[]
            t = Jexpresso.JX_LAST_SOLVE_TIME[]
            setup = solver === :ps  ? best(() -> Jexpresso.FourierCollocationPoissonSolver((N, N), (2π, 2π))) :
                    solver === :fft ? best(() -> Jexpresso.FFTPoissonSolver((N, N), (2π, 2π); flags = Jexpresso.FFTW.MEASURE)) :
                    NaN
            push!(rows, "$solver,$n,$(N^2),$(e.linf),$(e.l2rel),$t,$setup")
            println("SWEEP ", rows[end])
        end
    end
finally
    write(DECK, orig)                                              # restore the deck
end
write(OUT, join(rows, "\n") * "\n")
println("wrote ", OUT)
