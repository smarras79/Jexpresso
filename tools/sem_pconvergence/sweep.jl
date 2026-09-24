#=============================================================================
 tools/sem_pconvergence/sweep.jl — p-convergence of the SEM Poisson solve

 Runs Elliptic/poisson_periodic_sem (-∇²u = f on [0,2π]², doubly periodic,
 16×16 elements, direct SEM solve) for polynomial orders 2..8 and writes the
 L∞ and relative L² errors against the exact solution to a CSV:

     julia --project=. tools/sem_pconvergence/sweep.jl [out.csv]

 then draw the README figure with
     python3 tools/sem_pconvergence/plot.py
=============================================================================#
using Jexpresso

const OUT  = isempty(ARGS) ? joinpath(@__DIR__, "sem_periodic_poisson_pconv.csv") : ARGS[1]
const DECK = joinpath(dirname(dirname(@__DIR__)), "problems", "Elliptic", "poisson_periodic_sem", "user_inputs.jl")
const NOPS = 2:8
const NEL  = 16                 # elements per direction in the deck's mesh

orig = read(DECK, String)
occursin(":nop                  => 4,", orig) || error("sweep.jl: could not find the :nop line in $DECK")
rows = ["nop,dofs,linf,l2rel"]
try
    for n in NOPS
        write(DECK, replace(orig, ":nop                  => 4," => ":nop                  => $n,"))
        Jexpresso.run_case("Elliptic", "poisson_periodic_sem")   # (overwrites ARGS)
        e = Jexpresso.JX_LAST_SOLVE_ERR[]
        push!(rows, "$n,$((NEL*n)^2),$(e.linf),$(e.l2rel)")
        println("SWEEP ", rows[end])
    end
finally
    write(DECK, orig)                                              # restore the deck
end
write(OUT, join(rows, "\n") * "\n")
println("wrote ", OUT)
