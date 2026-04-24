# NOTE: Needs to copy-pasted into the REPL from the home directory of the repo

using LinearAlgebra
using Test


push!(empty!(ARGS), "CompEuler", "kelvinHelmholtzChan2022_test")

include("./src/Jexpresso.jl")

solution = Jexpresso.with_mpi() do distribute

    solution = Jexpresso.driver(Jexpresso.nparts,
           distribute,
           Jexpresso.inputs,
           Jexpresso.OUTPUT_DIR,
           Jexpresso.TFloat)
    return solution
end

using LinearAlgebra
l2norm = LinearAlgebra.norm(solution.u[2])
ref_l2norm = 396.70649542713255

@test isapprox(l2norm, ref_l2norm; atol=1e-10)
