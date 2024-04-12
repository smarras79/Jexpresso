module JExpressoRunTests

using Test

@time @testset "CompEuler" begin include("../problems/equations/CI-runs/CompEuler/runtests.jl") end

end # module
