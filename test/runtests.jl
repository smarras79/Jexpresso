module JExpressoRunTests

using Test

@time @testset "CompEuler" begin include("CI-runs/CompEuler/runtests.jl") end

end # module
