module CompEulerthetaTests

include(joinpath(pwd(), "test", "solnCompare.jl"))

using .solnCompare
using Test

@testset "JEXPRESSO Examples" begin run_example("CompEuler", "theta") end
end # module
