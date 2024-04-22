module CompEulerthetaTests

project_root = dirname(Base.current_project())
include(joinpath(project_root, "test", "solnCompare.jl"))

using .solnCompare
using Test

@testset "JEXPRESSO Examples" begin run_example("CompEuler", "wave1d") end
end # module
