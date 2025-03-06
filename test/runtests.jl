module JExpressoRunTests

project_root = dirname(Base.current_project())
include(joinpath(project_root, "test", "solnCompare.jl"))

using .solnCompare
using Test


@testset "JEXPRESSO Examples: CompEuler, $alg_case" for alg_case in ("theta", "thetaTracers")
    run_example("CompEuler", "$alg_case")
end

end # module
