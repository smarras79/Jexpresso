module JExpressoRunTests

project_root = dirname(Base.current_project())
include(joinpath(project_root, "test", "solnCompare.jl"))

using .solnCompare
using Test


@testset "JEXPRESSO Examples: CompEuler, $alg_case" for alg_case in ("theta", "thetaTracers")
    run_example("CompEuler", "$alg_case")
end

# IMEX/JACC portable linear-algebra kernels (CPU default backend). Self-
# contained: it only needs JACC + SparseArrays + LinearAlgebra and does not
# build the full Jexpresso module or require a mesh, so it always runs in CI.
# The GPU equivalent lives in test/test_imex_jacc_gpu.jl (run on a GPU box).
include(joinpath(project_root, "test", "test_imex_jacc.jl"))

end # module
