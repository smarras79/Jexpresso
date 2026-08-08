#==============================================================================
 test/runtests.jl — entry point of `Pkg.test()` and of the CI workflow.

 Every case listed in test/ci_cases.jl is run with CI_MODE=true and its HDF5
 output is compared against the reference solution committed in
 test/CI-ref/. Adding a test to this suite means adding one line to
 test/ci_cases.jl — this file never changes.

 Run the whole suite locally with

     julia --project=. -e 'using Pkg; Pkg.test()'

 or a single case with

     julia --project=. test/runtests.jl CompEuler/theta
==============================================================================#
module JexpressoRunTests

using Test

include(joinpath(@__DIR__, "ci_compare.jl"))
using .CICompare
using .CICompare.CICases

const SELECTION = isempty(ARGS) ? "all" : join(ARGS, ",")
const CASES     = select_cases(SELECTION)

problems = validate()
isempty(problems) || error("test/ci_cases.jl is inconsistent with the " *
                           "repository:\n  " * join(problems, "\n  "))

@testset "Jexpresso CI suite" begin
    if isempty(CASES)
        @warn "no CI cases selected — nothing to test"
    end
    for c in CASES
        @testset "$(case_name(c))" begin
            run_ci_case(c)
            compare_case(c)
        end
    end
end

end # module JexpressoRunTests
