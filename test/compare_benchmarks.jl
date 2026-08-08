#==============================================================================
 test/compare_benchmarks.jl

 Compares the HDF5 output of one or more CI cases against the reference
 solutions committed in test/CI-ref/. Used by the `benchmarks` workflow,
 which runs the simulation in one step and this script in the next, so that
 the comparison does not have to load the whole Jexpresso package (only
 HDF5.jl is needed here).

     julia --project=. test/compare_benchmarks.jl                  # all cases
     julia --project=. test/compare_benchmarks.jl CompEuler/theta  # one case

 The list of known cases and their tolerances lives in test/ci_cases.jl.
 Exit status is non-zero when any comparison fails.
==============================================================================#
module BenchmarkComparison

using Test

include(joinpath(@__DIR__, "ci_compare.jl"))
using .CICompare
using .CICompare.CICases

const SELECTION = isempty(ARGS) ? "all" : join(ARGS, ",")
const CASES     = select_cases(SELECTION)

isempty(CASES) && @warn "no CI cases selected — nothing to compare"

@testset "Benchmark comparisons" begin
    compare_cases(CASES)
end

end # module BenchmarkComparison
