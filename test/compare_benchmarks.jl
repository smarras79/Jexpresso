#==============================================================================
 test/compare_benchmarks.jl

 Compares the HDF5 output of one or more CI cases against the reference
 solutions committed in test/CI-ref/. Used by the `benchmarks` workflow,
 which runs the simulation in one step and this script in the next, so that
 the comparison does not have to load the whole Jexpresso package (only
 HDF5.jl is needed here).

     julia --project=. test/compare_benchmarks.jl                  # all cases
     julia --project=. test/compare_benchmarks.jl CompEuler/theta  # one case

 Options:
   --ref DIR   compare against the reference tree in DIR (laid out as
               <eqs>/<case>/output/) instead of test/CI-ref. Use it to check
               a reference produced somewhere else — e.g. the ciref-* artifact
               of a generate-ci-ref run on Linux, against your local output —
               without committing it first.

 The list of known cases and their tolerances lives in test/ci_cases.jl.
 Exit status is non-zero when any comparison fails.
==============================================================================#
module BenchmarkComparison

using Test

include(joinpath(@__DIR__, "ci_compare.jl"))
using .CICompare
using .CICompare.CICases

function parse_args(args::AbstractVector{<:AbstractString})
    cases    = String[]
    ref_root = joinpath(CICases.project_root(), "test", "CI-ref")

    i = 1
    while i <= length(args)
        arg = String(args[i])
        if arg == "--ref"
            i == length(args) && error("--ref needs a directory")
            ref_root = String(args[i + 1])
            i += 1
        elseif startswith(arg, "--ref=")
            ref_root = String(split(arg, '=', limit = 2)[2])
        elseif startswith(arg, "-")
            error("unknown option \"$arg\" — usage: " *
                  "julia --project=. test/compare_benchmarks.jl " *
                  "[--ref DIR] [<eqs>/<case> ...]")
        else
            push!(cases, arg)
        end
        i += 1
    end

    return (isempty(cases) ? "all" : join(cases, ","), abspath(ref_root))
end

const SELECTION, REF_ROOT = parse_args(ARGS)
const CASES = select_cases(SELECTION)

isempty(CASES) && @warn "no CI cases selected — nothing to compare"
isdir(REF_ROOT) || @warn "reference tree $REF_ROOT does not exist"
println("Comparing against references in ", REF_ROOT)

@testset "Benchmark comparisons" begin
    compare_cases(CASES; ref_root = REF_ROOT)
end

end # module BenchmarkComparison
