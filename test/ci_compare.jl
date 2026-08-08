#==============================================================================
 test/ci_compare.jl

 Shared machinery for the CI test suite: run a registered case and compare
 its HDF5 output against the committed reference solution.

 Used by
   * test/runtests.jl            (`Pkg.test()`, the CI workflow)
   * test/compare_benchmarks.jl  (the benchmarks workflow)

 The list of cases lives in test/ci_cases.jl — this file never hardcodes a
 case name.

 Pass/fail criteria of `compare_case`:
   SKIP : no reference files in test/CI-ref/<eqs>/<case>/output yet
          (run the `generate-ci-ref` workflow to create them)
   FAIL : the run produced no HDF5 output (crash, or non-hdf5 :outformat)
   FAIL : a reference file has no counterpart in the run's output
   FAIL : any field differs from the reference by more than the case `atol`
   PASS : every field of every file matches within tolerance
==============================================================================#
module CICompare

using HDF5
using Test

include(joinpath(@__DIR__, "ci_cases.jl"))
using .CICases

export run_ci_case, compare_case, compare_cases

#------------------------------------------------------------------------------
# Helpers
#------------------------------------------------------------------------------
"Sorted list of the .h5 files in `dir` (empty when `dir` does not exist)."
function find_h5_files(dir::AbstractString)
    isdir(dir) || return String[]
    return sort(filter(f -> endswith(f, ".h5"), readdir(dir, join = true)))
end

read_h5(path::AbstractString) =
    h5open(path, "r") do file
        Dict{String,Any}(name => read(file[name]) for name in keys(file))
    end

"""
    compare_h5_files(ref_path, gen_path; atol) -> Bool

`true` when every dataset of the reference file is present in the generated
file and matches it within `atol` (floating point) or exactly (everything
else). Differences are reported with `@error` so the CI log points at the
offending field.
"""
function compare_h5_files(ref_path::AbstractString, gen_path::AbstractString;
                          atol::Real = 1e-5)
    ref_data = read_h5(ref_path)
    gen_data = read_h5(gen_path)

    ok = true
    for key in sort(collect(keys(ref_data)))
        if !haskey(gen_data, key)
            @error "field '$key' missing from $(basename(gen_path))"
            ok = false
            continue
        end
        vref, vgen = ref_data[key], gen_data[key]

        if vref isa AbstractArray{<:AbstractFloat} && vgen isa AbstractArray{<:AbstractFloat}
            if size(vref) != size(vgen)
                @error "field '$key' has size $(size(vgen)), reference is $(size(vref))"
                ok = false
            elseif !isapprox(vref, vgen; atol = atol)
                @error "field '$key' differs from the reference beyond atol=$atol " *
                       "(max |Δ| = $(maximum(abs.(vgen .- vref))))"
                ok = false
            end
        elseif vref != vgen
            @error "field '$key' differs from the reference"
            ok = false
        end
    end
    return ok
end

#------------------------------------------------------------------------------
# Run
#------------------------------------------------------------------------------
"""
    run_ci_case(c::CICase)

Run one registered case with `CI_MODE=true`, i.e. from
`test/CI-runs/<eqs>/<case>/`, writing its HDF5 output to the `output/`
subdirectory of that same folder. Wrapped in a `@testset` so a crash is a
test failure instead of an aborted suite.
"""
function run_ci_case(c::CICase)
    @testset "run $(case_name(c))" begin
        try
            # Load the package lazily (this file is also included by
            # compare_benchmarks.jl, which must stay HDF5-only and fast) and
            # call through invokelatest because the module is loaded at run
            # time, i.e. in a newer world age than this function.
            jexpresso = Base.require(Main, :Jexpresso)
            Base.invokelatest(jexpresso.run_case, c.eqs, c.case; CI_MODE = true)
            @test true
        catch err
            message = sprint(showerror, err)
            println("Error while running $(case_name(c)): ",
                    message[1:min(2000, end)])
            @test false
        end
    end
    return nothing
end

#------------------------------------------------------------------------------
# Compare
#------------------------------------------------------------------------------
"""
    compare_case(c::CICase; root = project_root())

Compare the output of case `c` against its reference solution inside a
`@testset` named after the case.
"""
function compare_case(c::CICase; root::AbstractString = CICases.project_root())
    @testset "compare $(case_name(c))" begin
        reference = ref_dir(c, root)
        generated = CICases.output_dir(c, root)

        ref_files = find_h5_files(reference)
        gen_files = find_h5_files(generated)

        if isempty(ref_files)
            @warn "$(case_name(c)): no reference HDF5 files in " *
                  "$(relpath(reference, root)) — run the 'generate-ci-ref' " *
                  "workflow to create them"
            @test_skip "no reference files"
        elseif isempty(gen_files)
            @error "$(case_name(c)): the run produced no HDF5 output in " *
                   "$(relpath(generated, root))"
            @test false
        else
            gen_by_name = Dict(basename(f) => f for f in gen_files)
            for ref_file in ref_files
                name = basename(ref_file)
                @testset "$name" begin
                    if !haskey(gen_by_name, name)
                        @error "$(case_name(c)): $name is missing from the run output"
                        @test false
                    else
                        @test compare_h5_files(ref_file, gen_by_name[name];
                                               atol = c.atol)
                    end
                end
            end
        end
    end
    return nothing
end

"""
    compare_cases(cases; root = project_root())

`compare_case` for every case in `cases`.
"""
function compare_cases(cases::AbstractVector{CICase};
                       root::AbstractString = CICases.project_root())
    for c in cases
        compare_case(c; root = root)
    end
    return nothing
end

end # module CICompare
