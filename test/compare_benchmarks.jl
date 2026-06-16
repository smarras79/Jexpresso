# Compares HDF5 output produced by the CI benchmark runs against committed
# reference solutions in test/CI-ref/.
#
# Pass/fail criteria:
#   - SKIP  : no reference files exist yet (run generate-ci-ref workflow first)
#   - FAIL  : simulation produced no HDF5 output (code crashed or wrong format)
#   - FAIL  : any field differs from the reference beyond atol=1e-5
#   - PASS  : all fields match within tolerance

module BenchmarkComparison

using HDF5
using Test

project_root = dirname(Base.current_project())

const BENCHMARKS = [
    ("CompEuler",    "theta"),
    ("CompEuler",    "3d"),
    ("CompEuler",    "wave1d"),
    ("Burgers",      "case1"),
    ("Burgers",      "case2d"),
    ("Elliptic",     "2dlaplace"),
    ("ShallowWater", "TC2"),
]

function find_h5_files(dir::String)
    isdir(dir) || return String[]
    sort(filter(f -> endswith(f, ".h5"), readdir(dir, join=true)))
end

function compare_h5_files(ref_path::String, gen_path::String)::Bool
    ref_data = h5open(ref_path, "r") do f
        Dict(k => read(f[k]) for k in keys(f))
    end
    gen_data = h5open(gen_path, "r") do f
        Dict(k => read(f[k]) for k in keys(f))
    end
    for key in keys(ref_data)
        if !haskey(gen_data, key)
            @error "Field '$key' missing from generated output"
            return false
        end
        vr, vg = ref_data[key], gen_data[key]
        ok = vr isa Array{Float64} ? isapprox(vr, vg; atol=1e-5) : vr == vg
        if !ok
            @error "Field '$key' differs from reference"
            return false
        end
    end
    return true
end

@testset "Benchmark comparisons" begin
    for (eqs, case) in BENCHMARKS
        @testset "$eqs/$case" begin
            ref_dir = joinpath(project_root, "test", "CI-ref",  eqs, case, "output")
            gen_dir = joinpath(project_root, "test", "CI-runs", eqs, case, "output")

            ref_files = find_h5_files(ref_dir)
            gen_files = find_h5_files(gen_dir)

            if isempty(ref_files)
                @warn "$eqs/$case: no reference HDF5 files in $ref_dir — " *
                      "run the 'generate-ci-ref' workflow to create them"
                @test_skip "no reference files"
            elseif isempty(gen_files)
                @error "$eqs/$case: simulation produced no HDF5 output in $gen_dir"
                @test false
            elseif length(ref_files) != length(gen_files)
                @error "$eqs/$case: file count mismatch " *
                       "(ref=$(length(ref_files)), gen=$(length(gen_files)))"
                @test false
            else
                for (ref_f, gen_f) in zip(ref_files, gen_files)
                    @testset "$(basename(gen_f))" begin
                        @test compare_h5_files(ref_f, gen_f)
                    end
                end
            end
        end
    end
end

end # module
