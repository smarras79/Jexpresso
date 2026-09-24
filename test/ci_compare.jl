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

export run_ci_case, run_in_ci_mode, compare_case, compare_cases,
       vtk_smoke_case

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
    run_in_ci_mode(c::CICase)

Call `Jexpresso.run_case(...; CI_MODE = true)` for `c`. Throws on failure —
the callers below decide what a failure means.

Two details it takes care of. Decks address their mesh relative to the
repository root (`:gmsh_filename => "./problems/…/x.msh"`) while `Pkg.test()`
runs with the working directory set to test/, so the run happens from the
root; and the package is loaded lazily through `Base.require` +
`invokelatest`, because this file is also included by compare_benchmarks.jl,
which must stay HDF5-only and fast, and because a module loaded at run time
lives in a newer world age than this function.

The caller is responsible for restoring the working directory.
"""
function run_in_ci_mode(c::CICase)
    cd(CICases.project_root())
    jexpresso = Base.require(Main, :Jexpresso)
    Base.invokelatest(jexpresso.run_case, c.eqs, c.case; CI_MODE = true)
    return nothing
end

"""
    run_in_fresh_process(c::CICase)

Run `c` with `CI_MODE = true` in a Julia process of its own, and return `true`
when it exited successfully.

WHY A SEPARATE PROCESS, and not the cheaper in-process call above.

Every case defines the same six function NAMES — `user_source!`,
`user_flux!`, `initialize`, … — and `src/run.jl` `include`s them into one
module. Run two cases in one process and both sets of methods are live at
once. They do not simply overwrite each other, because the signatures differ
in ways that look cosmetic:

    CompEuler/theta      user_source!(S, q, qe, npoin::Int64, ::CL, ::PERT; …)
    AdvDiff/Wave_Train   user_source!(S, q, qe, npoin,        ::CL, ::PERT; …)

`npoin` is an `Int64` at the call site, so Julia picks the more specific
method — theta's — no matter which case is running. theta has four equations
and writes `S[3]`; Wave_Train has two, so its 2-element source buffer raises
a BoundsError. Observed exactly that: Wave_Train passes 4/4 alone and dies
after theta in the same suite.

A BoundsError is the lucky outcome. Two cases whose source terms happen to
have compatible arities would silently run each OTHER's physics and produce
a plausible, wrong answer — which the reference comparison would then either
fail for reasons nobody can trace, or, if the reference was generated the same
polluted way, pass forever.

So case isolation is not a nicety here, it is what makes a multi-case result
mean anything: the suite's verdict must not depend on the order of
`CI_CASES`. A process boundary is the only isolation Julia offers against
method-table pollution — a module cannot be unloaded. It costs one Julia
startup and package load per case, which is what the benchmarks and
generate-ci-ref workflows already pay by giving each case its own job.

Set `JEXPRESSO_CI_INPROCESS=1` to force the old in-process call — faster, and
useful when debugging a single case interactively, but only sound one case at
a time.
"""
function run_in_fresh_process(c::CICase)
    script = "using Jexpresso; Jexpresso.run_case($(repr(c.eqs)), " *
             "$(repr(c.case)); CI_MODE = true)"

    # Base.julia_cmd() carries this process's own flags (--check-bounds,
    # --depwarn, -O, …), so the child runs in the configuration the CI
    # workflow chose; dir = project_root() is what lets a deck address its
    # mesh relative to the repository root.
    cmd = Cmd(`$(Base.julia_cmd()) --project=$(Base.active_project()) -e $script`;
              dir = CICases.project_root())

    return success(pipeline(cmd; stdout = stdout, stderr = stderr))
end

"Whether to run each case in its own process (the default; see `run_in_fresh_process`)."
use_fresh_process() =
    !(lowercase(strip(get(ENV, "JEXPRESSO_CI_INPROCESS", "0"))) in ("1", "true", "yes", "on"))

"""
    run_ci_case(c::CICase)

Run one registered case with `CI_MODE=true`, i.e. from
`test/CI-runs/<eqs>/<case>/`, writing its HDF5 output to the `output/`
subdirectory of that same folder. Wrapped in a `@testset` so a crash is a
test failure instead of an aborted suite.
"""
function run_ci_case(c::CICase)
    @testset "run $(case_name(c))" begin
        previous_dir = pwd()

        # Clear the previous run's output first, so what compare_case reads is
        # necessarily what THIS run wrote.
        #
        # test/CI-runs/**/output/ is gitignored, so locally it survives from
        # one invocation to the next. Without this, a run that writes nothing —
        # diagnostics that never fire, a solver that returns early, an output
        # path that moved — leaves the previous run's files in place and the
        # comparison happily passes them against the reference. The suite then
        # reports success for a case whose result it has not seen. A CI runner
        # starts from a fresh clone and never had the leftovers, which is
        # exactly what makes this a local-only false pass: the developer's
        # green run is the one that is lying.
        for stale in find_h5_files(CICases.output_dir(c))
            rm(stale; force = true)
        end

        try
            ran_ok = use_fresh_process() ? run_in_fresh_process(c) :
                                           (run_in_ci_mode(c); true)

            # Returning normally is not the same as having solved the case.
            # An unstable solve does not throw: the integrator reports
            # "Instability detected. Aborting", run_case returns, and the only
            # trace left is an empty output directory. Attributed to the
            # comparison ("no HDF5 output") that reads as a problem with the
            # reference machinery, when what actually happened is that the
            # solver blew up — which is exactly what a wrong RHS tends to do.
            # So the run is what fails, and it says why.
            if !ran_ok
                @error "$(case_name(c)): the run process exited with a " *
                       "failure — its output is above."
                @test false
            elseif isempty(find_h5_files(CICases.output_dir(c)))
                @error "$(case_name(c)): the solver returned without writing " *
                       "any HDF5 output. CI mode forces HDF5, one write at " *
                       ":tend, so this means the solve never reached :tend. " *
                       "Look above for \"Instability detected. Aborting\" " *
                       "(a diverging solution — check the RHS, :Δt and the " *
                       "viscosity settings) or for an early return."
                @test false
            else
                @test true
            end
        catch err
            message = sprint(showerror, err)
            println("Error while running $(case_name(c)): ",
                    message[1:min(2000, end)])
            @test false
        finally
            cd(previous_dir)
        end
    end
    return nothing
end

#------------------------------------------------------------------------------
# Compare
#------------------------------------------------------------------------------
"""
    compare_case(c::CICase; root = project_root(), ref_root = <root>/test/CI-ref)

Compare the output of case `c` against its reference solution inside a
`@testset` named after the case.

`ref_root` points the comparison at a different reference tree — that is how
a reference produced elsewhere (the `ciref-*` artifact of a `generate-ci-ref`
run, say) is checked against a local run without being committed first.
"""
function compare_case(c::CICase; root::AbstractString = CICases.project_root(),
                      ref_root::AbstractString = joinpath(root, "test", "CI-ref"))
    # verbose: one row per reference file, not just a total. Which of the
    # solution variables moved is the first thing you want to know when a
    # comparison fails, and on a green run it is the evidence that every
    # field really was checked rather than the set being empty.
    @testset verbose = true "compare $(case_name(c))" begin
        reference = joinpath(ref_root, c.eqs, c.case, "output")
        generated = CICases.output_dir(c, root)

        ref_files = find_h5_files(reference)
        gen_files = find_h5_files(generated)

        if isempty(ref_files)
            @warn "$(case_name(c)): no reference HDF5 files in " *
                  "$reference — generate them with test/generate_ci_ref.jl"
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
    compare_cases(cases; root = project_root(), ref_root = <root>/test/CI-ref)

`compare_case` for every case in `cases`.
"""
function compare_cases(cases::AbstractVector{CICase};
                       root::AbstractString = CICases.project_root(),
                       ref_root::AbstractString = joinpath(root, "test", "CI-ref"))
    for c in cases
        compare_case(c; root = root, ref_root = ref_root)
    end
    return nothing
end

#------------------------------------------------------------------------------
# VTK writer smoke test
#------------------------------------------------------------------------------
"Every .vtu/.pvtu file under `dir`, at any depth (pvtk_grid writes a subdirectory)."
function find_vtk_files(dir::AbstractString)
    found = String[]
    isdir(dir) || return found
    for (root, _, files) in walkdir(dir), f in files
        (endswith(f, ".vtu") || endswith(f, ".pvtu")) && push!(found, joinpath(root, f))
    end
    return sort(found)
end

"""
    vtk_smoke_case(c::CICase)

Run `c` a second time with VTK output and check that the writer produced
non-empty files.

This is deliberately NOT a numerical comparison. VTK is the format production
runs use, so a break in `write_vtk` — a bad variable list, a malformed .pvtu,
an exception on some `SOL_VARS_TYPE` — would otherwise sail through a suite
that only ever writes HDF5. What is asserted is only that the writer ran and
produced output; the numbers are the HDF5 comparison's job.

Enable per case with `vtk_smoke = true` in test/ci_cases.jl, or for one run
with `julia --project=. test/runtests.jl --vtk <eqs>/<case>`.
"""
function vtk_smoke_case(c::CICase)
    @testset "vtk writer $(case_name(c))" begin
        output       = CICases.output_dir(c)
        previous_dir = pwd()
        previous_fmt = get(ENV, "JEXPRESSO_CI_OUTFORMAT", nothing)

        # Clear any .vtu left by an earlier run, so what is checked is what
        # this run wrote and not a leftover.
        for stale in find_vtk_files(output)
            rm(stale; force = true)
        end

        try
            ENV["JEXPRESSO_CI_OUTFORMAT"] = "vtk"
            run_in_ci_mode(c)
            @test true
        catch err
            message = sprint(showerror, err)
            if err isa MethodError && occursin("write_vtk", message)
                # Jexpresso implements write_vtk for NSD_2D and NSD_3D only;
                # 1D cases write PNG/ASCII (and HDF5). There is nothing to
                # smoke-test, so say that instead of dumping a MethodError.
                println("$(case_name(c)): this case has no VTK writer — " *
                        "Jexpresso implements write_vtk for 2D and 3D only, " *
                        "1D cases write PNG/ASCII. Drop `vtk_smoke = true` " *
                        "(or --vtk) for this case.")
            else
                println("Error while running $(case_name(c)) with VTK output: ",
                        message[1:min(2000, end)])
            end
            @test false
        finally
            previous_fmt === nothing ? delete!(ENV, "JEXPRESSO_CI_OUTFORMAT") :
                                       (ENV["JEXPRESSO_CI_OUTFORMAT"] = previous_fmt)
            cd(previous_dir)
        end

        written = find_vtk_files(output)
        if isempty(written)
            @error "$(case_name(c)): the VTK run wrote no .vtu/.pvtu in $output"
            @test false
        else
            empties = filter(f -> filesize(f) == 0, written)
            isempty(empties) ||
                @error "$(case_name(c)): empty VTK file(s): " *
                       join(map(basename, empties), ", ")
            @test isempty(empties)
            println("   VTK writer produced $(length(written)) file(s): ",
                    join(map(basename, written[1:min(4, end)]), ", "),
                    length(written) > 4 ? ", …" : "")
        end
    end
    return nothing
end

end # module CICompare
