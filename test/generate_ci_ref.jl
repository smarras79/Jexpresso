#==============================================================================
 test/generate_ci_ref.jl

 (Re)generate the reference ("golden") solutions in test/CI-ref/ that the CI
 comparison checks every run against.

 This is the same job the `generate-ci-ref` GitHub workflow does — the
 workflow calls this very script — so you can do it locally and simply commit
 the result:

     julia --project=. test/generate_ci_ref.jl                  # all registered cases
     julia --project=. test/generate_ci_ref.jl CompEuler/theta  # one case
     git add test/CI-ref && git commit -m "Update CI reference solutions"

 What it does, per case:
   1. runs it with CI_MODE=true, i.e. from test/CI-runs/<eqs>/<case>/, which
      writes HDF5 output into test/CI-runs/<eqs>/<case>/output/;
   2. deletes the stale .h5 files in test/CI-ref/<eqs>/<case>/output/ (a
      reference set must match the run exactly — a leftover file from an
      older run would be reported as missing output forever);
   3. copies the fresh .h5 files there, plus the user_inputs.jl the run used,
      so it is always visible which deck produced the reference.

 Options:
   --copy-only, --no-run   skip the simulation and publish the output already
                           sitting in test/CI-runs/<eqs>/<case>/output/ (use
                           after running a case by hand). Needs no packages.
   --dest DIR              write the references somewhere other than
                           test/CI-ref (the workflow stages them elsewhere
                           before committing them from a single job).
   -h, --help              this message.

 Cases are named "<eqs>/<case>" and must be registered in test/ci_cases.jl;
 with no case given, every registered case is regenerated.

 Exit status is non-zero if any case failed to run or produced no output.
==============================================================================#
module GenerateCIRef

include(joinpath(@__DIR__, "ci_cases.jl"))
using .CICases

const DEFAULT_DEST = joinpath(CICases.project_root(), "test", "CI-ref")

function usage(io::IO = stdout)
    println(io, """
    usage: julia --project=. test/generate_ci_ref.jl [options] [<eqs>/<case> ...]

      --copy-only, --no-run   do not run the cases, publish existing output
      --dest DIR              destination tree (default: test/CI-ref)
      -h, --help              show this message

    With no case given, every case registered in test/ci_cases.jl is used.
    """)
end

struct Options
    selection::String
    dest::String
    run_sim::Bool
end

function parse_options(args::AbstractVector{<:AbstractString})
    cases   = String[]
    dest    = DEFAULT_DEST
    run_sim = true

    i = 1
    while i <= length(args)
        arg = String(args[i])
        if arg in ("-h", "--help")
            usage()
            exit(0)
        elseif arg in ("--copy-only", "--no-run")
            run_sim = false
        elseif arg == "--dest"
            i == length(args) && error("--dest needs a directory")
            dest = String(args[i + 1])
            i += 1
        elseif startswith(arg, "--dest=")
            dest = String(split(arg, '=', limit = 2)[2])
        elseif startswith(arg, "-")
            usage(stderr)
            error("unknown option \"$arg\"")
        else
            push!(cases, arg)
        end
        i += 1
    end

    return Options(isempty(cases) ? "all" : join(cases, ","),
                   abspath(dest), run_sim)
end

#------------------------------------------------------------------------------
"Run one case exactly as CI does. Returns `true` on success."
function run_case(c::CICase)
    println("\n", "="^78)
    println("Running $(case_name(c))  (CI_MODE = true)")
    println("="^78)
    try
        # Loaded lazily so that --copy-only works with a bare `julia`, i.e.
        # without instantiating the project.
        jexpresso = Base.require(Main, :Jexpresso)
        Base.invokelatest(jexpresso.run_case, c.eqs, c.case; CI_MODE = true)
        return true
    catch err
        println(stderr, "\nERROR while running $(case_name(c)):")
        showerror(stderr, err)
        println(stderr)
        return false
    end
end

"""
Publish the output of `c` into `dest`. Returns the number of .h5 files
published (0 means the run produced nothing to publish).
"""
function publish_case(c::CICase, dest::AbstractString)
    source = CICases.output_dir(c)
    h5 = isdir(source) ?
         sort(filter(f -> endswith(f, ".h5"), readdir(source))) : String[]

    if isempty(h5)
        println(stderr, "WARNING: $(case_name(c)): no .h5 files in $source — ",
                "nothing to publish (did the run write :outformat => \"hdf5\"?)")
        return 0
    end

    target = joinpath(dest, c.eqs, c.case, "output")
    mkpath(target)

    # A reference set must correspond to exactly one run.
    for stale in filter(f -> endswith(f, ".h5"), readdir(target))
        stale in h5 || println("  removing stale reference $stale")
        rm(joinpath(target, stale))
    end

    for f in h5
        cp(joinpath(source, f), joinpath(target, f); force = true)
    end
    # Keep the deck the reference was produced with next to it.
    inputs = joinpath(source, "user_inputs.jl")
    isfile(inputs) && cp(inputs, joinpath(target, "user_inputs.jl"); force = true)

    println("  published $(length(h5)) file(s) → $target")
    return length(h5)
end

#------------------------------------------------------------------------------
function main(args::AbstractVector{<:AbstractString})
    options = parse_options(args)
    cases   = select_cases(options.selection)

    problems = validate()
    isempty(problems) || error("test/ci_cases.jl is inconsistent with the " *
                               "repository:\n  " * join(problems, "\n  "))

    println("Reference generation")
    println("  cases       : ", join(map(case_name, cases), ", "))
    println("  destination : ", options.dest)
    println("  simulate    : ", options.run_sim ? "yes" : "no (--copy-only)")

    results = Tuple{String,Bool,Int}[]     # case, ran ok, files published
    for c in cases
        ran = options.run_sim ? run_case(c) : true
        published = publish_case(c, options.dest)
        push!(results, (case_name(c), ran, published))
    end

    println("\n", "="^78)
    println("Summary")
    println("="^78)
    failed = false
    for (name, ran, published) in results
        status = !ran           ? "FAILED (run)" :
                 published == 0 ? "FAILED (no output)" :
                                  "ok ($published file(s))"
        (ran && published > 0) || (failed = true)
        println("  ", rpad(name, 32), status)
    end

    if failed
        println(stderr, "\nSome cases did not produce a reference solution.")
        exit(1)
    end

    if options.dest == abspath(DEFAULT_DEST)
        println("""

        Done. Review and commit the new references:

            git status --short test/CI-ref
            git add test/CI-ref
            git commit -m "Update CI reference solutions"
            git push
        """)
    end
    return nothing
end

end # module GenerateCIRef

if abspath(PROGRAM_FILE) == @__FILE__
    GenerateCIRef.main(ARGS)
end
