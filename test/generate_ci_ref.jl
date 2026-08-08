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
   0. copies problems/<eqs>/<case>/ to test/CI-runs/<eqs>/<case>/ if that CI
      deck does not exist yet (output directories are not copied), so adding
      a case is just naming it here;
   1. runs it with CI_MODE=true, i.e. from test/CI-runs/<eqs>/<case>/, which
      writes HDF5 output into test/CI-runs/<eqs>/<case>/output/. CI mode
      forces the deck's output settings — HDF5, next to the case inputs, no
      timestamped directory, and a single write at t=:tend however many
      diagnostics the original deck asked for — so a deck copied straight out
      of problems/ needs no editing at all (see src/run.jl);
   2. deletes the stale .h5 files in test/CI-ref/<eqs>/<case>/output/ (a
      reference set must match the run exactly — a leftover file from an
      older run would be reported as missing output forever);
   3. copies the fresh .h5 files there, plus the user_inputs.jl the run used,
      so it is always visible which deck produced the reference;
   4. adds the case to CI_CASES in test/ci_cases.jl if it is not there yet,
      with a timeout derived from how long the run actually took, so CI
      starts running it.

 So adding a case to CI is one command — everything below is done for you:

     julia --project=. test/generate_ci_ref.jl MyEquations/mycase
     git add test/CI-ref test/CI-runs test/ci_cases.jl && git commit && git push

 The one thing the script cannot decide for you is how long the case should
 run: :tend in the copied deck is whatever problems/ said. Shorten it if the
 case is slow, and re-run this script.

 Options:
   --copy-only, --no-run   skip the simulation and publish the output already
                           sitting in test/CI-runs/<eqs>/<case>/output/ (use
                           after running a case by hand). Needs no packages.
   --refresh-deck          re-copy the deck from problems/ even when the CI
                           copy already exists (discards edits to that copy).
   --no-register           leave test/ci_cases.jl alone.
   --dest DIR              write the references somewhere other than
                           test/CI-ref (the workflow stages them elsewhere
                           before committing them from a single job).
   -h, --help              this message.

 Cases are named "<eqs>/<case>". With no case given, every case already
 registered in test/ci_cases.jl is regenerated.

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
      --refresh-deck          re-copy the deck from problems/ even if the one
                              in test/CI-runs/ already exists (discards edits
                              made to the CI copy)
      --no-register           do not add new cases to test/ci_cases.jl
      --dest DIR              destination tree (default: test/CI-ref)
      -h, --help              show this message

    With no case given, every case registered in test/ci_cases.jl is used.
    A case that is not registered yet still works: its deck is copied from
    problems/<eqs>/<case>/ and its reference is generated, and you are told
    which line to add to test/ci_cases.jl.
    """)
end

struct Options
    selection::String
    dest::String
    run_sim::Bool
    refresh_deck::Bool
    register::Bool
end

#------------------------------------------------------------------------------
# Case decks
#------------------------------------------------------------------------------
"`problems/<eqs>/<case>` — where the case is developed."
problems_dir(c::CICase) =
    joinpath(CICases.project_root(), "problems", c.eqs, c.case)

"""
Recursive copy that skips output directories, so a deck that has been run
locally does not drag its own results into test/CI-runs (they would be
mistaken for the CI run's output and published as a reference).
"""
function copy_deck(source::AbstractString, target::AbstractString)
    mkpath(target)
    for entry in readdir(source)
        startswith(entry, "output") && continue
        from, to = joinpath(source, entry), joinpath(target, entry)
        isdir(from) ? copy_deck(from, to) : cp(from, to; force = true)
    end
    return nothing
end

"""
    ensure_deck!(c; refresh = false) -> Bool

Make sure `test/CI-runs/<eqs>/<case>/` exists, copying it from
`problems/<eqs>/<case>/` when it does not (or when `refresh` is set). This is
step one of adding a case to CI, done for you. Returns `false` if neither
directory exists.
"""
function ensure_deck!(c::CICase; refresh::Bool = false)
    ci_deck  = runs_dir(c)
    original = problems_dir(c)
    relative = path -> relpath(path, CICases.project_root())

    if isdir(ci_deck) && !refresh
        return true
    end

    if !isdir(original)
        if isdir(ci_deck)
            return true    # --refresh-deck with nothing to refresh from
        end
        println(stderr, "ERROR: $(case_name(c)): neither $(relative(ci_deck)) ",
                "nor $(relative(original)) exists — nothing to run")
        return false
    end

    action = isdir(ci_deck) ? "refreshing" : "creating"
    println(" # $action $(relative(ci_deck)) from $(relative(original))")
    copy_deck(original, ci_deck)

    isdir(ci_deck) && !isfile(joinpath(ci_deck, "user_inputs.jl")) &&
        println(stderr, "WARNING: $(case_name(c)): the copied deck has no ",
                "user_inputs.jl — the run will fail")

    println("   the CI run shortens nothing but the output cadence: if this ",
            "case is slow,\n   reduce :tend in $(relative(joinpath(ci_deck, "user_inputs.jl")))")
    return true
end

function parse_options(args::AbstractVector{<:AbstractString})
    cases   = String[]
    dest    = DEFAULT_DEST
    run_sim  = true
    refresh  = false
    register = true

    i = 1
    while i <= length(args)
        arg = String(args[i])
        if arg in ("-h", "--help")
            usage()
            exit(0)
        elseif arg in ("--copy-only", "--no-run")
            run_sim = false
        elseif arg in ("--refresh-deck", "--recopy")
            refresh = true
        elseif arg == "--no-register"
            register = false
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
                   abspath(dest), run_sim, refresh, register)
end

"""
    resolve_cases(selection) -> Vector{CICase}

The registered cases named by `selection`, plus ad-hoc entries for names not
in `test/ci_cases.jl` yet — that is what lets a brand new case be generated
before it is registered (`register_case!` adds it once it has run).
"""
function resolve_cases(selection::AbstractString)
    (isempty(strip(selection)) || lowercase(strip(selection)) == "all") &&
        return select_cases("all")

    cases = CICase[]
    for token in split(selection, ',')
        name = String(strip(token))
        isempty(name) && continue
        known = findfirst(c -> case_name(c) == name, CI_CASES)
        if known !== nothing
            push!(cases, CI_CASES[known])
            continue
        end

        parts = split(name, '/')
        length(parts) == 2 && all(!isempty, parts) ||
            error("\"$name\" is not a case name — use the form <eqs>/<case>, " *
                  "e.g. CompEuler/theta")
        push!(cases, CICase(eqs = String(parts[1]), case = String(parts[2])))

        println(" # $name is new: it will be run, referenced, and added to " *
                "test/ci_cases.jl")
    end
    return cases
end

#------------------------------------------------------------------------------
"""
Run one case exactly as CI does. Returns `(ok, seconds)` — the wall time is
what the registry timeout is derived from.
"""
function run_case(c::CICase)
    println("\n", "="^78)
    println("Running $(case_name(c))  (CI_MODE = true)")
    println("="^78)
    started = time()
    try
        # Loaded lazily so that --copy-only works with a bare `julia`, i.e.
        # without instantiating the project.
        jexpresso = Base.require(Main, :Jexpresso)
        Base.invokelatest(jexpresso.run_case, c.eqs, c.case; CI_MODE = true)
        return true, time() - started
    catch err
        println(stderr, "\nERROR while running $(case_name(c)):")
        showerror(stderr, err)
        println(stderr)
        return false, time() - started
    end
end

#------------------------------------------------------------------------------
# Registration
#------------------------------------------------------------------------------
const REGISTRY_FILE   = joinpath(CICases.project_root(), "test", "ci_cases.jl")
const REGISTRY_ANCHOR = "#<< test/generate_ci_ref.jl inserts new cases above this line >>"

"""
    suggest_timeout(seconds) -> Int

Minutes to allow this case in CI: three times the wall time measured here,
rounded up, never less than 10. A GitHub runner is slower than a laptop and a
timeout is a safety net, not a target.
"""
suggest_timeout(seconds::Real) = max(10, ceil(Int, 3 * seconds / 60))

"""
    register_case!(c, timeout) -> Bool

Add `c` to `CI_CASES` in test/ci_cases.jl, in place, so CI actually runs it.
Already-registered cases are left alone. Returns `true` if the file was
modified.
"""
function register_case!(c::CICase, timeout::Int)
    any(registered -> case_name(registered) == case_name(c), CI_CASES) && return false

    # The anchor carries its own indentation in the file, so the entry is
    # inserted unindented and the anchor is pushed onto the next line.
    entry = "CICase(eqs = \"$(c.eqs)\", case = \"$(c.case)\", timeout = $timeout),"
    text  = read(REGISTRY_FILE, String)

    if !occursin(REGISTRY_ANCHOR, text)
        println(stderr, """
        WARNING: could not find the insertion anchor in test/ci_cases.jl.
                 Add this line to CI_CASES by hand so CI runs the case:

            $entry
        """)
        return false
    end

    write(REGISTRY_FILE,
          replace(text, REGISTRY_ANCHOR => entry * "\n    " * REGISTRY_ANCHOR, count = 1))
    println(" # registered in test/ci_cases.jl:\n    ", entry)
    return true
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
                "nothing to publish. CI mode forces HDF5 output, so this means ",
                "the run wrote nothing at all: it failed early, its diagnostics ",
                "settings produce no output, or it was run with ",
                "JEXPRESSO_CI_OUTPUT=0.")
        return 0
    end

    # t.h5 alone is not a solution: write_hdf5 writes the time stamp first and
    # one var_<i>_<rank>.h5 per variable after it, so "only t.h5" means the
    # per-variable writes failed. Publishing that would replace a good
    # reference with an empty one.
    if !any(startswith(f, "var_") for f in h5)
        println(stderr, "WARNING: $(case_name(c)): $source has no var_*.h5 ",
                "files, only $(join(h5, ", ")) — the per-variable writes did ",
                "not happen, so there is no solution to publish")
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
    cases   = resolve_cases(options.selection)

    println("Reference generation")
    println("  cases       : ", join(map(case_name, cases), ", "))
    println("  destination : ", options.dest)
    println("  simulate    : ", options.run_sim ? "yes" : "no (--copy-only)")

    # Before anything else: make sure every case has a deck in test/CI-runs,
    # copying it from problems/ when this is a brand new case. Done here and
    # not after `validate()` precisely because a new case is what makes the
    # registry look "inconsistent".
    for c in cases
        ensure_deck!(c; refresh = options.refresh_deck) ||
            error("$(case_name(c)): no case deck to run")
    end

    problems = validate()
    isempty(problems) || error("test/ci_cases.jl is inconsistent with the " *
                               "repository:\n  " * join(problems, "\n  "))

    results   = Tuple{String,Bool,Int}[]   # case, ran ok, files published
    registered = String[]
    for c in cases
        ran, seconds = options.run_sim ? run_case(c) : (true, 0.0)
        published    = publish_case(c, options.dest)
        push!(results, (case_name(c), ran, published))

        # A reference nobody compares against is pointless, so a case that
        # produced one gets added to the registry — the last manual step.
        if options.register && ran && published > 0 &&
           options.dest == abspath(DEFAULT_DEST)
            register_case!(c, suggest_timeout(seconds)) &&
                push!(registered, case_name(c))
        end
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
        if !isempty(registered)
            println("\nAdded to the CI suite: ", join(registered, ", "))
        end
        println("""

        Done. Review and commit — this is everything CI needs:

            git status --short test/CI-ref test/CI-runs test/ci_cases.jl
            git add test/CI-ref test/CI-runs test/ci_cases.jl
            git commit -m "Add CI reference solutions"
            git push
        """)
    end
    return nothing
end

end # module GenerateCIRef

if abspath(PROGRAM_FILE) == @__FILE__
    GenerateCIRef.main(ARGS)
end
