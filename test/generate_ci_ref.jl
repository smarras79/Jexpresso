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
      deck does not exist yet — meshes are symlinked rather than copied, and
      :gmsh_filename is retargeted at the link (output directories are not
      copied) — so adding a case is just naming it here;
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

    Cases may be written CompEuler/sod1d, or as a bare pair CompEuler sod1d,
    or as "CompEuler", "sod1d" — quotes and commas are ignored.

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
Recursive copy of a case deck, skipping

  * output directories — a deck that has been run locally would otherwise
    drag its own results into test/CI-runs, where they would be mistaken for
    the CI run's output and published as a reference;
  * meshes — `link_meshes!` symlinks those instead of duplicating them.
"""
function copy_deck(source::AbstractString, target::AbstractString)
    mkpath(target)
    for entry in readdir(source)
        startswith(entry, "output") && continue
        endswith(entry, ".msh")     && continue
        from, to = joinpath(source, entry), joinpath(target, entry)
        isdir(from) ? copy_deck(from, to) : cp(from, to; force = true)
    end
    return nothing
end

"""
    link_meshes!(c)

Symlink every mesh of `problems/<eqs>/<case>/` into the CI deck, so the CI
case has the mesh at its own path without a second copy of the bytes in git
(a mesh is the largest thing a case owns). The link is *relative*, so it
resolves in every clone and in the source archive a runner unpacks; git
stores it as a link (mode 120000), which is the same trick test/meshes uses.

Falls back to a real copy where symlinks cannot be created (Windows without
developer mode), which is correct if wasteful.
"""
function link_meshes!(c::CICase)
    source = problems_dir(c)
    target = runs_dir(c)
    isdir(source) || return nothing

    for entry in readdir(source)
        endswith(entry, ".msh") || continue

        link = joinpath(target, entry)
        # test/CI-runs/<eqs>/<case>/<mesh> → four levels up is the repo root.
        destination = joinpath("..", "..", "..", "..",
                               "problems", c.eqs, c.case, entry)

        (islink(link) || ispath(link)) && rm(link; force = true)
        try
            symlink(destination, link)
            println("   linked $entry → problems/$(c.eqs)/$(c.case)/$entry")
        catch err
            @warn "$(case_name(c)): cannot symlink $entry " *
                  "($(sprint(showerror, err))) — copying it instead"
            cp(joinpath(source, entry), link; force = true)
        end
    end
    return nothing
end

"""
    retarget_mesh!(c)

Point the copied deck's `:gmsh_filename` at the mesh link that now sits in
the CI deck, rather than at the original path under problems/.

`:gmsh_filename` is resolved from the repository root, so without this the
link beside the deck would be dead weight and the CI deck would depend on the
problems/ path staying exactly as it is. Only a live (uncommented) entry whose
file is actually present in the CI deck is rewritten; a deck pointing
somewhere else entirely is left alone for `validate()` to judge.
"""
function retarget_mesh!(c::CICase)
    deck = joinpath(runs_dir(c), "user_inputs.jl")
    isfile(deck) || return nothing

    lines   = readlines(deck)
    changed = false
    for (i, line) in enumerate(lines)
        startswith(lstrip(line), "#") && continue
        m = match(r"^(.*:gmsh_filename\s*=>\s*\")([^\"]+)(\".*)$", line)
        m === nothing && continue

        mesh = basename(m.captures[2])
        isfile(joinpath(runs_dir(c), mesh)) || continue

        retargeted = "./test/CI-runs/$(c.eqs)/$(c.case)/$mesh"
        m.captures[2] == retargeted && continue
        lines[i] = string(m.captures[1], retargeted, m.captures[3])
        changed  = true
        println("   :gmsh_filename → $retargeted (the mesh link in the CI deck)")
    end

    changed && write(deck, join(lines, "\n") * "\n")
    return nothing
end

"""
    deck_drift(c) -> Vector{String}

The case files that differ, byte for byte, between `problems/<eqs>/<case>/`
and the CI copy in `test/CI-runs/<eqs>/<case>/`.
"""
function deck_drift(c::CICase)
    original, ci_deck = problems_dir(c), runs_dir(c)
    (isdir(original) && isdir(ci_deck)) || return String[]

    files = (CICases.REQUIRED_CASE_FILES..., "user_analytic.jl", "user_exactGeo.jl")
    return filter(collect(files)) do f
        a, b = joinpath(original, f), joinpath(ci_deck, f)
        isfile(a) && isfile(b) ? read(a) != read(b) : (isfile(a) || isfile(b))
    end
end

"""
    deck_keys(user_inputs_file) -> Set{Symbol}

The `:key =>` entries a deck actually sets, read textually — commented-out
lines skipped, since decks routinely keep alternatives commented above the
live one. Nothing here evaluates `user_inputs()`, so it costs nothing and
works before any package is instantiated.
"""
function deck_keys(file::AbstractString)
    found = Set{Symbol}()
    isfile(file) || return found
    for line in eachline(file)
        code = lstrip(line)
        startswith(code, "#") && continue
        m = match(r"^:([A-Za-z_][A-Za-z0-9_!]*)\s*=>", code)
        m === nothing || push!(found, Symbol(m.captures[1]))
    end
    return found
end

"""
    report_deck_drift(c)

Say which files the CI deck and the development deck disagree on.

The CI copy is authoritative once it exists (see `ensure_deck!`), so a fix
made in `problems/<eqs>/<case>/` does not reach CI on its own — the case
carries on being tested from the copy taken the day it was added. That is the
right default (the CI deck is deliberately a reduced version, with a shorter
`:tend` and its mesh retargeted), and it is also how a case silently stops
testing what its author is developing.

`user_inputs.jl` differing is therefore expected and gets no comment beyond
the listing. Any other file differing means the CI case runs different
physics from the case in `problems/`, which is worth saying out loud.
"""
function report_deck_drift(c::CICase)
    drifted = deck_drift(c)
    isempty(drifted) && return nothing

    println(" # note: $(case_name(c)) is run from its CI deck, ",
            "test/CI-runs/$(c.eqs)/$(c.case)/, which differs from ",
            "problems/$(c.eqs)/$(c.case)/ in: ", join(drifted, ", "))

    # The dangerous half of a user_inputs.jl drift is not a value that was
    # deliberately changed (:tend is MEANT to be shorter here) — it is a key
    # the development deck has and the CI copy does not. That key does not
    # announce itself: mod_inputs_user_inputs! quietly fills in the default,
    # and the case runs with different numerics from the one being developed.
    # :visc_model is the cautionary example. Its default is AV(), while
    # CompEuler/sod1d needs DSGS() to keep the shock from ringing, so a CI
    # copy predating that line diverges and writes no output at all.
    if "user_inputs.jl" in drifted
        missing_keys = setdiff(deck_keys(joinpath(problems_dir(c), "user_inputs.jl")),
                               deck_keys(joinpath(runs_dir(c),     "user_inputs.jl")))
        if !isempty(missing_keys)
            println("   set in problems/ but NOT in the CI deck: ",
                    join(sort!(collect(missing_keys)), ", "))
            println("   Those fall back to their defaults, silently. If the CI ",
                    "deck is simply out of date, re-copy it with --refresh-deck.")
        end
    end

    physics = filter(!isequal("user_inputs.jl"), drifted)
    isempty(physics) ||
        println("   ", join(physics, ", "), " differ too, so CI is exercising ",
                "different physics from the case you develop in problems/. ",
                "--refresh-deck re-copies the deck (discarding edits to the CI copy).")
    return nothing
end

"""
    ensure_deck!(c; refresh = false) -> Bool

Make sure `test/CI-runs/<eqs>/<case>/` exists, copying it from
`problems/<eqs>/<case>/` when it does not (or when `refresh` is set). This is
step one of adding a case to CI, done for you. Returns `false` if neither
directory exists.

Note what this does NOT do: refresh an existing CI deck. Once the copy
exists it is what the case is run from, and `problems/` is not consulted
again — so a setting added to the development deck (`:visc_model => DSGS()`,
say) never reaches CI, and the case quietly runs with the default instead
(`AV()`). `report_deck_drift` names the files that have diverged;
`--refresh-deck` re-copies them.
"""
function ensure_deck!(c::CICase; refresh::Bool = false)
    ci_deck  = runs_dir(c)
    original = problems_dir(c)
    relative = path -> relpath(path, CICases.project_root())

    if isdir(ci_deck) && !refresh
        report_deck_drift(c)
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
    link_meshes!(c)
    retarget_mesh!(c)

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

    return Options(isempty(cases) ? "all" : join(normalize_case_names(cases), ","),
                   abspath(dest), run_sim, refresh, register)
end

"""
    normalize_case_names(tokens) -> Vector{String}

Accept a case list written any of the ways people actually write one:

    CompEuler/sod1d                     the canonical form
    CompEuler sod1d                     a bare eqs/case pair
    "CompEuler", "sod1d"                copied out of Julia or a run command
    CompEuler/sod1d,CompEuler/theta     comma separated

Quotes and stray commas are stripped, and two consecutive bare words (neither
containing a `/`) are joined into one `eqs/case` name.
"""
function normalize_case_names(tokens::AbstractVector{<:AbstractString})
    cleaned = String[]
    for t in tokens, piece in split(t, ',')
        word = strip(piece, [' ', '"', '\'', ',', '\t'])
        isempty(word) || push!(cleaned, String(word))
    end

    names = String[]
    i = 1
    while i <= length(cleaned)
        word = cleaned[i]
        if occursin('/', word)
            push!(names, word)
            i += 1
        elseif i < length(cleaned) && !occursin('/', cleaned[i + 1])
            # "CompEuler", "sod1d" → CompEuler/sod1d
            push!(names, string(word, "/", cleaned[i + 1]))
            i += 2
        else
            push!(names, word)   # let resolve_cases produce the error
            i += 1
        end
    end
    return names
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
    started      = time()
    previous_dir = pwd()
    try
        # Decks address their mesh relative to the repository root, so run
        # from there no matter where this script was invoked from.
        cd(CICases.project_root())

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
    finally
        cd(previous_dir)
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
    reenable_case!(c) -> Bool

Uncomment the `CICase(...)` line for `c` in test/ci_cases.jl, if there is a
commented-out one, and return `true`.

A case that was disabled in place ("commented out sod1d for now") is
re-enabled where it stands rather than appended a second time: the registry
would otherwise carry two lines for the same case, one live and one dead, and
the next reader could not tell which one CI obeys. It also preserves the
timeout and tolerance the disabled line was carrying.
"""
function reenable_case!(c::CICase)
    lines = readlines(REGISTRY_FILE)
    for (i, line) in enumerate(lines)
        code = lstrip(line)
        startswith(code, "#") || continue
        body = lstrip(code, ['#', ' ', '\t'])
        startswith(body, "CICase(") || continue

        m_eqs  = match(r"eqs\s*=\s*\"([^\"]+)\"",  body)
        m_case = match(r"case\s*=\s*\"([^\"]+)\"", body)
        (m_eqs === nothing || m_case === nothing) && continue
        (m_eqs.captures[1] == c.eqs && m_case.captures[1] == c.case) || continue

        indent   = line[1:something(findfirst(!isspace, line), 1) - 1]
        lines[i] = string(indent, body)
        write(REGISTRY_FILE, join(lines, "\n") * "\n")
        println(" # re-enabled in test/ci_cases.jl:\n    ", body)
        return true
    end
    return false
end

"""
    register_case!(c, timeout) -> Bool

Add `c` to `CI_CASES` in test/ci_cases.jl, in place, so CI actually runs it.
Already-registered cases are left alone; a commented-out one is uncommented
rather than duplicated. Returns `true` if the file was modified.
"""
function register_case!(c::CICase, timeout::Int)
    any(registered -> case_name(registered) == case_name(c), CI_CASES) && return false

    reenable_case!(c) && return true

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
        # The most common cause by far is a solve that diverged. It does not
        # throw: the integrator prints "Instability detected. Aborting",
        # run_case returns as if nothing were wrong, and the only symptom is
        # this empty directory. Say so, rather than leaving the reader to
        # scroll back through a few thousand lines of solver output looking
        # for something that does not announce itself as an error.
        println(stderr, "WARNING: $(case_name(c)): no .h5 files in $source — ",
                "nothing to publish. CI mode forces HDF5 output and one write ",
                "at :tend, so this means the solve never reached :tend.\n",
                "  Look above for \"Instability detected. Aborting\": the ",
                "solution diverged, which is what a wrong RHS, too large a ",
                ":Δt or too little viscosity does. If you have been editing ",
                "the solver, check that the working tree is clean:\n",
                "      git status --short src/\n",
                "  Other causes: the run failed early (look for a stacktrace), ",
                "or it was run with JEXPRESSO_CI_OUTPUT=0.")
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

    # orphans = false: this script is the fix for an unregistered reference
    # (register_case! below adds or re-enables the entry), so it must not be
    # stopped by one. Everything else — decks, meshes, timeouts — still counts.
    problems = validate(cases = cases, orphans = false)
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
