#==============================================================================
 test/ci_cases.jl

 SINGLE SOURCE OF TRUTH for the Jexpresso CI test suite.

 Everything that runs a case in CI reads this file:

   * test/runtests.jl            (`Pkg.test()` / the CI workflow)
   * test/compare_benchmarks.jl  (solution-vs-reference comparison)
   * .github/workflows/benchmarks.yml       (one job per case)
   * .github/workflows/generate-ci-ref.yml  (one job per case)

 To ADD a test to CI you add ONE line to `CI_CASES` below. Nothing else in
 the CI machinery has to change: the GitHub Actions matrices are generated
 from this list at run time.

 See test/CIdescription.md for the full checklist (the case also needs a
 test/CI-runs/<eqs>/<case>/ directory and, once generated, a reference
 solution in test/CI-ref/<eqs>/<case>/output/).

 Command line interface (used by the workflows, handy locally too):

   julia test/ci_cases.jl list                 # one "eqs/case" per line
   julia test/ci_cases.jl matrix [selection]   # GitHub Actions matrix JSON
   julia test/ci_cases.jl count  [selection]   # number of selected cases
   julia test/ci_cases.jl validate             # check the registry vs. disk

 `selection` is either "all" (default) or a comma-separated list of
 "eqs/case" entries, e.g. "CompEuler/theta,Burgers/case1".
==============================================================================#
module CICases

export CICase, CI_CASES, case_name, select_cases, matrix_json, validate,
       project_root, runs_dir, ref_dir, referenced_cases

"""
    CICase(eqs, case; timeout, atol)

One CI test case.

  * `eqs`     : equation-set directory name, e.g. `"CompEuler"`.
  * `case`    : case directory name inside it, e.g. `"theta"`.
                The case is loaded from `test/CI-runs/<eqs>/<case>/`.
  * `timeout` : per-case timeout in MINUTES used by the GitHub Actions jobs.
  * `atol`    : absolute tolerance used when comparing the HDF5 output
                against the reference solution in `test/CI-ref/`.
  * `vtk_smoke`: also run the case a second time writing VTK, and check that
                the writer produced non-empty files. Off by default: it is a
                second full run of the case, and it is a check on the writer,
                not on the solution (no reference, no tolerance). Turn it on
                for a case whose VTK output you care about keeping working.
"""
Base.@kwdef struct CICase
    eqs::String
    case::String
    timeout::Int     = 30
    atol::Float64    = 1e-5
    vtk_smoke::Bool  = false
end

CICase(eqs::AbstractString, case::AbstractString; kwargs...) =
    CICase(; eqs = String(eqs), case = String(case), kwargs...)

case_name(c::CICase) = string(c.eqs, "/", c.case)

#------------------------------------------------------------------------------
# THE REGISTRY
#
# Add one line per test. Keep it alphabetical-ish; the order only affects the
# order of the jobs in the GitHub Actions run.
#------------------------------------------------------------------------------
const CI_CASES = CICase[
    CICase(eqs = "CompEuler", case = "theta", timeout = 40, atol = 1e-5),
    CICase(eqs = "CompEuler", case = "sod1d", timeout = 10),
    CICase(eqs = "AdvDiff", case = "kopriva", timeout = 10),
    CICase(eqs = "Helmholtz", case = "case1_laguerre", timeout = 10),
    CICase(eqs = "Helmholtz", case = "case1", timeout = 10),
    CICase(eqs = "AdvDiff", case = "2d_Laguerre", timeout = 18),
    CICase(eqs = "ShallowWater", case = "SoliWaveIsland", timeout = 11),
    # The spherical shell. These two are the same Galewsky jet on the same
    # cubed sphere, differing only in the stabilization (modal filter + constant
    # ν vs. the residual-based DynSGS model), so running both is what keeps the
    # comparison between them honest.
    #
    # ~170 s and ~140 s of solve respectively on a developer machine, plus
    # process start-up; the timeout is set well above that because a GitHub
    # runner is the slower machine. Their references carry, besides the
    # conservative state, the RELATIVE VORTICITY — the field the Galewsky test
    # is judged on, since h barely moves while the instability develops — and,
    # for the DynSGS case, the four per-equation eddy viscosities the model
    # produced. Checking the state alone would leave the thing each case exists
    # to compute uncovered.
    CICase(eqs = "ShallowWater", case = "SWsphere",     timeout = 20),
    CICase(eqs = "ShallowWater", case = "SWsphereDSGS", timeout = 20),
    #<< test/generate_ci_ref.jl inserts new cases above this line >>

    #--------------------------------------------------------------------------
    # Cases that used to be in the suite and are waiting to be re-enabled.
    #
    # To bring one back: uncomment its line and run
    #     julia --project=. test/generate_ci_ref.jl <EQS>/<CASE>
    # which recreates the CI deck from problems/<EQS>/<CASE>/ if it is no
    # longer under test/CI-runs/, and generates the reference solution. Some
    # of these decks have been removed and some still have (stale) references
    # under test/CI-ref/ — either way the reference has to be regenerated,
    # since it predates the current HDF5 file naming.
    #
    CICase(eqs = "CompEuler",    case = "thetaTracers",   timeout = 40),
    CICase(eqs = "CompEuler",    case = "theta_laguerre", timeout = 40),
    # CICase(eqs = "CompEuler",    case = "3d",             timeout = 40),
    # CICase(eqs = "CompEuler",    case = "wave1d",         timeout = 20),
    CICase(eqs = "CompEuler",    case = "wave1d_lag",     timeout = 20),
    # CICase(eqs = "Burgers",      case = "case1",          timeout = 20),
    # CICase(eqs = "Burgers",      case = "case2d",         timeout = 30),
    # CICase(eqs = "Elliptic",     case = "2dlaplace",      timeout = 20),
    # CICase(eqs = "ShallowWater", case = "TC2",            timeout = 40),
    #--------------------------------------------------------------------------
]

#------------------------------------------------------------------------------
# Paths
#------------------------------------------------------------------------------
"Repository root, i.e. the directory that holds Project.toml."
project_root() = normpath(joinpath(@__DIR__, ".."))

"Directory holding the CI version of the case inputs."
runs_dir(c::CICase, root::AbstractString = project_root()) =
    joinpath(root, "test", "CI-runs", c.eqs, c.case)

"Directory holding the committed reference solution."
ref_dir(c::CICase, root::AbstractString = project_root()) =
    joinpath(root, "test", "CI-ref", c.eqs, c.case, "output")

"Directory a CI run writes its HDF5 output to (`CI_MODE=true`, `output_dir=none`)."
output_dir(c::CICase, root::AbstractString = project_root()) =
    joinpath(runs_dir(c, root), "output")

#------------------------------------------------------------------------------
# Selection
#------------------------------------------------------------------------------
"""
    select_cases(selection = "all")

Return the subset of `CI_CASES` named by `selection`: either `"all"`
(or an empty string) for the whole registry, or a comma-separated list of
`"eqs/case"` entries. Unknown entries throw.
"""
function select_cases(selection::AbstractString = "all")
    wanted = strip(selection)
    if isempty(wanted) || lowercase(wanted) == "all"
        return copy(CI_CASES)
    end

    selected = CICase[]
    for token in split(wanted, ',')
        name = strip(token)
        isempty(name) && continue
        idx = findfirst(c -> case_name(c) == name, CI_CASES)
        if idx === nothing
            known = join(map(case_name, CI_CASES), ", ")
            error("unknown CI case \"$name\". Registered cases: $known " *
                  "(add it to CI_CASES in test/ci_cases.jl first)")
        end
        push!(selected, CI_CASES[idx])
    end
    return selected
end

#------------------------------------------------------------------------------
# GitHub Actions matrix
#------------------------------------------------------------------------------
"Minutes added to a case timeout to cover checkout, setup-julia and buildpkg."
const JOB_TIMEOUT_HEADROOM = 30

json_escape(s::AbstractString) = replace(String(s), '\\' => "\\\\", '"' => "\\\"")

function matrix_entry(c::CICase)
    string("{",
           "\"name\":\"",    json_escape(case_name(c)), "\",",
           "\"eqs\":\"",     json_escape(c.eqs),        "\",",
           "\"case\":\"",    json_escape(c.case),       "\",",
           "\"slug\":\"",    json_escape(string(c.eqs, "-", c.case)), "\",",
           "\"timeout\":",   c.timeout, ",",
           # GitHub Actions expressions have no arithmetic, so the job-level
           # allowance (case timeout + setup/build headroom) is precomputed.
           "\"job_timeout\":", c.timeout + JOB_TIMEOUT_HEADROOM,
           "}")
end

"""
    matrix_json(selection = "all")

Single-line JSON object suitable for
`strategy.matrix: \${{ fromJSON(needs.<job>.outputs.matrix) }}`:

    {"include":[{"name":"CompEuler/theta","eqs":"CompEuler","case":"theta",
                 "slug":"CompEuler-theta","timeout":40}]}
"""
matrix_json(selection::AbstractString = "all") =
    string("{\"include\":[",
           join(map(matrix_entry, select_cases(selection)), ","),
           "]}")

#------------------------------------------------------------------------------
# Validation
#------------------------------------------------------------------------------
const REQUIRED_CASE_FILES = ("user_inputs.jl", "initialize.jl", "user_flux.jl",
                             "user_source.jl", "user_bc.jl", "user_primitives.jl")

"""
    deck_mesh(case_dir) -> String or nothing

The GMSH file a deck reads, or `nothing` when it builds its mesh itself
(`:lread_gmsh` absent or false — the solver then enforces a 1D problem).

Read textually rather than by evaluating `user_inputs()`, so the check costs
nothing and runs before any package is instantiated. Commented-out lines are
skipped: decks routinely keep alternative meshes commented above the live one.
"""
function deck_mesh(case_dir::AbstractString)
    file = joinpath(case_dir, "user_inputs.jl")
    isfile(file) || return nothing

    reads_gmsh = false
    mesh       = nothing
    for line in eachline(file)
        code = lstrip(line)
        startswith(code, "#") && continue
        m = match(r":lread_gmsh\s*=>\s*(true|false)", code)
        m === nothing || (reads_gmsh = (m.captures[1] == "true"))
        m = match(r":gmsh_filename\s*=>\s*\"([^\"]+)\"", code)
        m === nothing || (mesh = String(m.captures[1]))
    end
    return reads_gmsh ? mesh : nothing
end

"""
    tracked_references(root) -> Set{String} or nothing

The `test/CI-ref/…` paths git has under version control, or `nothing` when
that cannot be determined (git missing, not a repository, no worktree).

Only a *committed* reference is a claim about what the suite covers. A local
`test/CI-ref/<eqs>/<case>/` left behind by an experiment with
`generate_ci_ref.jl` is one developer's scratch space: it is invisible to
everyone else and to CI, so it must not fail their `Pkg.test()`. Paths come
back repository-relative with forward slashes, which is how they are built in
`referenced_cases`.
"""
function tracked_references(root::AbstractString)
    try
        listing = read(Cmd(`git -C $root ls-files -z -- test/CI-ref`), String)
        return Set(split(listing, '\0'; keepempty = false))
    catch
        # No git and no repository means an unpacked source archive, which
        # contains committed files and nothing else — the plain directory
        # scan is then already the right answer.
        return nothing
    end
end

"""
    referenced_cases(root = project_root()) -> Vector{String}

Every `"<eqs>/<case>"` that has a committed reference solution, i.e. at least
one version-controlled `.h5` file in `test/CI-ref/<eqs>/<case>/output/`.

This is the repository's view of what the suite covers. `CI_CASES` is the
registry's view. `validate` compares the two.
"""
function referenced_cases(root::AbstractString = project_root())
    ref_root = joinpath(root, "test", "CI-ref")
    found    = String[]
    isdir(ref_root) || return found

    tracked = tracked_references(root)
    is_reference = (eqs, case, file) ->
        endswith(file, ".h5") &&
        (tracked === nothing ||
         "test/CI-ref/$eqs/$case/output/$file" in tracked)

    for eqs in readdir(ref_root)
        eqs_dir = joinpath(ref_root, eqs)
        isdir(eqs_dir) || continue
        for case in readdir(eqs_dir)
            out = joinpath(eqs_dir, case, "output")
            isdir(out) || continue
            any(f -> is_reference(eqs, case, f), readdir(out)) || continue
            push!(found, string(eqs, "/", case))
        end
    end
    return sort(found)
end

"""
    disabled_entry(name; file = <this file>) -> String or nothing

The commented-out `CICase(...)` line for `name` in the registry, or `nothing`
when there is none.

A case disabled by commenting its line out is the failure mode this whole
check exists for: the deck, the reference and the documentation all still say
the case is tested, and the one line that decides whether it runs says it is
not. When that is what happened, `validate` can say so exactly instead of
reporting a generic "not registered".
"""
function disabled_entry(name::AbstractString; file::AbstractString = @__FILE__)
    isfile(file) || return nothing
    for line in eachline(file)
        code = lstrip(line)
        startswith(code, "#") || continue
        body = lstrip(code, ['#', ' ', '\t'])
        startswith(body, "CICase(") || continue

        m_eqs  = match(r"eqs\s*=\s*\"([^\"]+)\"",  body)
        m_case = match(r"case\s*=\s*\"([^\"]+)\"", body)
        (m_eqs === nothing || m_case === nothing) && continue
        string(m_eqs.captures[1], "/", m_case.captures[1]) == name || continue
        return body
    end
    return nothing
end

"""
    validate(; cases = CI_CASES, root = project_root(), orphans = true) -> Vector{String}

Check cases against the repository: the `test/CI-runs` directory must exist
and contain the six `user_*.jl` files the solver includes unconditionally,
and any mesh the deck reads must be present. A missing reference solution is
reported as a note, not an error — `test/generate_ci_ref.jl` creates it.

With `orphans = true` it also checks the other direction, which is the one
that fails silently: every case with a committed reference under
`test/CI-ref/` must be registered in `CI_CASES`. A reference nobody compares
against does not make CI fail — it makes CI *green*, however wrong the
solution gets, while the deck, the reference and `test/CIdescription.md` all
still claim the case is covered. That check is against the whole registry, so
it does not depend on `cases`; pass `orphans = false` in the tool that is
about to register the case (`test/generate_ci_ref.jl`), which would otherwise
be blocked by the very inconsistency it fixes.

`cases` defaults to the whole registry, which is what the CI job checks. Pass
a subset when acting on specific cases, so that an unrelated broken entry
does not block the work at hand.

Returns the list of problems found (empty when everything checks out).
"""
function validate(; cases::AbstractVector{CICase} = CI_CASES,
                    root::AbstractString = project_root(),
                    orphans::Bool = true)
    problems = String[]

    if isempty(cases)
        push!(problems, "no cases to validate — CI_CASES is empty")
    end

    seen = Set{String}()
    for c in cases
        name = case_name(c)
        name in seen && push!(problems, "$name: duplicate entry in CI_CASES")
        push!(seen, name)

        dir = runs_dir(c, root)
        if !isdir(dir)
            push!(problems, "$name: missing case directory $(relpath(dir, root))")
            continue
        end
        for f in REQUIRED_CASE_FILES
            isfile(joinpath(dir, f)) ||
                push!(problems, "$name: missing $(relpath(joinpath(dir, f), root))")
        end

        # The mesh has to be IN THIS REPOSITORY. meshes/ is the developer's
        # link to smarras79/JexpressoMeshes and does not exist on a runner, so
        # a case pointing there dies a quarter of an hour into the job with
        # "Msh file not found". A CI case keeps its mesh next to its deck.
        mesh = deck_mesh(dir)
        if mesh !== nothing && !isfile(normpath(joinpath(root, mesh)))
            deck = relpath(joinpath(dir, "user_inputs.jl"), root)
            hint = occursin(r"(^|/)meshes/", mesh) ?
                "That path is inside meshes/, the link to the " *
                "smarras79/JexpressoMeshes repository, which a CI runner does " *
                "not have. Commit the mesh into this repository next to the " *
                "case deck — problems/$(c.eqs)/$(c.case)/ — and point " *
                ":gmsh_filename at it, e.g. " *
                "\"./problems/$(c.eqs)/$(c.case)/$(basename(mesh))\"." :
                "Commit it, or point :gmsh_filename at a mesh that is tracked " *
                "in this repository."
            push!(problems,
                  "$name: mesh $mesh does not exist (:gmsh_filename in $deck). " *
                  hint * " Paths are resolved from the repository root.")
        end

        c.timeout > 0 || push!(problems, "$name: timeout must be positive")
        c.atol    > 0 || push!(problems, "$name: atol must be positive")
    end

    #--------------------------------------------------------------------------
    # References nothing runs against.
    #
    # runtests.jl iterates CI_CASES, so a case that is not in it is simply not
    # run — no failure, no warning, no line in the test summary. Committing its
    # reference solution does not change that: the reference just sits there,
    # and the suite reports "all passed" while the case it belongs to has not
    # been evaluated since the day it was disabled.
    #--------------------------------------------------------------------------
    if orphans
        registered = Set(map(case_name, CI_CASES))
        for name in referenced_cases(root)
            name in registered && continue

            reference = "test/CI-ref/$name/output"
            disabled  = disabled_entry(name)
            fix = disabled === nothing ?
                "Add it to CI_CASES (or run `julia --project=. " *
                "test/generate_ci_ref.jl $name`, which registers it for you)." :
                "The registry has it commented out — uncomment that line:\n      " *
                disabled
            push!(problems,
                  "$name: $reference/ holds a reference solution, but the case " *
                  "is not in CI_CASES, so nothing ever compares against it and " *
                  "the suite stays green however far the solution drifts. " *
                  fix * "\n      If the reference is obsolete instead, delete " *
                  "test/CI-ref/$name/.")
        end
    end

    return problems
end

#------------------------------------------------------------------------------
# Command line interface
#------------------------------------------------------------------------------
function main(args::AbstractVector{<:AbstractString})
    command   = isempty(args) ? "list" : String(args[1])
    selection = length(args) >= 2 ? String(args[2]) : "all"

    if command == "list"
        for c in select_cases(selection)
            println(case_name(c))
        end
    elseif command == "matrix"
        println(matrix_json(selection))
    elseif command == "count"
        println(length(select_cases(selection)))
    elseif command == "validate"
        problems = validate()

        # A registered case with no reference yet is not an error — the
        # comparison skips it and says so — but it IS a case CI is not
        # actually checking, so it is worth a line here rather than only in
        # the middle of a test log.
        referenced = Set(referenced_cases())
        for c in CI_CASES
            case_name(c) in referenced ||
                println("NOTE: $(case_name(c)) has no reference solution in " *
                        "test/CI-ref/$(case_name(c))/output/ — the comparison " *
                        "will be skipped. Generate it with `julia --project=. " *
                        "test/generate_ci_ref.jl $(case_name(c))`.")
        end

        if isempty(problems)
            println("CI registry OK — $(length(CI_CASES)) case(s): ",
                    join(map(case_name, CI_CASES), ", "))
        else
            for p in problems
                println(stderr, "ERROR: ", p)
            end
            exit(1)
        end
    else
        println(stderr, "usage: julia test/ci_cases.jl " *
                        "[list|matrix|count|validate] [selection]")
        exit(2)
    end
    return nothing
end

end # module CICases

if abspath(PROGRAM_FILE) == @__FILE__
    CICases.main(ARGS)
end
