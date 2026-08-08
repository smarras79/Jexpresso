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
       project_root, runs_dir, ref_dir

"""
    CICase(eqs, case; timeout, atol)

One CI test case.

  * `eqs`     : equation-set directory name, e.g. `"CompEuler"`.
  * `case`    : case directory name inside it, e.g. `"theta"`.
                The case is loaded from `test/CI-runs/<eqs>/<case>/`.
  * `timeout` : per-case timeout in MINUTES used by the GitHub Actions jobs.
  * `atol`    : absolute tolerance used when comparing the HDF5 output
                against the reference solution in `test/CI-ref/`.
"""
Base.@kwdef struct CICase
    eqs::String
    case::String
    timeout::Int     = 30
    atol::Float64    = 1e-5
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
    #<< test/generate_ci_ref.jl inserts new cases above this line >>

    #--------------------------------------------------------------------------
    # Cases that used to be in the suite and are waiting to be re-enabled.
    # Uncomment a line once the case has been verified and its reference
    # solution regenerated with the `generate-ci-ref` workflow.
    #
    # CICase(eqs = "CompEuler",    case = "thetaTracers",   timeout = 40),
    # CICase(eqs = "CompEuler",    case = "theta_laguerre", timeout = 40),
    # CICase(eqs = "CompEuler",    case = "3d",             timeout = 40),
    # CICase(eqs = "CompEuler",    case = "wave1d",         timeout = 20),
    # CICase(eqs = "CompEuler",    case = "wave1d_lag",     timeout = 20),
    # CICase(eqs = "CompEuler",    case = "sod1d",          timeout = 20),
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
    validate([root]) -> Vector{String}

Check every registered case against the repository: the `test/CI-runs`
directory must exist and contain the six `user_*.jl` files the solver
includes unconditionally. A missing reference solution is reported as a
note, not an error — the `generate-ci-ref` workflow creates it.

Returns the list of problems found (empty when the registry is healthy).
"""
function validate(root::AbstractString = project_root())
    problems = String[]

    if isempty(CI_CASES)
        push!(problems, "CI_CASES is empty — no test would run")
    end

    seen = Set{String}()
    for c in CI_CASES
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

        # The mesh has to be IN THE REPOSITORY, not merely on the machine the
        # deck was written on: meshes/ is gitignored, so a mesh that was never
        # force-added does not exist on a CI runner, and the case dies a
        # quarter of an hour into the job with "Msh file not found".
        mesh = deck_mesh(dir)
        if mesh !== nothing && !isfile(normpath(joinpath(root, mesh)))
            push!(problems,
                  "$name: mesh $mesh does not exist. It is required by " *
                  ":gmsh_filename in $(relpath(joinpath(dir, "user_inputs.jl"), root)). " *
                  "Note that meshes/ is listed in .gitignore, so a mesh that " *
                  "only lives on your machine is invisible to CI — commit it " *
                  "(git add -f $mesh) or point the case at a mesh that is tracked")
        end

        c.timeout > 0 || push!(problems, "$name: timeout must be positive")
        c.atol    > 0 || push!(problems, "$name: atol must be positive")
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
