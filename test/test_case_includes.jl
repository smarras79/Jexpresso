#---------------------------------------------------------------------------------
# test/test_case_includes.jl — no case's user_*.jl may include ITSELF.
#
#   julia --project=. test/test_case_includes.jl
#
# The bug this pins down: `problems/ShallowWater/SWsphere/user_flux.jl` (and its
# four siblings) ended up containing
#
#     include(joinpath(@__DIR__, "..", "SWsphere", "user_flux.jl"))
#
# The line is correct in a case that *delegates* to SWsphere — it was written
# for the now-removed `SWsphere_visc` — but when those stub files were moved
# INTO SWsphere/ the path resolved back to the file doing the including. Julia
# recursed until the stack ran out, and because the failure happens while the
# case is being loaded, the StackOverflowError could not even be reported:
# `run_case("ShallowWater","SWsphere")` printed a wall of
#
#     Internal error: during type inference of
#     kwcall(..., Tuple{StackOverflowError, ...})
#     Encountered stack overflow.
#
# with no file, no line and no mention of SWsphere anywhere in it.
#
# Delegating between cases is a legitimate pattern (it is how two decks stay
# byte-identical apart from user_inputs.jl), so this test does not ban it. It
# bans only the degenerate case: an include whose target IS the including file.
#
# Nothing here loads Jexpresso or evaluates any case file — it is a textual
# scan of problems/ and test/CI-runs/, so it costs milliseconds.
#---------------------------------------------------------------------------------
using Test

const ROOT      = dirname(@__DIR__)
const CASE_ROOTS = (joinpath(ROOT, "problems"), joinpath(ROOT, "test", "CI-runs"))

# The per-case files src/run.jl includes by name.
const CASE_FILES = ("user_inputs.jl", "user_flux.jl", "user_source.jl",
                    "user_bc.jl", "initialize.jl", "user_primitives.jl",
                    "user_analytic.jl")

# `include(joinpath(@__DIR__, "..", "Foo", "user_flux.jl"))` and
# `include(joinpath(@__DIR__, "user_flux.jl"))` — capture the joinpath
# arguments after @__DIR__ so the target can be resolved against the file's
# own directory. Anything more dynamic than a list of string literals is left
# alone: this test is a tripwire, not an include resolver.
const INCLUDE_RE = r"include\s*\(\s*joinpath\s*\(\s*@__DIR__\s*,\s*([^)]*)\)\s*\)"
const STRING_RE  = r"\"([^\"]*)\""

"""
    self_include_targets(path) -> Vector{String}

Every `include(joinpath(@__DIR__, ...))` in `path` whose target resolves to
`path` itself, as the raw source line.
"""
function self_include_targets(path::AbstractString)
    hits = String[]
    dir  = dirname(path)
    for m in eachmatch(INCLUDE_RE, read(path, String))
        parts = [x.captures[1] for x in eachmatch(STRING_RE, m.captures[1])]
        isempty(parts) && continue
        target = normpath(joinpath(dir, parts...))
        target == normpath(path) && push!(hits, strip(m.match))
    end
    return hits
end

"""
    case_dirs() -> Vector{String}

Every directory under problems/ or test/CI-runs/ that holds at least one of
the per-case user files.
"""
function case_dirs()
    dirs = String[]
    for root in CASE_ROOTS
        isdir(root) || continue
        for (dir, _, files) in walkdir(root)
            any(f -> f in CASE_FILES, files) && push!(dirs, dir)
        end
    end
    return sort(dirs)
end

@testset "case user_*.jl files do not include themselves" begin
    dirs = case_dirs()
    @test !isempty(dirs)   # the scan found nothing at all → the paths are wrong

    offenders = String[]
    for dir in dirs, name in CASE_FILES
        path = joinpath(dir, name)
        isfile(path) || continue
        for line in self_include_targets(path)
            push!(offenders, string(relpath(path, ROOT), ": ", line))
        end
    end

    if !isempty(offenders)
        @info """A case file includes itself — loading it recurses until the
                 stack overflows, and the error cannot be reported because the
                 failure happens during loading. Point the include at the case
                 it should delegate to, or inline the real content:
                 """ * join(offenders, "\n                 ")
    end
    @test isempty(offenders)
end

# SWsphere is the case that regressed; check its files carry real content
# rather than a stub, so a future "fold case B back into case A" cannot empty
# them out again unnoticed.
@testset "SWsphere carries its own case files" begin
    swdir = joinpath(ROOT, "problems", "ShallowWater", "SWsphere")
    @test isdir(swdir)
    for (name, fn) in ("user_flux.jl"       => "user_flux!",
                       "user_source.jl"     => "user_source!",
                       "user_bc.jl"         => "user_bc_dirichlet!",
                       "user_primitives.jl" => "user_primitives!",
                       "initialize.jl"      => "function initialize(",
                       "user_inputs.jl"     => "function user_inputs(")
        path = joinpath(swdir, name)
        @test isfile(path)
        @test occursin(fn, read(path, String))
    end
end
