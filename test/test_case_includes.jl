#---------------------------------------------------------------------------------
# test/test_case_includes.jl — no case file may include itself, directly or
# through a cycle.
#
#   julia test/test_case_includes.jl
#
# NO `using Jexpresso`: this only parses files, so it runs in a bare Julia and
# costs a fraction of a second.
#
# WHAT IT GUARDS. src/run.jl loads a case by `include`ing the six user_*.jl in
# its directory. If one of those includes a path that resolves back to itself,
# `include` recurses until the stack is gone — and because the failure happens
# inside Julia's own loader, there is no Jexpresso output and no Julia-level
# stack trace, just a wall of
#
#     Internal error: during type inference of ...
#     Encountered stack overflow.
#
# which points nowhere near the actual file. That is expensive to diagnose from
# cold, so it is worth a test.
#
# It really happened. ShallowWater/SWsphere_visc used to consist of six one-line
# shims of the form
#
#     include(joinpath(@__DIR__, "..", "SWsphere", "user_flux.jl"))
#
# which is correct FROM SWsphere_visc/. Commit 6225a0b folded the visc case into
# SWsphere and copied the shims over SWsphere's real implementations, at which
# point @__DIR__ was itself .../SWsphere and all five files pointed at
# themselves. Every run of the case died on load until 9c06e63 restored them.
#
# Cycles through two or more files fail exactly the same way, so they are
# checked too, along with includes that point at a file which is not there.
#---------------------------------------------------------------------------------

using Test

const CASE_ROOTS = [joinpath(@__DIR__, "..", "problems"),
                    joinpath(@__DIR__, "..", "test", "CI-runs")]

#
# Resolve the argument of an `include(...)` to a path, when it can be done
# statically. `dir` is the directory of the file doing the including, which is
# both what @__DIR__ expands to and what Julia resolves a bare relative string
# against. Returns nothing for anything computed at run time.
#
function _resolve_include(arg, dir::String)
    arg isa String && return normpath(joinpath(dir, arg))
    if arg isa Expr && arg.head === :call && arg.args[1] === :joinpath
        parts = String[]
        for a in arg.args[2:end]
            if a isa String
                push!(parts, a)
            elseif a === Symbol("@__DIR__") ||
                   (a isa Expr && a.head === :macrocall && a.args[1] === Symbol("@__DIR__"))
                push!(parts, dir)
            else
                return nothing              # runtime-computed: cannot check here
            end
        end
        return normpath(joinpath(parts...))
    end
    return nothing
end

# Every include target of one file, found by walking the parsed AST rather than
# by regex, so a commented-out or string-embedded `include` cannot fool it.
function _include_targets(path::String)
    dir = dirname(path)
    out = Tuple{String,Bool}[]          # (target, resolved?)
    ast = try
        Meta.parseall(read(path, String); filename = path)
    catch
        return out
    end
    walk(x) = begin
        if x isa Expr
            if x.head === :call && !isempty(x.args) && x.args[1] === :include && length(x.args) >= 2
                r = _resolve_include(x.args[2], dir)
                push!(out, r === nothing ? ("", false) : (r, true))
            end
            foreach(walk, x.args)
        end
    end
    walk(ast)
    return out
end

function _all_case_files()
    files = String[]
    for root in CASE_ROOTS
        isdir(root) || continue
        for (dp, _, fns) in walkdir(root), fn in fns
            endswith(fn, ".jl") && push!(files, normpath(joinpath(dp, fn)))
        end
    end
    return files
end


@testset "case files never include themselves" begin

    files = _all_case_files()
    @test !isempty(files)                       # the walk found something at all

    graph = Dict{String,Vector{String}}()
    unresolved = String[]
    for f in files
        tgts = String[]
        for (t, ok) in _include_targets(f)
            ok ? push!(tgts, t) : push!(unresolved, f)
        end
        isempty(tgts) || (graph[f] = tgts)
    end

    @testset "no direct self-include" begin
        for (f, tgts) in graph, t in tgts
            @test t != f
            t == f && @error "case file includes ITSELF — every run of this case will die with a StackOverflowError on load" file=f
        end
    end

    @testset "no include cycle" begin
        # DFS over the include graph; a target already on the current stack is a
        # cycle and recurses just like a self-include.
        function findcycle(node, stack)
            for nxt in get(graph, node, String[])
                if nxt in stack
                    return vcat(stack[findfirst(==(nxt), stack):end], nxt)
                elseif haskey(graph, nxt)
                    c = findcycle(nxt, vcat(stack, nxt))
                    c === nothing || return c
                end
            end
            return nothing
        end
        for f in keys(graph)
            c = findcycle(f, [f])
            @test c === nothing
            c === nothing || @error "include cycle" cycle=join(c, " -> ")
        end
    end

    @testset "include targets exist" begin
        for (f, tgts) in graph, t in tgts
            @test isfile(t)
            isfile(t) || @error "include points at a missing file" file=f target=t
        end
    end

    # Not a failure — just surfaced, since a runtime-computed include path is
    # the one shape the checks above cannot see through.
    if !isempty(unresolved)
        @info "include paths that could not be resolved statically (review by hand)" files=unique(unresolved)
    end
end
