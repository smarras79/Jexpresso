#!/usr/bin/env julia
#
# scope_check.jl -- find names a function READS but never BINDS.
#
# WHY THIS EXISTS, AND WHY syntax_check.jl CANNOT DO IT
# -----------------------------------------------------
# In Julia a bare name inside a function that is not a local is not an error at
# parse time or at lowering time -- it is a global lookup, resolved when the
# line RUNS. So a value computed in one function and used in another, without
# being passed, is invisible to every static check that stops at lowering:
#
#     function build(...)          function report(..., checks)
#         sinf = ...                   vfast, sres = checks
#         report(..., (vfast, sres))   @printf("...", sinf)   # <- UndefVarError
#     end                          end                        #    at run time
#
# That shipped, and it threw only after a 256-rank job had done its mesh read,
# its operator build and its self-check -- the most expensive possible place to
# discover a typo. This catches it in a second, without loading the package.
#
# HOW IT WORKS. For each top-level definition, walk the AST and collect the
# names it BINDS (arguments, keyword arguments, assignments, for-loop targets,
# destructuring, `local`, `do` and `catch` variables, inner definitions) and
# the names it REFERENCES. Subtract the bindings, then subtract everything the
# package defines at column 0, everything Base and Core export, and the modules
# it uses. What is left is either a genuine missing binding or a name from a
# dependency, so the output is a REVIEW LIST, not a verdict -- run it on the
# files you changed and read the handful of lines it prints.
#
# An earlier version of this file worked on LOWERED code instead. Two blind
# spots killed it, and both reported "clean" on the bug it was written for:
# the method body's CodeInfo was not reachable by walking Expr args, and
# `"""doc""" function f() end` parses as Expr(:macrocall, @doc, ...), so
# matching Expr(:function) skipped nearly every function in the codebase.
#
#   julia tools/scope_check.jl src/kernel/solvers/hevi/imex3d.jl

# Collect names BOUND and names REFERENCED inside one function definition.
sym(x) = x isa Symbol ? x : nothing

function argnames!(out, a)
    a isa Symbol && (push!(out, a); return)
    a isa Expr || return
    if a.head === :(::) || a.head === :(...) || a.head === :kw
        argnames!(out, a.args[1])
    elseif a.head === :parameters || a.head === :tuple
        for x in a.args; argnames!(out, x); end
    end
end

function lhsnames!(out, a)
    a isa Symbol && (push!(out, a); return)
    a isa Expr || return
    if a.head === :tuple || a.head === :(::) || a.head === :(...)
        for x in a.args; lhsnames!(out, x); end
    elseif a.head === :ref || a.head === :(.)
        # a[i] = ... and a.b = ... bind nothing new
    end
end

function collect_binds!(b, e)
    e isa Expr || return
    if e.head === :(=)
        # `z() = zeros(...)` and `_r(dt) = ...` inside a function body bind BOTH
        # the name and its arguments. Treating them as plain assignments left
        # `z`, `_r` and `dt` on the review list on every run, which is the kind
        # of noise that gets a checker ignored.
        if e.args[1] isa Expr && e.args[1].head === :call
            f = e.args[1].args[1]
            f isa Symbol && push!(b, f)
            for a in e.args[1].args[2:end]; argnames!(b, a); end
        else
            lhsnames!(b, e.args[1])
        end
    elseif e.head === :for
        spec = e.args[1]
        for s in (spec isa Expr && spec.head === :block ? spec.args : [spec])
            s isa Expr && s.head === :(=) && lhsnames!(b, s.args[1])
        end
    elseif e.head === :local
        for x in e.args; lhsnames!(b, x); end
    elseif e.head === :function || e.head === :(->)
        argnames!(b, e.args[1])
        if e.args[1] isa Expr && e.args[1].head === :call
            for a in e.args[1].args[2:end]; argnames!(b, a); end
            s = sym(e.args[1].args[1]); s === nothing || push!(b, s)
        end
    elseif e.head === :do
        length(e.args) >= 2 && e.args[2] isa Expr && argnames!(b, e.args[2].args[1])
    elseif e.head === :try
        length(e.args) >= 2 && (s = sym(e.args[2])) !== nothing && push!(b, s)
    elseif e.head in (:generator, :comprehension)
        for a in e.args; collect_binds!(b, a); end
    end
    for a in e.args; collect_binds!(b, a); end
end

function collect_refs!(r, e)
    if e isa Symbol; push!(r, e); return; end
    e isa Expr || return
    if e.head === :(.)                       # a.b -> only `a` is a name
        collect_refs!(r, e.args[1]); return
    elseif e.head === :quote || e.head === :meta
        return
    elseif e.head === :kw                    # f(key = val) -> only val
        collect_refs!(r, e.args[2]); return
    elseif e.head === :macrocall
        for a in e.args[2:end]; collect_refs!(r, a); end; return
    end
    for a in e.args; collect_refs!(r, a); end
end

function fn_signature_binds(ex)
    b = Set{Symbol}()
    sig = ex.args[1]
    while sig isa Expr && sig.head === :where; sig = sig.args[1]; end
    if sig isa Expr && sig.head === :call
        for a in sig.args[2:end]; argnames!(b, a); end
    end
    return b
end

const FILES = isempty(ARGS) ? String[] : ARGS
isempty(FILES) && (println("usage: julia tools/scope_check.jl <file.jl> ..."); exit(2))

# ---- names the package defines at top level, harvested from source -----------
const DEFINED = Set{Symbol}()
# ANCHORED AT COLUMN 0, DELIBERATELY. An earlier version allowed leading
# whitespace on the bare-assignment pattern, so `    sinf = NaN` INSIDE a
# function was harvested as a package global -- and the checker then declared
# the very bug it was written for to be clean. A top-level definition in this
# codebase starts at column 0; anything indented is a local, and treating it
# otherwise is exactly the blind spot that matters.
# A LEADING MACRO IS STILL A TOP-LEVEL DEFINITION. `@inline function f`,
# `@with_kw struct S` and friends are how a good deal of this codebase defines
# things, and requiring `function` at column 0 missed all of them -- which put
# axpy_field!, hevi_trace and PhysicalConst on the review list as if they did
# not exist. `M` below is that optional prefix.
const M = raw"(?:@[\w.]+\s+)*"
defpat = [Regex("^" * M * raw"function\s+([A-Za-z_][A-Za-z0-9_!]*)", "m"),
          Regex("^" * M * raw"([A-Za-z_][A-Za-z0-9_!]*)\s*\([^)]*\)\s*=", "m"),
          Regex("^" * M * raw"const\s+([A-Za-z_][A-Za-z0-9_!]*)", "m"),
          Regex("^" * M * raw"(?:mutable\s+)?struct\s+([A-Za-z_][A-Za-z0-9_]*)", "m"),
          Regex("^" * M * raw"abstract\s+type\s+([A-Za-z_][A-Za-z0-9_]*)", "m"),
          Regex("^" * M * raw"macro\s+([A-Za-z_][A-Za-z0-9_!]*)", "m"),
          Regex("^" * M * raw"([A-Za-z_][A-Za-z0-9_!]*)\s*=\s*", "m")]
for (r, _, fs) in walkdir("src"), f in fs
    endswith(f, ".jl") || continue
    src = try read(joinpath(r, f), String) catch; continue end
    for p in defpat, m in eachmatch(p, src)
        push!(DEFINED, Symbol(m.captures[1]))
    end
    for m in eachmatch(r"^\s*(?:using|import)\s+([A-Za-z_][\w\.]*)"m, src)
        push!(DEFINED, Symbol(split(m.captures[1], '.')[1]))
    end
end
for mod in (Base, Core)
    for n in names(mod; all = false, imported = true); push!(DEFINED, n); end
end
for n in names(Main; all = true); push!(DEFINED, n); end

globalrefs(x, out) = nothing
function globalrefs(e::Expr, out)
    for a in e.args; globalrefs(a, out); end
end
globalrefs(g::GlobalRef, out) = push!(out, g.name)
function globalrefs(ci::Core.CodeInfo, out)
    for st in ci.code; globalrefs(st, out); end
end


bad = 0
for f in FILES
    src = read(f, String)
    pos = 1
    while pos <= lastindex(src)
        local ex
        try; ex, pos = Meta.parse(src, pos; raise = true); catch; break; end
        ex === nothing && break
        ex isa Expr || continue
        # UNWRAP THE DOCSTRING FIRST. `"""doc""" function f() end` parses as
        # Expr(:macrocall, @doc, ...), NOT Expr(:function).
        while ex isa Expr && ex.head === :macrocall &&
              (ex.args[1] === Symbol("@doc") ||
               (ex.args[1] isa GlobalRef && ex.args[1].name === Symbol("@doc")))
            ex = ex.args[end]
        end
        (ex isa Expr && ex.head in (:function, :(=))) || continue
        (ex.args[1] isa Expr && ex.args[1].head in (:call, :where)) || continue
        b = fn_signature_binds(ex)
        collect_binds!(b, ex.args[2])
        r = Set{Symbol}(); collect_refs!(r, ex.args[2])
        sig = ex.args[1]
        while sig isa Expr && sig.head === :where; sig = sig.args[1]; end
        fname = string(sig.args[1])
        unknown = sort(collect(setdiff(r, b, DEFINED)))
        filter!(n -> occursin(r"^[A-Za-z_][A-Za-z0-9_!]*$", string(n)), unknown)
        if !isempty(unknown)
            println("$f  in $fname:")
            for n in unknown; println("    reads `$n` -- not bound here and not defined at top level in src/"); end
            global bad += length(unknown)
        end
    end
end
println(bad == 0 ? "scope_check: no unresolved globals" : "\nscope_check: $bad name(s) to review")
exit(0)
