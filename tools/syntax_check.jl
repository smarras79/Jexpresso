#!/usr/bin/env julia
#
# syntax_check.jl -- parse AND lower every .jl file, without loading the package.
#
# WHY LOWERING AND NOT JUST PARSING
# ---------------------------------
# `Meta.parse` is not enough. A stray `$` outside a string parses cleanly -- it
# becomes Expr(:$, ...) -- and only fails when the expression is LOWERED:
#
#     julia> Meta.parse("f(m) = g((\$m + 4) * 2)")     # no error
#     julia> Meta.lower(Main, ans)
#     :($(Expr(:error, "\"\$\" expression outside quote")))
#
# That exact typo shipped in krylov.jl and broke precompilation on every rank
# of a 256-rank job, after a 133-second precompile, for a one-character error a
# check like this catches in under a second. A parse-only check reported the
# file as fine.
#
# This is deliberately NOT `using Jexpresso`: it needs no instantiated
# environment and no MPI, so it runs anywhere, including a laptop or a login
# node with nothing set up.
#
#   julia tools/syntax_check.jl [paths...]      (default: src test problems tools)
#
# Exit status is 1 if anything fails, so it can gate a commit or a job script.

const ROOTS = isempty(ARGS) ? ["src", "test", "problems", "tools"] : ARGS

# Some files in this repo are not Julia at all -- src/arrays/Maps.jl is a saved
# GitHub HTML page, committed long ago and included by nothing. Reporting it as
# a syntax failure would make this check always exit 1 and therefore useless as
# a gate, so it is called out separately instead.
isnotjulia(src) = occursin(r"^\s*<(!DOCTYPE|html)"i, first(src, 200))

function check(path)
    src = read(path, String)
    isnotjulia(src) && return :notjulia
    pos, n = 1, 0
    while pos <= lastindex(src)
        local ex
        try
            ex, pos = Meta.parse(src, pos; raise = true)
        catch e
            return "parse: " * sprint(showerror, e)
        end
        ex === nothing && break
        n += 1
        # Lowering expands macros, so it THROWS on any macro this process has
        # not loaded (@add_arg_table, @kernel, ...). That is not a defect in the
        # file -- skip it. Only a returned Expr(:error) is a real syntax fault;
        # that is the form `$` outside a quote takes, which is the case this
        # whole check exists for.
        lo = try
            Meta.lower(Main, ex)
        catch
            continue
        end
        if lo isa Expr && lo.head === :error
            return "lower: " * string(lo.args[1]) * "  (top-level expression $n)"
        end
    end
    return nothing
end

bad = Tuple{String,String}[]
skipped = String[]
nfiles = 0
for root in ROOTS
    if !isdir(root)
        if isfile(root) && endswith(root, ".jl")
            global nfiles += 1
            e = check(root)
            e === nothing || (e === :notjulia ? push!(skipped, root) : push!(bad, (root, e)))
        end
        continue
    end
    for (dir, _, files) in walkdir(root), f in files
        endswith(f, ".jl") || continue
        p = joinpath(dir, f)
        global nfiles += 1
        e = check(p)
        e === nothing && continue
        e === :notjulia ? push!(skipped, p) : push!(bad, (p, e))
    end
end
for p in skipped
    println("syntax_check: SKIP $p (not a Julia source file)")
end

if isempty(bad)
    println("syntax_check: $nfiles files OK")
    exit(0)
else
    for (p, e) in bad
        println("syntax_check: FAIL $p\n    $e")
    end
    println("\nsyntax_check: $(length(bad)) of $nfiles files failed")
    exit(1)
end
