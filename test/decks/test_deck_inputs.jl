#=============================================================================
 test_deck_inputs.jl -- CALL every deck's user_inputs() and check the Dict it
                        actually builds.

 WHY THIS EXISTS, AND WHY tools/syntax_check.jl IS NOT ENOUGH
 ------------------------------------------------------------
 syntax_check.jl parses AND lowers every file, which is more than a parse -- but
 it never CALLS anything. A whole class of deck bug is valid Julia that only
 fails when the Dict is built, and the worst of them is this one:

     :statistics_time => (10.0:10:100.0..., 9000.0:10:10800.0...)   # fine
     :statistics_time =>  10.0:10:100.0..., 9000.0:10:10800.0...    # 191 args

 Lose the parentheses and the splats land in the Pair call itself, so the deck
 asks for `Pair(::Symbol, ::Float64, ::Float64, ... x191)`. That parses, it
 lowers, and it dies at startup on every rank of a 256-rank job with a
 MethodError whose stack trace points at the `inputs = Dict(` line -- Julia
 attributes every Pair in a multi-line Dict to the line the Dict opens on --
 so it does not even tell you which key it was.

 This file evaluates the Dict and, when the MethodError is that one, NAMES THE
 KEY, which is the whole difference between a two-minute fix and an afternoon.

 Three decks in this repo carry splatted ranges (:statistics_time,
 :diagnostics_at_times), and their comments already record having been bitten
 twice, once per range. That is the test.

 STUBS, AND WHY A MISSING ONE IS A SKIP AND NOT A FAILURE
 --------------------------------------------------------
 Decks name solver and closure types from Jexpresso. deck_probe_worker.jl stubs
 the common ones, which costs nothing and needs no instantiated environment and
 no MPI -- the same reason syntax_check.jl avoids `using Jexpresso`. A deck that
 reaches for something not stubbed there -- a GPU backend, a Gridap model --
 reports SKIP rather than failing: this test is about the SHAPE of the Dict, and
 a missing stub is a fact about the test, not about the deck. The skip list is
 printed, and at least half the decks must really evaluate, so it cannot quietly
 grow until it covers everything.

 One worker process per deck (test/decks/deck_probe_worker.jl): a deck defines
 top-level names, and two decks in the same session would clobber each other and
 const-redefine their way into false failures.

     julia test/decks/test_deck_inputs.jl              # all decks
     julia test/decks/test_deck_inputs.jl <paths...>   # just these
=============================================================================#

using Test, Printf

say(s) = println(s)

function deck_paths()
    isempty(ARGS) || return ARGS
    ps = String[]
    for (root, _, fs) in walkdir("problems"), f in fs
        f == "user_inputs.jl" && push!(ps, joinpath(root, f))
    end
    return sort(ps)
end

const WORKER = joinpath(@__DIR__, "deck_probe_worker.jl")
const JULIA  = joinpath(Sys.BINDIR, Base.julia_exename())

paths = deck_paths()
say(@sprintf("\n=== evaluating user_inputs() for %d decks ===", length(paths)))

ok = String[]; skipped = Tuple{String,String}[]; bad = Tuple{String,String}[]

for p in paths
    # `ignorestatus`: the worker exits non-zero on every verdict except OK, and
    # that verdict is the whole point -- without it the failure is swallowed and
    # every deck reports the same useless "worker exited" line.
    out = try
        readchomp(pipeline(ignorestatus(`$JULIA --startup-file=no $WORKER $p`),
                           stderr = devnull))
    catch e
        "ERR could not run the worker: " * sprint(showerror, e)
    end
    line = isempty(out) ? "ERR the worker printed no verdict" : last(split(out, '\n'))

    if startswith(line, "OK")
        push!(ok, p)
    elseif startswith(line, "SKIP")
        push!(skipped, (p, strip(line[5:end])))
    elseif startswith(line, "PAIRARITY")
        f = split(line)
        push!(bad, (p, "key :$(f[2]) built a Pair with $(f[3]) arguments -- a `...` " *
                       "splat reached the Pair call itself. Put the value in " *
                       "parentheses: `:$(f[2]) => (a..., b...)`."))
    else
        push!(bad, (p, strip(line)))
    end
end

say(@sprintf("  %d evaluated, %d skipped (missing stub), %d failed",
             length(ok), length(skipped), length(bad)))
for (p, w) in skipped; say("  SKIP  $p\n          $w is not stubbed in deck_probe_worker.jl"); end
for (p, w) in bad;     say("  FAIL  $p\n          $w"); end

@testset "every deck's user_inputs() builds a Dict" begin
    @test isempty(bad)
    # A deck that cannot be evaluated is not a pass. This keeps the skip list
    # from quietly growing until it covers everything and the suite means
    # nothing.
    @test length(ok) >= 0.5 * length(paths)
end

say("\n=== done ===\n")
