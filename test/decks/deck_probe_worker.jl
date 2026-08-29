#=============================================================================
 deck_probe_worker.jl -- evaluate ONE deck's user_inputs() and print a verdict.

 Driven by test/decks/test_deck_inputs.jl, one worker process per deck: a deck
 defines top-level names, and two decks in the same session would clobber each
 other's definitions and const-redefine their way into false failures.

     julia test/decks/deck_probe_worker.jl <path/to/user_inputs.jl>

 One line on stdout, and the exit status matches it:
     OK <nkeys>                  0
     SKIP <name>                 4   deck needs a type this file does not stub
     PAIRARITY <key> <nargs>     4   a `...` splat reached a Pair call
     ERR <message>               4
     NOTADICT <type>             3
=============================================================================#

module DeckStubs
struct SMAG end
struct VREM end
struct WALE end
struct DSMAG end
struct DSGS end
struct DSGS_MHD end
struct AV end
struct TOTAL end
struct PERT end
struct CL end
struct NCL end
struct FD end
struct CarpenterKennedy2N54 end
struct SSPRK33 end
struct SSPRK53 end
struct SSPRK54 end
struct ORK256 end
struct MSRK5 end
struct SPLIT_EXPLICIT end
struct GMRES; a; end
GMRES() = GMRES(nothing)
struct IMEX_ARK; tableau; end
IMEX_ARK() = IMEX_ARK(:ARS343)
struct HEVI_ARK; tableau; end
HEVI_ARK() = HEVI_ARK(:ARS232)
struct ARK; tableau; end
ARK() = ARK(:default)
struct RDPK3SpFSAL49 end          # OrdinaryDiffEq
struct ranocha end                # two-point volume flux; lowercase in the deck
end

for n in names(DeckStubs; all = true)
    isdefined(DeckStubs, n)            || continue
    startswith(String(n), "#")         && continue
    getfield(DeckStubs, n) isa Type    || continue
    Core.eval(Main, :(const $n = DeckStubs.$n))
end

try
    Base.include(Main, abspath(ARGS[1]))
catch e
    # A deck that will not even load is worth reporting as itself, not as a
    # missing verdict.
    if e isa UndefVarError
        println("SKIP ", e.var); exit(4)
    end
    println("ERR could not include the deck: ",
            sprint(showerror, e)[1:min(end, 400)])
    exit(4)
end

try
    d = Main.user_inputs()
    if !(d isa AbstractDict)
        println("NOTADICT ", typeof(d)); exit(3)
    end
    println("OK ", length(d))
catch e
    # THE ONE THIS FILE EXISTS FOR. `:k => a..., b...` (parentheses lost) asks
    # for Pair(::Symbol, ::T, ::T, ...). Julia attributes every Pair in a
    # multi-line Dict to the line the Dict OPENS on, so the run-time stack
    # trace never names the key. e.args[1] does.
    if e isa MethodError && e.f === Pair && length(e.args) > 2
        println("PAIRARITY ", e.args[1], " ", length(e.args))
    elseif e isa UndefVarError
        println("SKIP ", e.var)
    else
        println("ERR ", sprint(showerror, e)[1:min(end, 400)])
    end
    exit(4)
end
