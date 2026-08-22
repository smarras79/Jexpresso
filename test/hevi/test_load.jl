#=============================================================================
 test_load.jl -- catch method-table conflicts before precompilation does.

 Julia treats a method that overwrites another as a WARNING at the REPL and as
 a hard ERROR during module precompilation:

     WARNING: Method definition (::Type{HEVI_ARK})() ... overwritten at ...
     ERROR: Method overwriting is not permitted during Module precompilation.

 So a clash of this kind is invisible to anything that loads the sources with
 `include` -- which is what every other test here does -- and fatal the moment
 the package is precompiled. This test closes that gap: it loads the HEVI
 sources into a module in a fresh subprocess and fails on any overwrite
 warning, without needing the Jexpresso dependency stack to precompile.

 The clash it was written for: keyword-argument and default-argument
 constructors collide, because BOTH generate the same zero-argument positional
 method.

     HEVI_ARK(name::Symbol = :ARS343) = ...   # defines HEVI_ARK()
     HEVI_ARK(; tableau::Symbol = :ARS343)    # ALSO defines HEVI_ARK()

 Keyword arguments do not participate in dispatch, so these are the same
 signature and the second silently replaces the first. This repo has been bitten
 by exactly this before -- see the long comment above assemble_mpi! in
 src/kernel/mpi/mpi_communications.jl.
=============================================================================#

using Test

const HEVI_SRC = joinpath(@__DIR__, "..", "..", "src", "kernel", "solvers", "hevi")
const HEVI_FILES = ("columns.jl", "operator.jl", "factorize.jl", "ark.jl",
                    "hevi.jl", "cfl_diagnostics.jl")

incl(p) = string("include(raw\"", abspath(p), "\")")

code = string("module HeviLoadProbe\n",
              "using MPI, LinearAlgebra, Printf, OrdinaryDiffEq\n",
              incl(joinpath(@__DIR__, "mock_sem.jl")), "\n",
              join((incl(joinpath(HEVI_SRC, f)) for f in HEVI_FILES), "\n"),
              "\nend\nprintln(\"HEVI sources loaded\")\n")

@testset "HEVI sources load without method conflicts" begin
    err = IOBuffer()
    out = IOBuffer()
    ok = try
        # --project is NOT part of julia_cmd(), so the subprocess would
        # otherwise start in the default environment and fail on `using MPI`
        # for reasons that have nothing to do with what is being tested.
        # --warn-overwrite=yes is the whole point: Julia does NOT warn about
        # method overwriting by default, which is why a clash like this sails
        # through every include-based test and only surfaces when the package
        # is precompiled, where the same condition is a hard error.
        #
        # --project is not part of julia_cmd(), so without it the subprocess
        # starts in the default environment and fails on `using MPI` for
        # reasons that have nothing to do with what is being tested.
        run(pipeline(`$(Base.julia_cmd()) --startup-file=no --warn-overwrite=yes
                      --project=$(Base.active_project()) -e $code`;
                     stdout = out, stderr = err))
        true
    catch
        false
    end
    e = String(take!(err))
    o = String(take!(out))

    if !isempty(e)
        println("---- stderr from the load probe ----")
        println(e)
        println("------------------------------------")
    end

    # The specific failure: two methods with the same signature.
    @test !occursin("overwritten", e)
    # Anything else that made the load fail is worth failing on too.
    @test ok
    @test occursin("HEVI sources loaded", o)
end
