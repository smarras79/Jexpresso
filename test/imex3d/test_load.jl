#=============================================================================
 test_load.jl -- catch method-table conflicts before precompilation does.

 Julia treats a method that overwrites another as a WARNING at the REPL and as
 a hard ERROR during module precompilation:

     WARNING: Method definition (::Type{IMEX_ARK})() ... overwritten at ...
     ERROR: Method overwriting is not permitted during Module precompilation.

 So a clash of this kind is invisible to anything that loads the sources with
 `include` -- which is what every other test here does -- and fatal the moment
 the package is precompiled. This test closes that gap: it loads the sources
 into a module in a fresh subprocess with --warn-overwrite=yes and fails on
 any overwrite warning, without needing the Jexpresso dependency stack.

 The clash it exists for, verbatim from the HEVI version: keyword-argument and
 default-argument constructors collide, because BOTH generate the same
 zero-argument positional method.

     IMEX_ARK(name::Symbol = :ARS343) = ...   # defines IMEX_ARK()
     IMEX_ARK(; tableau::Symbol = :ARS343)    # ALSO defines IMEX_ARK()

 There is a second reason to have this file rather than to rely on
 test/hevi/test_load.jl. IMEX3D shares ark.jl, operator.jl and hevi.jl with
 HEVI, and it adds methods to functions those files declare -- `ark_hooks`,
 `ark_relinearize!`. A new method on `ark_relinearize!(params, ::IMEX3DCache,
 u)` is fine; one whose signature happens to match HEVI's would silently
 replace it, and every HEVI run would then relinearise through the 3D path.
 Nothing else in either suite would notice.
=============================================================================#

using Test

const IMEX_SRC = joinpath(@__DIR__, "..", "..", "src", "kernel", "solvers", "hevi")
# The full include order Jexpresso.jl uses, so the probe fails on the same
# things a real precompile would.
const IMEX_FILES = ("columns.jl", "vdiffusion.jl", "operator.jl", "factorize.jl", "acoustic.jl",
                    "ark.jl", "hevi.jl", "cfl_diagnostics.jl", "substep.jl",
                    "krylov.jl", "precond_api.jl", "imex3d.jl")

incl(p) = string("include(raw\"", abspath(p), "\")")

code = string("module Imex3DLoadProbe\n",
              "using MPI, LinearAlgebra, Printf, OrdinaryDiffEq\n",
              incl(joinpath(@__DIR__, "..", "hevi", "mock_sem.jl")), "\n",
              join((incl(joinpath(IMEX_SRC, f)) for f in IMEX_FILES), "\n"),
              "\nend\nprintln(\"IMEX3D sources loaded\")\n")

@testset "IMEX3D sources load without method conflicts" begin
    err = IOBuffer()
    out = IOBuffer()
    ok = try
        # --warn-overwrite=yes is the whole point: Julia does NOT warn about
        # method overwriting by default, which is why a clash sails through
        # every include-based test and only surfaces when the package is
        # precompiled, where the same condition is a hard error.
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

    @test !occursin("overwritten", e)
    @test ok
    @test occursin("IMEX3D sources loaded", o)
end
