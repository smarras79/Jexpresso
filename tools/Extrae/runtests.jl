#==============================================================================
# runtests.jl  --  test the Extrae.jl integration end-to-end
#
# This is the test you can run on your MacBook Air.  It verifies that:
#
#   1. the ExtraeShim API (init/register/emit/user_function/@user_function/
#      finish) is callable and never throws, on ANY platform;
#   2. the serial axpy example runs to completion;
#   3. the MPI example runs to completion on N ranks (default 4 — one per
#      core on a 4-core MacBook Air; override with JEXPRESSO_EXTRAE_NRANKS).
#
# On Linux with Extrae installed these same runs additionally produce real
# Paraver traces; on macOS the instrumentation is a no-op (see ExtraeShim.jl)
# but the test still passes, proving the instrumented code paths are correct.
#
# Run from the Jexpresso project root:
#
#     julia --project=. tools/Extrae/runtests.jl
#==============================================================================
module ExtraeIntegrationTests

using Test

const HERE         = @__DIR__
const PROJECT_ROOT = dirname(dirname(HERE))
const NRANKS       = parse(Int, get(ENV, "JEXPRESSO_EXTRAE_NRANKS", "4"))

include(joinpath(HERE, "ExtraeShim.jl"))
using .ExtraeShim

# A subprocess helper: run a Julia script with this project active and assert
# it exits 0.  We forward stdout/stderr so failures are easy to read.
function run_julia_script(script::String; nranks::Int = 0)
    julia = Base.julia_cmd()
    if nranks == 0
        cmd = `$julia --project=$PROJECT_ROOT $script`
    else
        # Use MPI.jl's bundled launcher so we honour the project's MPI config.
        mpiexec_cmd = readchomp(`$julia --project=$PROJECT_ROOT -e "using MPI; print(mpiexec())"`)
        cmd = `$mpiexec_cmd -n $nranks $julia --project=$PROJECT_ROOT $script`
    end
    @info "running" cmd
    proc = run(ignorestatus(setenv(cmd; dir = PROJECT_ROOT)))
    return proc.exitcode
end

@testset "Extrae.jl integration" begin

    @testset "shim API is always callable" begin
        # These must never throw, regardless of platform.
        @test init() isa Bool
        @test is_active() isa Bool
        @test register(84210, "Vector length") === nothing
        @test register(84211, "Element type", (16, 32, 64),
                       ("Float16", "Float32", "Float64")) === nothing
        @test emit(84210, 1_000_000) === nothing
        result = @user_function begin
            s = 0.0
            for i in 1:1000
                s += sqrt(i)
            end
            s
        end
        @test result > 0.0           # @user_function returns the body value
        @test user_function(true)  === nothing
        @test user_function(false) === nothing
        @test finish() === nothing
        @test is_active() == false   # finish() always clears the active flag
    end

    @testset "serial axpy example runs" begin
        rc = run_julia_script(joinpath(HERE, "extrae_axpy.jl"))
        @test rc == 0
    end

    @testset "MPI example runs on $NRANKS ranks" begin
        rc = run_julia_script(joinpath(HERE, "extrae_mpi_jexpresso_pattern.jl");
                              nranks = NRANKS)
        @test rc == 0
    end
end

end # module
