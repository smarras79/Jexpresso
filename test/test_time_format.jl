# time_unit: the readable time of the run log, and that it allocates nothing.
# Run: julia --startup-file=no test/test_time_format.jl
using Test, Printf

module TimeFormatTest
    using Printf
    include(joinpath(@__DIR__, "..", "src", "auxiliary", "time_format.jl"))
end
const TF = TimeFormatTest

alloc_of(t) = @allocated TF.time_unit(t)
line_of(t)  = @sprintf("t=%.6g s (%.2f %s)", t, TF.time_unit(t)...)

@testset "time_unit" begin
    @test TF.time_unit(12.5)    == (12.5, "s")
    @test TF.time_unit(150.0)   == (2.5, "min")
    @test TF.time_unit(5400.0)  == (1.5, "h")
    @test TF.time_unit(7.776e6) == (90.0, "d")
    @test TF.time_unit(-7200.0) == (-2.0, "h")
    @test TF.time_unit(600f0)   == (10.0, "min")      # Float32 model time
    @test TF.time_unit(3)       == (3.0, "s")          # integer
    @test line_of(7.776e6)  == "t=7.776e+06 s (90.00 d)"
    @test line_of(1.2708e8) == "t=1.2708e+08 s (1470.83 d)"
    for t in (0.0, 12.5, 150.0, 5400.0, 7.776e6, 1.0e300, NaN, Inf, 600f0)
        alloc_of(t)                                    # compile first
        @test alloc_of(t) == 0
    end
end
