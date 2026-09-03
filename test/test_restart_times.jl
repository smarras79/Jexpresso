#---------------------------------------------------------------------------------
# test/test_restart_times.jl — the two forms of :restart_time.
#
#   julia --project=. test/test_restart_times.jl
#
# The HDF5 restart callback in time_loop! accepts :restart_time either as a
# period (a Number, 0.0 = never) or as a list of explicit times in the same
# form :diagnostics_at_times takes. The list form used to reach
# `rem(t, restart_time)` untouched and kill the run on the first step with
#
#     MethodError: no method matching rem(::Float64, ::NTuple{209, Float64})
#
# This pins down restart_times/restart_due on both forms, with the exact deck
# value that produced the 209-tuple. Nothing here touches the filesystem, a
# mesh, or an MPI collective.
#---------------------------------------------------------------------------------
using Test
using Jexpresso: restart_times, restart_due, flatten_times

@testset ":restart_time forms" begin
    @testset "period" begin
        @test restart_times(9000.0) === nothing
        @test restart_times(9000)   === nothing
        @test restart_times(0.0)    === nothing

        @test  restart_due(9000.0,  9000.0, nothing)
        @test  restart_due(18000.0, 9000.0, nothing)
        @test !restart_due(9001.0,  9000.0, nothing)
        # 0.0 disables restarts rather than dividing by zero.
        @test !restart_due(0.0,   0.0, nothing)
        @test !restart_due(100.0, 0.0, nothing)
    end

    @testset "list" begin
        tend = 10800.0
        spec = (0.0:100.0:1000.0..., 1000.0:500.0:9000.0..., 9000.0:10.0:tend...)
        @test spec isa NTuple{209, Float64}   # the deck value that hit rem()

        times = restart_times(spec)
        @test times isa Vector{Float64}
        @test times == flatten_times(spec)
        @test length(times) == 209

        @test  restart_due(0.0,     spec, times)
        @test  restart_due(9000.0,  spec, times)
        @test  restart_due(10800.0, spec, times)
        @test !restart_due(9005.0,  spec, times)
        @test !restart_due(50.0,    spec, times)

        # Every other list shape a deck can write.
        @test restart_times(9000.0:10.0:9030.0) == [9000.0, 9010.0, 9020.0, 9030.0]
        @test restart_times([1.0, 2.0])         == [1.0, 2.0]
        @test restart_times((100, 200.0))       == [100.0, 200.0]
        # A tuple that mixes numbers and an unsplatted range still flattens.
        @test restart_times((1.0, 2.0:1.0:3.0)) == [1.0, 2.0, 3.0]
    end
end
