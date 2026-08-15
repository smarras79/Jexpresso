#---------------------------------------------------------------------------------
# test/test_sphere_output.jl — WHEN the spherical-shell time loop writes
# (sphere_output_times in src/kernel/solvers/sphere_time_loop.jl).
#
#   julia --project=. test/test_sphere_output.jl
#
# This is a pure function of the inputs dict, so the test needs no grid and runs
# in seconds. It exists because the shell loop silently wrote NOTHING under
# CI_MODE: it understood :ndiagnostics_outputs and not :diagnostics_at_times,
# and the two are coupled — mod_inputs_user_inputs! zeroes the first whenever
# the deck (or CI) sets the second. ShallowWater/SWsphere therefore solved all
# the way to :tend and produced not one file, which surfaced as
# test/ci_compare.jl's "the solver returned without writing any HDF5 output" and
# reads like a blow-up rather than an empty schedule.
#
# An empty schedule is a legitimate request, so it cannot be an error. What it
# must not be is a SURPRISE, hence this test and the banner line next to it.
#---------------------------------------------------------------------------------

using Test
using Jexpresso
using Jexpresso: sphere_output_times

@testset verbose = true "spherical shell: output schedule" begin

    #-----------------------------------------------------------------------
    # The decks in problems/: n equally spaced writes, the last exactly at
    # :tend. This is the branch SWsphere/SWsphereDSGS take, and it has to keep
    # landing on :tend to the last bit — the final file IS the result.
    #-----------------------------------------------------------------------
    @testset ":ndiagnostics_outputs" begin
        touts = sphere_output_times(Dict{Symbol,Any}(:ndiagnostics_outputs => 24),
                                    0.0, 864000.0)
        @test length(touts) == 24
        @test touts[1]   ≈ 36000.0
        @test touts[end] == 864000.0          # exactly, not ≈
        @test issorted(touts)

        # a non-zero :tinit shifts the whole schedule rather than rescaling it
        touts2 = sphere_output_times(Dict{Symbol,Any}(:ndiagnostics_outputs => 4),
                                     100.0, 500.0)
        @test touts2 == [200.0, 300.0, 400.0, 500.0]
    end

    #-----------------------------------------------------------------------
    # THE CI PATH. run.jl sets :diagnostics_at_times => [tend]; mod_inputs then
    # sets :ndiagnostics_outputs = 0 (its else branch). Before the fix that
    # combination produced no writes at all.
    #-----------------------------------------------------------------------
    @testset ":diagnostics_at_times (the CI_MODE combination)" begin
        ci = Dict{Symbol,Any}(:ndiagnostics_outputs => 0,
                              :diagnostics_at_times => [864000.0])
        @test sphere_output_times(ci, 0.0, 864000.0) == [864000.0]

        # the deck's own :ndiagnostics_outputs must NOT win over an explicit
        # list once mod_inputs has zeroed it — that zero is the signal
        @test sphere_output_times(Dict{Symbol,Any}(:ndiagnostics_outputs => 0,
                                                   :diagnostics_at_times => [1.0, 2.0]),
                                  0.0, 2.0) == [1.0, 2.0]

        # scalar form: mod_inputs defaults :diagnostics_at_times to :tend
        # itself (a Number, not a vector) when the deck sets neither key
        @test sphere_output_times(Dict{Symbol,Any}(:diagnostics_at_times => 500.0),
                                  0.0, 500.0) == [500.0]

        # ranges are accepted, as TimeIntegrators.jl accepts them
        @test sphere_output_times(Dict{Symbol,Any}(:diagnostics_at_times => 0:10:30),
                                  0.0, 30.0) == [10.0, 20.0, 30.0]
    end

    #-----------------------------------------------------------------------
    # Normalisation: sorted, de-duplicated, and nothing at or before :tinit.
    # The monitor walks the list with one index and writes while
    # t >= touts[iout+1], so an unsorted or duplicated list would either skip
    # entries or burn several file numbers on one step.
    #-----------------------------------------------------------------------
    @testset "normalisation" begin
        touts = sphere_output_times(
            Dict{Symbol,Any}(:diagnostics_at_times => [30.0, 10.0, 20.0, 10.0, 0.0, -5.0]),
            0.0, 30.0)
        @test touts == [10.0, 20.0, 30.0]     # sorted, unique, > tinit

        # :tinit itself is dropped: the initial state is written by
        # :lwrite_initial, not by the schedule
        @test sphere_output_times(Dict{Symbol,Any}(:diagnostics_at_times => [5.0, 7.0]),
                                  5.0, 7.0) == [7.0]
    end

    #-----------------------------------------------------------------------
    # An empty schedule is legal and must be reported as empty, not as one
    # write at some default time.
    #-----------------------------------------------------------------------
    @testset "no schedule at all" begin
        @test isempty(sphere_output_times(Dict{Symbol,Any}(), 0.0, 100.0))
        @test isempty(sphere_output_times(Dict{Symbol,Any}(:ndiagnostics_outputs => 0),
                                          0.0, 100.0))
        # everything requested is in the past
        @test isempty(sphere_output_times(
            Dict{Symbol,Any}(:diagnostics_at_times => [1.0, 2.0]), 10.0, 100.0))
    end
end
