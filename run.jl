#!/usr/bin/env julia
#
# Top-level entry point for Jexpresso problems.
#
# Usage:
#   julia --project=. ./run.jl <eqs> [<eqs_case>]
#
# Shortcuts (single-arg) for problems that are still stubs:
#   julia --project=. ./run.jl swe_sphere
#
# For full problems use the standard two-arg form, e.g.:
#   julia --project=. ./run.jl CompEuler theta
# which forwards to ./src/Jexpresso.jl.

const _SHORTCUTS = Dict(
    "swe_sphere" => ("ShallowWater", "swe_sphere"),
)

const _STUBS = Set(["swe_sphere"])

if length(ARGS) == 1 && haskey(_SHORTCUTS, ARGS[1])
    name = ARGS[1]
    if name in _STUBS
        println("$(name): not implemented yet")
        exit(0)
    end
    eqs, eqs_case = _SHORTCUTS[name]
    push!(empty!(ARGS), eqs, eqs_case)
end

include(joinpath(@__DIR__, "src", "Jexpresso.jl"))
