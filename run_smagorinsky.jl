#!/usr/bin/env julia

"""
Run Jexpresso with Smagorinsky Turbulence Model

This script demonstrates how to run the Smagorinsky implementation.

Usage:
  julia run_smagorinsky.jl

Or from Julia REPL:
  julia> include("run_smagorinsky.jl")
"""

println("="^70)
println("  Jexpresso - Spectral Element with Smagorinsky Turbulence Model")
println("="^70)
println()

# Set up the case to run
push!(empty!(ARGS), "CompEuler", "theta_smagorinsky")

# Print configuration
println("Running Configuration:")
println("  Equations: CompEuler")
println("  Case:      theta_smagorinsky")
println("  Features:  Smagorinsky LES model enabled")
println()

# Run Jexpresso
include("src/Jexpresso.jl")
