#!/usr/bin/env julia

"""
Direct dependency installer for EquationGenerator.jl

This bypasses Project.toml issues and directly installs packages.
"""

println("=" ^ 70)
println("Installing Dependencies for EquationGenerator.jl")
println("=" ^ 70)

using Pkg

# First, ensure we're in the right directory
cd(@__DIR__)

println("\n📦 Creating fresh project environment...")

# Remove any problematic Manifest/Project artifacts
for file in ["Manifest.toml"]
    if isfile(file)
        println("  Removing old $file...")
        rm(file)
    end
end

println("\n📥 Installing packages directly...")

# Install each package explicitly
packages = ["HTTP", "JSON3", "Mustache", "ArgParse"]

for pkg in packages
    println("\n  Installing $pkg...")
    try
        Pkg.add(pkg)
        println("    ✓ $pkg installed")
    catch e
        println("    ✗ Failed: $e")
        println("\nTrying with explicit version...")
        try
            Pkg.add(name=pkg)
            println("    ✓ $pkg installed on retry")
        catch e2
            println("    ✗ Still failed: $e2")
        end
    end
end

println("\n✅ Installation attempt complete!")
println("\nCheck installed packages:")
Pkg.status()

println("\n" * "=" ^ 70)
println("Next: Test with julia -i init.jl")
println("=" ^ 70)
