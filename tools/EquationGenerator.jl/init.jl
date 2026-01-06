#!/usr/bin/env julia

"""
Load EquationGenerator from source files directly
This bypasses package loading issues
"""

println("=" ^ 70)
println("Loading EquationGenerator.jl")
println("=" ^ 70)

# Change to package directory
cd(@__DIR__)

# Ensure packages are available
println("\n📦 Checking dependencies...")
required_packages = ["HTTP", "JSON3", "Mustache", "ArgParse"]

for pkg_name in required_packages
    try
        eval(Meta.parse("using $pkg_name"))
        println("  ✓ $pkg_name")
    catch
        println("  ✗ $pkg_name not found!")
        println("\n⚠ Please run: julia install_deps.jl")
        exit(1)
    end
end

# Load source files directly (bypassing module system)
println("\n📂 Loading source files...")

include("src/pdf_parser.jl")
println("  ✓ pdf_parser.jl")

include("src/claude_client.jl")
println("  ✓ claude_client.jl")

include("src/equation_analyzer.jl")
println("  ✓ equation_analyzer.jl")

include("src/templates.jl")
println("  ✓ templates.jl")

include("src/code_generator.jl")
println("  ✓ code_generator.jl")

# Load main module functions
println("\n📚 Loading EquationGenerator functions...")
include("src/EquationGenerator.jl")
using .EquationGenerator

println("\n✅ EquationGenerator loaded successfully!")
println("\n" * "=" ^ 70)
println("Available functions:")
println("  - generate_problem(pdf_path; kwargs...)")
println("  - parse_pdf(pdf_path)")
println("  - analyze_equations(pdf_text, api_key)")
println("  - generate_code(eq_info, output_dir, category)")
println("\nExample:")
println("  generate_problem(\"./examples/euler3d.pdf\", category=\"CompEuler\")")
println("=" ^ 70)
