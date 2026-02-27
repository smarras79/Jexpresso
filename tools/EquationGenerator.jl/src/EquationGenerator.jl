"""
# EquationGenerator

AI-driven tool to generate Jexpresso problem directories from PDF equation specifications.

This Julia package provides functionality to:
- Parse PDF files containing equation specifications
- Analyze equations using Claude AI
- Generate complete Jexpresso problem directories with Julia code

## Usage

```julia
using EquationGenerator

# Generate from PDF
generate_problem("path/to/equations.pdf", output_dir="problems")

# Or use the CLI
# julia --project=. src/cli.jl equations.pdf
```
"""
module EquationGenerator

using HTTP
using JSON3
using Mustache
using PDFIO

export generate_problem, parse_pdf, analyze_equations, generate_code

include("pdf_parser.jl")
include("claude_client.jl")
include("equation_analyzer.jl")
include("code_generator.jl")
include("templates.jl")

"""
    generate_problem(pdf_path::String; kwargs...)

Main function to generate a complete problem directory from a PDF.

# Arguments
- `pdf_path::String`: Path to the PDF file containing equations
- `output_dir::String="../../problems"`: Output directory
- `category::Union{String,Nothing}=nothing`: Problem category (e.g., "CompEuler")
- `api_key::Union{String,Nothing}=nothing`: Anthropic API key (defaults to ENV["ANTHROPIC_API_KEY"])
- `save_json::Bool=false`: Save extracted equation info to JSON
- `dry_run::Bool=false`: Analyze but don't generate code

# Returns
- Path to created problem directory (or nothing if dry_run=true)

# Example
```julia
prob_dir = generate_problem("my_equations.pdf", category="CompEuler")
```
"""
function generate_problem(
    pdf_path::String;
    output_dir::String="../../problems",
    category::Union{String,Nothing}=nothing,
    api_key::Union{String,Nothing}=nothing,
    save_json::Bool=false,
    dry_run::Bool=false
)
    println("=" ^ 70)
    println("Jexpresso Equation Generator (Julia)")
    println("=" ^ 70)

    # Step 1: Parse PDF
    println("\n📄 Parsing PDF: $pdf_path")
    pdf_text = parse_pdf(pdf_path)
    println("✓ Extracted text from PDF")

    # Step 2: Analyze equations
    println("\n🤖 Analyzing equations with Claude AI...")
    eq_info = analyze_equations(pdf_text, api_key)
    validate_equation_info(eq_info)

    println("✓ Detected equation system: $(eq_info.problem_name)")
    println("  - $(eq_info.num_equations) equations")
    println("  - $(get(eq_info, :spatial_dimensions, 2))D problem")
    println("  - Variables: $(join(eq_info.variables, ", "))")

    # Save JSON if requested
    if save_json
        json_file = replace(pdf_path, r"\.pdf$" => ".json")
        open(json_file, "w") do f
            JSON3.pretty(f, eq_info)
        end
        println("\n💾 Saved equation info to: $json_file")
    end

    # Stop if dry run
    if dry_run
        println("\n🏁 Dry run complete (no files generated)")
        print_equation_info(eq_info)
        return nothing
    end

    # Step 3: Generate code
    println("\n📝 Generating Julia code...")
    problem_dir = generate_code(eq_info, output_dir, category)

    println("\n✅ Problem directory created: $problem_dir")
    println("\n📋 Generated files:")
    for file in sort(readdir(problem_dir))
        if endswith(file, ".jl") || file == "README.md"
            println("  - $file")
        end
    end

    print_summary(eq_info, problem_dir)

    return problem_dir
end

"""
    validate_equation_info(eq_info)

Validate that equation information is consistent.
"""
function validate_equation_info(eq_info)
    num_eqs = eq_info.num_equations
    num_vars = length(eq_info.variables)
    num_flux_x = length(eq_info.flux_x)
    num_source = length(eq_info.source_terms)

    num_vars == num_eqs || error("Number of variables ($num_vars) doesn't match number of equations ($num_eqs)")
    num_flux_x == num_eqs || error("Number of x-flux terms ($num_flux_x) doesn't match number of equations ($num_eqs)")
    num_source == num_eqs || error("Number of source terms ($num_source) doesn't match number of equations ($num_eqs)")

    spatial_dim = get(eq_info, :spatial_dimensions, 2)
    if spatial_dim >= 2
        flux_y = get(eq_info, :flux_y, [])
        num_flux_y = length(flux_y)
        num_flux_y == num_eqs || error("Number of y-flux terms ($num_flux_y) doesn't match number of equations ($num_eqs)")
    end

    return true
end

"""Print detailed equation information."""
function print_equation_info(eq_info)
    println("\n" * "=" ^ 70)
    println("Equation Information")
    println("=" ^ 70)
    println("\nProblem: $(eq_info.problem_name)")
    println("Description: $(get(eq_info, :description, "N/A"))")
    println("\nEquations: $(eq_info.num_equations)")
    println("Dimensions: $(get(eq_info, :spatial_dimensions, 2))D")

    println("\nVariables:")
    for (i, var) in enumerate(eq_info.variables)
        desc = get(get(eq_info, :variable_descriptions, Dict()), Symbol(var), "")
        println("  $i. $var" * (isempty(desc) ? "" : " - $desc"))
    end

    println("\nFlux (X-direction):")
    for (i, flux) in enumerate(eq_info.flux_x)
        println("  $i. $flux")
    end

    if haskey(eq_info, :flux_y) && !isempty(eq_info.flux_y)
        println("\nFlux (Y-direction):")
        for (i, flux) in enumerate(eq_info.flux_y)
            println("  $i. $flux")
        end
    end

    println("\nSource Terms:")
    for (i, source) in enumerate(eq_info.source_terms)
        println("  $i. $source")
    end
end

"""Print summary and next steps."""
function print_summary(eq_info, problem_dir)
    println("\n" * "=" ^ 70)
    println("Next Steps")
    println("=" ^ 70)
    println("\n1. Review the generated files, especially:")
    println("   - $problem_dir/initialize.jl (initial conditions)")
    println("   - $problem_dir/user_bc.jl (boundary conditions)")
    println("   - $problem_dir/user_flux.jl (flux calculations)")
    println("\n2. Update user_inputs.jl with problem-specific parameters")
    println("\n3. Test with your mesh and verify results")
    println("\n4. See README.md for detailed TODO list")
    println("\n" * "=" ^ 70)
end

end # module
