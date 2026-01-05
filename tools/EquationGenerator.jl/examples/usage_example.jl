"""
Example usage of EquationGenerator.jl

This file demonstrates different ways to use the package.
"""

using EquationGenerator

# Example 1: Basic usage
println("=" ^ 70)
println("Example 1: Basic PDF Processing")
println("=" ^ 70)

# Generate problem from PDF (replace with actual PDF path)
pdf_path = "my_equations.pdf"

if isfile(pdf_path)
    problem_dir = generate_problem(pdf_path)
    println("✓ Created problem directory: $problem_dir")
else
    println("⚠ Create a PDF file at: $pdf_path")
    println("  See examples/euler_theta_tracers_example.txt for format")
end

# Example 2: Dry run to preview
println("\n" * "=" ^ 70)
println("Example 2: Dry Run (Preview)")
println("=" ^ 70)

if isfile(pdf_path)
    generate_problem(pdf_path, dry_run=true)
else
    println("⚠ Need PDF file for dry run")
end

# Example 3: With all options
println("\n" * "=" ^ 70)
println("Example 3: Full Options")
println("=" ^ 70)

if isfile(pdf_path)
    problem_dir = generate_problem(
        pdf_path,
        output_dir = "../../problems",
        category = "CompEuler",
        save_json = true,
        dry_run = false
    )
    println("✓ Generated with all options")
else
    println("⚠ Need PDF file")
end

# Example 4: Step-by-step process
println("\n" * "=" ^ 70)
println("Example 4: Step-by-Step")
println("=" ^ 70)

if isfile(pdf_path)
    # Step 1: Parse PDF
    println("Step 1: Parsing PDF...")
    pdf_text = parse_pdf(pdf_path)
    println("  ✓ Extracted $(length(pdf_text)) characters")

    # Step 2: Analyze equations
    println("\nStep 2: Analyzing equations with Claude AI...")
    eq_info = analyze_equations(pdf_text)
    println("  ✓ Problem: $(eq_info.problem_name)")
    println("  ✓ Equations: $(eq_info.num_equations)")
    println("  ✓ Variables: $(join(eq_info.variables, ", "))")

    # Step 3: Validate
    println("\nStep 3: Validating...")
    validate_equation_info(eq_info)
    println("  ✓ Validation passed")

    # Step 4: Generate code
    println("\nStep 4: Generating code...")
    output_dir = generate_code(eq_info, "../../problems", "CompEuler")
    println("  ✓ Created: $output_dir")
else
    println("⚠ Need PDF file for step-by-step demo")
end

# Example 5: Error handling
println("\n" * "=" ^ 70)
println("Example 5: Error Handling")
println("=" ^ 70)

try
    generate_problem("nonexistent.pdf")
catch e
    println("✓ Properly caught error: ", typeof(e))
    println("  Message: ", e.msg)
end

# Example 6: Testing API key
println("\n" * "=" ^ 70)
println("Example 6: API Key Test")
println("=" ^ 70)

if test_api_key()
    println("✓ API key is valid and working")
else
    println("✗ API key issue - check ANTHROPIC_API_KEY environment variable")
end

println("\n" * "=" ^ 70)
println("Examples Complete!")
println("=" ^ 70)
println("\nNext steps:")
println("1. Create a PDF with your equation system")
println("2. Run: generate_problem(\"your_equations.pdf\")")
println("3. Review and customize generated files")
println("\nSee QUICKSTART.md for detailed guide.")
