```julia
# EquationGenerator.jl

An AI-driven Julia package that automatically generates Jexpresso problem directories from PDF equation specifications using Claude AI.

## Overview

This native Julia implementation provides:
1. PDF parsing to extract equation text
2. Claude AI analysis to understand equation structure
3. Automatic Julia code generation for complete problem directories

The generated directories follow the same structure as `problems/CompEuler/thetaTracers`.

## Installation

**Recommended Method:**

```bash
cd tools/EquationGenerator.jl
julia install_deps.jl
```

This installs all required dependencies to your global Julia environment:
- `HTTP.jl` - API communication
- `JSON3.jl` - JSON parsing
- `Mustache.jl` - Template rendering
- `ArgParse.jl` - CLI argument parsing

**Why this approach?**
- Bypasses Project.toml/UUID conflicts
- Installs to global environment (works everywhere)
- Simpler and more reliable

**Alternative (Manual):**
```bash
julia -e 'using Pkg; Pkg.add(["HTTP", "JSON3", "Mustache", "ArgParse"])'
```

## Setup

Set your Anthropic API key:

```bash
export ANTHROPIC_API_KEY="your-api-key-here"
```

Or in Julia:
```julia
ENV["ANTHROPIC_API_KEY"] = "your-api-key"
```

## Usage

### As a Julia Package

**Option 1: Using init.jl (Easiest - Recommended)**
```bash
cd tools/EquationGenerator.jl
julia -i init.jl
```

This automatically:
- Checks all dependencies are installed
- Loads all source files
- Drops you into Julia REPL with functions ready to use

```julia
# Now you can use the package
julia> generate_problem("equations.pdf", category="CompEuler")

# Dry run to preview
julia> generate_problem("equations.pdf", dry_run=true)

# Save equation analysis to JSON
julia> generate_problem("equations.pdf", save_json=true)

# Generate with all options
julia> problem_dir = generate_problem(
           "path/to/equations.pdf",
           output_dir = "../../problems",
           category = "CompEuler",
           save_json = true
       )
```

### As a Command-Line Tool

```bash
# Basic usage
julia --project=. src/cli.jl equations.pdf

# With options
julia --project=. src/cli.jl equations.pdf \
    --output-dir ../../problems \
    --category CompEuler \
    --save-json

# Dry run
julia --project=. src/cli.jl equations.pdf --dry-run
```

### Command-Line Options

```
positional arguments:
  pdf_file              Path to PDF containing equations

optional arguments:
  -o, --output-dir      Output directory (default: ../../problems)
  -c, --category        Problem category (e.g., CompEuler)
  -k, --api-key         Anthropic API key
  -s, --save-json       Save equation info to JSON
  -d, --dry-run         Analyze but don't generate code
  -h, --help            Show help message
```

## PDF Format Requirements

Your PDF should clearly describe the equation system:

```
∂_t(q) + ∂_x(F) + ∂_y(G) = S

where:
- q = conserved variables [ρ, ρu, ρv, ρθ, c1, c2]
- F = x-direction flux
- G = y-direction flux
- S = source terms
```

See `examples/` for example PDF content.

## Generated Files

For each equation system, generates:

```
problems/Category/ProblemName/
├── initialize.jl          # Variable definitions & initial conditions
├── user_flux.jl          # Flux functions (F, G)
├── user_source.jl        # Source terms (S)
├── user_bc.jl            # Boundary conditions
├── user_primitives.jl    # Variable transformations
├── user_inputs.jl        # Configuration
└── README.md             # Documentation
```

All files include:
- ✅ CPU and GPU implementations
- ✅ TOTAL and PERT formulations
- ✅ NCL (non-conservative) variants
- ✅ TODO markers for customization

## API Reference

### Main Functions

```julia
generate_problem(pdf_path; kwargs...)
```
Generate complete problem directory from PDF.

**Arguments:**
- `pdf_path::String` - Path to PDF file
- `output_dir::String` - Output directory (default: "../../problems")
- `category::Union{String,Nothing}` - Problem category
- `api_key::Union{String,Nothing}` - API key
- `save_json::Bool` - Save JSON (default: false)
- `dry_run::Bool` - Preview mode (default: false)

**Returns:** Path to created directory

```julia
parse_pdf(pdf_path::String) -> String
```
Extract text from PDF file.

```julia
analyze_equations(pdf_text::String, api_key=nothing) -> NamedTuple
```
Analyze equations using Claude AI.

**Returns:** NamedTuple with:
- `problem_name::String`
- `num_equations::Int`
- `variables::Vector{String}`
- `flux_x::Vector{String}`
- `flux_y::Vector{String}`
- `source_terms::Vector{String}`
- And more...

```julia
generate_code(eq_info::NamedTuple, output_dir, category=nothing) -> String
```
Generate Julia files from equation info.

## Examples

### Example 1: Shallow Water Equations

```julia
using EquationGenerator

# Assume shallow_water.pdf contains:
# ∂h/∂t + ∂(hu)/∂x + ∂(hv)/∂y = 0
# ∂(hu)/∂t + ∂(hu² + 0.5gh²)/∂x + ∂(huv)/∂y = 0
# ∂(hv)/∂t + ∂(huv)/∂x + ∂(hv² + 0.5gh²)/∂y = 0

dir = generate_problem("shallow_water.pdf", category="ShallowWater")
# Creates: problems/ShallowWater/shallow_water/
```

### Example 2: Custom Advection-Diffusion

```julia
# With dry run first
generate_problem("advection.pdf", dry_run=true)

# Review output, then generate
dir = generate_problem("advection.pdf", category="AdvDiff")
```

### Example 3: Programmatic Usage

```julia
# Extract and analyze separately
pdf_text = parse_pdf("equations.pdf")
eq_info = analyze_equations(pdf_text)

# Inspect
println("Problem: ", eq_info.problem_name)
println("Variables: ", eq_info.variables)

# Generate if looks good
generate_code(eq_info, "../../problems", "CompEuler")
```

## Troubleshooting

### API Key Errors

```
Error: API key required
```

**Solution:**
```bash
export ANTHROPIC_API_KEY="your-key"
```

### PDF Parsing Errors

If PDFIO.jl fails to parse your PDF:
1. Try a different PDF generator
2. Simplify PDF formatting
3. Use plain text instead (convert manually)

### Package Installation Issues

```bash
# Clean and reinstall
julia --project=. -e 'using Pkg; Pkg.rm("PDFIO"); Pkg.add("PDFIO")'
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Template Rendering Errors

If Mustache.jl fails:
- Check that equation info is complete
- Verify all required fields are present
- Use dry-run to inspect extracted data

## Development

### Running Tests

```julia
using Pkg
Pkg.test("EquationGenerator")
```

### Adding Custom Templates

Edit `src/templates.jl` to customize generated code:

```julia
const MY_TEMPLATE = """
# Your custom template using Mustache syntax
# {{variable}} - interpolation
# {{#list}} ... {{/list}} - iteration
"""
```

### Extending Functionality

To add new features:
1. Add functions to appropriate module file
2. Export new functions in `EquationGenerator.jl`
3. Update documentation
4. Add tests

## Comparison with Python Version

| Feature | Python | Julia |
|---------|--------|-------|
| Speed | Fast | **Faster** |
| Integration | External | **Native** |
| Dependencies | pip | Pkg |
| PDF Parsing | pdfplumber | PDFIO |
| Templates | Jinja2 | Mustache |
| CLI | click | ArgParse |
| Type Safety | Runtime | **Compile-time** |

## Performance

Typical performance on standard equations PDF:
- PDF parsing: < 1s
- AI analysis: ~3-5s (depends on Claude API)
- Code generation: < 0.5s
- **Total: ~5-7 seconds**

## Contributing

This package is part of the Jexpresso project. To contribute:
1. Test with various equation systems
2. Report issues
3. Suggest improvements to templates
4. Add example PDFs

## License

Part of the Jexpresso project.

## Support

For issues:
1. Check this README
2. Review generated README.md in problem directory
3. Examine `src/templates.jl` for customization
4. File issue on Jexpresso repository

## See Also

- Python version: `tools/equation_generator/`
- Example reference: `problems/CompEuler/thetaTracers/`
- Jexpresso documentation: `docs/`
```
