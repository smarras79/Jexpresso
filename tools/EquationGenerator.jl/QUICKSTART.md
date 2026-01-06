# Quick Start Guide - EquationGenerator.jl

Get started in 5 minutes with the Julia version!

## Step 1: Install Dependencies

```bash
cd tools/EquationGenerator.jl
julia install_deps.jl
```

This installs all dependencies (HTTP, JSON3, Mustache, ArgParse) to your global Julia environment.

## Step 2: Set API Key

```bash
export ANTHROPIC_API_KEY="your-api-key-here"
```

Get your key from: https://console.anthropic.com/

## Step 3: Prepare Your PDF

Create a PDF describing your equation system. Must include:
- Equations in form: ∂_t(q) + ∂_x(F) + ∂_y(G) = S
- Variable definitions
- Clear flux and source terms

See `examples/` for format.

## Step 4: Run Generator

### Option A: Using init.jl (Easiest - Recommended)

```bash
julia -i init.jl

# In Julia REPL (after loading completes):
julia> generate_problem("equations.pdf")
```

**What happens:**
- Checks all dependencies are installed
- Loads source files directly
- Drops you into Julia REPL ready to use

### Option B: Command Line

```bash
julia --project=. src/cli.jl equations.pdf
```

### Option C: With Options

```bash
julia --project=. src/cli.jl equations.pdf \
    --output-dir ../../problems \
    --category CompEuler \
    --save-json
```

## Step 5: Review and Customize

Generated files will be in:
```
problems/Category/YourProblem/
├── initialize.jl
├── user_bc.jl
├── user_flux.jl
├── user_inputs.jl
├── user_primitives.jl
├── user_source.jl
└── README.md
```

**Customize:**
1. ✅ `initialize.jl` - Set initial conditions
2. ✅ `user_bc.jl` - Define boundary conditions
3. ✅ `user_flux.jl` - Verify flux calculations

## Example: Complete Workflow

```bash
# 1. Install dependencies
julia install_deps.jl

# 2. Load EquationGenerator
julia -i init.jl
```

```julia
# 3. First, do a dry run
julia> generate_problem("my_equations.pdf", dry_run=true)

# 4. Review output, then generate
julia> dir = generate_problem(
           "my_equations.pdf",
           output_dir="../../problems",
           category="CompEuler",
           save_json=true
       )

# 5. Directory created!
julia> println("Created: ", dir)
```

## Common Commands

```bash
# Dry run (preview only)
julia --project=. src/cli.jl equations.pdf --dry-run

# Generate with category
julia --project=. src/cli.jl equations.pdf -c CompEuler

# Save analysis to JSON
julia --project=. src/cli.jl equations.pdf --save-json

# Custom output directory
julia --project=. src/cli.jl equations.pdf -o /path/to/output

# Help
julia --project=. src/cli.jl --help
```

## Tips

1. **Test first:** Always use `--dry-run` to preview
2. **Save JSON:** Use `--save-json` to inspect extracted data
3. **Check API key:** Make sure `ANTHROPIC_API_KEY` is set
4. **Read README:** Generated README.md has problem-specific TODOs

## Troubleshooting

### "API key required"
```bash
export ANTHROPIC_API_KEY="your-key"
```

### "Package not found" during init.jl
```bash
# Re-run dependency installation
julia install_deps.jl
```

### "PDF parsing failed"
- Ensure PDF is text-based (not scanned image)
- Try simplifying PDF content
- Check PDF file permissions

### Installation issues
If `install_deps.jl` fails:
```bash
# Install packages manually
julia -e 'using Pkg; Pkg.add(["HTTP", "JSON3", "Mustache", "ArgParse"])'
```

## Next Steps

- Read full [README.md](README.md)
- Check [examples/](examples/) for reference PDFs
- Review generated code in created directory
- Test with your mesh files

## Performance Note

Julia version is typically **2-3x faster** than Python version after first run (due to JIT compilation).

First run: ~10-15 seconds (includes compilation)
Subsequent runs: ~5-7 seconds

## Support

For help:
1. Check error messages
2. Try `--dry-run` to debug
3. Review [README.md](README.md)
4. Check Jexpresso documentation
