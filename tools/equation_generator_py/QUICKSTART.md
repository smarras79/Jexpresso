# Quick Start Guide

Get started with the Jexpresso Equation Generator in 5 minutes!

## Step 1: Install Dependencies

```bash
cd tools/equation_generator
pip install -r requirements.txt
```

## Step 2: Set API Key

Get your Anthropic API key from https://console.anthropic.com/

```bash
export ANTHROPIC_API_KEY="your-api-key-here"
```

Or add to your `~/.bashrc` or `~/.zshrc`:
```bash
echo 'export ANTHROPIC_API_KEY="your-api-key-here"' >> ~/.bashrc
source ~/.bashrc
```

## Step 3: Prepare Your PDF

Create a PDF that describes your equation system. See `examples/euler_theta_tracers_example.txt` for the required format.

Your PDF should include:
- The equation system in the form: ∂_t(q) + ∂_x(F) + ∂_y(G) = S
- Variable definitions
- Clear flux terms
- Source terms

## Step 4: Run the Generator

```bash
# Make the script executable (Linux/Mac)
chmod +x main.py

# Run with your PDF
python main.py path/to/your/equations.pdf
```

Example with all options:
```bash
python main.py my_equations.pdf \
  --output-dir ../../problems \
  --category CompEuler \
  --save-json
```

## Step 5: Review and Customize

The tool will create a directory like:
```
problems/CompEuler/YourProblemName/
├── initialize.jl
├── user_bc.jl
├── user_flux.jl
├── user_inputs.jl
├── user_primitives.jl
├── user_source.jl
└── README.md
```

**Important:** Review and customize:
1. ✅ `initialize.jl` - Set initial conditions
2. ✅ `user_bc.jl` - Define boundary conditions
3. ✅ `user_flux.jl` - Verify flux calculations
4. ✅ `user_primitives.jl` - Check variable conversions

## Common Issues

### "API key required"
Set your `ANTHROPIC_API_KEY` environment variable

### "PDF file not found"
Check the path to your PDF file

### "Failed to parse equations"
Your PDF might not be in the expected format. Try:
- Using `--dry-run` to see what was extracted
- Checking the example format in `examples/`
- Simplifying your PDF content

## Example: Test Run

Try a dry run first to see what the tool will generate:

```bash
# Create a simple text file with equation description
cat > test_equations.txt << 'EOF'
Shallow Water Equations

∂_t(h) + ∂_x(hu) + ∂_y(hv) = 0
∂_t(hu) + ∂_x(hu² + 0.5gh²) + ∂_y(huv) = 0
∂_t(hv) + ∂_x(huv) + ∂_y(hv² + 0.5gh²) = 0

Variables:
h - water height
u - velocity x
v - velocity y
g - gravity (9.81)
EOF

# Convert to PDF (requires tools like pandoc or wkhtmltopdf)
# Or just use a PDF editor to create a simple PDF with this content

# Run dry run
python main.py your_equations.pdf --dry-run
```

## Next Steps

- Read the full [README.md](README.md) for detailed documentation
- Check [examples/](examples/) for reference equation formats
- Customize the templates in [templates/](templates/) for your needs

## Support

For help:
1. Check error messages carefully
2. Try `--dry-run` to debug
3. Use `--save-json` to inspect extracted data
4. Review the generated `README.md` in your problem directory
