# Jexpresso Equation Generator

An AI-driven tool that automatically generates Jexpresso problem directories from PDF equation specifications using Claude AI.

## Overview

This tool allows you to:
1. Provide a PDF file containing a system of partial differential equations (PDEs)
2. Let Claude AI analyze and extract the equation structure
3. Automatically generate a complete problem directory with all necessary Julia files

The generated directory follows the same structure as `problems/CompEuler/thetaTracers` and includes:
- `initialize.jl` - Variable definitions and initial conditions
- `user_flux.jl` - Flux function definitions (F, G)
- `user_source.jl` - Source term definitions
- `user_bc.jl` - Boundary conditions
- `user_primitives.jl` - Conservative to primitive variable transformations
- `user_inputs.jl` - Configuration parameters
- `README.md` - Documentation and TODO list

## Installation

1. Install required Python packages:
```bash
cd tools/equation_generator
pip install -r requirements.txt
```

2. Set up your Anthropic API key:
```bash
export ANTHROPIC_API_KEY="your-api-key-here"
```

## Usage

### Basic Usage

```bash
python main.py path/to/equations.pdf
```

This will:
1. Parse the PDF
2. Analyze equations with Claude AI
3. Generate a problem directory in `../../problems/`

### Command Line Options

```bash
python main.py [OPTIONS] PDF_FILE

Options:
  -o, --output-dir PATH   Output directory (default: ../../problems)
  -c, --category TEXT     Problem category (e.g., CompEuler, AdvDiff)
  -k, --api-key TEXT      Anthropic API key
  -s, --save-json         Save equation info to JSON file
  -d, --dry-run           Analyze but don't generate code
  --help                  Show help message
```

### Examples

**Example 1: Generate with specific category**
```bash
python main.py my_equations.pdf -c CompEuler -o ../../problems
```

**Example 2: Dry run to see what would be generated**
```bash
python main.py my_equations.pdf --dry-run
```

**Example 3: Save equation analysis to JSON**
```bash
python main.py my_equations.pdf --save-json
```

## PDF Format Requirements

Your PDF should clearly describe the equation system in the form:

```
∂_t(q) + ∂_x(F) + ∂_y(G) = S
```

Where:
- **q** = vector of conserved variables (e.g., [ρ, ρu, ρv, ρθ, c1, c2])
- **F** = flux in x-direction
- **G** = flux in y-direction (for 2D problems)
- **S** = source terms

### Example PDF Content

```
Compressible Euler Equations with Theta and Tracers

System of equations:

∂_t(ρ)   + ∂_x(ρu)      + ∂_y(ρv)      = 0
∂_t(ρu)  + ∂_x(ρu²+p)   + ∂_y(ρvu)     = 0
∂_t(ρv)  + ∂_x(ρuv)     + ∂_y(ρv²+p)   = -ρg
∂_t(ρθ)  + ∂_x(ρθu)     + ∂_y(ρθv)     = 0
∂_t(c1)  + ∂_x(c1u)     + ∂_y(c1v)     = 0
∂_t(c2)  + ∂_x(c2u)     + ∂_y(c2v)     = 0

Variables:
- ρ: density
- u: x-velocity
- v: y-velocity
- θ: potential temperature
- c1: tracer 1
- c2: tracer 2
- p: pressure
- g: gravitational acceleration
```

## How It Works

### 1. PDF Parser (`pdf_parser.py`)
- Extracts text from PDF using `pdfplumber`
- Handles multi-page documents
- Can extract images if needed

### 2. Equation Analyzer (`equation_analyzer.py`)
- Sends PDF text to Claude AI
- Uses structured prompts to extract:
  - Variable names
  - Flux terms (F, G)
  - Source terms (S)
  - Problem metadata
- Validates consistency (number of equations matches variables, fluxes, etc.)

### 3. Code Generator (`code_generator.py`)
- Uses Jinja2 templates to generate Julia files
- Fills in equation-specific information
- Creates problem directory structure
- Generates README with TODO items

## Generated Files

All generated files include:
- ✅ Correct equation structure
- ✅ Placeholder comments for customization
- ✅ Both CPU and GPU implementations
- ✅ TOTAL and PERT formulations
- ⚠️ TODO markers for user review

### Important: Post-Generation Steps

After generation, you **must** review and customize:

1. **initialize.jl**
   - Define initial conditions for your specific case
   - Set appropriate values for each variable

2. **user_bc.jl**
   - Implement boundary conditions
   - Currently contains placeholder logic

3. **user_flux.jl**
   - Review flux calculations
   - Add pressure/equation of state if needed
   - Verify all terms are correct

4. **user_primitives.jl**
   - Define conservative → primitive conversions
   - Important for post-processing

5. **user_inputs.jl**
   - Update with problem-specific parameters
   - Set mesh file, time stepping, etc.

## Troubleshooting

### API Key Not Found
```
Error: API key required. Set ANTHROPIC_API_KEY environment variable
```
**Solution:** Export your API key:
```bash
export ANTHROPIC_API_KEY="your-key"
```

### Equation Extraction Failed
If Claude AI fails to extract equations correctly:
1. Check your PDF formatting
2. Try the `--dry-run` flag to see what was extracted
3. Use `--save-json` to inspect the extracted data
4. Manually create/edit the JSON and use it as input

### Validation Errors
```
Error: Number of flux terms doesn't match number of equations
```
**Solution:** This usually means the PDF equations are ambiguous. Clarify the equation system in your PDF.

## Advanced Usage

### Using Saved JSON

If you've saved equation info to JSON:
```bash
# First, extract and save
python main.py equations.pdf --save-json --dry-run

# Edit the JSON file if needed
nano equations.json

# TODO: Add feature to generate from JSON
```

### Custom Templates

To customize generated code:
1. Edit templates in `templates/` directory
2. Templates use Jinja2 syntax
3. Available variables: see `code_generator.py::_prepare_context()`

## Development

### Adding New Templates

1. Create `templates/your_file.jl.j2`
2. Use Jinja2 syntax with available context variables
3. Add to `files_to_generate` in `code_generator.py`

### Testing

```bash
# Test with example PDF
python main.py examples/euler_theta_tracers.pdf --dry-run

# Validate generation
python main.py examples/euler_theta_tracers.pdf -o /tmp/test
```

## License

Part of the Jexpresso project.

## Support

For issues or questions:
1. Check this README
2. Review generated README.md in problem directory
3. Examine template files for customization options
4. File an issue on the Jexpresso repository
