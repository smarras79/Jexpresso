# Jexpresso Tools

This directory contains utilities and tools for the Jexpresso framework.

## Available Tools

### Equation Generator (`equation_generator/`)

An AI-driven tool that automatically generates problem directories from PDF equation specifications.

**What it does:**
- Reads a PDF containing a system of PDEs
- Uses Claude AI to extract equation structure
- Generates a complete problem directory with all necessary Julia files

**Quick Start:**
```bash
cd equation_generator
pip install -r requirements.txt
export ANTHROPIC_API_KEY="your-key"
python main.py your_equations.pdf
```

**Documentation:**
- [README.md](equation_generator/README.md) - Full documentation
- [QUICKSTART.md](equation_generator/QUICKSTART.md) - Quick start guide
- [examples/](equation_generator/examples/) - Example equation specifications

**Use Cases:**
- Quickly prototype new equation systems
- Generate boilerplate code for new problems
- Standardize problem directory structure
- Reduce manual coding errors

## Adding New Tools

To add a new tool to this directory:

1. Create a new subdirectory: `tools/your_tool_name/`
2. Add a README.md explaining what it does
3. Update this file to list the new tool
4. Follow Python best practices (requirements.txt, proper structure)

## Tool Development Guidelines

- Each tool should be self-contained
- Include comprehensive documentation
- Add examples and test cases
- Use requirements.txt for dependencies
- Keep tools focused on a single purpose
