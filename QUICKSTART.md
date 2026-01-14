# Jexpresso Quick Start Guide

Get up and running with Jexpresso in minutes!

## Installation (One Command)

```bash
./install.sh
```

or

```bash
julia install_dependencies.jl
```

That's it! The installation script will:
- ✓ Check your Julia version
- ✓ Install all dependencies with correct versions
- ✓ Verify the installation
- ✓ Precompile everything

## What Gets Installed

The following packages will be installed with specific versions for maximum compatibility:

- MPI 0.20.22
- MPIPreferences 0.1.11
- PackageCompiler 2.2.1
- Thermodynamics 0.12.7
- PrettyTables 2.4.0
- Crayons 4.1.1
- UnicodePlots 3.7.2
- Gridap v0.18.12
- GridapDistributed v0.4.7
- GridapGmsh v0.7.2
- GridapP4est v0.3.11

## Running Your First Simulation

After installation:

```bash
# Activate the project
julia --project=.

# Run a test case
julia> include("examples/your_example.jl")
```

## Common Issues

### "Julia version too old"
**Solution:** Install Julia 1.11.2 or higher from https://julialang.org/downloads/

### "Package version conflict"
**Solution:** The installer automatically handles this, but if issues persist:
```bash
rm Manifest.toml
./install.sh
```

### "Precompilation failed"
**Solution:** Usually not critical. Try:
```julia
using Pkg
Pkg.precompile()
```

## Need Help?

- 📖 Full documentation: [INSTALLATION.md](INSTALLATION.md)
- 🌐 Online docs: https://smarras79.github.io/Jexpresso/dev/
- 📧 Email: smarras@njit.edu

## Next Steps

1. ✅ Installation complete? Great!
2. 📚 Read the [documentation](https://smarras79.github.io/Jexpresso/dev/)
3. 🚀 Try the example cases in `examples/`
4. 🎓 Learn about spectral element methods
5. 🤝 Join the Jexpresso community!

---

**Pro Tip:** For faster startup times on subsequent runs, consider using PackageCompiler to create a custom system image. See INSTALLATION.md for details.
