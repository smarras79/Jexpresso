# Jexpresso Installation Guide

This guide provides detailed instructions for installing Jexpresso and its dependencies with the correct package versions.

## Prerequisites

- **Julia 1.11.2 or higher** (required)
- **Git** (for cloning the repository)
- **MPI** (for parallel computing support)

## Quick Start

### 1. Clone the Repository

```bash
git clone https://github.com/smarras79/Jexpresso.git
cd Jexpresso
```

### 2. Checkout the yt/wallmodel branch if you are not already there

The `yt/wallmodel` branch is the most developed and will soon replace current master:

```bash
git checkout yt/wallmodel
```

### 3. Run the Installation Script

We provide a robust installation script that automatically installs all dependencies with the correct versions:

```bash
julia install_dependencies.jl
```

The script will:
- ✓ Check your Julia version
- ✓ Activate the project environment
- ✓ Create a backup of your current Manifest.toml
- ✓ Install all required packages with specific versions
- ✓ Resolve dependencies
- ✓ Verify installations
- ✓ Precompile packages

### 4. Verify Installation

After installation, you can verify that everything is correctly installed by checking the package versions:

```julia
using Pkg
Pkg.status()
```

## Required Package Versions

The following packages must be installed with specific versions to ensure compatibility:

| Package | Version |
|---------|---------|
| MPI | 0.20.22 |
| MPIPreferences | 0.1.11 |
| PackageCompiler | 2.2.1 |
| Thermodynamics | 0.12.7 |
| PrettyTables | 2.4.0 |
| Crayons | 4.1.1 |
| UnicodePlots | 3.7.2 |
| Gridap | 0.18.12 |
| GridapDistributed | 0.4.7 |
| GridapGmsh | 0.7.2 |
| GridapP4est | 0.3.11 |

## Troubleshooting

### Installation Script Fails

If the installation script fails, try the following:

1. **Update Julia**: Ensure you have Julia 1.11.2 or higher
   ```bash
   julia --version
   ```

2. **Manual Installation**: Install packages manually if the script fails:
   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.add(name="MPI", version="0.20.22")
   Pkg.add(name="MPIPreferences", version="0.1.11")
   # ... continue for other packages
   ```

3. **Clear Package Cache**: If you encounter persistent errors, try clearing the package cache:
   ```julia
   using Pkg
   Pkg.gc()
   ```

4. **Remove Manifest.toml**: If there are conflicts, remove the Manifest.toml and run the installation script again:
   ```bash
   rm Manifest.toml
   julia install_dependencies.jl
   ```

### Version Conflicts

If you experience version conflicts:

1. The installation script will attempt to pin packages to the correct versions automatically
2. If automatic pinning fails, you can manually pin packages:
   ```julia
   using Pkg
   Pkg.pin(name="PackageName", version="x.y.z")
   ```

### Precompilation Errors

If precompilation fails:

1. Try precompiling individual packages:
   ```julia
   using Pkg
   Pkg.precompile("PackageName")
   ```

2. Some packages may have compatibility issues with your system. Check the error messages for specific guidance.

## Advanced Installation

### Installing from Scratch

If you want a completely fresh installation:

```bash
# Remove existing Julia environment files
rm -f Manifest.toml

# Run installation script
julia install_dependencies.jl
```

### Installing on HPC Systems

For HPC systems with module systems:

```bash
# Load Julia module
module load julia/1.11.2

# Set Julia depot path (optional, for custom package location)
export JULIA_DEPOT_PATH=/path/to/your/depot

# Run installation
julia install_dependencies.jl
```

### Installing with Different MPI Implementations

If you need to use a specific MPI implementation:

```julia
using MPIPreferences
MPIPreferences.use_system_binary()  # Use system MPI
# or
MPIPreferences.use_jll_binary()     # Use Julia-provided MPI
```

Then run the installation script.

## Verification

After installation, verify that Jexpresso works correctly:

```bash
# Run a simple test case
julia --project=. test/run_tests.jl
```

## Backup and Recovery

The installation script automatically creates a backup of your `Manifest.toml` file as `Manifest.toml.backup`. If something goes wrong, you can restore it:

```bash
mv Manifest.toml.backup Manifest.toml
```

## Getting Help

If you encounter issues not covered in this guide:

1. Check the [Jexpresso documentation](https://smarras79.github.io/Jexpresso/dev/)
2. Contact the developers:
   - Simone Marras: smarras@njit.edu
   - Yassine Tissaoui: tissaoui@wisc.edu
   - Hang Wang: hang.wang@njit.edu
3. Open an issue on GitHub: https://github.com/smarras79/Jexpresso/issues

## Updating Dependencies

If package versions need to be updated in the future:

1. Edit the `REQUIRED_PACKAGES` dictionary in `install_dependencies.jl`
2. Run the installation script again
3. Test thoroughly before committing changes

## Performance Tips

After installation:

1. **Precompile packages**: This is done automatically by the script, but can be repeated:
   ```julia
   using Pkg
   Pkg.precompile()
   ```

2. **Use PackageCompiler**: For frequently-run scripts, consider using PackageCompiler to create a custom system image:
   ```julia
   using PackageCompiler
   create_sysimage(["Jexpresso"], sysimage_path="jexpresso.so")
   ```

3. **Enable threading**: Run Julia with multiple threads for better performance:
   ```bash
   julia --project=. --threads=auto
   ```


## Setup and Run with MPI

JEXPRESSO supports parallel execution using either OpenMPI or MPICH. Follow these steps to configure and run with your preferred MPI implementation.

### 1. Install MPI Implementation

Choose either OpenMPI or MPICH:

#### OpenMPI Installation
```bash
# Ubuntu/Debian
sudo apt install libopenmpi-dev openmpi-bin

# macOS (Homebrew)
brew install open-mpi

# Verify installation
mpiexec --version
```

#### MPICH Installation
```bash
# Ubuntu/Debian
sudo apt install mpich libmpich-dev

# macOS (Homebrew) 
brew install mpich

# Verify installation
mpiexec --version
```

### 2. Configure MPI Preferences

#### Automatic Configuration (Default Path)
Use this command when MPI (OpenMPI/MPICH) is installed in standard system paths (`/usr/bin`, `/usr/local/bin`, etc.):
```bash
julia --project=. -e 'using Pkg; Pkg.add("MPIPreferences"); using MPIPreferences; MPIPreferences.use_system_binary()'
```

#### Manual Configuration (For Multiple MPI Installations or MPI not in Default Path)
For MPI installations in non-standard locations (e.g., /opt/openmpi, or custom paths):
```bash
julia --project=. -e 'using Pkg; Pkg.add("MPIPreferences"); using MPIPreferences; MPIPreferences.use_system_binary(;extra_paths = ["/where/your/mpi/lib"])'
```
If MPI is installed via homebrew on macOS, the MPI lib path is:
```bash
/opt/homebrew/lib
```

### 3. Running with MPI

#### Basic Execution
```bash
mpiexec -n <NPROCS> julia --project=. -e 'push!(empty!(ARGS), "<EQUATIONS>", "<CASE_NAME>"); include("./src/Jexpresso.jl")'
```

#### Implementation-Specific Examples
```bash
mpiexec -n 4 julia --project=. -e 'push!(empty!(ARGS), "CompEuler", "3d"); include("./src/Jexpresso.jl")'
```
#### Script
You can simplify the run steps with a `runjexpresso` script like this:
```bash
#!/bin/bash

MPIRUN=/YOUR/PATH/TO/mpirun
JULIA=/YOUR/PATH/TO/julia

$MPIRUN -np $1 $JULIA --project=. -e 'push!(empty!(ARGS), "'"$2"'", "'"$3"'"); include("./src/Jexpresso.jl")' "$@"
```
and run it like this:
```bash
./runjexpresso 4 CompEuler theta
```

### Troubleshooting

- **Library conflicts:** Clear existing preferences:
  ```bash
  rm -f LocalPreferences.toml
  ```
- **Path issues:** Verify paths with:
  ```bash
  which mpiexec
  which mpirun
  ```
  You may have to use the full aboslute path to mpiexec or mpirun and to julia like this if necessary:
  ```
  /opt/homebrew/Cellar/open-mpi/5.0.6/bin/mpirun -n 4 /Applications/Julia-1.11.app/Contents/Resources/julia/bin/julia --project=. -e 'push!(empty!(ARGS), "CompEuler", "theta"); include("./src/Jexpresso.jl")'
  ```


- **Version mismatches:** Ensure consistent versions:
  ```bash
  mpicc --version
  mpif90 --version
  ```

<!-- <img src="assets/mpi_performance_comparison.png" width="700" alt="MPI Performance Comparison"> -->

## Plotting
Files can be written to VTK (recommended) or png (png is now only used for 1D results). For the png plots, we use [Makie](https://github.com/MakieOrg/Makie.jl). If you want to use a different package,
modify ./src/io/plotting/jplots.jl accordinly.

## Contacts
[Simone Marras](mailto:smarras@njit.edu), [Yassine Tissaoui](mailto:tissaoui@wisc.edu), [Hang Wang](mailto:hang.wang@njit.edu)

## License

Jexpresso is open-source software. Please see the LICENSE file for details.

## Citation

If you use Jexpresso, please cite:

```bibtex
@article{tissaoui2024,
  author = {Y. Tissaoui and J. F. Kelly and S. Marras},
  title = {Efficient Spectral Element Method for the Euler Equations on Unbounded Domains},
  volume = {487},
  pages = {129080},
  year = {2024},
  journal = {App. Math. Comput.},
}
```

