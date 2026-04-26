# Getting Started with Jexpresso

Jexpresso is a Julia-based solver for partial differential equations (PDEs) in atmospheric science and beyond.

---

## Requirements

Before you begin, make sure you have:

- **Julia 1.12.5** (or 1.11.2+) — [download here](https://julialang.org/downloads/)
- **Git** — [download here](https://git-scm.com/)
- **MPI** (for parallel runs) — OpenMPI or MPICH

### Install MPI (if needed)

**Linux (Ubuntu/Debian):**
```bash
sudo apt install libopenmpi-dev openmpi-bin
```

**macOS (Homebrew):**
```bash
brew install open-mpi
```

---

## Step 1 — Download

Clone the repository:

```bash
git clone https://github.com/smarras79/Jexpresso.git
cd Jexpresso
```

Then switch to the recommended branch:

```bash
git checkout sm/newmaster
```

---

## Step 2 — Install

Run these three commands **in order**:

```bash
# 1. Download all dependencies (without precompiling yet)
julia --project=. -e 'ENV["JULIA_PKG_PRECOMPILE_AUTO"]=0; using Pkg; Pkg.instantiate()'

# 2. Configure MPI to use your system's installation
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary()'

# 3. Precompile everything
julia --project=. -e 'using Pkg; Pkg.precompile()'
```

> **Note:** Precompilation can take several minutes the first time.

---

## Step 3 — Run

Jexpresso runs are defined by two arguments: the **equation set** and the **case name**.

### Single-process run

Open a Julia session from the `Jexpresso/` folder:

```bash
julia --project=.
```

Then inside Julia:

```julia
push!(empty!(ARGS), "CompEuler", "theta")
include("./src/Jexpresso.jl")
```

### Parallel run (MPI)

```bash
mpiexec -n 4 julia --project=. -e 'push!(empty!(ARGS), "CompEuler", "theta"); include("./src/Jexpresso.jl")'
```

Replace `4` with the number of processes you want to use.

---

## Available Equation Sets and Cases

| Equation Set | Example Cases |
|---|---|
| `CompEuler` | `theta`, `thetaTracers`, `3d` |
| `AdvDiff` | `Wave_Train` |
| `ShallowWater` | (see `problems/ShallowWater/`) |
| `Burgers` | (see `problems/Burgers/`) |
| `Helmholtz` | `case1_laguerre` |
| `Elliptic` | (see `problems/Elliptic/`) |

All available cases live under the `problems/` directory.

---

## Verify Your Installation

Run the built-in test suite:

```bash
julia --project=. test/runtests.jl
```

---

## Output

Results are written in **VTK format** (recommended, open with [ParaView](https://www.paraview.org/)). 1D results can also be written as PNG images.

---

## Troubleshooting

**Installation fails or produces conflicts:**
```bash
rm Manifest.toml
julia --project=. -e 'ENV["JULIA_PKG_PRECOMPILE_AUTO"]=0; using Pkg; Pkg.instantiate()'
```

**MPI not found at a standard path (e.g., macOS Homebrew):**
```bash
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary(extra_paths=["/opt/homebrew/lib"])'
```

**Clear package cache:**
```julia
using Pkg
Pkg.gc()
```

---

## Help and Documentation

- Full docs: https://smarras79.github.io/Jexpresso/dev/
- Issues: https://github.com/smarras79/Jexpresso/issues
- Contact: [Simone Marras](mailto:smarras@njit.edu), [Yassine Tissaoui](mailto:tissaoui@wisc.edu), [Hang Wang](mailto:hang.wang@njit.edu)
