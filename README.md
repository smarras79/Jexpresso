# JexpressoMeshes

This repository contains:
1. A collection of gmsh `*.geo` and `*.msh` meshes often used in Jexpresso benchmarks
2. **NEW:** `JexpressoMeshGenerator.jl` - A standalone Julia module for programmatic mesh generation

The meshes directory used to be part of the main Jexpresso repo but was growing in size, so we
decided to remove it from the Jexpresso history and have it stand alone.

## Using Pre-Generated Meshes

After cloning the JexpressoMeshes repo, follow the instructions below to use the meshes with Jexpresso:

```bash
cd Jexpresso
ln -s JexpressoMeshes/meshes .
```

## Using the Mesh Generator (NEW!)

**JexpressoMeshGenerator.jl** allows you to generate meshes directly from within Julia using GMSH.jl,
eliminating the need for external GMSH scripting.

### Quick Start

```julia
using Pkg
Pkg.add("Gmsh")

push!(LOAD_PATH, "/path/to/JexpressoMeshes/src")
using JexpressoMeshGenerator

# Generate a 3D periodic LES mesh
params = MeshParams3D(64, 64, 36, 0.0, 10240.0, 0.0, 10240.0, 0.0, 3500.0)
generate_3d_periodic_mesh(params, periodic_x=true, periodic_y=true,
                          output_file="LESICP_64x64x36.msh")
```

### Documentation

- **Quick Start:** See [QUICKSTART.md](QUICKSTART.md)
- **Full Documentation:** See [MESH_GENERATOR_README.md](MESH_GENERATOR_README.md)
- **Examples:** See `examples/basic_usage.jl` and `examples/advanced_usage.jl`
- **Tests:** Run `julia test/runtests.jl`

### Features

- 2D and 3D structured mesh generation
- Periodic boundary conditions
- Vertical stretching for boundary layer resolution
- Transfinite meshing for hexahedral elements
- Full compatibility with existing Jexpresso `.msh` format
- Programmatic mesh generation for parametric studies
