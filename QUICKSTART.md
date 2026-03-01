# Quick Start Guide - JexpressoMeshGenerator

Get started with generating meshes for Jexpresso in 5 minutes!

## Step 1: Install Dependencies

```julia
using Pkg
Pkg.add("Gmsh")
```

## Step 2: Load the Module

```julia
# Add the src directory to your Julia path
push!(LOAD_PATH, "/path/to/JexpressoMeshes/src")
using JexpressoMeshGenerator
```

Or simply include it directly:

```julia
include("/path/to/JexpressoMeshes/src/JexpressoMeshGenerator.jl")
using .JexpressoMeshGenerator
```

## Step 3: Generate Your First Mesh

### Example 1: Simple 2D Mesh

```julia
# Create a 20x20 mesh on a unit square
params = MeshParams2D(20, 20, 0.0, 1.0, 0.0, 1.0)
generate_2d_structured_mesh(params, output_file="my_first_mesh_2d.msh")
```

### Example 2: Simple 3D Mesh

```julia
# Create a 10x10x10 mesh on a unit cube
params = MeshParams3D(10, 10, 10, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
generate_3d_structured_mesh(params, output_file="my_first_mesh_3d.msh")
```

### Example 3: Periodic LES Mesh

```julia
# Create a Large Eddy Simulation mesh with periodic boundaries
# 10 km x 10 km x 3.5 km domain, 64x64x36 elements
params = MeshParams3D(
    64,      # elements in x
    64,      # elements in y
    36,      # elements in z
    0.0,     # xmin (meters)
    10240.0, # xmax (10.24 km)
    0.0,     # ymin
    10240.0, # ymax (10.24 km)
    0.0,     # zmin (surface)
    3500.0   # zmax (3.5 km height)
)

generate_3d_periodic_mesh(
    params,
    periodic_x = true,
    periodic_y = true,
    periodic_z = false,
    output_file = "LESICP_64x64x36.msh"
)
```

## Step 4: Use the Mesh in Jexpresso

Now use your generated `.msh` file in your Jexpresso simulations:

```julia
# In your Jexpresso simulation
mesh = read_mesh("LESICP_64x64x36.msh")
# ... continue with your simulation setup
```

## Common Mesh Types

### Atmospheric Boundary Layer
```julia
params = MeshParams3D(48, 48, 50, 0.0, 5000.0, 0.0, 5000.0, 0.0, 2000.0)
generate_3d_periodic_mesh(params, periodic_x=true, periodic_y=true,
                          output_file="ABL_mesh.msh")
```

### Channel Flow
```julia
params = MeshParams3D(100, 40, 80, 0.0, 12.56, -1.0, 1.0, 0.0, 6.28)
create_stretched_mesh_3d(params, 1.08, output_file="channel.msh")
```

### Convection Cell
```julia
params = MeshParams2D(40, 40, 0.0, 1.0, 0.0, 1.0)
generate_2d_periodic_mesh(params, periodic_x=true, periodic_y=false,
                          output_file="convection.msh")
```

## Next Steps

- Read the full documentation: `MESH_GENERATOR_README.md`
- Explore examples: `examples/basic_usage.jl` and `examples/advanced_usage.jl`
- Customize boundary conditions and mesh parameters for your specific case

## Troubleshooting

**"Gmsh not found"**
```julia
using Pkg
Pkg.add("Gmsh")
```

**"Module not found"**
Make sure the path in `push!(LOAD_PATH, ...)` is correct and points to the `src` directory.

**"Mesh file not created"**
Check that you have write permissions in the output directory.

## Need Help?

- Check `MESH_GENERATOR_README.md` for detailed documentation
- Review examples in `examples/` directory
- Consult Jexpresso.jl documentation for simulation setup

Happy meshing! 🎉
