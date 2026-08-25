#!/usr/bin/env julia
#
# generate_mesh.jl -- build rtb3d_20x20x80.msh from the .geo beside it, using
#                     the gmsh SDK that GridapGmsh already ships.
#
#   julia --project=. problems/CompEuler/rtb3d_schur/generate_mesh.jl [out.msh]
#
# WHY THIS EXISTS RATHER THAN "just run gmsh"
# -------------------------------------------
# The mesh is 2.1M gridpoints and 32000 elements, which is too big to keep in
# the repository and trivial to rebuild -- but rebuilding it needs a gmsh, and
# `gmsh` is often not on PATH even on a machine that can run Jexpresso. It does
# not have to be: GridapGmsh is already a dependency and it ships the SDK with
# its Julia API, so the mesher is present in any environment where the code
# runs at all. This drives that one.
#
# It writes format 4.1, which is what every other .msh in problems/ is and what
# `tools/inspect_msh.jl` and `tools/reorder_msh_columns.jl` parse. 2.2 also
# loads fine through GridapGmsh and is far easier to read by eye -- a flat node
# list and a flat element list rather than per-entity blocks -- but those two
# tools reject it (2.2's `4 -5000 -5000 10000` node line trips the entity-block
# parser), and keeping the repository's own tooling working on this mesh is
# worth more than being able to `head` it.
using GridapGmsh: gmsh

const HERE = @__DIR__
const GEO  = joinpath(HERE, "rtb3d_20x20x80.geo")
const OUT  = length(ARGS) >= 1 ? ARGS[1] : joinpath(HERE, "rtb3d_20x20x80.msh")

isfile(GEO) || error("generate_mesh.jl: cannot find $GEO")

println("gmsh: reading  ", GEO)
gmsh.initialize()
try
    # Quiet, but not silent: warnings about the geometry are worth seeing.
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("General.Verbosity", 3)
    # 4.1, for the reason in the header.
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    # The .geo already asks for Transfinite + Recombine everywhere, so this is
    # a structured hex mesh and the 3D algorithm has nothing to decide. Saying
    # so explicitly keeps a future gmsh from picking a tet algorithm and
    # silently producing a mesh Jexpresso's hex-only reader cannot use.
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)

    gmsh.open(GEO)
    gmsh.model.mesh.generate(3)

    println("gmsh: writing  ", OUT)
    gmsh.write(OUT)

    # Report what came out, so a mesh that is quietly wrong is visible here
    # rather than three minutes into a 25-rank job.
    ntag, = gmsh.model.mesh.getElementsByType(5)      # 8-node hexahedra
    nodes, = gmsh.model.mesh.getNodes()
    println("gmsh: ", length(ntag), " hexahedra, ", length(nodes), " vertices")
    length(ntag) == 20*20*80 ||
        @warn "expected $(20*20*80) hexahedra, got $(length(ntag)) -- the .geo " *
              "and this check disagree about the element counts"
finally
    gmsh.finalize()
end

println("done. Point :gmsh_filename at ", OUT)
