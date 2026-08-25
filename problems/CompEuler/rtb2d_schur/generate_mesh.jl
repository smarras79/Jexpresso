#!/usr/bin/env julia
#
# generate_mesh.jl -- build rtb2d_20x1x80.msh from the .geo beside it, using
#                     the gmsh SDK that GridapGmsh already ships.
#
#   julia --project=. problems/CompEuler/rtb2d_schur/generate_mesh.jl [out.msh]
#
# WHY THIS EXISTS RATHER THAN "just run gmsh"
# -------------------------------------------
# Both meshes this case uses are small enough to commit and are committed. This
# exists so a VARIANT is one command rather than an installed mesher -- but rebuilding it needs a gmsh, and
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

#
# SIZE. With no arguments this builds the 20 x 1 x 80 the deck defaults to.
# Pass three integers for a different one:
#
#   julia --project=. problems/CompEuler/rtb2d_schur/generate_mesh.jl 10 1 40
#
# and point DBG_MESH at the result.
#
# KEEP nelemy = 1. More than one element in y is no longer the 2D analogue --
# it is a thin 3D box, the bubble is still a cylinder so nothing new happens
# physically, and it costs a factor of nelemy in points for that nothing.
#
# KEEP nelemz = 4*nelemx WHILE THE DOMAIN IS 10 km x 10 km: that is what makes
# h_x/h_z = 4:1, which is the whole reason this case exists (see the .geo).
# 10 x 1 x 40 is the useful small one -- 33,005 points against 130,005, same
# aspect ratio, and it runs on a single core.
#
# The three element counts are substituted into the ONE .geo rather than kept
# in a second copy of it, so a variant cannot drift from the geometry,
# the boundary tags or the transfinite structure it is derived from.
const HERE = @__DIR__
const GEO  = joinpath(HERE, "rtb2d_20x1x80.geo")

isfile(GEO) || error("generate_mesh.jl: cannot find $GEO")

const NX, NY, NZ = length(ARGS) >= 3 ?
    (parse(Int, ARGS[1]), parse(Int, ARGS[2]), parse(Int, ARGS[3])) : (20, 1, 80)
const OUT = joinpath(HERE, "rtb2d_$(NX)x$(NY)x$(NZ).msh")

# Derive the .geo for this size. Only the three counts change; everything
# structural is whatever the committed .geo says.
const GEOSRC = let src = read(GEO, String)
    src = replace(src, r"(?m)^nelemx = \d+;" => "nelemx = $NX;")
    src = replace(src, r"(?m)^nelemy = \s*\d+;" => "nelemy = $NY;")
    src = replace(src, r"(?m)^nelemz = \d+;" => "nelemz = $NZ;")
    src
end
occursin("nelemx = $NX;", GEOSRC) && occursin("nelemy = $NY;", GEOSRC) &&
    occursin("nelemz = $NZ;", GEOSRC) ||
    error("generate_mesh.jl: could not substitute the element counts into $GEO -- ",
          "its `nelemx/y/z = N;` lines must have changed shape.")
const GEOTMP = (NX, NY, NZ) == (20, 1, 80) ? GEO :
               (t = joinpath(tempdir(), "rtb2d_$(NX)x$(NY)x$(NZ).geo");
                write(t, GEOSRC); t)

println("gmsh: reading  ", GEOTMP, "  (", NX, " x ", NY, " x ", NZ, ")")
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

    gmsh.open(GEOTMP)
    gmsh.model.mesh.generate(3)

    println("gmsh: writing  ", OUT)
    gmsh.write(OUT)

    # Report what came out, so a mesh that is quietly wrong is visible here
    # rather than three minutes into a 25-rank job.
    ntag, = gmsh.model.mesh.getElementsByType(5)      # 8-node hexahedra
    nodes, = gmsh.model.mesh.getNodes()
    println("gmsh: ", length(ntag), " hexahedra, ", length(nodes), " vertices")
    length(ntag) == NX*NY*NZ ||
        @warn "expected $(NX*NY*NZ) hexahedra, got $(length(ntag)) -- the .geo " *
              "and this check disagree about the element counts"
finally
    gmsh.finalize()
end

println("done. Point :gmsh_filename at ", OUT)
