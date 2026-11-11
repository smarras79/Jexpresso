#!/usr/bin/env python3
"""
Regenerate the gmsh meshes used by the Elasticity cases.

    pip install gmsh && python3 problems/Elasticity/make_meshes.py

Both meshes are structured quadrilateral grids (transfinite + recombine),
which is what Jexpresso's spectral elements want: one gmsh quad becomes one
spectral element carrying (nop+1)^2 LGL nodes.

Physical group names are the strings the solver dispatches on:

  * the surface is always "domain";
  * "periodicx" / "periodicy" are the tags mod_mesh_read_gmsh! recognises for
    periodic directions (in 2D "periodicy" is remapped to "periodicz"
    internally, and BCs.jl then skips those edges entirely);
  * any other curve name is handed to user_bc_dirichlet! as `tag`.
"""
import os
import gmsh

HERE = os.path.dirname(os.path.abspath(__file__))


def structured_rectangle(x0, x1, y0, y1, nx, ny, curve_groups, path, name):
    """One rectangle meshed into nx*ny quads.

    curve_groups maps a physical name to the sides it covers, sides being
    any of "bottom", "right", "top", "left".
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(name)

    p = [gmsh.model.geo.addPoint(x, y, 0.0)
         for x, y in ((x0, y0), (x1, y0), (x1, y1), (x0, y1))]
    sides = {
        "bottom": gmsh.model.geo.addLine(p[0], p[1]),
        "right":  gmsh.model.geo.addLine(p[1], p[2]),
        "top":    gmsh.model.geo.addLine(p[2], p[3]),
        "left":   gmsh.model.geo.addLine(p[3], p[0]),
    }
    loop = gmsh.model.geo.addCurveLoop([sides["bottom"], sides["right"],
                                        sides["top"], sides["left"]])
    surf = gmsh.model.geo.addPlaneSurface([loop])

    # Structured: nx+1 nodes along the horizontal sides, ny+1 along the
    # vertical ones, then recombine the triangles into quads.
    for s in ("bottom", "top"):
        gmsh.model.geo.mesh.setTransfiniteCurve(sides[s], nx + 1)
    for s in ("right", "left"):
        gmsh.model.geo.mesh.setTransfiniteCurve(sides[s], ny + 1)
    gmsh.model.geo.mesh.setTransfiniteSurface(surf)
    gmsh.model.geo.mesh.setRecombine(2, surf)

    gmsh.model.geo.synchronize()

    for gname, members in curve_groups.items():
        tag = gmsh.model.addPhysicalGroup(1, [sides[s] for s in members])
        gmsh.model.setPhysicalName(1, tag, gname)
    tag = gmsh.model.addPhysicalGroup(2, [surf])
    gmsh.model.setPhysicalName(2, tag, "domain")

    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(path)

    nnode = len(gmsh.model.mesh.getNodes()[0])
    gmsh.finalize()
    print(f"  {os.path.relpath(path, HERE)}: {nx}x{ny} quads, {nnode} vertices")


def structured_quarter_annulus(a, b, nr, nt, path, name):
    """The quarter annulus a <= r <= b, 0 <= theta <= pi/2, meshed into
    nr (radial) x nt (circumferential) quads.

    The two arcs are gmsh `Circle`s, but gmsh still writes a STRAIGHT-SIDED
    grid: the high-order nodes Jexpresso adds would sit on the chords. The case
    curves them back onto the circles with :exact_geometry, which is why the
    two arcs carry their own physical names.

    Boundary groups:

      "inner"   r = a          the loaded surface (internal pressure)
      "outer"   r = b          traction-free
      "xaxis"   y = 0, a<=x<=b symmetry: u_y = 0
      "yaxis"   x = 0, a<=y<=b symmetry: u_x = 0
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(name)

    geo = gmsh.model.geo
    o = geo.addPoint(0.0, 0.0, 0.0)
    pa = geo.addPoint(a, 0.0, 0.0)      # A: inner radius on the x axis
    pb = geo.addPoint(b, 0.0, 0.0)      # B: outer radius on the x axis
    pc = geo.addPoint(0.0, b, 0.0)
    pd = geo.addPoint(0.0, a, 0.0)

    xaxis = geo.addLine(pa, pb)
    outer = geo.addCircleArc(pb, o, pc)
    yaxis = geo.addLine(pc, pd)
    inner = geo.addCircleArc(pd, o, pa)

    loop = geo.addCurveLoop([xaxis, outer, yaxis, inner])
    surf = geo.addPlaneSurface([loop])

    for c in (xaxis, yaxis):
        geo.mesh.setTransfiniteCurve(c, nr + 1)
    for c in (outer, inner):
        geo.mesh.setTransfiniteCurve(c, nt + 1)
    geo.mesh.setTransfiniteSurface(surf)
    geo.mesh.setRecombine(2, surf)

    geo.synchronize()

    for gname, curve in (("xaxis", xaxis), ("outer", outer),
                         ("yaxis", yaxis), ("inner", inner)):
        tag = gmsh.model.addPhysicalGroup(1, [curve])
        gmsh.model.setPhysicalName(1, tag, gname)
    tag = gmsh.model.addPhysicalGroup(2, [surf])
    gmsh.model.setPhysicalName(2, tag, "domain")

    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(path)

    nnode = len(gmsh.model.mesh.getNodes()[0])
    gmsh.finalize()
    print(f"  {os.path.relpath(path, HERE)}: {nr}x{nt} quads, {nnode} vertices")


if __name__ == "__main__":
    print("Elasticity meshes:")

    # plane_wave: doubly periodic unit square.
    structured_rectangle(
        0.0, 1.0, 0.0, 1.0, 16, 16,
        {"periodicx": ("left", "right"), "periodicy": ("bottom", "top")},
        os.path.join(HERE, "plane_wave2d", "square_periodic_16x16.msh"),
        "square_periodic")

    # beam2d: simply supported beam, span 1 and depth 0.2 (L/h = 5), loaded
    # on a patch of its top surface at midspan. Three boundary groups, because
    # each gets a different condition:
    #
    #   "support"  the two end faces      hinged: v = 0, σxx = 0
    #   "top"      the loaded surface     σyy = -p(x,t), σxy = 0
    #   "bottom"   the free underside     σyy = 0, σxy = 0
    #
    # 32 x 6 elements, near-square at 0.03125 x 0.03333.
    structured_rectangle(
        0.0, 1.0, -0.1, 0.1, 32, 6,
        {"support": ("left", "right"), "top": ("top",), "bottom": ("bottom",)},
        os.path.join(HERE, "beam2d", "beam_32x6.msh"),
        "beam")

    # viscoelastic_cylinder2d: Artioli et al. (2017) Sec. 4.2, a quarter of the
    # thick-walled cylinder a = 2, b = 4. 4 radial x 8 circumferential
    # elements: at nop = 4 that is 17 x 33 LGL points, and the spatial error is
    # far below the time-integration error the case is about (see its README).
    structured_quarter_annulus(
        2.0, 4.0, 4, 8,
        os.path.join(HERE, "viscoelastic_cylinder2d", "quarter_annulus_4x8.msh"),
        "quarter_annulus")
