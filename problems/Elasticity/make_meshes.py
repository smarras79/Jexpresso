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


def structured_half_annulus(a, b, nr, ntq, path, name):
    """The half annulus a <= r <= b, 0 <= theta <= pi (an arch), as two
    structured quarters (a gmsh arc must be shorter than pi): nr quads across
    the wall, ntq along each quarter. The arcs are curved back onto their
    circles by the case's :exact_geometry.

    Boundary groups:

      "load"    theta = 0:  y = 0,  a <= x <= b      the free end, loaded
      "clamp"   theta = pi: y = 0, -b <= x <= -a     clamped
      "inner"   r = a                                traction-free
      "outer"   r = b                                traction-free
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(name)

    geo = gmsh.model.geo
    o   = geo.addPoint(0.0, 0.0, 0.0)
    pa0 = geo.addPoint(a, 0.0, 0.0)
    pb0 = geo.addPoint(b, 0.0, 0.0)
    pb1 = geo.addPoint(0.0, b, 0.0)
    pa1 = geo.addPoint(0.0, a, 0.0)
    pb2 = geo.addPoint(-b, 0.0, 0.0)
    pa2 = geo.addPoint(-a, 0.0, 0.0)

    load  = geo.addLine(pa0, pb0)
    mid   = geo.addLine(pb1, pa1)
    clamp = geo.addLine(pb2, pa2)
    o1 = geo.addCircleArc(pb0, o, pb1); o2 = geo.addCircleArc(pb1, o, pb2)
    i1 = geo.addCircleArc(pa1, o, pa0); i2 = geo.addCircleArc(pa2, o, pa1)

    s1 = geo.addPlaneSurface([geo.addCurveLoop([load, o1, mid, i1])])
    s2 = geo.addPlaneSurface([geo.addCurveLoop([-mid, o2, clamp, i2])])

    for c in (load, mid, clamp):
        geo.mesh.setTransfiniteCurve(c, nr + 1)
    for c in (o1, o2, i1, i2):
        geo.mesh.setTransfiniteCurve(c, ntq + 1)
    for srf in (s1, s2):
        geo.mesh.setTransfiniteSurface(srf)
        geo.mesh.setRecombine(2, srf)

    geo.synchronize()

    for gname, curves in (("load", [load]), ("clamp", [clamp]),
                          ("inner", [i1, i2]), ("outer", [o1, o2])):
        tag = gmsh.model.addPhysicalGroup(1, curves)
        gmsh.model.setPhysicalName(1, tag, gname)
    tag = gmsh.model.addPhysicalGroup(2, [s1, s2])
    gmsh.model.setPhysicalName(2, tag, "domain")

    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(path)

    nnode = len(gmsh.model.mesh.getNodes()[0])
    gmsh.finalize()
    print(f"  {os.path.relpath(path, HERE)}: {nr}x{2*ntq} quads, {nnode} vertices")


def structured_perforated_plate(L, H, R, nt, nr, ny, grading, path, name):
    """A quarter of a strip 2L x 2H with a central hole of radius R, as three
    structured blocks: two fans around the hole, split by the 45-degree
    diagonal, and the rectangle above them,

        block I    A=(R,0) -> P1=(L,0) -> Q=(L,L) -> M (hole at 45 deg)
        block II   M -> Q -> D=(0,L) -> B=(0,R)
        block III  D -> Q -> P2=(L,H) -> P3=(0,H)          (needs H > L)

    nt elements along each 45-degree arc (and across the strip), nr along the
    radial lines — graded towards the hole by `grading` (> 1), where the
    plastic zone starts — and ny up the rectangle.

    Boundary groups:

      "hole"     the quarter circle r = R              traction-free
      "bottom"   y = 0, R <= x <= L                    symmetry: u_y = 0
      "left"     x = 0, R <= y <= H                    symmetry: u_x = 0
      "right"    x = L                                 lateral restraint
      "top"      y = H                                 the prescribed displacement
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(name)

    geo = gmsh.model.geo
    c45 = R/2**0.5
    o  = geo.addPoint(0.0, 0.0, 0.0)
    pa = geo.addPoint(R, 0.0, 0.0)          # A
    p1 = geo.addPoint(L, 0.0, 0.0)
    q  = geo.addPoint(L, L, 0.0)
    m  = geo.addPoint(c45, c45, 0.0)
    d  = geo.addPoint(0.0, L, 0.0)
    pb = geo.addPoint(0.0, R, 0.0)          # B
    p2 = geo.addPoint(L, H, 0.0)
    p3 = geo.addPoint(0.0, H, 0.0)

    # radial lines run from the hole outwards, so the grading is consistent
    c_bot  = geo.addLine(pa, p1)
    c_diag = geo.addLine(m, q)
    c_left = geo.addLine(pb, d)
    arc1   = geo.addCircleArc(pa, o, m)
    arc2   = geo.addCircleArc(m, o, pb)
    c_r1   = geo.addLine(p1, q)
    c_t1   = geo.addLine(q, d)
    c_r2   = geo.addLine(q, p2)
    c_top  = geo.addLine(p2, p3)
    c_l2   = geo.addLine(d, p3)

    s1 = geo.addPlaneSurface([geo.addCurveLoop([c_bot, c_r1, -c_diag, -arc1])])
    s2 = geo.addPlaneSurface([geo.addCurveLoop([c_diag, c_t1, -c_left, -arc2])])
    s3 = geo.addPlaneSurface([geo.addCurveLoop([-c_t1, c_r2, c_top, -c_l2])])

    for c in (c_bot, c_diag, c_left):
        geo.mesh.setTransfiniteCurve(c, nr + 1, "Progression", grading)
    for c in (arc1, arc2, c_r1, c_t1, c_top):
        geo.mesh.setTransfiniteCurve(c, nt + 1)
    for c in (c_r2, c_l2):
        geo.mesh.setTransfiniteCurve(c, ny + 1)
    for s in (s1, s2, s3):
        geo.mesh.setTransfiniteSurface(s)
        geo.mesh.setRecombine(2, s)

    geo.synchronize()

    for gname, curves in (("hole", [arc1, arc2]), ("bottom", [c_bot]), ("left", [c_left, c_l2]),
                          ("right", [c_r1, c_r2]), ("top", [c_top])):
        tag = gmsh.model.addPhysicalGroup(1, curves)
        gmsh.model.setPhysicalName(1, tag, gname)
    tag = gmsh.model.addPhysicalGroup(2, [s1, s2, s3])
    gmsh.model.setPhysicalName(2, tag, "domain")

    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(path)

    nquad = 2*nt*nr + nt*ny
    nnode = len(gmsh.model.mesh.getNodes()[0])
    gmsh.finalize()
    print(f"  {os.path.relpath(path, HERE)}: {nquad} quads, {nnode} vertices")


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

    # plastic_cylinder2d: the same quarter cylinder, 8 radial elements so the
    # elastic-plastic front crossing the wall is resolved.
    structured_quarter_annulus(
        2.0, 4.0, 8, 8,
        os.path.join(HERE, "plastic_cylinder2d", "quarter_annulus_8x8.msh"),
        "quarter_annulus")

    # perforated_plate2d: Artioli et al. (2017) Sec. 4.3, a quarter of the
    # 200 x 360 mm strip with its 100 mm hole. 132 quads, graded towards the
    # hole, where the plastic zone starts.
    structured_perforated_plate(
        100.0, 180.0, 50.0, 6, 8, 6, 1.12,
        os.path.join(HERE, "perforated_plate2d", "quarter_plate_hole_132.msh"),
        "quarter_plate_hole")

    # sma_arch2d: Artioli et al. (2017) Sec. 4.4, the semicircular SMA arch,
    # 3.5 <= r <= 4.5 mm: 2 elements across the wall, 24 along the arch.
    structured_half_annulus(
        3.5, 4.5, 2, 12,
        os.path.join(HERE, "sma_arch2d", "half_annulus_2x24.msh"),
        "half_annulus")
