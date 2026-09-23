#!/usr/bin/env python3
"""
Generate the structured quad meshes of problems/MHD/astroJetWuShu2018 WITHOUT a
gmsh installation.

    tools/astro_jet_mesh.py                  # 40x60 and 100x150, the two the case ships
    tools/astro_jet_mesh.py 60 90            # any nx ny, written as AJ_<nx>x<ny>.msh

The output is a gmsh 4.1 ASCII .msh with the same entity / node-block /
element-block layout gmsh itself produces for the transfinite, recombined
rectangle of problems/MHD/astroJetWuShu2018/AJ.geo, so GridapGmsh reads it
exactly as it reads a gmsh-generated file. `gmsh -2 AJ.geo -o AJ_40x60.msh`
remains the reference path; this script exists so that the case runs out of the
box on a machine without gmsh, and so that the shipped .msh files are
reproducible.

THE ONE CONSTRAINT ON nx. The nozzle of the jet is {y = 0, |x| <= 0.05} and the
inflow is imposed at NODES (user_bc.jl), so the lip |x| = 0.05 must fall on an
element boundary: 0.05*nx must be an integer, i.e. nx must be a multiple of 20.
The script refuses anything else.

Keep the elements square (ny = 1.5*nx): the DynSGS element length scale is
Delta = min(edge)/(N+1), so a stretched element gets the dissipation of its
short side over its long one.

Physical curves: "bottom" (y=0), "right" (x=+0.5), "top" (y=1.5), "left"
(x=-0.5); physical surface "domain". Those four strings are the tags that reach
user_bc_dirichlet!.
"""

import os
import sys

# The domain of Wu & Shu (2018), Example 5.6.
XMIN, XMAX = -0.5, 0.5
YMIN, YMAX = 0.0, 1.5
XNOZZLE = 0.05

CURVE_NAMES = ("bottom", "right", "top", "left")


def _fmt(v):
    """Short decimal, integers as integers — the way gmsh prints coordinates."""
    return str(int(v)) if v == int(v) else repr(float(v))


def write_msh(path, nx, ny, xmin=XMIN, xmax=XMAX, ymin=YMIN, ymax=YMAX):
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny
    X = [xmin + i * dx for i in range(nx + 1)]
    Y = [ymin + j * dy for j in range(ny + 1)]
    X[0], X[nx] = xmin, xmax          # exact endpoints, not xmin + nx*dx
    Y[0], Y[ny] = ymin, ymax

    # ---- node tags ---------------------------------------------------------
    # Nodes are attached to the entity that owns them, as gmsh does it: the four
    # corner points, then the interior of each curve in the direction that curve
    # is oriented, then the interior of the surface.
    nid = [[0] * (ny + 1) for _ in range(nx + 1)]
    tag = 0
    pt_blocks = []
    for (i, j) in ((0, 0), (nx, 0), (nx, ny), (0, ny)):
        tag += 1
        nid[i][j] = tag
        pt_blocks.append((tag, X[i], Y[j]))

    curves = (
        [(i, 0) for i in range(1, nx)],            # 1 bottom  P1 -> P2
        [(nx, j) for j in range(1, ny)],           # 2 right   P2 -> P3
        [(i, ny) for i in range(nx - 1, 0, -1)],   # 3 top     P3 -> P4
        [(0, j) for j in range(ny - 1, 0, -1)],    # 4 left    P4 -> P1
    )
    cv_blocks = []
    for ci, lst in enumerate(curves, start=1):
        blk = []
        for (i, j) in lst:
            tag += 1
            nid[i][j] = tag
            blk.append((tag, X[i], Y[j]))
        cv_blocks.append((ci, blk))

    sf_blk = []
    for j in range(1, ny):
        for i in range(1, nx):
            tag += 1
            nid[i][j] = tag
            sf_blk.append((tag, X[i], Y[j]))
    nnodes = tag

    # ---- elements: the boundary lines, then the quads ----------------------
    etag = 0
    line_blocks = []
    for ci, lst in ((1, [(i, 0) for i in range(nx + 1)]),
                    (2, [(nx, j) for j in range(ny + 1)]),
                    (3, [(i, ny) for i in range(nx, -1, -1)]),
                    (4, [(0, j) for j in range(ny, -1, -1)])):
        blk = []
        for k in range(len(lst) - 1):
            etag += 1
            blk.append((etag,
                        nid[lst[k][0]][lst[k][1]],
                        nid[lst[k + 1][0]][lst[k + 1][1]]))
        line_blocks.append((ci, blk))

    quads = []
    for j in range(ny):
        for i in range(nx):
            etag += 1
            # counter-clockwise, so the Jacobian is positive
            quads.append((etag, nid[i][j], nid[i + 1][j],
                          nid[i + 1][j + 1], nid[i][j + 1]))
    nelem = etag

    # ---- write -------------------------------------------------------------
    o = ["$MeshFormat\n4.1 0 8\n$EndMeshFormat\n", "$PhysicalNames\n5\n"]
    for d, t, nm in ((1, 2, CURVE_NAMES[0]), (1, 3, CURVE_NAMES[1]),
                     (1, 4, CURVE_NAMES[2]), (1, 5, CURVE_NAMES[3]),
                     (2, 1, "domain")):
        o.append('%d %d "%s"\n' % (d, t, nm))
    o.append("$EndPhysicalNames\n")

    o.append("$Entities\n4 4 1 0\n")
    for (t, x, y) in pt_blocks:
        o.append("%d %s %s 0 0 \n" % (t, _fmt(x), _fmt(y)))
    # curveTag  bbox(min) bbox(max)  numPhys phys  numBoundingPoints pts
    cbox = {1: (xmin, ymin, xmax, ymin, 1, -2),
            2: (xmax, ymin, xmax, ymax, 2, -3),
            3: (xmin, ymax, xmax, ymax, 3, -4),
            4: (xmin, ymin, xmin, ymax, 4, -1)}
    cphys = {1: 2, 2: 3, 3: 4, 4: 5}
    for ci in (1, 2, 3, 4):
        ax, ay, bx, by, p0, p1 = cbox[ci]
        o.append("%d %s %s 0 %s %s 0 1 %d 2 %d %d \n"
                 % (ci, _fmt(ax), _fmt(ay), _fmt(bx), _fmt(by), cphys[ci], p0, p1))
    o.append("1 %s %s 0 %s %s 0 1 1 4 1 2 3 4 \n"
             % (_fmt(xmin), _fmt(ymin), _fmt(xmax), _fmt(ymax)))
    o.append("$EndEntities\n")

    o.append("$Nodes\n9 %d 1 %d\n" % (nnodes, nnodes))
    for (t, x, y) in pt_blocks:
        o.append("0 %d 0 1\n%d\n%s %s 0\n" % (t, t, _fmt(x), _fmt(y)))
    for ci, blk in cv_blocks:
        o.append("1 %d 0 %d\n" % (ci, len(blk)))
        o.extend("%d\n" % b[0] for b in blk)
        o.extend("%s %s 0\n" % (_fmt(b[1]), _fmt(b[2])) for b in blk)
    o.append("2 1 0 %d\n" % len(sf_blk))
    o.extend("%d\n" % b[0] for b in sf_blk)
    o.extend("%s %s 0\n" % (_fmt(b[1]), _fmt(b[2])) for b in sf_blk)
    o.append("$EndNodes\n")

    o.append("$Elements\n5 %d 1 %d\n" % (nelem, nelem))
    for ci, blk in line_blocks:
        o.append("1 %d 1 %d\n" % (ci, len(blk)))
        o.extend("%d %d %d \n" % e for e in blk)
    o.append("2 1 3 %d\n" % len(quads))
    o.extend("%d %d %d %d %d \n" % q for q in quads)
    o.append("$EndElements\n")

    with open(path, "w") as f:
        f.write("".join(o))
    return nnodes, nelem


def main(argv):
    here = os.path.dirname(os.path.abspath(__file__))
    outdir = os.path.join(here, os.pardir, "problems", "MHD", "astroJetWuShu2018")

    if len(argv) == 3:
        sizes = [(int(argv[1]), int(argv[2]))]
    elif len(argv) == 1:
        sizes = [(40, 60), (100, 150)]
    else:
        sys.exit("usage: astro_jet_mesh.py [nx ny]")

    for nx, ny in sizes:
        if abs(XNOZZLE * nx - round(XNOZZLE * nx)) > 1e-12:
            sys.exit("nx = %d puts the nozzle lip |x| = %g inside an element; "
                     "nx must be a multiple of 20." % (nx, XNOZZLE))
        path = os.path.join(outdir, "AJ_%dx%d.msh" % (nx, ny))
        nn, ne = write_msh(path, nx, ny)
        h = (XMAX - XMIN) / nx
        print("%s  %d nodes  %d elements  h = %g x %g  (%d elements across the nozzle)"
              % (os.path.relpath(path), nn, ne, h, (YMAX - YMIN) / ny,
                 round(2 * XNOZZLE / h)))


if __name__ == "__main__":
    main(sys.argv)
