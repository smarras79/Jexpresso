#!/usr/bin/env python3
"""
generate_mesh.py -- write ramp15.msh (gmsh 4.1 ASCII) for the 2D hypersonic
compression ramp of Cao, Hao, Klioutchnikov, Olivier & Wen, JFM 912 A3 (2021).

This is the same grid ramp15.geo describes; it is written directly so the
case can be regenerated without a gmsh installation.  Every number below is
read from ramp15.geo -- keep the two in step.

Topology.  The upper boundary is the wall contour shifted VERTICALLY by H,
so every line of constant streamwise index is vertical and the three
transfinite blocks are one logical (NX+1) x (NY+1) structured grid:

    node(I, J) = ( x[I],  y_wall[I] + f[J]*H )

with f the wall-normal Progression and x / y_wall the streamwise
distribution.  Nodes on the geometric corners belong to the point entities,
nodes on the ten boundary/interface curves to the curve entities, the rest
to the three surface entities -- which is what gmsh itself writes.

    python3 generate_mesh.py [-o ramp15.msh]
"""

import argparse
import math

# ---------------------------------------------------------------- geometry
L        = 0.1                    # flat plate length [m]
ALPHA    = math.radians(15.0)     # ramp deflection
RAMP_LEN = 0.1                    # ramp length along the ramp surface [m]
H        = 0.06                   # domain height above the wall contour [m]

X_END = L + RAMP_LEN*math.cos(ALPHA)
Y_END =     RAMP_LEN*math.sin(ALPHA)

# -------------------------------------------------------------- resolution
NX_PLATE, NX_RAMP, NY = 137, 132, 60
PX, PY = 1.006, 1.0808            # Progression ratios (plate, wall-normal)

NX = NX_PLATE + NX_RAMP           # 269
I_CORNER = NX_PLATE               # streamwise index of the compression corner

# physical tags, matching ramp15.geo
PH_DOMAIN, PH_INFLOW, PH_OUTFLOW, PH_WALL, PH_TOP = 1, 2, 3, 4, 5


def progression(n, p):
    """gmsh `Using Progression p` over n intervals: fractions 0 .. 1."""
    if abs(p - 1.0) < 1.0e-12:
        return [i/n for i in range(n + 1)]
    den = p**n - 1.0
    return [(p**i - 1.0)/den for i in range(n + 1)]


def build_grid():
    """Return xs[0..NX], yws[0..NX] (the wall contour) and fy[0..NY]."""
    xs, yws = [], []

    # A: flat plate, clustered at the leading edge, y = 0
    fp = progression(NX_PLATE, PX)
    for i in range(NX_PLATE):
        xs.append(L*fp[i])
        yws.append(0.0)

    # B: ramp, uniform in arclength
    for i in range(NX_RAMP + 1):
        s = RAMP_LEN*i/NX_RAMP
        xs.append(L + s*math.cos(ALPHA))
        yws.append(s*math.sin(ALPHA))

    assert len(xs) == NX + 1
    return xs, yws, progression(NY, PY)


# ------------------------------------------------------------------ entities
# Corner (I, J) of each of the eight geometric points, in gmsh point order.
def corner_ij():
    return {1: (0, 0),   2: (I_CORNER, 0),  3: (NX, 0),
            4: (NX, NY), 5: (I_CORNER, NY), 6: (0, NY)}


# Each curve: (list of (I,J) along it INCLUDING both end points, physical tag
#              or 0, the two bounding point tags in curve orientation)
def curves():
    horiz = lambda j, i0, i1: [(i, j) for i in range(i0, i1 + 1)]
    vert  = lambda i, j0, j1: [(i, j) for j in range(j0, j1 + 1)]
    return {
        1: (horiz(0,  0,        I_CORNER), PH_WALL,    (1, 2)),
        2: (horiz(0,  I_CORNER, NX),       PH_WALL,    (2, 3)),
        3: (vert(NX,  0,        NY),       PH_OUTFLOW, (3, 4)),
        4: (horiz(NY, I_CORNER, NX),       PH_TOP,     (5, 4)),
        5: (horiz(NY, 0,        I_CORNER), PH_TOP,     (6, 5)),
        6: (vert(0,   0,        NY),       PH_INFLOW,  (1, 6)),
        7: (vert(I_CORNER, 0,   NY),       0,          (2, 5)),
    }


# Each surface: (I range of its elements, bounding curve tags with sign)
def surfaces():
    return {
        1: (range(0, I_CORNER),  [1, 7, -5, -6]),
        2: (range(I_CORNER, NX), [2, 3, -4, -7]),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--output", default="ramp15.msh")
    args = ap.parse_args()

    xs, yws, fy = build_grid()
    node = lambda i, j: j*(NX + 1) + i + 1          # 1-based global node tag
    coord = lambda i, j: (xs[i], yws[i] + fy[j]*H)

    cij, crv, srf = corner_ij(), curves(), surfaces()
    corner_nodes = {node(*ij) for ij in cij.values()}

    # interior nodes of each curve (the ends belong to the point entities)
    curve_nodes = {t: [node(i, j) for (i, j) in pts[1:-1]] for t, (pts, _, _) in crv.items()}
    on_curve = {n for ns in curve_nodes.values() for n in ns}

    # ------------------------------------------------------------ assemble
    out = []
    w = out.append
    w("$MeshFormat")
    w("4.1 0 8")
    w("$EndMeshFormat")

    w("$PhysicalNames")
    w("5")
    for tag, name in ((PH_INFLOW, "inflow"), (PH_OUTFLOW, "outflow"),
                      (PH_WALL, "wall"), (PH_TOP, "top")):
        w('1 %d "%s"' % (tag, name))
    w('2 %d "domain"' % PH_DOMAIN)
    w("$EndPhysicalNames")

    w("$Entities")
    w("6 7 2 0")
    for t in range(1, 7):
        x, y = coord(*cij[t])
        w("%d %.16g %.16g 0 0" % (t, x, y))
    for t in range(1, 8):
        pts, ph, (pa, pb) = crv[t]
        cs = [coord(i, j) for (i, j) in pts]
        bb = (min(c[0] for c in cs), min(c[1] for c in cs),
              max(c[0] for c in cs), max(c[1] for c in cs))
        phs = "1 %d" % ph if ph else "0"
        w("%d %.16g %.16g 0 %.16g %.16g 0 %s 2 %d -%d"
          % (t, bb[0], bb[1], bb[2], bb[3], phs, pa, pb))
    for t in range(1, 3):
        irange, bnd = srf[t]
        i0, i1 = irange.start, irange.stop
        cs = [coord(i, j) for i in range(i0, i1 + 1) for j in (0, NY)]
        w("%d %.16g %.16g 0 %.16g %.16g 0 1 %d 4 %s"
          % (t, min(c[0] for c in cs), min(c[1] for c in cs),
             max(c[0] for c in cs), max(c[1] for c in cs), PH_DOMAIN,
             " ".join(str(c) for c in bnd)))
    w("$EndEntities")

    # ------------------------------------------------------------- nodes
    nnodes = (NX + 1)*(NY + 1)
    w("$Nodes")
    w("15 %d 1 %d" % (nnodes, nnodes))
    for t in range(1, 7):
        i, j = cij[t]
        x, y = coord(i, j)
        w("0 %d 0 1" % t)
        w("%d" % node(i, j))
        w("%.16g %.16g 0" % (x, y))
    for t in range(1, 8):
        pts = crv[t][0][1:-1]
        w("1 %d 0 %d" % (t, len(pts)))
        for (i, j) in pts:
            w("%d" % node(i, j))
        for (i, j) in pts:
            x, y = coord(i, j)
            w("%.16g %.16g 0" % (x, y))
    for t in range(1, 3):
        irange = srf[t][0]
        # a surface owns only its interior nodes: the corners and the four
        # bounding curves (the two vertical interfaces included) are already
        # claimed by the point and curve entities above.
        pts = [(i, j) for i in irange for j in range(1, NY)
               if node(i, j) not in corner_nodes and node(i, j) not in on_curve]
        w("2 %d 0 %d" % (t, len(pts)))
        for (i, j) in pts:
            w("%d" % node(i, j))
        for (i, j) in pts:
            x, y = coord(i, j)
            w("%.16g %.16g 0" % (x, y))
    w("$EndNodes")

    # ---------------------------------------------------------- elements
    nlines = sum(len(crv[t][0]) - 1 for t in crv)
    nquads = NX*NY
    nelem = 6 + nlines + nquads
    w("$Elements")
    w("15 %d 1 %d" % (nelem, nelem))
    e = 1
    for t in range(1, 7):
        w("0 %d 15 1" % t)
        w("%d %d" % (e, node(*cij[t])))
        e += 1
    for t in range(1, 8):
        pts = crv[t][0]
        w("1 %d 1 %d" % (t, len(pts) - 1))
        for k in range(len(pts) - 1):
            w("%d %d %d" % (e, node(*pts[k]), node(*pts[k + 1])))
            e += 1
    for t in range(1, 3):
        irange = srf[t][0]
        w("2 %d 3 %d" % (t, len(irange)*NY))
        for i in irange:
            for j in range(NY):
                w("%d %d %d %d %d" % (e, node(i, j), node(i + 1, j),
                                      node(i + 1, j + 1), node(i, j + 1)))
                e += 1
    w("$EndElements")
    w("")

    with open(args.output, "w") as fh:
        fh.write("\n".join(out))

    # --------------------------------------------------------- report
    dy1 = (fy[1] - fy[0])*H
    dx_plate1 = xs[1] - xs[0]
    lgl1 = 0.5*(1.0 - math.sqrt(3.0/7.0))          # first LGL interval, nop=4
    print("%s: %d nodes, %d quads (%d x %d elements)"
          % (args.output, nnodes, nquads, NX, NY))
    print("  x nodes at nop=4 : %d   (paper G1: 1080)" % (4*NX + 1))
    print("  y nodes at nop=4 : %d   (paper G1:  240)" % (4*NY + 1))
    print("  first y element  : %.4e m -> Dy_wall = %.3e m (paper 8e-6)"
          % (dy1, lgl1*dy1))
    print("  first x element  : %.4e m at the leading edge" % dx_plate1)
    print("  domain           : x in [%.4f, %.4f], H = %.3f m" % (xs[0], xs[-1], H))


if __name__ == "__main__":
    main()
