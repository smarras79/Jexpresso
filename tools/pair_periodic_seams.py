#!/usr/bin/env python3
"""Make the nodes paired across the periodic seams of a 2D gmsh mesh coincide, in place.

Jexpresso pairs the nodes of opposite periodic boundaries by position, within a tolerance, and
each element keeps its own coordinates. gmsh places the nodes of an edge with a small smooth error
(~1e-12 of the edge length) that follows the edge's direction, so when opposite edges of the .geo
run opposite ways the paired nodes are up to a few 1e-12 apart: the faces on the two sides of the
seam then differ in length, the flux through them does not cancel, and mass, momentum and energy
leak at ~1e-13 wherever the flow varies along the seam (problems/MHD/alfvenWaveBabbarRanocha2026).

The cure at the source is a .geo whose opposite edges run the same way. For a mesh without its
.geo this tool does the same to the file: every node on the right boundary takes the y of its
partner on the left, every node on the top the x of its partner on the bottom (partners in sorted
order). Only a direction whose two boundaries carry a physical name containing "periodic" is
touched, so walls are left alone. Only those coordinate tokens change; the node numbering, the
elements, the tags and every other byte stay as they are.

    python3 tools/pair_periodic_seams.py mesh.msh [more.msh ...]     # fix in place
    python3 tools/pair_periodic_seams.py --check mesh.msh            # report only

For an axis-aligned rectangular domain in gmsh 4.1 ASCII format.
"""
import sys


def read(lines):
    """Physical names, the side of the domain each periodic curve lies on, and the node rows."""
    names = {}
    if "$PhysicalNames" in lines:
        i = lines.index("$PhysicalNames") + 1
        for k in range(int(lines[i])):
            t = lines[i + 1 + k].split(maxsplit=2)
            names[(int(t[0]), int(t[1]))] = t[2].strip('"')
    i = lines.index("$Entities") + 1
    npts, ncrv = map(int, lines[i].split()[:2])
    curves = {}                                    # curve tag -> (bbox, is periodic)
    for k in range(ncrv):
        t = lines[i + 1 + npts + k].split()
        tag, box, nphys = int(t[0]), tuple(map(float, t[1:7])), int(t[7])
        phys = [int(p) for p in t[8:8 + nphys]]
        curves[tag] = (box, any("periodic" in names.get((1, p), "").lower() for p in phys))
    i = lines.index("$Nodes") + 1
    nblocks = int(lines[i].split()[0]); i += 1
    rows = []                                      # (line index, x, y)
    for _ in range(nblocks):
        cnt = int(lines[i].split()[3]); i += 1 + cnt
        for k in range(cnt):
            t = lines[i + k].split()
            rows.append((i + k, float(t[0]), float(t[1])))
        i += cnt
    return curves, rows


def pair(path, check=False):
    lines = open(path).read().split("\n")
    if not lines[1].startswith("4"):
        sys.exit(f"{path}: gmsh format {lines[1].split()[0]}, only 4.x is handled")
    curves, rows = read(lines)
    x0, x1 = min(r[1] for r in rows), max(r[1] for r in rows)
    y0, y1 = min(r[2] for r in rows), max(r[2] for r in rows)
    tol = 1e-6*max(x1 - x0, y1 - y0)

    def periodic_side(coord, val):                 # a periodic curve lying on x = val (coord 0) or y = val (1)
        return any(p and abs(b[coord] - val) < tol and abs(b[coord + 3] - val) < tol for b, p in curves.values())

    report, moved, worst = [], 0, 0.0
    # (coordinate fixed on the boundary, master value, slave value, coordinate copied): x-seam, then y-seam
    for fix, vm, vs, along, lab in ((1, x0, x1, 2, "left/right"), (2, y0, y1, 1, "bottom/top")):
        if not (periodic_side(fix - 1, vm) and periodic_side(fix - 1, vs)):
            report.append(f"{lab} not periodic, left alone")
            continue
        m = sorted((r for r in rows if abs(r[fix] - vm) < tol), key=lambda r: r[along])
        s = sorted((r for r in rows if abs(r[fix] - vs) < tol), key=lambda r: r[along])
        if len(m) != len(s):
            sys.exit(f"{path}: {lab}: {len(m)} nodes against {len(s)}, not a structured seam")
        gap = max(abs(a[along] - b[along]) for a, b in zip(m, s))
        report.append(f"{lab} paired nodes {gap:.2e} apart")
        for rm, rs in zip(m, s):
            if rs[along] != rm[along]:
                t = lines[rs[0]].split(" ")
                t[along - 1] = repr(rm[along])
                lines[rs[0]] = " ".join(t)
                moved += 1
                worst = max(worst, abs(rs[along] - rm[along]))
    if not check and moved:
        open(path, "w").write("\n".join(lines))
    action = "checked" if check else f"{moved} seam nodes moved onto their partners (largest move {worst:.2e})"
    print(f"{path}: {'; '.join(report)}; {action}")


if __name__ == "__main__":
    args = sys.argv[1:]
    check = "--check" in args
    files = [a for a in args if a != "--check"]
    if not files:
        sys.exit(__doc__)
    for f in files:
        pair(f, check)
