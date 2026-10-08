#!/usr/bin/env python3
"""log10(rho/rho0) and magnetic field lines from Jexpresso MPI output (iter_N.pvtu -> iter_N/iter_N_*.vtu), on the native grid.
Field lines: isolines A_n = k/(N+1), k = 1..N, of the flux function A (Bx = dA/dz, Bz = -dA/dx), A_n = (A - Amin)/(Amax - Amin).
  python3 tools/plot_fluxemergence_vtu.py output/MHD/<case>/output-<date> [--steps 5 10] [--vectors] [--outdir figs] [--format pdf]
Needs numpy and matplotlib; scipy only for --a-method lsq."""
import argparse, base64, glob, os, re, sys
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.tri import Triangulation

# ParaView "Rainbow Desaturated" preset, RGB interpolation
_RD = [(0.000, (0.278431372549, 0.278431372549, 0.858823529412)),
       (0.143, (0.000000000000, 0.000000000000, 0.360784313725)),
       (0.286, (0.000000000000, 1.000000000000, 1.000000000000)),
       (0.429, (0.000000000000, 0.501960784314, 0.000000000000)),
       (0.571, (1.000000000000, 1.000000000000, 0.000000000000)),
       (0.714, (1.000000000000, 0.380392156863, 0.000000000000)),
       (0.857, (0.419607843137, 0.000000000000, 0.000000000000)),
       (1.000, (0.878431372549, 0.301960784314, 0.301960784314))]
RAINBOW_DESATURATED = LinearSegmentedColormap.from_list("rainbow_desaturated", _RD, N=1024)

_NP = {"Float64": "f8", "Float32": "f4", "Int64": "i8", "Int32": "i4", "Int8": "i1",
       "UInt64": "u8", "UInt32": "u4", "UInt8": "u1"}


def read_vtu(path):
    """{(section, name): array} of one serial VTU piece (appended raw, inline base64 or ascii)."""
    raw = open(path, "rb").read()
    k = raw.find(b"<AppendedData")
    xml = (raw if k < 0 else raw[:k]).decode("utf-8", errors="replace")
    blob = None if k < 0 else raw.index(b"_", raw.index(b">", k)) + 1
    vf = re.search(r"<VTKFile\b[^>]*>", xml).group(0)
    if "compressor=" in vf:
        sys.exit(f"{path}: compressed VTU is not supported (Jexpresso writes compress=false)")
    bo = ">" if 'byte_order="BigEndian"' in vf else "<"
    m = re.search(r'header_type="(\w+)"', vf)
    htype = np.dtype(bo + _NP[m.group(1) if m else "UInt32"])

    def block(buf, start, dt):
        n = int(np.frombuffer(buf, htype, 1, start)[0])
        return np.frombuffer(buf, dt, n // dt.itemsize, start + htype.itemsize)

    out = {}
    for sec in ("Points", "Cells", "PointData", "FieldData"):
        s = re.search(rf"<{sec}\b[^>]*>(.*?)</{sec}>", xml, re.S)
        if s is None:
            continue
        for a in re.finditer(r"<DataArray\b([^>]*?)(?:/>|>(.*?)</DataArray>)", s.group(1), re.S):
            at = dict(re.findall(r'([\w:.-]+)="([^"]*)"', a.group(1)))
            dt = np.dtype(bo + _NP[at["type"]])
            fmt = at.get("format", "appended")
            if fmt == "appended":
                v = block(raw, blob + int(at["offset"]), dt)
            elif fmt == "binary":
                v = block(base64.b64decode(a.group(2).strip()), 0, dt)
            else:
                v = np.array(a.group(2).split(), dtype=dt)
            nc = int(at.get("NumberOfComponents", "1"))
            out[(sec, at.get("Name", ""))] = v.reshape(-1, nc) if nc > 1 else v
    return out


def _cluster(v, tol):
    """Integer id per value; values closer than tol share an id. Also returns the id centers."""
    o = np.argsort(v, kind="stable")
    ids = np.empty(len(v), np.int64)
    ids[o] = np.concatenate([[0], np.cumsum(np.diff(v[o]) > tol)])
    return ids, np.bincount(ids, weights=v) / np.bincount(ids)


def read_pvtu(path):
    """Merge the rank pieces of iter_N.pvtu (all point fields); nodes shared by ranks are merged."""
    pieces = [os.path.join(os.path.dirname(path), p.get("Source")) for p in ET.parse(path).getroot().iter("Piece")]
    xs, ys, quads, vals, t, n0 = [], [], [], {}, None, 0
    for f in pieces:
        d = read_vtu(f)
        P = next(v for (s, _), v in d.items() if s == "Points").reshape(-1, 3)
        if not np.all(d[("Cells", "types")] == 9):
            sys.exit(f"{f}: only VTK_QUAD cells are supported")
        quads.append(d[("Cells", "connectivity")].reshape(-1, 4).astype(np.int64) + n0)
        xs.append(P[:, 0]); ys.append(P[:, 1])
        for (sec, nm), v in d.items():
            if sec == "PointData" and v.ndim == 1:
                vals.setdefault(nm, []).append(np.asarray(v, float))
        if ("FieldData", "TimeValue") in d:
            t = float(d[("FieldData", "TimeValue")][0])
        n0 += len(P)
    x, y, q = np.concatenate(xs), np.concatenate(ys), np.concatenate(quads)
    tol = 1e-8 * max(np.ptp(x), np.ptp(y))
    ix, xc = _cluster(x, tol)
    iy, yc = _cluster(y, tol)
    _, first, inv = np.unique(ix * len(yc) + iy, return_index=True, return_inverse=True)
    fields = {nm: np.concatenate(v)[first] for nm, v in vals.items() if len(v) == len(pieces)}
    if "log10_ρ" not in fields and "ρ" in fields:
        fields["log10_ρ"] = np.log10(np.maximum(fields["ρ"], 1e-300))
    return dict(x=x[first], y=y[first], quads=inv.ravel()[q], ix=ix[first], iy=iy[first], xc=xc, yc=yc,
                t=t, fields=fields)


def _field(g, name, path):
    if name not in g["fields"]:
        sys.exit(f"{path}: no point field '{name}'. Available: " + ", ".join(sorted(g["fields"])))
    return g["fields"][name]


def velocity_vectors(ax, g, a, path):
    """White arrows at the native nodes nearest an nx-by-ny lattice; reference arrow bottom-left (as the PNG writer)."""
    U, V = _field(g, a.vec_fields[0], path), _field(g, a.vec_fields[1], path)
    x, y = g["x"], g["y"]
    Lx, Ly = np.ptp(x), np.ptp(y)
    xs = x.min() + Lx * np.arange(1, a.vec_n[0] + 1) / (a.vec_n[0] + 1)
    ys = y.min() + Ly * np.arange(1, a.vec_n[1] + 1) / (a.vec_n[1] + 1)
    ip = np.array([np.argmin((x - px) ** 2 + (y - py) ** 2) for py in ys for px in xs])
    ref = a.vec_ref if a.vec_ref else max(np.max(np.hypot(U, V)), 1e-300)
    ip = ip[np.hypot(U[ip], V[ip]) >= 0.01 * ref]
    s = ref / (0.09 * Lx)  # the reference speed is drawn 0.09 Lx long
    kw = dict(angles="xy", scale_units="xy", scale=s, color=a.vec_color, width=a.vec_width,
              headwidth=4, headlength=4, headaxislength=3.5, zorder=3)
    ax.quiver(x[ip], y[ip], U[ip], V[ip], **kw)
    x0, y0 = x.min() + 0.02 * Lx, y.min() + 0.05 * Ly
    ax.quiver([x0], [y0], [ref], [0.0], **kw)
    ax.text(x0 + 0.09 * Lx + 0.01 * Lx, y0, f"= {ref:g}", color=a.vec_color, fontsize=8, va="center", zorder=3)


def triangles(g):
    """Two triangles per native SEM sub-cell; sub-cells wrapping across a periodic boundary are dropped."""
    q = g["quads"]
    keep = np.ptp(g["x"][q], axis=1) < 0.5 * np.ptp(g["x"])
    if not keep.all():
        print(f"   note: {np.count_nonzero(~keep)} periodic wrap-around cells not drawn")
    q = q[keep]
    return np.vstack([q[:, [0, 1, 2]], q[:, [0, 2, 3]]])


def flux_function_integrate(g, Bx, Bz):
    """Trapezoidal A = -int Bz(x, zmin) dx + int Bx dz over the tensor-product LGL nodes."""
    nx, ny = len(g["xc"]), len(g["yc"])
    if nx * ny != len(g["x"]):
        return None
    X, Z = g["xc"], g["yc"]
    bx = np.empty((ny, nx)); bz = np.empty((ny, nx))
    bx[g["iy"], g["ix"]] = Bx; bz[g["iy"], g["ix"]] = Bz
    A = np.zeros((ny, nx))
    A[0, 1:] = -np.cumsum(0.5 * (bz[0, 1:] + bz[0, :-1]) * np.diff(X))
    A[1:, :] = A[0, :] + np.cumsum(0.5 * (bx[1:, :] + bx[:-1, :]) * np.diff(Z)[:, None], axis=0)
    return A[g["iy"], g["ix"]]


def flux_function_lsq(g, tri, Bx, Bz):
    """Least-squares P1 fit of grad A = (-Bz, Bx) on the native triangulation (any mesh)."""
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
    x, y = g["x"], g["y"]
    X, Y = x[tri], y[tri]
    area2 = (X[:, 1] - X[:, 0]) * (Y[:, 2] - Y[:, 0]) - (X[:, 2] - X[:, 0]) * (Y[:, 1] - Y[:, 0])
    b = np.stack([Y[:, 1] - Y[:, 2], Y[:, 2] - Y[:, 0], Y[:, 0] - Y[:, 1]], 1) / area2[:, None]
    c = np.stack([X[:, 2] - X[:, 1], X[:, 0] - X[:, 2], X[:, 1] - X[:, 0]], 1) / area2[:, None]
    w = 0.5 * np.abs(area2)
    Gx, Gy = -Bz[tri].mean(1), Bx[tri].mean(1)
    Kl = w[:, None, None] * (b[:, :, None] * b[:, None, :] + c[:, :, None] * c[:, None, :])
    n = len(x)
    K = sp.coo_matrix((Kl.ravel(), (np.repeat(tri, 3, 1).ravel(), np.tile(tri, 3).ravel())), (n, n)).tocsr()
    f = np.bincount(tri.ravel(), (w[:, None] * (b * Gx[:, None] + c * Gy[:, None])).ravel(), n)
    A = np.zeros(n)
    A[1:] = spla.spsolve(K[1:, 1:].tocsc(), f[1:])
    return A


def plot_one(path, a):
    g = read_pvtu(path)
    field, Bx, Bz = _field(g, a.var, path), _field(g, a.bx, path), _field(g, a.bz, path)
    tri = triangles(g)
    T = Triangulation(g["x"], g["y"], tri)
    A = flux_function_integrate(g, Bx, Bz) if a.a_method in ("auto", "integrate") else None
    if A is None:
        if a.a_method == "integrate":
            sys.exit(f"{path}: the nodes are not a tensor-product grid; use --a-method lsq")
        A = flux_function_lsq(g, tri, Bx, Bz)
    An = (A - A.min()) / max(np.ptp(A), 1e-300)

    vmin, vmax = a.clim if a.clim else (np.nanmin(field), np.nanmax(field))
    Lx, Ly = np.ptp(g["x"]), np.ptp(g["y"])
    fig, ax = plt.subplots(figsize=(a.width, a.width * Ly / Lx + 1.2))
    pc = ax.tripcolor(T, np.clip(field, vmin, vmax), shading="gouraud", cmap=a.cmap,
                      vmin=vmin, vmax=vmax, rasterized=True)
    if np.ptp(A) > 0:
        ax.tricontour(T, An, levels=np.arange(1, a.nlevels + 1) / (a.nlevels + 1),
                      colors="k", linewidths=a.linewidth)
    if a.vectors:
        velocity_vectors(ax, g, a, path)
    ax.set_aspect("equal")
    ax.set_xlim(g["x"].min(), g["x"].max()); ax.set_ylim(g["y"].min(), g["y"].max())
    ax.set_xlabel(a.xlabel); ax.set_ylabel(a.ylabel)
    label = r"$\log_{10}(\rho/\rho_0)$" if a.var == "log10_ρ" else a.var
    ax.set_title(label + (f"   t = {g['t']:g}{a.time_unit}" if g["t"] is not None else ""))
    fig.colorbar(pc, ax=ax, shrink=0.8, pad=0.02)
    stem = re.sub(r"\.pvtu$", "", os.path.basename(path))
    out = os.path.join(a.outdir or os.path.dirname(path),
                       f"{a.var.replace('ρ', 'rho')}-fieldlines-{stem}.{a.format}")
    fig.savefig(out, dpi=a.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f" {path} -> {out}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("paths", nargs="+", help="output directories (all iter_*.pvtu) and/or iter_N.pvtu files")
    p.add_argument("--steps", type=int, nargs="*", help="only these N of iter_N.pvtu")
    p.add_argument("--var", default="log10_ρ", help="colored point field (default log10_ρ)")
    p.add_argument("--clim", type=float, nargs=2, default=None, help="color range (default -8.1 0 for log10_ρ)")
    p.add_argument("--cmap", default="rainbow_desaturated", help="matplotlib colormap name or rainbow_desaturated")
    p.add_argument("--bx", default="Bx", help="horizontal field component (default Bx)")
    p.add_argument("--bz", default="By", help="vertical field component (Jexpresso's By, default)")
    p.add_argument("--nlevels", type=int, default=29, help="number of field lines N (default 29)")
    p.add_argument("--linewidth", type=float, default=0.8)
    p.add_argument("--a-method", choices=("auto", "integrate", "lsq"), default="auto",
                   help="integrate: trapezoids on the tensor LGL nodes; lsq: P1 least squares (any mesh)")
    p.add_argument("--vectors", action="store_true", help="draw velocity vectors (white, as the PNG writer)")
    p.add_argument("--vec-fields", nargs=2, default=("u", "v"), help="velocity components (default u v)")
    p.add_argument("--vec-n", type=int, nargs=2, default=(30, 13), help="arrows in x and z (default 30 13)")
    p.add_argument("--vec-ref", type=float, default=5.0, help="reference speed, drawn 0.09 Lx long (default 5; 0 = max)")
    p.add_argument("--vec-color", default="white")
    p.add_argument("--vec-width", type=float, default=0.0015, help="shaft width, fraction of the axes width")
    p.add_argument("--outdir", default=None, help="default: next to each pvtu")
    p.add_argument("--format", default="png", help="png, pdf, svg, ...")
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--width", type=float, default=12.0, help="figure width in inches")
    p.add_argument("--xlabel", default=r"$X/H_0$")
    p.add_argument("--ylabel", default=r"$Z/H_0$")
    p.add_argument("--time-unit", default=r" $\tau_0$")
    a = p.parse_args()
    if a.var == "log10_rho":
        a.var = "log10_ρ"
    if a.clim is None and a.var == "log10_ρ":
        a.clim = (-8.1, 0.0)
    a.cmap = RAINBOW_DESATURATED if a.cmap == "rainbow_desaturated" else plt.get_cmap(a.cmap)
    if a.outdir:
        os.makedirs(a.outdir, exist_ok=True)

    step = lambda f: int(re.search(r"iter_(\d+)\.pvtu$", f).group(1))
    files = []
    for pth in a.paths:
        files += glob.glob(os.path.join(pth, "iter_*.pvtu")) if os.path.isdir(pth) else [pth]
    files = sorted(set(files), key=step)
    if a.steps:
        files = [f for f in files if step(f) in a.steps]
    if not files:
        sys.exit("no iter_N.pvtu files found")
    for f in files:
        plot_one(f, a)


if __name__ == "__main__":
    main()
