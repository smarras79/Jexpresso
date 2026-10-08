#!/usr/bin/env python3
"""log10(rho/rho0) and magnetic field lines from Jexpresso MPI output (iter_N.pvtu -> iter_N/iter_N_*.vtu), on the native grid.
Field lines: isolines A_n = k/(N+1), k = 1..N, of the flux function A (Bx = dA/dz, Bz = -dA/dx), A_n = (A - Amin)/(Amax - Amin).
  python3 tools/plot_fluxemergence_vtu.py output/MHD/<case>/output-<date> [--steps 5,10,20-40] [--vectors] [--outdir figs] [--format pdf]
Needs numpy and matplotlib; scipy only for --a-method lsq."""
import argparse, base64, gc, os, re, sys
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.tri import Triangulation

# ParaView "Rainbow Desaturated" preset, RGB interpolation
_RD = [(0.000, (0.278431372549, 0.278431372549, 0.858823529412)),
       (0.143, (0.000000000000, 0.000000000000, 0.360784313725)),
       (0.285, (0.000000000000, 1.000000000000, 1.000000000000)),
       (0.429, (0.000000000000, 0.501960784314, 0.000000000000)),
       (0.571, (1.000000000000, 1.000000000000, 0.000000000000)),
       (0.714, (1.000000000000, 0.380392156863, 0.000000000000)),
       (0.857, (0.419607843137, 0.000000000000, 0.000000000000)),
       (1.000, (0.878431372549, 0.301960784314, 0.301960784314))]
RAINBOW_DESATURATED = LinearSegmentedColormap.from_list("rainbow_desaturated", _RD, N=1024)

_NP = {"Float64": "f8", "Float32": "f4", "Int64": "i8", "Int32": "i4", "Int16": "i2", "Int8": "i1",
       "UInt64": "u8", "UInt32": "u4", "UInt16": "u2", "UInt8": "u1"}
_GREEK = {"rho": "ρ", "beta": "β", "psi": "ψ"}


def greek(name):
    """ASCII aliases for the Greek field names: rho, log10_rho, beta, psi."""
    return re.sub(r"(?<![A-Za-z])(rho|beta|psi)(?![A-Za-z])", lambda m: _GREEK[m.group(1)], name)


def ascii_name(name):
    return name.translate(str.maketrans({v: k for k, v in _GREEK.items()}))


def read_vtu(path):
    """{(section, name): array} of one serial VTU piece (appended raw, inline base64 or ascii)."""
    raw = open(path, "rb").read()
    k = raw.find(b"<AppendedData")
    xml = (raw if k < 0 else raw[:k]).decode("utf-8", errors="replace")
    blob = None if k < 0 else raw.find(b"_", raw.find(b">", k)) + 1
    vf = re.search(r"<VTKFile\b[^>]*>", xml)
    if vf is None or (k >= 0 and blob <= 0):
        raise ValueError(f"{path}: not a complete VTU file")
    vf = vf.group(0)
    if "compressor=" in vf:
        raise ValueError(f"{path}: compressed VTU is not supported (Jexpresso writes compress=false)")
    bo = ">" if 'byte_order="BigEndian"' in vf else "<"
    m = re.search(r'header_type="(\w+)"', vf)
    htype = np.dtype(bo + _NP[m.group(1) if m else "UInt32"])

    def block(buf, start, dt):
        if start + htype.itemsize > len(buf):
            raise ValueError(f"{path}: truncated data (piece still being written?)")
        n = int(np.frombuffer(buf, htype, 1, start)[0])
        if start + htype.itemsize + n > len(buf):
            raise ValueError(f"{path}: truncated data (piece still being written?)")
        return np.frombuffer(buf, dt, n // dt.itemsize, start + htype.itemsize)

    def b64block(text, dt):
        s, hc = "".join(text.split()), 4 * -(-htype.itemsize // 3)
        if s[hc - 1] == "=":  # WriteVTK encodes header and data as two base64 streams
            return block(base64.b64decode(s[:hc]) + base64.b64decode(s[hc:]), 0, dt)
        return block(base64.b64decode(s), 0, dt)

    out = {}
    for sec in ("Points", "Cells", "PointData", "FieldData"):
        s = re.search(rf"<{sec}\b[^>]*>(.*?)</{sec}>", xml, re.S)
        if s is None:
            continue
        for a in re.finditer(r"<DataArray\b([^>]*?)(?:/>|>(.*?)</DataArray>)", s.group(1), re.S):
            at = dict(re.findall(r'([\w:.-]+)="([^"]*)"', a.group(1)))
            if at.get("type") not in _NP:
                continue
            dt = np.dtype(bo + _NP[at["type"]])
            fmt = at.get("format", "appended")
            if fmt == "appended":
                if blob is None:
                    raise ValueError(f"{path}: truncated data (piece still being written?)")
                v = block(raw, blob + int(at["offset"]), dt)
            elif fmt == "binary":
                v = b64block(a.group(2), dt)
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
            raise ValueError(f"{f}: only VTK_QUAD cells are supported")
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
        raise KeyError(f"no point field '{name}'. Available: " + ", ".join(sorted(g["fields"])))
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


def plot_one(path, a, out):
    g = read_pvtu(path)
    field, Bx, Bz = _field(g, a.var, path), _field(g, a.bx, path), _field(g, a.bz, path)
    tri = triangles(g)
    T = Triangulation(g["x"], g["y"], tri)
    A = flux_function_integrate(g, Bx, Bz) if a.a_method in ("auto", "integrate") else None
    if A is None:
        if a.a_method == "integrate":
            raise ValueError("the nodes are not a tensor-product grid; use --a-method lsq")
        A = flux_function_lsq(g, tri, Bx, Bz)
    An = (A - A.min()) / max(np.ptp(A), 1e-300)

    vmin, vmax = a.clim if a.clim else (np.nanmin(field), np.nanmax(field))
    vmax = vmax if vmax > vmin else vmin + 1.0
    Lx, Ly = np.ptp(g["x"]), np.ptp(g["y"])
    fig, ax = plt.subplots(figsize=(a.width, a.width * Ly / Lx + 1.2))
    fc = np.clip(field, vmin, vmax)
    if a.shading == "gouraud":
        ax.tripcolor(T, fc, shading="gouraud", cmap=a.cmap, vmin=vmin, vmax=vmax, rasterized=True)
    else:  # linear in the scalar on each triangle, 256 bands like ParaView's lookup table
        pc = ax.tricontourf(T, fc, levels=np.linspace(vmin, vmax, 257), cmap=a.cmap, vmin=vmin, vmax=vmax,
                            antialiased=False)
        pc.set_rasterized(True)
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
    fig.colorbar(plt.cm.ScalarMappable(Normalize(vmin, vmax), a.cmap), ax=ax, shrink=0.8, pad=0.02)
    fig.savefig(out, dpi=a.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f" {path} -> {out}")


def resolve_cmap(name):
    """RainbowDesaturated (any spelling, _r reverses) or a matplotlib colormap name."""
    key = name.lower().replace(" ", "").replace("_", "").replace("-", "")
    if key in ("rainbowdesaturated", "rainbowdesaturatedr"):
        return RAINBOW_DESATURATED.reversed() if key.endswith("dr") else RAINBOW_DESATURATED
    return plt.get_cmap(name)


def pvtu_time(path):
    """TimeValue of iter_N.pvtu (from its FieldData, else from its first piece); None if absent."""
    root = ET.parse(path).getroot()
    for da in root.iter("DataArray"):
        if da.get("Name") == "TimeValue" and (da.text or "").strip():
            return float(da.text.split()[0])
    piece = next(root.iter("Piece"), None)
    if piece is None:
        return None
    d = read_vtu(os.path.join(os.path.dirname(path), piece.get("Source")))
    return float(d[("FieldData", "TimeValue")][0]) if ("FieldData", "TimeValue") in d else None


def parse_steps(text):
    """'5,10,20-40' -> {5, 10, 20, ..., 40}"""
    out = set()
    for tok in text.replace(" ", ",").split(","):
        if tok:
            lo, _, hi = tok.partition("-")
            out.update(range(int(lo), int(hi or lo) + 1))
    return out


def step_of(f):
    m = re.search(r"iter_(\d+)\.pvtu$", os.path.basename(f))
    return int(m.group(1)) if m else -1


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("paths", nargs="+", help="output directories (all iter_N.pvtu) and/or .pvtu files")
    p.add_argument("--steps", type=parse_steps, help="only these N of iter_N.pvtu, e.g. 5,10,20-40")
    p.add_argument("--var", default="log10_ρ", help="colored point field (default log10_ρ; rho/beta/psi accepted)")
    p.add_argument("--clim", type=float, nargs=2, default=None, help="color range (default -8.1 0 for log10_ρ)")
    p.add_argument("--cmap", default="RainbowDesaturated", help="RainbowDesaturated (default) or a matplotlib name; _r reverses")
    p.add_argument("--shading", choices=("contourf", "gouraud"), default="contourf",
                   help="contourf: 256 bands, scalar-linear per triangle (default); gouraud: faster, blends colors")
    p.add_argument("--bx", default="Bx", help="horizontal field component (default Bx)")
    p.add_argument("--bz", default="By", help="vertical field component (Jexpresso's By, default)")
    p.add_argument("--nlevels", type=int, default=29, help="number of field lines N (default 29)")
    p.add_argument("--linewidth", type=float, default=0.8)
    p.add_argument("--a-method", choices=("auto", "integrate", "lsq"), default="auto",
                   help="integrate: trapezoids on the tensor LGL nodes; lsq: P1 least squares (any mesh, needs scipy)")
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
    a.var, a.bx, a.bz = greek(a.var), greek(a.bx), greek(a.bz)
    a.vec_fields = [greek(v) for v in a.vec_fields]
    if a.clim is None and a.var == "log10_ρ":
        a.clim = (-8.1, 0.0)
    try:
        a.cmap = resolve_cmap(a.cmap)
    except ValueError as e:
        p.error(str(e))
    if a.outdir:
        os.makedirs(a.outdir, exist_ok=True)

    files = []
    for pth in a.paths:
        if os.path.isdir(pth):
            files += [os.path.join(pth, n) for n in os.listdir(pth) if re.fullmatch(r"iter_\d+\.pvtu", n)]
        elif os.path.isfile(pth):
            files.append(pth)
        else:
            p.error(f"{pth}: no such file or directory")
    rundir = lambda f: os.path.dirname(os.path.abspath(f))
    files = sorted(set(map(os.path.normpath, files)), key=lambda f: (rundir(f), step_of(f), f))
    if a.steps:
        missing = sorted(a.steps - {step_of(f) for f in files})
        if missing:
            print(f"warning: no iter_N.pvtu for steps {missing}", file=sys.stderr)
        files = [f for f in files if step_of(f) in a.steps]
    if not files:
        sys.exit("no iter_N.pvtu files found")
    multi = a.outdir and len({rundir(f) for f in files}) > 1

    failed = []
    for f in files:
        stem = re.sub(r"\.pvtu$", "", os.path.basename(f))
        prefix = os.path.basename(rundir(f)) + "-" if multi else ""
        out = os.path.join(a.outdir or rundir(f), f"{prefix}{ascii_name(a.var)}-fieldlines-{stem}.{a.format}")
        try:
            plot_one(f, a, out)
        except Exception as e:
            plt.close("all")
            print(f" {f}: skipped ({type(e).__name__}: {e})", file=sys.stderr)
            failed.append(f)
        gc.collect()
    if failed:
        sys.exit(f"{len(failed)} of {len(files)} snapshots failed")


if __name__ == "__main__":
    main()
