#!/usr/bin/env python3
"""Orszag-Tang vortex: Bx (left) and By (right) at several resolutions (rows), one shared horizontal colorbar.
Reads the MPI output (iter_N.pvtu -> iter_N/iter_N_*.vtu) of each run on its native grid, at one common time.
  python3 tools/plot_orszagtang_vtu.py [--base output/MHD/orszagTangBormanis2024] [--res 128 256 512] [--time 1.0]
                                       [--symmetric] [--single [--no-grid]] [--no-colorbar] [--hpad 0.01]

   python3 tools/plot_orszagtang_vtu.py --res 128 256 512 --time 1.0 --symmetric --format pdf --out figs/OT_B.pdf

Time grid (one 2x2 figure per resolution; columns = times, rows = variables, one colorbar per variable at the bottom):
  python3 tools/plot_orszagtang_vtu.py --timegrid [--res 64 128 256 512] [--times 0.5 1.0] [--rowvars rho mu_rho_Bx_By_Bz]
Needs numpy and matplotlib (reuses tools/plot_fluxemergence_vtu.py).

   python3 tools/plot_orszagtang_vtu.py --timegrid --res 256 512 --rowvars rho mu_dsgs_rho_Bx_By_Bz_psi
"""
import argparse, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_fluxemergence_vtu import (read_pvtu, triangles, resolve_cmap, pvtu_time, greek, ascii_name, step_of, fill)
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, ListedColormap

# ParaView "Cool to Warm (Extended)" (Remoting/Views/ColorMaps.json): (x, r, g, b), interpolated in CIELAB
_COOL_TO_WARM_EXT = [
    (0, 0, 0, 0.34902), (0.03125, 0.039216, 0.062745, 0.380392), (0.0625, 0.062745, 0.117647, 0.411765),
    (0.09375, 0.090196, 0.184314, 0.45098), (0.125, 0.12549, 0.262745, 0.501961), (0.15625, 0.160784, 0.337255, 0.541176),
    (0.1875, 0.2, 0.396078, 0.568627), (0.21875, 0.239216, 0.454902, 0.6), (0.25, 0.286275, 0.521569, 0.65098),
    (0.28125, 0.337255, 0.592157, 0.701961), (0.3125, 0.388235, 0.654902, 0.74902), (0.34375, 0.466667, 0.737255, 0.819608),
    (0.375, 0.572549, 0.819608, 0.878431), (0.40625, 0.654902, 0.866667, 0.909804), (0.4375, 0.752941, 0.917647, 0.941176),
    (0.46875, 0.823529, 0.956863, 0.968627), (0.5, 0.988235, 0.960784, 0.901961), (0.5, 0.941176, 0.984314, 0.988235),
    (0.52, 0.988235, 0.945098, 0.85098), (0.54, 0.980392, 0.898039, 0.784314), (0.5625, 0.968627, 0.835294, 0.698039),
    (0.59375, 0.94902, 0.733333, 0.588235), (0.625, 0.929412, 0.65098, 0.509804), (0.65625, 0.909804, 0.564706, 0.435294),
    (0.6875, 0.878431, 0.458824, 0.352941), (0.71875, 0.839216, 0.388235, 0.286275), (0.75, 0.760784, 0.294118, 0.211765),
    (0.78125, 0.701961, 0.211765, 0.168627), (0.8125, 0.65098, 0.156863, 0.129412), (0.84375, 0.6, 0.094118, 0.094118),
    (0.875, 0.54902, 0.066667, 0.098039), (0.90625, 0.501961, 0.05098, 0.12549), (0.9375, 0.45098, 0.054902, 0.172549),
    (0.96875, 0.4, 0.054902, 0.192157), (1, 0.34902, 0.070588, 0.211765)]
_PARAVIEW_LAB = {"cooltowarmextended": ("CoolToWarmExtended", _COOL_TO_WARM_EXT)}
_D65 = np.array([0.95047, 1.0, 1.08883])
_M = np.array([[0.4124, 0.3576, 0.1805], [0.2126, 0.7152, 0.0722], [0.0193, 0.1192, 0.9505]])


def _rgb2lab(c):
    c = np.where(c > 0.04045, ((c + 0.055) / 1.055) ** 2.4, c / 12.92)
    f = (_M @ c) / _D65
    f = np.where(f > 0.008856, np.cbrt(f), 7.787 * f + 16 / 116)
    return np.array([116 * f[1] - 16, 500 * (f[0] - f[1]), 200 * (f[1] - f[2])])


def _lab2rgb(L):
    fy = (L[0] + 16) / 116
    f = np.array([fy + L[1] / 500, fy, fy - L[2] / 200])
    xyz = np.where(f ** 3 > 0.008856, f ** 3, (f - 16 / 116) / 7.787) * _D65
    c = np.linalg.solve(_M, xyz)
    c = np.where(c > 0.0031308, 1.055 * np.abs(c) ** (1 / 2.4) - 0.055, 12.92 * c)
    return np.clip(c, 0, 1)


def _lab_cmap(name, pts, n=256):
    """ListedColormap sampled from (x, r, g, b) control points interpolated in CIELAB, as ParaView does."""
    x = np.array([p[0] for p in pts])
    lab = [_rgb2lab(np.array(p[1:])) for p in pts]
    cols = []
    for s in np.linspace(0, 1, n):
        k = min(max(np.searchsorted(x, s, side="right") - 1, 0), len(x) - 2)
        w = 0.0 if x[k + 1] == x[k] else (s - x[k]) / (x[k + 1] - x[k])
        cols.append(_lab2rgb((1 - w) * lab[k] + w * lab[k + 1]))
    return ListedColormap(cols, name=name)


def get_cmap(name):
    """ParaView Lab maps defined here (e.g. CoolToWarmExtended, "Cool to Warm (Extended)"), else resolve_cmap(); _r reverses."""
    rev = name.endswith("_r")
    key = re.sub(r"[^a-z]", "", (name[:-2] if rev else name).lower())
    if key in _PARAVIEW_LAB:
        cm = _lab_cmap(*_PARAVIEW_LAB[key])
        return cm.reversed() if rev else cm
    return resolve_cmap(name)


def snapshots(d):
    """{time: path} of the iter_N.pvtu files in run directory d."""
    if not os.path.isdir(d):
        sys.exit(f"{d}: no such directory")
    out = {}
    for n in sorted(os.listdir(d), key=lambda n: step_of(n)):
        if re.fullmatch(r"iter_\d+\.pvtu", n):
            t = pvtu_time(os.path.join(d, n))
            if t is not None:
                out[t] = os.path.join(d, n)
    if not out:
        sys.exit(f"{d}: no iter_N.pvtu with a TimeValue")
    return out


def pick(snaps, t, tol):
    k = min(snaps, key=lambda s: abs(s - t))
    return snaps[k] if abs(k - t) <= tol else None


def _title(name):
    return "$" + name[0] + "_" + name[1:] + "$" if name in ("Bx", "By", "Bz") else name


def save(fig, out, dpi):
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    fig.savefig(out, dpi=dpi, bbox_inches="tight")


def tighten(fig, axs, hpad, iters=6):
    """Minimal vertical spacing: tiny constrained-layout pads, then shrink the figure height by the
    slack that equal-aspect axes leave in their layout slots, until the rows (and colorbar) sit tight."""
    fig.get_layout_engine().set(h_pad=hpad, hspace=0.0)
    col = np.atleast_2d(axs)[:, 0]
    for _ in range(iters):
        fig.draw_without_rendering()
        W, H = fig.get_size_inches()
        slack = sum(ax.get_position(original=True).height - ax.get_position().height for ax in col) * H
        if slack < 1e-3:
            break
        fig.set_size_inches(W, H - slack)
    fig.draw_without_rendering()


def _label(name):
    """Pretty label for a field name (ρ -> $\\rho$, Bx -> $B_x$, anything else as is)."""
    if name in ("ρ", "rho"):
        return r"$\rho$"
    return _title(name)


def _find(g, raw, f):
    """Field array for variable raw, accepting the raw name or its greek() form."""
    for n in (greek(raw), raw):
        if n in g["fields"]:
            return g["fields"][n]
    sys.exit(f"{f}: no point field {raw!r}. Available: " + ", ".join(g["names"]))


def timegrid(a, dirs, tags, labels, cmap):
    """One figure per run: rows = a.times, columns = a.rowvars, one colorbar per variable at the bottom."""
    nv, nt = len(a.rowvars), len(a.times)
    if a.rowcmaps and len(a.rowcmaps) != nv:
        sys.exit("--gridcmaps needs one colormap per --gridvars entry")
    names = a.rowcmaps or [a.mu_cmap if v.startswith("mu") else None for v in a.rowvars]
    try:
        cmaps = [get_cmap(c) if c else cmap for c in names]
    except ValueError as e:
        sys.exit(str(e))
    want = set(a.rowvars) | {greek(v) for v in a.rowvars}
    for d, tag, lab in zip(dirs, tags, labels):
        snaps = snapshots(d)
        tol = 1e-6 * max(1.0, max(snaps))
        grid = []  # grid[i][j] = (g, field, extent) for variable i at time j
        cols = []
        for t in a.times:
            f = pick(snaps, t, tol)
            if f is None:
                sys.exit(f"{d}: no snapshot at t = {t:g} (has t = {', '.join(f'{s:g}' for s in sorted(snaps))})")
            g = read_pvtu(f, want)
            cols.append((g, f, (g["x"].min(), g["x"].max(), g["y"].min(), g["y"].max())))
            print(f" t = {t:g}: {f} ({len(g['x'])} nodes)")
        for v in a.rowvars:
            grid.append([(g, _find(g, v, f), ext) for g, f, ext in cols])

        # one color range per variable, shared across the times
        norms = []
        for row in grid:
            vmin = min(np.nanmin(fld) for _, fld, _ in row)
            vmax = max(np.nanmax(fld) for _, fld, _ in row)
            if a.symmetric and vmin < 0 < vmax:
                vmax = max(abs(vmin), abs(vmax)); vmin = -vmax
            vmax = vmax if vmax > vmin else vmin + 1.0
            norms.append(Normalize(vmin, vmax))

        x0, x1, y0, y1 = cols[0][2]
        aspect = (y1 - y0) / (x1 - x0)
        # rows = times (top to bottom), columns = variables (left to right)
        fig, axs = plt.subplots(nt, nv, figsize=(a.width, a.width / nv * nt * aspect + 1.4),
                                sharex=True, sharey=True, squeeze=False, layout="constrained")
        for j, col in enumerate(grid):          # variable j -> column j
            for i, (g, fld, ext) in enumerate(col):   # time i -> row i
                ax = axs[i, j]
                fill(ax, g, lambda g=g: triangles(g), fld, cmaps[j], norms[j].vmin, norms[j].vmax, a.shading)
                ax.set_aspect("equal")
                ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
                if i == 0:
                    ax.set_title(_label(a.rowvars[j]))
                if i == nt - 1:
                    ax.set_xlabel(a.xlabel)
                if j == 0:
                    ax.set_ylabel(f"t = {a.times[i]:g}\n" + a.ylabel)
        fig.suptitle(lab)
        # all colorbars at the very bottom: one under each variable's column, or stacked full width
        cbkw = dict(orientation="horizontal", location="bottom", fraction=0.04, pad=0.0, aspect=20 if not a.stack_colorbars else 40)
        for j in range(nv):
            host = axs if a.stack_colorbars else axs[:, j]
            cb = fig.colorbar(plt.cm.ScalarMappable(norms[j], cmaps[j]), ax=host, **cbkw)
            cb.set_label(_label(a.rowvars[j]))
        tighten(fig, axs, a.hpad)
        out = f"{os.path.join(a.base if not a.dirs else '.', 'OT')}_timegrid_{tag}_" + \
              "_".join(f"t{t:g}" for t in a.times) + f".{a.format}"
        save(fig, out, a.dpi)
        plt.close(fig)
        print(f" -> {out}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--base", default="output/MHD/orszagTangBormanis2024", help="holds output-<R>x<R>/ run directories")
    p.add_argument("--res", nargs="+", type=int, default=None,
                   help="resolutions R: one row each (default 128 256 512), or one figure each with --timegrid (default 64 128 256 512)")
    p.add_argument("--timegrid", action="store_true", help="one figure per resolution: rows = --times, columns = --gridvars")
    p.add_argument("--times", type=float, nargs="+", default=[0.5, 1.0], help="--timegrid rows, top to bottom (default 0.5 1.0)")
    p.add_argument("--gridvars", "--rowvars", dest="rowvars", nargs="+", default=["rho", "mu_rho_Bx_By_Bz"],
                   help="--timegrid columns, left to right (default rho mu_rho_Bx_By_Bz)")
    p.add_argument("--mu-cmap", default="CoolToWarmExtended",
                   help="--timegrid: colormap for variables starting with 'mu' (default ParaView's \"Cool to Warm (Extended)\")")
    p.add_argument("--gridcmaps", "--rowcmaps", dest="rowcmaps", nargs="+", default=None,
                   help="--timegrid: one colormap per --gridvars entry (overrides --cmap and --mu-cmap)")
    p.add_argument("--stack-colorbars", action="store_true",
                   help="--timegrid: stack the colorbars full width instead of one under each column")
    p.add_argument("--dirs", nargs="+", default=None, help="run directories instead of --base/--res (rows, top to bottom)")
    p.add_argument("--labels", nargs="+", default=None, help="row labels (default R×R or the directory names)")
    p.add_argument("--time", type=float, default=None, help="time to plot (default: latest time present in every run)")
    p.add_argument("--vars", nargs=2, default=("Bx", "By"), help="left and right column fields (default Bx By)")
    p.add_argument("--clim", type=float, nargs=2, default=None, help="common color range (default: over all panels)")
    p.add_argument("--symmetric", action="store_true", help="color range symmetric about zero: [-m, m], m = max |range|")
    p.add_argument("--colorbar", action=argparse.BooleanOptionalAction, default=True, help="draw the colorbar (default on)")
    p.add_argument("--single", action="store_true", help="also write one figure per panel (same color range)")
    p.add_argument("--grid", action=argparse.BooleanOptionalAction, default=True, help="write the combined figure (default on)")
    p.add_argument("--cmap", default="inferno", help="inferno (default, ParaView's \"Inferno (matplotlib)\"), RainbowDesaturated, "
                   "CoolToWarmExtended or any matplotlib name; _r reverses")
    p.add_argument("--shading", choices=("contourf", "gouraud"), default="contourf")
    p.add_argument("--out", default=None, help="output file (default <base>/OT_<vars>_t<time>.<format>)")
    p.add_argument("--format", default="png", help="png, pdf, svg, ...")
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--width", type=float, default=8.0, help="figure width in inches")
    p.add_argument("--hpad", type=float, default=0.01, help="vertical gap between rows / above the colorbar, inches")
    p.add_argument("--xlabel", default=r"$x$")
    p.add_argument("--ylabel", default=r"$y$")
    a = p.parse_args()
    try:
        cmap = get_cmap(a.cmap)
    except ValueError as e:
        p.error(str(e))
    names = [greek(v) for v in a.vars]
    if a.res is None:
        a.res = [64, 128, 256, 512] if a.timegrid else [128, 256, 512]

    dirs = a.dirs or [os.path.join(a.base, f"output-{r}x{r}") for r in a.res]
    tags = [f"{r}x{r}" for r in a.res] if not a.dirs else \
           [re.sub(r"^output-", "", os.path.basename(os.path.normpath(d))) for d in dirs]
    labels = a.labels or ([f"{r}" + r"$\times$" + f"{r}" for r in a.res] if not a.dirs else tags)
    if len(labels) != len(dirs):
        p.error("--labels needs one label per run")
    if a.timegrid:
        timegrid(a, dirs, tags, labels, cmap)
        return
    if not (a.grid or a.single):
        p.error("--no-grid needs --single")
    snaps = [snapshots(d) for d in dirs]
    tol = 1e-6 * max(1.0, max(max(s) for s in snaps))
    if a.time is None:
        common = [t for t in snaps[0] if all(pick(s, t, tol) for s in snaps[1:])]
        if not common:
            sys.exit("no output time is present in every run; pass --time")
        a.time = max(common)
    files = [pick(s, a.time, tol) for s in snaps]
    for d, f, s in zip(dirs, files, snaps):
        if f is None:
            sys.exit(f"{d}: no snapshot at t = {a.time:g} (has t = {', '.join(f'{t:g}' for t in sorted(s))})")

    panels = []  # (node grid g, [field per column], extent)
    for f in files:
        g = read_pvtu(f, set(names) | {"ρ"})
        miss = [n for n in names if n not in g["fields"]]
        if miss:
            sys.exit(f"{f}: no point field {miss}. Available: " + ", ".join(g["names"]))
        panels.append((g, [g["fields"][n] for n in names], (g["x"].min(), g["x"].max(), g["y"].min(), g["y"].max())))
        print(f" t = {a.time:g}: {f} ({len(g['x'])} nodes)")
    if a.clim:
        vmin, vmax = a.clim
    else:
        vmin = min(np.nanmin(v) for _, vs, _ in panels for v in vs)
        vmax = max(np.nanmax(v) for _, vs, _ in panels for v in vs)
    if a.symmetric:
        vmax = max(abs(vmin), abs(vmax)); vmin = -vmax
    vmax = vmax if vmax > vmin else vmin + 1.0
    norm = Normalize(vmin, vmax)

    nr = len(panels)
    x0, x1, y0, y1 = panels[0][2]
    aspect = (y1 - y0) / (x1 - x0)
    stem = os.path.join(a.base if not a.dirs else ".", "OT")
    cbkw = dict(orientation="horizontal", location="bottom", fraction=0.04, pad=0.0, aspect=40)
    if a.grid:
        # start deliberately tall; tighten() then removes the slack
        fig, axs = plt.subplots(nr, 2, figsize=(a.width, a.width / 2 * nr * aspect + (1.0 if a.colorbar else 0.6)),
                                sharex=True, sharey=True, squeeze=False, layout="constrained")
        for i, (g, vs, ext) in enumerate(panels):
            for j, v in enumerate(vs):
                ax = axs[i, j]
                fill(ax, g, lambda g=g: triangles(g), v, cmap, vmin, vmax, a.shading)
                ax.set_aspect("equal")
                ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
                if i == 0:
                    ax.set_title(_title(names[j]))
                if i == nr - 1:
                    ax.set_xlabel(a.xlabel)
                if j == 0:
                    ax.set_ylabel(labels[i] + "\n" + a.ylabel)
        fig.suptitle(f"t = {a.time:g}")
        if a.colorbar:
            fig.colorbar(plt.cm.ScalarMappable(norm, cmap), ax=axs, **cbkw)
        tighten(fig, axs, a.hpad)
        out = a.out or f"{stem}_{'_'.join(ascii_name(n) for n in names)}_t{a.time:g}.{a.format}"
        save(fig, out, a.dpi)
        plt.close(fig)
        print(f" -> {out}")
    if a.single:
        root = os.path.splitext(a.out)[0] if a.out else stem
        for i, (g, vs, ext) in enumerate(panels):
            for j, v in enumerate(vs):
                w = a.width / 2
                fig, ax = plt.subplots(figsize=(w, w * aspect + (1.0 if a.colorbar else 0.5)), layout="constrained")
                fill(ax, g, lambda g=g: triangles(g), v, cmap, vmin, vmax, a.shading)
                ax.set_aspect("equal")
                ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
                ax.set_xlabel(a.xlabel); ax.set_ylabel(a.ylabel)
                ax.set_title(f"{_title(names[j])}   {labels[i]}   t = {a.time:g}")
                if a.colorbar:
                    fig.colorbar(plt.cm.ScalarMappable(norm, cmap), ax=ax, **cbkw)
                tighten(fig, np.array([[ax]]), a.hpad)
                out = f"{root}_{ascii_name(names[j])}_{tags[i]}_t{a.time:g}.{a.format}"
                save(fig, out, a.dpi)
                plt.close(fig)
                print(f" -> {out}")


if __name__ == "__main__":
    main()
