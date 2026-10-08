#!/usr/bin/env python3
"""Orszag-Tang vortex: Bx (left) and By (right) at several resolutions (rows), one shared horizontal colorbar.
Reads the MPI output (iter_N.pvtu -> iter_N/iter_N_*.vtu) of each run on its native grid, at one common time.
  python3 tools/plot_orszagtang_vtu.py [--base output/MHD/orszagTangBormanis2024] [--res 128 256 512] [--time 1.0]
                                       [--symmetric] [--single [--no-grid]] [--no-colorbar]
Needs numpy and matplotlib (reuses tools/plot_fluxemergence_vtu.py)."""
import argparse, gc, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_fluxemergence_vtu import (read_pvtu, triangles, resolve_cmap, pvtu_time, greek, ascii_name, step_of)
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.tri import Triangulation


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


def draw(ax, T, v, cmap, norm, shading):
    """One field on the native triangulation; values outside the color range take the end colors."""
    if shading == "gouraud":
        ax.tripcolor(T, np.clip(v, norm.vmin, norm.vmax), shading="gouraud", cmap=cmap, norm=norm, rasterized=True)
    else:
        pc = ax.tricontourf(T, v, levels=np.linspace(norm.vmin, norm.vmax, 257), cmap=cmap, norm=norm,
                            extend="both", antialiased=False)
        pc.set_rasterized(True)
    ax.set_aspect("equal")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--base", default="output/MHD/orszagTangBormanis2024", help="holds output-<R>x<R>/ run directories")
    p.add_argument("--res", nargs="+", type=int, default=[128, 256, 512], help="resolutions R, one row each (top to bottom)")
    p.add_argument("--dirs", nargs="+", default=None, help="run directories instead of --base/--res (rows, top to bottom)")
    p.add_argument("--labels", nargs="+", default=None, help="row labels (default R×R or the directory names)")
    p.add_argument("--time", type=float, default=None, help="time to plot (default: latest time present in every run)")
    p.add_argument("--vars", nargs=2, default=("Bx", "By"), help="left and right column fields (default Bx By)")
    p.add_argument("--clim", type=float, nargs=2, default=None, help="common color range (default: over all panels)")
    p.add_argument("--symmetric", action="store_true", help="color range symmetric about zero: [-m, m], m = max |range|")
    p.add_argument("--colorbar", action=argparse.BooleanOptionalAction, default=True, help="draw the colorbar (default on)")
    p.add_argument("--single", action="store_true", help="also write one figure per panel (same color range)")
    p.add_argument("--grid", action=argparse.BooleanOptionalAction, default=True, help="write the combined figure (default on)")
    p.add_argument("--cmap", default="RainbowDesaturated", help="RainbowDesaturated (default) or a matplotlib name; _r reverses")
    p.add_argument("--shading", choices=("contourf", "gouraud"), default="contourf")
    p.add_argument("--out", default=None, help="output file (default <base>/OT_<vars>_t<time>.<format>)")
    p.add_argument("--format", default="png", help="png, pdf, svg, ...")
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--width", type=float, default=8.0, help="figure width in inches")
    p.add_argument("--xlabel", default=r"$x$")
    p.add_argument("--ylabel", default=r"$y$")
    a = p.parse_args()
    try:
        cmap = resolve_cmap(a.cmap)
    except ValueError as e:
        p.error(str(e))
    names = [greek(v) for v in a.vars]

    dirs = a.dirs or [os.path.join(a.base, f"output-{r}x{r}") for r in a.res]
    tags = [f"{r}x{r}" for r in a.res] if not a.dirs else \
           [re.sub(r"^output-", "", os.path.basename(os.path.normpath(d))) for d in dirs]
    labels = a.labels or ([f"{r}" + r"$\times$" + f"{r}" for r in a.res] if not a.dirs else tags)
    if len(labels) != len(dirs):
        p.error("--labels needs one label per run")
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

    panels = []  # (Triangulation, [field per column], extent)
    for f in files:
        g = read_pvtu(f)
        miss = [n for n in names if n not in g["fields"]]
        if miss:
            sys.exit(f"{f}: no point field {miss}. Available: " + ", ".join(sorted(g["fields"])))
        panels.append((Triangulation(g["x"], g["y"], triangles(g)), [g["fields"][n] for n in names],
                       (g["x"].min(), g["x"].max(), g["y"].min(), g["y"].max())))
        print(f" t = {a.time:g}: {f} ({len(g['x'])} nodes)")
        del g; gc.collect()
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
    cbkw = dict(orientation="horizontal", location="bottom", fraction=0.04, pad=0.02, aspect=40)
    if a.grid:
        fig, axs = plt.subplots(nr, 2, figsize=(a.width, a.width / 2 * nr * aspect + (1.0 if a.colorbar else 0.6)),
                                sharex=True, sharey=True, squeeze=False, layout="constrained")
        for i, (T, vs, ext) in enumerate(panels):
            for j, v in enumerate(vs):
                ax = axs[i, j]
                draw(ax, T, v, cmap, norm, a.shading)
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
        out = a.out or f"{stem}_{'_'.join(ascii_name(n) for n in names)}_t{a.time:g}.{a.format}"
        save(fig, out, a.dpi)
        plt.close(fig)
        print(f" -> {out}")
    if a.single:
        root = os.path.splitext(a.out)[0] if a.out else stem
        for i, (T, vs, ext) in enumerate(panels):
            for j, v in enumerate(vs):
                w = a.width / 2
                fig, ax = plt.subplots(figsize=(w, w * aspect + (1.0 if a.colorbar else 0.5)), layout="constrained")
                draw(ax, T, v, cmap, norm, a.shading)
                ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
                ax.set_xlabel(a.xlabel); ax.set_ylabel(a.ylabel)
                ax.set_title(f"{_title(names[j])}   {labels[i]}   t = {a.time:g}")
                if a.colorbar:
                    fig.colorbar(plt.cm.ScalarMappable(norm, cmap), ax=ax, **cbkw)
                out = f"{root}_{ascii_name(names[j])}_{tags[i]}_t{a.time:g}.{a.format}"
                save(fig, out, a.dpi)
                plt.close(fig)
                print(f" -> {out}")


if __name__ == "__main__":
    main()
