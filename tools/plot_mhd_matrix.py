#!/usr/bin/env python3
"""Matrix of 2D MHD solution fields: rows = output times, columns = variables (DynSGS coefficient included),
one color range and colorbar per variable at the bottom. Reads iter_N.pvtu -> iter_N/iter_N_*.vtu on the native grid.
  python3 tools/plot_mhd_matrix.py output/MHD/blastBalsaraSpicer1999/output-<date> --times 0.002 0.006 0.01 \
          --vars rho p pmag Mach mu_dsgs --log p mu_dsgs --fieldlines
  python3 tools/plot_mhd_matrix.py output/MHD/blastBalsaraSpicer1999/output-<date> --times 0.01 \
          --vars rho p speed Bmag --log rho p                    # Balsara (2004) Fig. 6
  python3 tools/plot_mhd_matrix.py output/MHD/rotorDaoNazarov2022/output-<date> --times 0.05 0.1 0.15 --vars rho p pmag Mach mu_dsgs
Default times: the outputs nearest tend/3, 2 tend/3, tend. Needs numpy and matplotlib.
"""
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_fluxemergence_vtu import read_pvtu, triangles, fill, greek, flux_function_integrate, flux_function_lsq, grid2d
from plot_orszagtang_vtu import get_cmap, snapshots, pick, tighten, save, _find
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.tri import Triangulation

LABELS = {"ρ": r"$\rho$", "p": r"$p$", "pmag": r"$|\mathbf{B}|^2/2$", "Mach": r"$M$", "mu_dsgs": r"$\nu_\mathrm{DSGS}$",
          "β": r"$\beta$", "ψ": r"$\psi$", "T": r"$T$", "u": r"$u$", "v": r"$v$", "w": r"$w$",
          "Bx": r"$B_x$", "By": r"$B_y$", "Bz": r"$B_z$", "vA": r"$|\mathbf{B}|/\sqrt{\rho}$", "speed": r"$|\mathbf{v}|$",
          "Bmag": r"$|\mathbf{B}|$"}
DERIVED = {"vA": ("ρ", "Bx", "By", "Bz"), "speed": ("u", "v", "w"), "Bmag": ("Bx", "By", "Bz")}


def label(v, log):
    s = LABELS.get(greek(v), LABELS.get(v, v))
    return r"$\log_{10}$ " + s if log else s


def field(g, v, f):
    """Point field v of g (ASCII aliases accepted), or a derived one: vA = |B|/sqrt(rho), speed = |v|, Bmag = |B|."""
    F = g["fields"]
    if v == "Bmag":
        return np.sqrt(F["Bx"]**2 + F["By"]**2 + F.get("Bz", 0.0)**2)
    if v == "vA":
        return np.sqrt((F["Bx"]**2 + F["By"]**2 + F.get("Bz", 0.0)**2) / np.maximum(F["ρ"], 1e-300))
    if v == "speed":
        return np.sqrt(F["u"]**2 + F["v"]**2 + F.get("w", 0.0)**2)
    return _find(g, v, f)


def flux_lines(ax, g, n, lw):
    """Isolines of the flux function A (Bx = dA/dy, By = -dA/dx): magnetic field lines."""
    Bx, By = g["fields"]["Bx"], g["fields"]["By"]
    A = flux_function_integrate(g, Bx, By)
    tri = None
    if A is None:
        tri = triangles(g)
        A = flux_function_lsq(g, tri, Bx, By)
    if np.ptp(A) <= 0:
        return
    An, lev = (A - A.min()) / np.ptp(A), np.arange(1, n + 1) / (n + 1)
    Z = grid2d(g, An)
    if Z is not None:
        ax.contour(g["xc"], g["yc"], Z, levels=lev, colors="k", linewidths=lw)
    else:
        ax.tricontour(Triangulation(g["x"], g["y"], tri), An, levels=lev, colors="k", linewidths=lw)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("outdir", help="run directory with iter_N.pvtu")
    p.add_argument("--times", type=float, nargs="+", help="rows, top to bottom (default: nearest tend/3, 2 tend/3, tend)")
    p.add_argument("--vars", nargs="+", default=["rho", "p", "pmag", "Mach", "mu_dsgs"],
                   help="columns: any point field (rho, psi, beta accepted), vA = |B|/sqrt(rho), speed = |v|, Bmag = |B|")
    p.add_argument("--log", nargs="*", default=[], help="variables drawn as log10 (values below max·1e-6 clipped)")
    p.add_argument("--clim", nargs=3, action="append", default=[], metavar=("VAR", "MIN", "MAX"),
                   help="color range of one variable (after log10 if --log), repeatable")
    p.add_argument("--fieldlines", action="store_true", help="magnetic field lines on every column except mu_*")
    p.add_argument("--nlevels", type=int, default=25)
    p.add_argument("--linewidth", type=float, default=0.4)
    p.add_argument("--cmap", default="inferno", help="solution fields (default inferno; RainbowDesaturated, any matplotlib name)")
    p.add_argument("--mu-cmap", default="CoolToWarmExtended", help="mu_* columns (default ParaView's Cool to Warm (Extended))")
    p.add_argument("--shading", choices=("contourf", "gouraud"), default="contourf")
    p.add_argument("--title", default=None)
    p.add_argument("--xlabel", default=r"$x$")
    p.add_argument("--ylabel", default=r"$y$")
    p.add_argument("--width", type=float, default=None, help="inches (default 2.2 per column)")
    p.add_argument("--hpad", type=float, default=0.01)
    p.add_argument("--out", default=None, help="default <outdir>/matrix_<vars>.<format>")
    p.add_argument("--format", default="png")
    p.add_argument("--dpi", type=int, default=200)
    a = p.parse_args()

    snaps = snapshots(a.outdir)
    ts = sorted(snaps)
    times = a.times or sorted({min(ts, key=lambda s: abs(s - f * ts[-1])) for f in (1 / 3, 2 / 3, 1)})
    tol = 1e-6 * max(1.0, ts[-1])
    logs = {greek(v) for v in a.log}
    clims = {greek(v): (float(lo), float(hi)) for v, lo, hi in a.clim}
    want = {greek(v) for v in a.vars} | {f for v in a.vars for f in DERIVED.get(v, ())}
    if a.fieldlines:
        want |= {"Bx", "By"}

    rows = []
    for t in times:
        f = pick(snaps, t, tol)
        if f is None:
            sys.exit(f"{a.outdir}: no snapshot at t = {t:g} (has t = {', '.join(f'{s:g}' for s in ts)})")
        g = read_pvtu(f, want)
        print(f" t = {t:g}: {f} ({len(g['x'])} nodes)")
        cols = []
        for v in a.vars:
            q = np.asarray(field(g, v, f), dtype=float)
            if greek(v) in logs:
                q = np.log10(np.maximum(q, max(np.nanmax(q), 1e-300) * 1e-6))
            cols.append(q)
        rows.append((t, g, cols))

    nv, nt = len(a.vars), len(rows)
    norms, cmaps = [], []
    for j, v in enumerate(a.vars):
        lo = min(np.nanmin(r[2][j]) for r in rows)
        hi = max(np.nanmax(r[2][j]) for r in rows)
        lo, hi = clims.get(greek(v), (lo, hi))
        norms.append(Normalize(lo, hi if hi > lo else lo + 1.0))
        cmaps.append(get_cmap(a.mu_cmap if v.startswith("mu") else a.cmap))

    g0 = rows[0][1]
    x0, x1, y0, y1 = g0["x"].min(), g0["x"].max(), g0["y"].min(), g0["y"].max()
    W = a.width or 2.2 * nv
    fig, axs = plt.subplots(nt, nv, figsize=(W, W / nv * nt * (y1 - y0) / (x1 - x0) + 1.2),
                            sharex=True, sharey=True, squeeze=False, layout="constrained")
    for i, (t, g, cols) in enumerate(rows):
        cache = {}
        for j, q in enumerate(cols):
            ax = axs[i, j]
            fill(ax, g, lambda: cache.setdefault("tri", triangles(g)), q, cmaps[j], norms[j].vmin, norms[j].vmax, a.shading)
            if a.fieldlines and not a.vars[j].startswith("mu"):
                flux_lines(ax, g, a.nlevels, a.linewidth)
            ax.set_aspect("equal")
            ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
            if i == 0:
                ax.set_title(label(a.vars[j], greek(a.vars[j]) in logs))
            if i == nt - 1:
                ax.set_xlabel(a.xlabel)
            if j == 0:
                ax.set_ylabel(f"$t = {t:g}$\n" + a.ylabel)
    if a.title:
        fig.suptitle(a.title)
    for j in range(nv):
        cb = fig.colorbar(plt.cm.ScalarMappable(norms[j], cmaps[j]), ax=axs[:, j], orientation="horizontal",
                          location="bottom", fraction=0.04, pad=0.0, aspect=12)
        cb.ax.tick_params(labelsize=7)
        cb.ax.xaxis.get_major_locator().set_params(nbins=4)
    tighten(fig, axs, a.hpad)
    out = a.out or os.path.join(a.outdir, "matrix_" + "_".join(greek(v) for v in a.vars).replace("ρ", "rho") + "." + a.format)
    save(fig, out, a.dpi)
    print(f" -> {out}")


if __name__ == "__main__":
    main()
