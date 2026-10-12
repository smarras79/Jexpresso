#!/usr/bin/env python3
"""Matrix of 2D MHD solution fields: rows = output times (one run) or runs at one time (several runs, e.g. a
parameter sweep), columns = variables (DynSGS coefficient included), one color range and colorbar per variable.
Reads iter_N.pvtu -> iter_N/iter_N_*.vtu on the native grid.
  python3 tools/plot_mhd_matrix.py output/MHD/blastBalsaraSpicer1999/output-<date> --times 0.002 0.006 0.01 \
          --vars rho p pmag Mach mu_dsgs --log p mu_dsgs --fieldlines
  python3 tools/plot_mhd_matrix.py output/MHD/blastBalsaraSpicer1999/output-<date> --times 0.01 \
          --vars rho p speed Bmag --log rho p                    # Balsara (2004) Fig. 6
  python3 tools/plot_mhd_matrix.py output/MHD/rotorDaoNazarov2022/output-<date> --times 0.05 0.1 0.15 --vars rho p pmag Mach mu_dsgs
  python3 tools/plot_mhd_matrix.py <run1> <run2> <run3> --labels '$\beta = 2.5$' '$\beta = 0.25$' '$\beta = 0.025$' \
          --times 0.01 --vars rho p mu_dsgs --log p mu_dsgs  # several runs: one row each, at one time
  python3 tools/plot_mhd_matrix.py output/MHD/astroJetWuShu2018/output-<date> \
          --vars rho p Mach_f mu_dsgs --log rho p           # a half-domain (AJ_half_*) run is drawn mirrored; --no-mirror to keep the half
  python3 tools/plot_mhd_matrix.py <run> --vars rho p pmag Mach mu_dsgs --cmaps pmag managua --cmaps Mach gist_ncar
  python3 tools/plot_mhd_matrix.py <run> --vars rho log10rho p log10p pmag Mach Bx By mu_dsgs --format pdf
          # log10<field>: computed when the file lacks it; --log VAR instead turns that column itself to log10
Colormaps: per-variable defaults in VAR_CMAPS (rho turbo, p inferno, pmag managua, Mach and Mach_f gist_ncar,
schlieren Yellow - Gray - Blue, schlieren_grad_rho the same reversed (both: shocks dark red), mu_dsgs Cool to Warm Extended); --cmaps VAR NAME overrides one column (ParaView's
XRay, CoolToWarmExtended and "Yellow - Gray - Blue", RainbowDesaturated, or any matplotlib name; _r reverses);
--cmap / --mu-cmap cover all other variables.
Default times: the outputs nearest tend/3, 2 tend/3, tend. Needs numpy and matplotlib.
"""
import argparse, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_fluxemergence_vtu import read_pvtu, triangles, fill, greek, flux_function_integrate, flux_function_lsq, grid2d, _cluster
from plot_orszagtang_vtu import get_cmap, snapshots, pick, tighten, save, _find
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.tri import Triangulation

LABELS = {"ρ": r"$\rho$", "p": r"$p$", "pmag": r"$|\mathbf{B}|^2/2$", "Mach": r"$M$", "mu_dsgs": r"$\nu_\mathrm{DSGS}$",
          "β": r"$\beta$", "ψ": r"$\psi$", "T": r"$T$", "u": r"$u$", "v": r"$v$", "w": r"$w$",
          "Bx": r"$B_x$", "By": r"$B_y$", "Bz": r"$B_z$", "vA": r"$|\mathbf{B}|/\sqrt{\rho}$", "speed": r"$|\mathbf{v}|$",
          "Bmag": r"$|\mathbf{B}|$", "Mach_f": r"$M_f$", "schlieren": "schlieren", "schlieren_grad_rho": r"$|\nabla\rho|$"}
DERIVED = {"vA": ("ρ", "Bx", "By", "Bz"), "speed": ("u", "v", "w"), "Bmag": ("Bx", "By", "Bz")}
ODD = ("u", "Bx")   # fields that change sign under the reflection x -> -x


def half_mesh(d):
    """True when the run's user_inputs.jl (copied into its output directory) names a half-domain mesh (AJ_half_*)."""
    try:
        m = re.search(r':gmsh_filename\s*=>\s*"([^"]+)"', open(os.path.join(d, "user_inputs.jl")).read())
    except OSError:
        return False
    return bool(m) and "half" in os.path.basename(m.group(1))


def mirror(g):
    """Reflect a half-domain solution about its left edge x = x_min (u and Bx change sign); edge nodes are not duplicated."""
    x, n = g["x"], len(g["x"])
    tol = 1e-8 * np.ptp(x)
    new = np.flatnonzero(x > x.min() + tol)
    m = np.arange(n); m[new] = n + np.arange(len(new))
    X = np.concatenate([x, 2 * x.min() - x[new]]); Y = np.concatenate([g["y"], g["y"][new]])
    F = {k: np.concatenate([v, -v[new] if k in ODD else v[new]]) for k, v in g["fields"].items()}
    ix, xc = _cluster(X, tol); iy, yc = _cluster(Y, tol)
    return dict(g, x=X, y=Y, quads=np.vstack([g["quads"], m[g["quads"]][:, ::-1]]), ix=ix, iy=iy, xc=xc, yc=yc, fields=F)


# per-variable colormaps (rotor figure style), keyed by greek(name); --cmaps VAR NAME overrides,
# --cmap / --mu-cmap cover variables not listed here
VAR_CMAPS = {"ρ": "turbo", "rho": "turbo", "p": "inferno", "pmag": "managua", "Mach": "gist_ncar",
             "Mach_f": "gist_ncar", "schlieren": "YellowGreyBlue", "schlieren_grad_rho": "YellowGreyBlue_r",
             "mu_dsgs": "CoolToWarmExtended"}


def resolve_cmap(name):
    """Matplotlib colormap by name if matplotlib has it, else the ParaView presets of plot_orszagtang_vtu.get_cmap."""
    import matplotlib
    if name in matplotlib.colormaps:
        return matplotlib.colormaps[name]
    if name == "managua":   # needs matplotlib >= 3.10
        print(" warning: colormap 'managua' needs matplotlib >= 3.10; using 'RdYlBu_r' for pmag")
        return matplotlib.colormaps["RdYlBu_r"]
    return get_cmap(name)


def log_base(v):
    """'log10p', 'log10rho', 'log10_mu_dsgs' -> the field under the log10; None for any other name."""
    return v[5:].lstrip("_") if v.startswith("log10") and len(v) > 5 else None


def label(v, log):
    if log_base(v):
        return label(log_base(v), True)
    s = LABELS.get(greek(v), LABELS.get(v, v))
    return r"$\log_{10}$ " + s if log else s


def field(g, v, f):
    """Point field v of g (ASCII aliases accepted), or a derived one: vA = |B|/sqrt(rho), speed = |v|, Bmag = |B|."""
    F = g["fields"]
    b = log_base(v)
    if b and greek(v) not in F and v not in F:   # log10 of any field, clipped like --log
        q = np.asarray(field(g, b, f), dtype=float)
        return np.log10(np.maximum(q, max(np.nanmax(q), 1e-300) * 1e-6))
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
    p.add_argument("outdirs", nargs="+", help="run directory with iter_N.pvtu; several = one row per run at one time")
    p.add_argument("--labels", nargs="+", help="row labels of several runs (default: the directory names)")
    p.add_argument("--times", type=float, nargs="+", help="rows, top to bottom (default: nearest tend/3, 2 tend/3, tend); "
                   "with several runs the first value only (default: the latest time all runs have)")
    p.add_argument("--vars", nargs="+", default=["rho", "p", "pmag", "Mach", "mu_dsgs"],
                   help="columns: any point field (rho, psi, beta accepted), vA = |B|/sqrt(rho), speed = |v|, Bmag = |B|, "
                        "and log10<field> (e.g. log10rho, log10p) when the file does not have it")
    p.add_argument("--log", nargs="*", default=[], help="variables drawn as log10 (values below max·1e-6 clipped)")
    p.add_argument("--clim", nargs=3, action="append", default=[], metavar=("VAR", "MIN", "MAX"),
                   help="color range of one variable (after log10 if --log), repeatable")
    p.add_argument("--fieldlines", action="store_true", help="magnetic field lines on every column except mu_*")
    p.add_argument("--mirror", action="store_true",
                   help="reflect the run about its left edge (the symmetry axis), for the plot only; automatic when the "
                        "run's mesh is a half-domain one (AJ_half_*)")
    p.add_argument("--no-mirror", action="store_true", help="never mirror, even a half-domain run")
    p.add_argument("--nlevels", type=int, default=25)
    p.add_argument("--linewidth", type=float, default=0.4)
    p.add_argument("--cmap", default="inferno",
                   help="solution fields not in VAR_CMAPS (default inferno; RainbowDesaturated, any matplotlib name)")
    p.add_argument("--mu-cmap", default="CoolToWarmExtended",
                   help="mu_* columns not in VAR_CMAPS (default ParaView's Cool to Warm (Extended))")
    p.add_argument("--cmaps", nargs=2, action="append", default=[], metavar=("VAR", "CMAP"),
                   help="colormap of one variable, repeatable (overrides VAR_CMAPS, --cmap and --mu-cmap)")
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

    runs = [(d, snapshots(d)) for d in a.outdirs]
    if len(runs) == 1:
        d, snaps = runs[0]
        ts = sorted(snaps)
        times = a.times or sorted({min(ts, key=lambda s: abs(s - f * ts[-1])) for f in (1 / 3, 2 / 3, 1)})
        items = [(d, snaps, t, f"$t = {t:g}$") for t in times]
    else:
        t = a.times[0] if a.times else min(max(s) for _, s in runs)
        labs = a.labels or [os.path.basename(os.path.normpath(d)) for d, _ in runs]
        if len(labs) != len(runs):
            sys.exit("--labels needs one label per run directory")
        items = [(d, s, t, lab) for (d, s), lab in zip(runs, labs)]
    logs = {greek(v) for v in a.log}
    clims = {greek(v): (float(lo), float(hi)) for v, lo, hi in a.clim}
    cmap_over = {greek(v): c for v, c in a.cmaps}
    cols = {greek(v) for v in a.vars} | {greek(log_base(v)) for v in a.vars if log_base(v)}
    for opt, names in (("--log", a.log), ("--clim", [c[0] for c in a.clim]), ("--cmaps", [c[0] for c in a.cmaps])):
        stray = [v for v in names if greek(v) not in cols]
        if stray:   # e.g. "--log rho p u v": u, v were meant as columns
            print(f" warning: {opt} {' '.join(stray)}: not a column (--vars {' '.join(a.vars)}); put columns after --vars")
    base = [log_base(v) or v for v in a.vars]
    want = ({greek(v) for v in a.vars + base} | set(a.vars) |
            {f for v in base for f in DERIVED.get(v, ())})
    if a.fieldlines:
        want |= {"Bx", "By"}

    rows = []
    for d, snaps, t, rlab in items:
        f = pick(snaps, t, 1e-6 * max(1.0, max(snaps)))
        if f is None:
            sys.exit(f"{d}: no snapshot at t = {t:g} (has t = {', '.join(f'{s:g}' for s in sorted(snaps))})")
        g = read_pvtu(f, want)
        mir = a.mirror or (half_mesh(d) and not a.no_mirror)
        xaxis = g["x"].min()
        g = mirror(g) if mir else g
        print(f" t = {t:g}: {f} ({len(g['x'])} nodes" + (f", mirrored about x = {xaxis:g} for the plot)" if mir else ")"))
        cols = []
        for v in a.vars:
            q = np.asarray(field(g, v, f), dtype=float)
            if greek(v) in logs:
                q = np.log10(np.maximum(q, max(np.nanmax(q), 1e-300) * 1e-6))
            cols.append(q)
        rows.append((rlab, g, cols))

    nv, nt = len(a.vars), len(rows)
    norms, cmaps = [], []
    for j, v in enumerate(a.vars):
        lo = min(np.nanmin(r[2][j]) for r in rows)
        hi = max(np.nanmax(r[2][j]) for r in rows)
        lo, hi = clims.get(greek(v), (lo, hi))
        norms.append(Normalize(lo, hi if hi > lo else lo + 1.0))
        b = log_base(v) or v
        gv = greek(b)
        cname = (cmap_over.get(greek(v)) or cmap_over.get(gv) or VAR_CMAPS.get(gv) or VAR_CMAPS.get(b)
                 or (a.mu_cmap if b.startswith("mu") else a.cmap))
        print(f" colormap {v}: {cname}")
        cmaps.append(resolve_cmap(cname))

    g0 = rows[0][1]
    x0, x1, y0, y1 = g0["x"].min(), g0["x"].max(), g0["y"].min(), g0["y"].max()
    W = a.width or 2.2 * nv
    fig, axs = plt.subplots(nt, nv, figsize=(W, W / nv * nt * (y1 - y0) / (x1 - x0) + 1.2),
                            sharex=True, sharey=True, squeeze=False, layout="constrained")
    for i, (rlab, g, cols) in enumerate(rows):
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
                ax.set_ylabel(rlab + "\n" + a.ylabel)
    if a.title or len(runs) > 1:
        fig.suptitle(a.title or f"$t = {items[0][2]:g}$")
    for j in range(nv):
        cb = fig.colorbar(plt.cm.ScalarMappable(norms[j], cmaps[j]), ax=axs[:, j], orientation="horizontal",
                          location="bottom", fraction=0.04, pad=0.0, aspect=12)
        cb.ax.tick_params(labelsize=7)
        cb.ax.xaxis.get_major_locator().set_params(nbins=4)
    tighten(fig, axs, a.hpad)
    base = a.outdirs[0] if len(runs) == 1 else os.path.dirname(os.path.normpath(a.outdirs[0]))
    out = a.out or os.path.join(base, "matrix_" + "_".join(greek(v) for v in a.vars).replace("ρ", "rho") + "." + a.format)
    save(fig, out, a.dpi)
    print(f" -> {out}")


if __name__ == "__main__":
    main()
