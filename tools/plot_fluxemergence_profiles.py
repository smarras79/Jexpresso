#!/usr/bin/env python3
"""Alfvén speed V_A/C_s = |B|/sqrt(rho) (or --var vz: V_z/C_s) along x = X_max/2 at several times, as in Son et al. (2025)
Figs. 5-6, with Shibata et al. (1989a)'s law a (z - z0)/H0 (a2 = 0.3 for V_A, a1 = 0.062 for V_z) and z_cor marked.
    python3 tools/plot_fluxemergence_profiles.py output/MHD/fluxEmergenceSon2025DSGS/output-<date> --times 33 40 47 51 54
"""
import argparse
import glob
import os
import re
import sys
import xml.etree.ElementTree as ET

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_fluxemergence_vtu import read_pvtu, grid2d, pvtu_time, step_of  # noqa: E402

RAMP = {1: ["#2a78d6"], 2: ["#86b6ef", "#0d366b"], 3: ["#86b6ef", "#2a78d6", "#0d366b"],
        4: ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"],
        5: ["#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b"]}          # ordinal blue ramp, validated
INK, INK2, THEORY, GRID = "#0b0b0b", "#52514e", "#8c8b87", "#e2e1dd"
VARS = {"vA": (r"$V_A/C_s$", 0.3, r"$V_A/C_s = a_2\,\Delta z$"),
        "vz": (r"$V_z/C_s$", 0.062, r"$V_z/C_s = a_1\,\Delta z$")}


def outputs(d):
    """[(t, path)] of the iter_N.pvtu in d, times from simulation.pvd when present."""
    files = sorted(glob.glob(os.path.join(d, "iter_*.pvtu")), key=step_of)
    tpvd = {}
    pvd = os.path.join(d, "simulation.pvd")
    if os.path.isfile(pvd):
        for ds in ET.parse(pvd).getroot().iter("DataSet"):
            tpvd[os.path.basename(ds.get("file", ""))] = float(ds.get("timestep"))
    out = []
    for f in files:
        t = tpvd.get(os.path.basename(f))
        out.append((pvtu_time(f) if t is None else t, f))
    return [(t, f) for t, f in out if t is not None]


def profile(path, var, x0):
    keep = {"ρ", "Bx", "By", "Bz"} if var == "vA" else {var if var != "vz" else "v"}
    g = read_pvtu(path, keep)
    F = g["fields"]
    if var == "vA":
        B2 = F["Bx"]**2 + F["By"]**2 + (F["Bz"]**2 if "Bz" in F else 0.0)
        q = np.sqrt(B2 / np.maximum(F["ρ"], 1e-300))
    else:
        q = F["v" if var == "vz" else var]
    Z = grid2d(g, q)
    if Z is None:
        sys.exit(f"{path}: the nodes do not form a tensor-product grid")
    xc = g["xc"]
    x0 = 0.5 * (xc[0] + xc[-1]) if x0 is None else x0
    i = int(np.clip(np.searchsorted(xc, x0), 1, len(xc) - 1))
    w = 1.0 if abs(xc[i] - x0) <= 1e-9 * (1 + abs(x0)) else (x0 - xc[i - 1]) / (xc[i] - xc[i - 1])
    return g["yc"], (1 - w) * Z[:, i - 1] + w * Z[:, i], x0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir", help="output directory with iter_N.pvtu (and simulation.pvd)")
    ap.add_argument("--var", default="vA", help="vA (default, |B|/sqrt(rho)), vz, or any field of the VTU (no law)")
    ap.add_argument("--times", nargs="+", type=float, default=[33, 40, 47, 51, 54],
                    help="output times in tau0 (nearest output used; default the paper's 33 40 47 51 54)")
    ap.add_argument("--x", type=float, help="x of the vertical cut in H0 (default: the middle of the domain, X_max/2)")
    ap.add_argument("--z0", type=float, default=4.0, help="base of the flux sheet, Delta z = (z - z0)/H0 (default 4)")
    ap.add_argument("--a", type=float, help="slope of the law (default a2 = 0.3 for vA, a1 = 0.062 for vz)")
    ap.add_argument("--zcor", type=float, default=18.0, help="height of the transition region (default 18)")
    ap.add_argument("--ylim", nargs=2, type=float)
    ap.add_argument("--width", type=float, default=3.4, help="inches (default 3.4, one journal column)")
    ap.add_argument("--height", type=float, default=2.6)
    ap.add_argument("--out", help="output file (default: profile_<var>.<format> in outdir)")
    ap.add_argument("--format", default="pdf")
    ap.add_argument("--dpi", type=int, default=300)
    a = ap.parse_args()

    outs = outputs(a.outdir)
    if not outs:
        sys.exit(f"{a.outdir}: no iter_N.pvtu with a time")
    tsel = []
    for T in a.times:
        t, f = min(outs, key=lambda o: abs(o[0] - T))
        if (t, f) not in tsel:
            tsel.append((t, f))
            if abs(t - T) > 0.5:
                print(f"t = {T:g}: nearest output is t = {t:g} ({os.path.basename(f)})")
    tsel.sort()
    if len(tsel) > len(RAMP):
        sys.exit(f"at most {len(RAMP)} times per figure (an ordered ramp of more steps is not readable)")
    label, slope, law = VARS.get(a.var, (a.var, None, None))
    slope = a.a if a.a is not None else slope

    plt.rcParams.update({
        "font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 9, "axes.labelsize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7.5, "axes.linewidth": 0.6,
        "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK, "xtick.color": INK2,
        "ytick.color": INK2, "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True,
        "ytick.right": True, "xtick.minor.visible": True, "ytick.minor.visible": True, "lines.linewidth": 1.2,
        "legend.frameon": False, "savefig.bbox": "tight", "savefig.pad_inches": 0.02, "pdf.fonttype": 42,
    })
    fig, ax = plt.subplots(figsize=(a.width, a.height))
    colors = RAMP[len(tsel)]
    ymin, ymax, zmax, x0 = 0.0, 0.0, 0.0, a.x
    for (t, f), c in zip(tsel, colors):
        z, q, x0 = profile(f, a.var, x0)
        ax.plot(z, q, color=c, label=rf"$t = {t:.4g}\,\tau_0$", zorder=3)
        k = int(np.nanargmax(q))
        print(f"t = {t:7.3f}: max {a.var} = {q[k]:.4g} at z = {z[k]:.3g} H0"
              + (f"  (law: {slope * (z[k] - a.z0):.4g})" if slope is not None else ""))
        ymin, ymax, zmax = min(ymin, float(np.nanmin(q))), max(ymax, float(np.nanmax(q))), max(zmax, float(z[-1]))
    ytop = a.ylim[1] if a.ylim else 1.12 * ymax
    if slope is not None:
        zl = np.array([a.z0, min(zmax, a.z0 + ytop / slope)])
        ax.plot(zl, slope * (zl - a.z0), color=THEORY, ls="--", lw=0.9, zorder=2,
                label=(law or rf"${a.a:g}\,\Delta z$") + rf", $z_0 = {a.z0:g}\,H_0$")
    ax.axvline(a.zcor, color=INK, ls="-.", lw=0.8, zorder=1)
    ax.annotate(rf"$z_\mathrm{{cor}} = {a.zcor:g}\,H_0$", xy=(a.zcor, 1), xycoords=("data", "axes fraction"),
                xytext=(3, -3), textcoords="offset points", ha="left", va="top", fontsize=7.5, color=INK2)
    ax.grid(True, which="major", color=GRID, lw=0.4, zorder=0)
    ax.set_xlim(0, zmax)
    ax.set_ylim(*(a.ylim if a.ylim else (0.0 if ymin >= 0 else ax.get_ylim()[0], ytop)))
    ax.set_xlabel(r"$z/H_0$")
    ax.set_ylabel(label)
    ax.legend(loc="best", handlelength=1.8,
              title=rf"$x = {x0:g}\,H_0$", title_fontsize=7.5)
    out = a.out or os.path.join(a.outdir, f"profile_{re.sub(r'[^A-Za-z0-9]+', '', a.var) or 'var'}.{a.format}")
    fig.savefig(out, dpi=a.dpi)
    print("wrote", out)


if __name__ == "__main__":
    main()
