#!/usr/bin/env python3
"""Mass and energy conservation from Jexpresso's conservation.dat (:conservation_every): |Q(t)-Q(0)|/|Q(0)| on log axes.

All runs are overlaid on the same panels (by default one for mass, one for total energy), one curve per run.
Legend labels are read from the element count in each file name, e.g. blast_conservation_256x256.dat -> $256^2$.

    python3 tools/plot_conservation.py blast_conservation_64x64.dat blast_conservation_128x128.dat \
        blast_conservation_256x256.dat --out blast_conservation.pdf

    python3 tools/plot_conservation.py output/MHD/blastBalsaraSpicer1999/blast_conservation_*.dat
"""
import argparse
import os
import re
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]  # fixed order
DASHES = ["-", "--", "-.", (0, (1, 1.2)), (0, (5, 1.5)), (0, (3, 1, 1, 1, 1, 1)), (0, (8, 2)), (0, (2, 2))]
MARKS = ["o", "s", "^", "D", "v", "P", "X", "*"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e2e1dd"
NAMES = {"ρ": ("mass", "M"), "rho": ("mass", "M"), "ρE": ("total energy", "E"), "rhoE": ("total energy", "E"),
         "E": ("total energy", "E"), "ρθ": (r"$\rho\theta$", r"\Theta"), "ρu": ("x-momentum", "P_x"),
         "ρv": ("y-momentum", "P_y"), "ρw": ("z-momentum", "P_z")}
MASS = ["ρ", "rho"]
ENERGY = ["ρE", "rhoE", "E"]
RES = re.compile(r"(\d+)x(\d+)(?:x(\d+))?")


def read(path):
    if os.path.isdir(path):
        path = os.path.join(path, "conservation.dat")
    with open(path, encoding="utf-8") as f:
        head = f.readline().lstrip("#").split()
    d = np.loadtxt(path, comments="#", ndmin=2)
    if head[:2] != ["t", "step"] or d.shape[1] != len(head):
        sys.exit(f"{path}: not a Jexpresso conservation.dat (header {head})")
    return path, head[2:], d[:, 0], d[:, 2:]


def resolution(path):
    """Element counts parsed from the last NxM(xK) in the file name (or its directory name), else None."""
    for name in (os.path.basename(path), os.path.basename(os.path.dirname(os.path.abspath(path)))):
        m = RES.findall(name)
        if m:
            return tuple(int(n) for n in m[-1] if n)
    return None


def label(path):
    n = resolution(path)
    if n is None:
        return os.path.splitext(os.path.basename(path))[0].replace("_", r"\_")
    if len(set(n)) == 1:
        return rf"${n[0]}^{len(n)}$"
    return "$" + r"\times".join(str(k) for k in n) + "$"


def default_vars(cols):
    """Mass and total energy if present, otherwise every tracked column."""
    pick = [next((c for c in group if c in cols), None) for group in (MASS, ENERGY)]
    pick = [c for c in pick if c]
    return pick or cols


def style():
    plt.rcParams.update({
        "font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 9, "axes.labelsize": 9,
        "axes.titlesize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
        "axes.linewidth": 0.6, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": INK2, "ytick.color": INK2, "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True, "xtick.minor.visible": True,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
        "lines.linewidth": 1.1, "legend.frameon": False, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="conservation files, e.g. blast_conservation_256x256.dat (or output dirs)")
    ap.add_argument("--labels", nargs="+", help="legend label per run (default: element count from the file name)")
    ap.add_argument("--vars", nargs="+", help="tracked columns to plot, by name (default: mass and total energy)")
    ap.add_argument("--no-sort", action="store_true", help="keep runs in command-line order instead of by resolution")
    ap.add_argument("--legend", choices=["each", "top"], default="each",
                    help="a legend inside every panel, or one shared legend above the figure")
    ap.add_argument("--legend-title", default="", help="legend title (default: none)")
    ap.add_argument("--absolute", action="store_true", help="plot the totals Q(t) instead of |Q(t)-Q(0)|/|Q(0)|")
    ap.add_argument("--no-eps", action="store_true", help="omit the machine-epsilon reference line")
    ap.add_argument("--markers", action="store_true", help="points only (clearer when the totals move by single ulps)")
    ap.add_argument("--floor", type=float, help="draw values below this (exact zeros included) at it, instead of omitting them")
    ap.add_argument("--layout", choices=["row", "column"], default="row", help="panels side by side or stacked")
    ap.add_argument("--width", type=float, help="figure width in inches (default 6.5 row, 3.4 column)")
    ap.add_argument("--height", type=float, help="panel height in inches (default 2.4)")
    ap.add_argument("--xlabel", default="$t$")
    ap.add_argument("--ylim", nargs=2, type=float, help="common y limits")
    ap.add_argument("--out", help="output file (default: <common prefix>.<format> next to the first run)")
    ap.add_argument("--format", default="pdf")
    ap.add_argument("--dpi", type=int, default=300)
    a = ap.parse_args()

    runs = [read(p) for p in a.runs]
    labels = a.labels or [label(p) for p, *_ in runs]
    if len(labels) != len(runs):
        sys.exit("--labels needs one label per run")
    if len(runs) > len(SERIES):
        sys.exit(f"at most {len(SERIES)} runs per figure: facet the rest into a second figure")
    if not a.no_sort and not a.labels:
        key = [np.prod(resolution(p) or (0,)) for p, *_ in runs]
        order = sorted(range(len(runs)), key=lambda i: key[i])
        runs, labels = [runs[i] for i in order], [labels[i] for i in order]
    names = a.vars or default_vars(runs[0][1])
    for p, nm, *_ in runs:
        miss = [v for v in names if v not in nm]
        if miss:
            sys.exit(f"{p}: no column {miss} (has {nm})")

    style()
    n = len(names)
    w = a.width or (6.5 if a.layout == "row" else 3.4)
    h = a.height or 2.4
    shape = (1, n) if a.layout == "row" else (n, 1)
    fig, axes = plt.subplots(*shape, figsize=(w, h * shape[0]), squeeze=False, sharex=True)
    axes = axes.ravel()
    eps = np.finfo(float).eps
    for ax, v in zip(axes, names):
        title, sym = NAMES.get(v, (v, "Q"))
        exact = []
        for k, ((p, nm, t, Q), lab) in enumerate(zip(runs, labels)):
            q = Q[:, nm.index(v)]
            if a.absolute:
                y = q
            else:
                q0 = q[0] if q[0] != 0 else 1.0
                r = np.abs(q - q[0]) / abs(q0)
                nz = int((r[1:] == 0).sum())
                print(f"{lab:>24s}  {title:<14s} max {r.max():.2e}  final {r[-1]:.2e}  exact zeros {nz}/{len(r) - 1}"
                      f"  ({len(t)} samples, t = {t[0]:g} .. {t[-1]:g})")
                if r.max() == 0 and a.floor is None:
                    exact.append(lab)
                y = np.maximum(r, a.floor) if a.floor is not None else np.where(r > 0, r, np.nan)
            hollow = a.markers and k > 0
            ax.semilogy(t, y, color=SERIES[k], ls="none" if a.markers else DASHES[k], marker=MARKS[k] if a.markers else "o",
                        ms=(3.2 if hollow else 2.4) if a.markers else 1.6, mew=0.6 if hollow else 0,
                        mfc="none" if hollow else SERIES[k], label=lab, zorder=2 + len(runs) - k)
        if exact:
            ax.text(0.03, 0.92, "exactly conserved" + ("" if len(runs) == 1 else ": " + ", ".join(exact)),
                    transform=ax.transAxes, ha="left", va="top", fontsize=8, color=INK2)
            if len(exact) == len(runs) and not a.ylim:
                ax.set_ylim(eps / 10, eps * 100)
        if not a.absolute and not a.no_eps:
            ax.axhline(eps, color=INK2, lw=0.6, ls=":", zorder=1)
            ax.annotate(r"$\varepsilon_{\mathrm{mach}}$", xy=(0, eps), xycoords=("axes fraction", "data"),
                        xytext=(3, 2), textcoords="offset points", ha="left", va="bottom", fontsize=7, color=INK2)
        ax.grid(True, which="major", color=GRID, lw=0.4, zorder=0)
        ax.set_title(title, loc="left", color=INK)
        ax.set_ylabel(rf"${sym}(t)$" if a.absolute else rf"$|{sym}(t)-{sym}(0)|\,/\,|{sym}(0)|$")
        ax.set_xlim(min(r[2][0] for r in runs), max(r[2][-1] for r in runs))
        if a.ylim:
            ax.set_ylim(*a.ylim)
        elif not a.absolute and not exact:
            ys = [np.asarray(l.get_ydata(), float) for l in ax.get_lines()]
            pos = np.concatenate([y[np.isfinite(y) & (y > 0)] for y in ys] + [[eps]])
            lo, hi = np.floor(np.log10(pos.min())), np.ceil(np.log10(pos.max()))
            ax.set_ylim(10.0**lo, 10.0**max(hi, lo + 2))
        if a.legend == "each" and len(runs) >= 2:
            if not a.ylim:  # headroom on the log axis so the in-panel legend sits above the curves
                lo, hi = ax.get_ylim()
                ax.set_ylim(lo, hi * (hi / lo) ** 0.25)
            ax.legend(loc="upper right", ncol=min(len(runs), 4), title=a.legend_title or None, title_fontsize=8,
                      handlelength=2.2, columnspacing=1.0, borderaxespad=0.3,
                      markerscale=2.0 if a.markers else 1.0)
    for ax in (axes if a.layout == "row" else axes[-1:]):
        ax.set_xlabel(a.xlabel)
    fig.tight_layout(w_pad=1.5, h_pad=0.6)
    if a.legend == "top":
        hs, ls = axes[0].get_legend_handles_labels()
        fig.legend(hs, ls, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=min(len(runs), 4), handlelength=2.6,
                   title=a.legend_title or None, title_fontsize=8, markerscale=2.0 if a.markers else 1.0)
    if a.out:
        out = a.out
    else:
        stems = [os.path.splitext(os.path.basename(p))[0] for p, *_ in runs]
        stem = RES.sub("", os.path.commonprefix(stems)).rstrip("_-") if len(runs) > 1 else stems[0]
        out = os.path.join(os.path.dirname(os.path.abspath(runs[0][0])), (stem or "conservation") + "." + a.format)
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig.savefig(out, dpi=a.dpi)
    print("wrote", out)


if __name__ == "__main__":
    main()
