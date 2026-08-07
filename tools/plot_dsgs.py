#!/usr/bin/env python3
"""
Comparative filled-contour figures for Jexpresso DynSGS runs.

Produces two kinds of PNG from Jexpresso VTK output, for any case in 1D or 2D:

  compare : the same fields from two runs side by side -- typically DynSGS on
            vs DynSGS off (or vs another stabilization) -- on a SHARED colour
            scale, plus a difference panel, so the effect of the model is the
            thing you actually see.

  coeff   : the per-equation DynSGS eddy viscosity mu_dsgs_* written by
            src/io/write_output.jl, one panel per equation.

Nothing here is case-specific: fields are discovered from the VTU, and 1D vs
2D is detected from the point coordinates. It works for sod1d, theta_dsgs,
orszagTangBormanis2024 or anything else that writes Jexpresso VTK.

--------------------------------------------------------------------------
Usage
--------------------------------------------------------------------------

  # DynSGS vs no DynSGS, last snapshot, auto-selected fields
  tools/plot_dsgs.py compare \\
      --a output/CompEuler/theta_dsgs/output      --label-a "DynSGS" \\
      --b output/CompEuler/theta_nodsgs/output    --label-b "no SGS" \\
      --out dsgs_theta.png

  # pick fields and a snapshot
  tools/plot_dsgs.py compare --a RUN_A --b RUN_B --fields rho,p,Bx --iter 12

  # the DynSGS coefficients themselves, log colour scale
  tools/plot_dsgs.py coeff --a output/MHD/orszagTangBormanis2024/output-... \\
      --log --out dsgs_mu_mhd.png

Field names are matched case- and Unicode-insensitively, so `--fields rho,theta`
finds the fields Jexpresso writes as `ρ` and `θ`.

Requires matplotlib and numpy only -- it parses the VTU appended-binary blocks
directly, so no VTK/pyvista install is needed.
"""

import argparse
import glob
import os
import re
import struct
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
from matplotlib import ticker


# --------------------------------------------------------------------------
# VTU reading
# --------------------------------------------------------------------------

def read_vtu(path):
    """Parse one serial .vtu with raw appended Float64 data -> {name: ndarray}."""
    with open(path, "rb") as fh:
        data = fh.read()
    i = data.find(b"<AppendedData")
    if i < 0:
        raise ValueError(f"{path}: no <AppendedData> block (not a raw-appended VTU?)")
    header = data[:i].decode("utf8", errors="replace")
    start = data.find(b"_", i) + 1

    out = {}
    pat = (r'<DataArray type="(\w+)" Name="([^"]+)" '
           r'NumberOfComponents="(\d+)" format="appended" offset="(\d+)"')
    for m in re.finditer(pat, header):
        typ, name, ncomp, off = m.group(1), m.group(2), int(m.group(3)), int(m.group(4))
        if typ != "Float64":
            continue
        base = start + off
        nbytes = struct.unpack("<Q", data[base:base + 8])[0]
        arr = np.frombuffer(data[base + 8:base + 8 + nbytes], dtype="<f8")
        out[name] = arr.reshape(-1, ncomp) if ncomp > 1 else arr
    return out


def iter_indices(run_dir):
    """Snapshot indices available in a run directory, sorted."""
    idx = []
    for p in glob.glob(os.path.join(run_dir, "iter_*.pvtu")):
        m = re.search(r"iter_(\d+)\.pvtu$", p)
        if m:
            idx.append(int(m.group(1)))
    if not idx:                      # serial runs may only have the piece dirs
        for p in glob.glob(os.path.join(run_dir, "iter_*")):
            m = re.search(r"iter_(\d+)$", p)
            if m and os.path.isdir(p):
                idx.append(int(m.group(1)))
    return sorted(idx)


def read_snapshot(run_dir, it):
    """Read snapshot `it`, concatenating all MPI rank pieces."""
    pieces = sorted(glob.glob(os.path.join(run_dir, f"iter_{it}", f"iter_{it}_*.vtu")))
    if not pieces:
        raise FileNotFoundError(f"no VTU pieces for iter {it} under {run_dir}")
    parts = [read_vtu(p) for p in pieces]
    if len(parts) == 1:
        return parts[0]
    merged = {}
    for k in parts[0]:
        if all(k in p for p in parts):
            merged[k] = np.concatenate([p[k] for p in parts], axis=0)
    return merged


# --------------------------------------------------------------------------
# Field selection
# --------------------------------------------------------------------------

# VTK arrays that are geometry/bookkeeping, not solution fields
_SKIP = {"Points", "connectivity", "offsets", "types", "part"}

# so `--fields rho` matches the field Jexpresso writes as `ρ`
_ALIAS = {
    "rho": "ρ", "theta": "θ", "psi": "ψ", "mu": "μ", "kappa": "κ",
    "rhou": "ρu", "rhov": "ρv", "rhow": "ρw", "rhoe": "ρE", "rhotheta": "ρθ",
}


def _norm(name):
    n = name.lower()
    for k, v in _ALIAS.items():
        n = n.replace(k, v)
    return n


def solution_fields(arrays):
    return [k for k in arrays
            if k not in _SKIP and not k.startswith("mu_dsgs")
            and arrays[k].ndim == 1]


def coeff_fields(arrays):
    return sorted(k for k in arrays if k.startswith("mu_dsgs"))


def resolve_fields(requested, available):
    """Map user-requested names onto what the file actually contains."""
    if not requested:
        return available
    chosen = []
    for want in requested:
        hit = next((a for a in available if _norm(a) == _norm(want)), None)
        if hit is None:
            hit = next((a for a in available if _norm(want) in _norm(a)), None)
        if hit is None:
            print(f"  ! field '{want}' not found; have: {', '.join(available)}",
                  file=sys.stderr)
        elif hit not in chosen:
            chosen.append(hit)
    return chosen


# --------------------------------------------------------------------------
# Plot helpers
# --------------------------------------------------------------------------

def geometry(arrays):
    """Return (x, y, is_1d). 1D is detected from a degenerate y-extent."""
    pts = arrays["Points"]
    x, y = pts[:, 0], pts[:, 1]
    span_x = x.max() - x.min()
    span_y = y.max() - y.min()
    return x, y, (span_y <= 1e-12 * max(span_x, 1.0))


def _contour(ax, x, y, v, vmin, vmax, cmap, nlev, log=False, diverging=False):
    if log:
        pos = v[v > 0]
        if pos.size == 0:
            ax.text(0.5, 0.5, "identically zero", ha="center", va="center",
                    transform=ax.transAxes)
            ax.set_xticks([]); ax.set_yticks([])
            return None
        lo = max(pos.min(), pos.max() * 1e-6)
        levels = np.logspace(np.log10(lo), np.log10(pos.max()), nlev)
        return ax.tricontourf(x, y, np.clip(v, lo, None), levels=levels,
                              norm=LogNorm(vmin=lo, vmax=pos.max()),
                              cmap=cmap, extend="both")
    if vmin == vmax:
        vmin, vmax = vmin - 0.5, vmax + 0.5
    if diverging:
        m = max(abs(vmin), abs(vmax))
        levels = np.linspace(-m, m, nlev)
        return ax.tricontourf(x, y, v, levels=levels,
                              norm=TwoSlopeNorm(0.0, -m, m), cmap=cmap, extend="both")
    levels = np.linspace(vmin, vmax, nlev)
    return ax.tricontourf(x, y, v, levels=levels, cmap=cmap, extend="both")


def _finish_2d(ax, x, y, title):
    ax.set_title(title, fontsize=10)
    ax.set_aspect("equal")
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())
    ax.tick_params(labelsize=7)


# --------------------------------------------------------------------------
# compare
# --------------------------------------------------------------------------

def cmd_compare(args):
    ia = args.iter if args.iter is not None else iter_indices(args.a)[-1]
    ib = args.iter if args.iter is not None else iter_indices(args.b)[-1]
    A, B = read_snapshot(args.a, ia), read_snapshot(args.b, ib)

    fields = resolve_fields(args.fields, [f for f in solution_fields(A)
                                          if f in solution_fields(B)])
    if not fields:
        sys.exit("no common solution fields to compare")

    x, y, is_1d = geometry(A)
    ncol = 2 if (is_1d or args.no_diff) else 3
    nrow = len(fields)
    fig, axes = plt.subplots(nrow, ncol,
                             figsize=(4.2 * ncol, (2.6 if is_1d else 3.6) * nrow),
                             squeeze=False)

    for r, f in enumerate(fields):
        va, vb = A[f], B[f]

        if is_1d:
            o = np.argsort(x)
            ax = axes[r][0]
            ax.plot(x[o], va[o], "-",  lw=1.4, label=args.label_a)
            ax.plot(x[o], vb[o], "--", lw=1.4, label=args.label_b)
            ax.set_ylabel(f, fontsize=10)
            ax.tick_params(labelsize=7)
            if r == 0:
                ax.legend(fontsize=8)
            ax2 = axes[r][1]
            ax2.plot(x[o], (va - vb)[o], "-", lw=1.2, color="crimson")
            ax2.axhline(0.0, color="k", lw=0.5)
            ax2.set_title(f"{args.label_a} − {args.label_b}", fontsize=9)
            ax2.tick_params(labelsize=7)
            continue

        # shared colour scale so the two panels are actually comparable
        vmin = min(va.min(), vb.min())
        vmax = max(va.max(), vb.max())
        for c, (v, lab) in enumerate(((va, args.label_a), (vb, args.label_b))):
            ax = axes[r][c]
            cf = _contour(ax, x, y, v, vmin, vmax, args.cmap, args.levels)
            _finish_2d(ax, x, y, f"{f} — {lab}")
            if cf is not None and c == 1:
                fig.colorbar(cf, ax=axes[r][:2].tolist(), fraction=0.025, pad=0.02)

        if ncol == 3:
            d = va - vb
            ax = axes[r][2]
            cf = _contour(ax, x, y, d, d.min(), d.max(), "RdBu_r",
                          args.levels, diverging=True)
            _finish_2d(ax, x, y, f"{f} — difference")
            if cf is not None:
                fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.02)

    fig.suptitle(args.title or
                 f"{args.label_a} vs {args.label_b}   (snapshot {ia})",
                 fontsize=12)
    if is_1d:
        fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
    print(f"wrote {args.out}  ({nrow} fields, snapshot {ia})")


# --------------------------------------------------------------------------
# coeff
# --------------------------------------------------------------------------

def cmd_coeff(args):
    it = args.iter if args.iter is not None else iter_indices(args.a)[-1]
    A = read_snapshot(args.a, it)

    fields = resolve_fields(args.fields, coeff_fields(A))
    if not fields:
        sys.exit("no mu_dsgs_* fields found — was this run made with a DSGS "
                 "visc_model, and is the Jexpresso build recent enough to "
                 "write them?")

    x, y, is_1d = geometry(A)

    if is_1d:
        fig, ax = plt.subplots(figsize=(7.5, 4.2))
        o = np.argsort(x)
        for f in fields:
            if A[f].max() <= 0:
                continue
            ax.plot(x[o], A[f][o], lw=1.3, label=f)
        if args.log:
            ax.set_yscale("log")
        ax.set_xlabel("x"); ax.set_ylabel(r"$\mu_{DSGS}$")
        ax.legend(fontsize=8); ax.tick_params(labelsize=8)
        ax.grid(alpha=0.3)
    else:
        ncol = min(3, len(fields))
        nrow = (len(fields) + ncol - 1) // ncol
        fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.7 * nrow),
                                 squeeze=False)
        for k, f in enumerate(fields):
            ax = axes[k // ncol][k % ncol]
            v = A[f]
            cf = _contour(ax, x, y, v, v.min(), v.max(), args.cmap,
                          args.levels, log=args.log)
            _finish_2d(ax, x, y, f"{f}\nmax={v.max():.3e}  mean={v.mean():.3e}")
            if cf is not None:
                cb = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.02)
                if args.log:
                    cb.locator = ticker.LogLocator(numticks=6)
                    cb.formatter = ticker.LogFormatterSciNotation()
                    cb.update_ticks()
                cb.ax.tick_params(labelsize=7)
        for k in range(len(fields), nrow * ncol):
            axes[k // ncol][k % ncol].axis("off")

    fig.suptitle(args.title or
                 f"DynSGS eddy viscosity per equation   (snapshot {it})",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
    print(f"wrote {args.out}  ({len(fields)} coefficients, snapshot {it})")


# --------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp):
        sp.add_argument("--a", required=True, metavar="DIR",
                        help="run directory containing iter_*.pvtu")
        sp.add_argument("--iter", type=int, default=None,
                        help="snapshot index (default: last)")
        sp.add_argument("--fields", type=lambda s: s.split(","), default=None,
                        help="comma-separated field names (default: all)")
        sp.add_argument("--out", default=None, help="output PNG")
        sp.add_argument("--title", default=None)
        sp.add_argument("--levels", type=int, default=32,
                        help="number of filled-contour levels (default 32)")
        sp.add_argument("--dpi", type=int, default=130)

    c = sub.add_parser("compare", help="two runs side by side")
    common(c)
    c.add_argument("--b", required=True, metavar="DIR", help="second run directory")
    c.add_argument("--label-a", default="A")
    c.add_argument("--label-b", default="B")
    c.add_argument("--cmap", default="viridis")
    c.add_argument("--no-diff", action="store_true",
                   help="omit the difference column")
    c.set_defaults(func=cmd_compare, out_default="dsgs_compare.png")

    k = sub.add_parser("coeff", help="per-equation DynSGS coefficients")
    common(k)
    k.add_argument("--cmap", default="magma")
    k.add_argument("--log", action="store_true",
                   help="logarithmic colour scale (mu spans orders of magnitude)")
    k.set_defaults(func=cmd_coeff, out_default="dsgs_coeff.png")

    args = p.parse_args()
    if args.out is None:
        args.out = args.out_default
    args.func(args)


if __name__ == "__main__":
    main()
