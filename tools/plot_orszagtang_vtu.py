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

Pressure cuts at y = 0.3125 and y = 0.4277, t = 0.5, against Athena (Babbar & Ranocha, JCP 563 (2026) 115060, Fig. 15):
  python3 tools/plot_orszagtang_vtu.py --cuts [--res 128 256 512 | --dirs RUN ...] [--time 0.5] [--nop 4]
"""
import argparse, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_fluxemergence_vtu import (read_pvtu, triangles, resolve_cmap, pvtu_time, greek, ascii_name, step_of, fill,
                                     grid2d)
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
# ParaView's "Yellow - Gray - Blue" (F. Samsel; Lab), from ParaView's Remoting/Views/ColorMaps.json
_YELLOW_GRAY_BLUE = [(0, 0.301961, 0.047059, 0.090196),
    (0.01639024615287782, 0.396078431372549, 0.0392156862745098, 0.058823529411764705),
    (0.03278049230575558, 0.49411764705882355, 0.054901960784313725, 0.03529411764705882),
    (0.049170738458633456, 0.5882352941176471, 0.11372549019607843, 0.023529411764705882),
    (0.06556098461151122, 0.6627450980392157, 0.16862745098039217, 0.01568627450980392),
    (0.08195123076438904, 0.7411764705882353, 0.22745098039215686, 0.00392156862745098),
    (0.09834147691726686, 0.788235294117647, 0.2901960784313726, 0),
    (0.11473172307014465, 0.8627450980392157, 0.3803921568627451, 0.011764705882352941),
    (0.13112196922302247, 0.9019607843137255, 0.4588235294117647, 0.027450980392156862),
    (0.14751221537590026, 0.9176470588235294, 0.5215686274509804, 0.047058823529411764),
    (0.16390246152877808, 0.9254901960784314, 0.5803921568627451, 0.0784313725490196),
    (0.18029270768165587, 0.9372549019607843, 0.6431372549019608, 0.12156862745098039),
    (0.19668295383453369, 0.9450980392156862, 0.7098039215686275, 0.1843137254901961),
    (0.2130731999874115, 0.9529411764705882, 0.7686274509803922, 0.24705882352941178),
    (0.2294634461402893, 0.9647058823529412, 0.8274509803921568, 0.3254901960784314),
    (0.24585369229316711, 0.9686274509803922, 0.8784313725490196, 0.4235294117647059),
    (0.26224393844604493, 0.9725490196078431, 0.9176470588235294, 0.5137254901960784),
    (0.27863418459892275, 0.9803921568627451, 0.9490196078431372, 0.596078431372549),
    (0.2950244307518005, 0.9803921568627451, 0.9725490196078431, 0.6705882352941176),
    (0.31141467690467833, 0.9882352941176471, 0.9882352941176471, 0.7568627450980392),
    (0.32714931321144103, 0.984313725490196, 0.9882352941176471, 0.8549019607843137),
    (0.32780492305755615, 0.9882352941176471, 0.9882352941176471, 0.8588235294117647),
    (0.32781492305755616, 0.9529411764705882, 0.9529411764705882, 0.8941176470588236),
    (0.32781492305755616, 0.9529411764705882, 0.9529411764705882, 0.8941176470588236),
    (0.34483569985938073, 0.8901960784313725, 0.8901960784313725, 0.807843137254902),
    (0.3618564766612053, 0.8274509803921568, 0.8235294117647058, 0.7372549019607844),
    (0.37887725346302986, 0.7764705882352941, 0.7647058823529411, 0.6784313725490196),
    (0.39589803026485443, 0.7254901960784313, 0.7137254901960784, 0.6274509803921569),
    (0.412918807066679, 0.6784313725490196, 0.6627450980392157, 0.5803921568627451),
    (0.42993958386850356, 0.6313725490196078, 0.6078431372549019, 0.5333333333333333),
    (0.4469603606703281, 0.5803921568627451, 0.5568627450980392, 0.48627450980392156),
    (0.46398113747215275, 0.5372549019607843, 0.5058823529411764, 0.44313725490196076),
    (0.4810019142739773, 0.4980392156862745, 0.4588235294117647, 0.40784313725490196),
    (0.4980226910758019, 0.4627450980392157, 0.4196078431372549, 0.37254901960784315),
    (0.5150434678776264, 0.43137254901960786, 0.38823529411764707, 0.34509803921568627),
    (0.532064244679451, 0.403921568627451, 0.3568627450980392, 0.3176470588235294),
    (0.5490850214812756, 0.37254901960784315, 0.3215686274509804, 0.29411764705882354),
    (0.5661057982831001, 0.34509803921568627, 0.29411764705882354, 0.26666666666666666),
    (0.5831265750849247, 0.3176470588235294, 0.2627450980392157, 0.23921568627450981),
    (0.6001473518867493, 0.28627450980392155, 0.23137254901960785, 0.21176470588235294),
    (0.6171681286885738, 0.2549019607843137, 0.2, 0.1843137254901961),
    (0.6341889054903984, 0.23137254901960785, 0.17254901960784313, 0.16470588235294117),
    (0.651209682292223, 0.2, 0.1450980392156863, 0.13725490196078433),
    (0.6682404590940476, 0.14902, 0.196078, 0.278431),
    (0.7014164131846428, 0.2, 0.2549019607843137, 0.34509803921568627),
    (0.734592367275238, 0.24705882352941178, 0.3176470588235294, 0.41568627450980394),
    (0.7677683213658333, 0.3058823529411765, 0.38823529411764707, 0.49411764705882355),
    (0.8009442754564285, 0.37254901960784315, 0.4588235294117647, 0.5686274509803921),
    (0.8341202295470238, 0.44313725490196076, 0.5333333333333333, 0.6431372549019608),
    (0.867296183637619, 0.5176470588235295, 0.615686274509804, 0.7254901960784313),
    (0.9004721377282143, 0.6, 0.6980392156862745, 0.8),
    (0.9336480918188095, 0.6862745098039216, 0.7843137254901961, 0.8705882352941177),
    (0.9668240459094047, 0.7607843137254902, 0.8588235294117647, 0.9294117647058824),
    (0.9834120229547024, 0.807843137254902, 0.9019607843137255, 0.9607843137254902),
    (1, 0.8901960784313725, 0.9568627450980393, 0.984313725490196)]
_PARAVIEW_LAB = {"cooltowarmextended": ("CoolToWarmExtended", _COOL_TO_WARM_EXT),
                 "yellowgrayblue": ("YellowGrayBlue", _YELLOW_GRAY_BLUE), "yellowgreyblue": ("YellowGrayBlue", _YELLOW_GRAY_BLUE)}
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
    """ParaView maps defined here (CoolToWarmExtended, "Yellow - Gray - Blue", "X Ray", any spelling), else resolve_cmap(); _r reverses."""
    rev = name.endswith("_r")
    key = re.sub(r"[^a-z]", "", (name[:-2] if rev else name).lower())
    if key in _PARAVIEW_LAB:
        cm = _lab_cmap(*_PARAVIEW_LAB[key])
        return cm.reversed() if rev else cm
    if key == "xray":   # ParaView's "X Ray": white at the minimum, black at the maximum
        cm = plt.get_cmap("gray_r")
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


# ---------------------------------------------------------------------------------------------------------------------
# --cuts: pressure along y = 0.3125 and y = 0.4277 at t = 0.5 against the Athena reference of Babbar & Ranocha,
# J. Comput. Phys. 563 (2026) 115060, Fig. 15.
#
# THE REFERENCE is problems/MHD/orszagTangBormanis2024/athena/orszag_tang_athena_0_{3125,4277}.csv: the red points of
# Fig. 15, copied from the paper's reproducibility repository (see the README there). They are the Athena solution on
# 1024^2 points of Rueda-Ramirez et al., J. Comput. Phys. 444 (2021) 110580, digitized from a plot.
#
# THE NORMALIZATION. Eq. 82 of the paper prints rho = 1, p = 1/gamma, B = (-sin 2 pi y, sin 4 pi x)/gamma. The run behind
# Figs. 14-15 (run/tang/run_mhd_tang.jl) and the Athena reference, however, use rho = 25/(36 pi), p = 5/(12 pi),
# B0 = 1/sqrt(4 pi), the normalization of orszagTangBormanis2024 (rho = gamma^2/(4 pi), p = gamma/(4 pi), B0 = 1/sqrt(4 pi)).
# So the pressures are compared as they are. A run in the normalization of eq. 82 has a pressure 4 pi/gamma^2 ~= 4.52
# times larger at the same time; --pscale 0.2210 (gamma^2/(4 pi)) brings it back.
#
# THE CUT is the numerical solution itself. The output nodes lie on a tensor-product lattice (the LGL nodes of a uniform
# quad mesh). At every node column x_i, p(x_i, y0) is evaluated with the degree-N polynomial in y of the element row
# that contains y0. Between node columns the curve is evaluated with the degree-N polynomial in x of each element
# (--refine points per node interval). --nop is the polynomial order of the runs (default 4, the deck's). If the lattice
# does not have the LGL structure of that order, the cut falls back to linear interpolation between the two bracketing
# node rows (what ParaView's Plot Over Line does, which made the paper's curves) and says so.
# ---------------------------------------------------------------------------------------------------------------------
ATHENA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "problems", "MHD", "orszagTangBormanis2024",
                          "athena")
CUT_COLORS = ("#2a78d6", "#1baf7a", "#7a5bd1", "#c2452d")   # the runs; the Athena points are ATHENA_COLOR
ATHENA_COLOR = "#eb6834"
ATHENA_TIME = 0.5                                            # the time of the reference cuts


def athena_cut(y0):
    """The Athena points of the cut at y0 as an (n, 2) array (x, p) sorted in x; None if there is no file for y0."""
    f = os.path.join(ATHENA_DIR, "orszag_tang_athena_0_" + f"{y0:.4f}".split(".")[1] + ".csv")
    if not os.path.isfile(f):
        return None
    with open(f) as io:
        xp = np.array([ln.split() for ln in io.read().splitlines()[1:] if ln.strip()], dtype=float)
    return xp[np.argsort(xp[:, 0])]


def lgl_nodes(N):
    """Legendre-Gauss-Lobatto nodes on [-1, 1]: -1, the zeros of P_N' (eigenvalues of the Jacobi matrix of P^(1,1)_(N-1)), 1."""
    k = np.arange(1, N - 1)
    b = np.sqrt(k * (k + 2.0) / ((2.0 * k + 1.0) * (2.0 * k + 3.0)))
    inner = np.linalg.eigvalsh(np.diag(b, 1) + np.diag(b, -1)) if N >= 2 else np.array([])
    return np.concatenate(([-1.0], np.sort(inner), [1.0]))


def lagrange(nodes, s):
    """l_k(s) for the Lagrange basis on nodes, shape (len(s), len(nodes))."""
    s = np.atleast_1d(np.asarray(s, float))
    L = np.ones((len(s), len(nodes)))
    for k, xk in enumerate(nodes):
        for m, xm in enumerate(nodes):
            if m != k:
                L[:, k] *= (s - xm) / (xk - xm)
    return L


def is_lgl_lattice(c, N, tol=1e-6):
    """True if the sorted node coordinates c are the LGL nodes of order N on consecutive elements."""
    n = len(c) - 1
    if N < 1 or n < N or n % N:
        return False
    ref = (lgl_nodes(N) + 1.0) / 2.0
    for e in range(n // N):
        lo, hi = c[e * N], c[(e + 1) * N]
        if np.max(np.abs(c[e * N:(e + 1) * N + 1] - (lo + ref * (hi - lo)))) > tol * (hi - lo):
            return False
    return True


def line_cut(g, v, y0, N, refine):
    """(x, v(x, y0), how) along the line y = y0 of the node grid g."""
    xc, yc = g["xc"], g["yc"]
    Z = grid2d(g, v)
    if Z is None:                                            # not a tensor lattice: linear on the triangles
        from matplotlib.tri import Triangulation, LinearTriInterpolator
        xs = np.linspace(xc[0], xc[-1], 4097)
        f = LinearTriInterpolator(Triangulation(g["x"], g["y"], triangles(g)), v)
        return xs, np.asarray(f(xs, np.full_like(xs, y0))), "linear on the triangles (no node lattice)"
    if not yc[0] <= y0 <= yc[-1]:
        sys.exit(f"y = {y0} is outside the data, y in [{yc[0]}, {yc[-1]}]")
    if is_lgl_lattice(yc, N) and is_lgl_lattice(xc, N):
        ref = lgl_nodes(N)
        e = min(int(np.searchsorted(yc[::N], y0, side="right")) - 1, (len(yc) - 1) // N - 1)
        lo, hi = yc[e * N], yc[(e + 1) * N]
        vc = lagrange(ref, 2.0 * (y0 - lo) / (hi - lo) - 1.0)[0] @ Z[e * N:(e + 1) * N + 1, :]   # at every node column
        s = np.linspace(-1.0, 1.0, refine * N + 1)
        Ls = lagrange(ref, s)
        xs, vs = [], []
        for ex in range((len(xc) - 1) // N):                 # between the columns: the element's polynomial in x
            lo, hi = xc[ex * N], xc[(ex + 1) * N]
            keep = slice(0 if ex == 0 else 1, None)          # the shared end point once
            xs.append((lo + (s + 1.0) / 2.0 * (hi - lo))[keep])
            vs.append((Ls @ vc[ex * N:(ex + 1) * N + 1])[keep])
        return np.concatenate(xs), np.concatenate(vs), f"element polynomials of order {N}"
    j = min(int(np.searchsorted(yc, y0, side="right")) - 1, len(yc) - 2)
    th = (y0 - yc[j]) / (yc[j + 1] - yc[j])
    return xc, (1.0 - th) * Z[j, :] + th * Z[j + 1, :], f"linear in y (not the LGL lattice of order {N}: check --nop)"


def cuts(a, dirs, tags, labels):
    """--cuts: one panel per a.cut_y, one curve per run, the Athena points where the paper has them."""
    t = 0.5 if a.time is None else a.time
    nops = a.nop * len(dirs) if len(a.nop) == 1 else a.nop
    if len(nops) != len(dirs):
        sys.exit("--nop: one value, or one per run")
    files = []
    for d in dirs:
        snaps = snapshots(d)
        f = pick(snaps, t, 1e-6 * max(1.0, t))
        if f is None:
            sys.exit(f"{d}: no snapshot at t = {t:g} (has t = {', '.join(f'{s:g}' for s in sorted(snaps))})")
        files.append(f)
    grids = []
    for f in files:
        g = read_pvtu(f, {a.pfield})
        if a.pfield not in g["fields"]:
            sys.exit(f"{f}: no point field {a.pfield!r}. Available: " + ", ".join(g["names"]))
        grids.append(g)
        print(f" t = {t:g}: {f} ({len(g['x'])} nodes)")
    stem = os.path.splitext(a.out)[0] if a.out else os.path.join(a.base if not a.dirs else ".", f"OT_p_cuts_t{t:g}")
    fig, axs = plt.subplots(1, len(a.cut_y), figsize=(a.width / 2 * len(a.cut_y), a.width * 0.42),
                            squeeze=False, layout="constrained")
    at_ref = abs(t - ATHENA_TIME) <= 1e-6 * max(1.0, t)
    if not at_ref:
        print(f" the Athena cuts are at t = {ATHENA_TIME:g}: not drawn at t = {t:g}")
    for ax, y0 in zip(axs[0], a.cut_y):
        ref = athena_cut(y0) if at_ref else None
        if ref is not None:
            ax.plot(ref[:, 0], ref[:, 1], "o", ms=3.0, mfc="none", mec=ATHENA_COLOR, mew=0.9, zorder=3,
                    label=r"Athena, $1024^2$ (Babbar & Ranocha 2026)")
        print(f" y = {y0:g}: " + (f"{len(ref)} Athena points" if ref is not None else "no Athena reference drawn"))
        for i, (g, tag, lab, N) in enumerate(zip(grids, tags, labels, nops)):
            x, v, how = line_cut(g, a.pscale * g["fields"][a.pfield], y0, N, a.refine)
            msg = f"   {tag}: {how}"
            if ref is not None:
                dv = np.interp(ref[:, 0], x, v) - ref[:, 1]
                k = int(np.argmax(np.abs(dv)))
                msg += (f"; against Athena: mean |dp| = {np.mean(np.abs(dv)):.4f} "
                        f"({100 * np.mean(np.abs(dv)) / np.ptp(ref[:, 1]):.1f}% of its range), "
                        f"max |dp| = {abs(dv[k]):.4f} at x = {ref[k, 0]:.3f}")
            print(msg)
            ax.plot(x, v, "-", lw=1.3, color=CUT_COLORS[i % len(CUT_COLORS)], zorder=2, label=f"Jexpresso, {lab}")
            with open(f"{stem}_y{y0:g}_{tag}.csv", "w") as io:
                io.write("x p\n")
                io.writelines(f"{xi:.10g} {vi:.10g}\n" for xi, vi in zip(x, v))
        ax.set_title(f"$y = {y0:g}$,  $t = {t:g}$")
        ax.set_xlabel(a.xlabel)
        ax.set_ylabel(r"$p$")
        ax.set_xlim(min(g["xc"][0] for g in grids), max(g["xc"][-1] for g in grids))
        ax.ticklabel_format(axis="y", useOffset=False)
        ax.grid(True, ls="--", lw=0.5, alpha=0.5)
        ax.legend(fontsize=8, frameon=False)
    out = a.out or f"{stem}.{a.format}"
    save(fig, out, a.dpi)
    plt.close(fig)
    print(f" -> {out} (the cuts: {stem}_y<y>_<run>.csv)")


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
    p.add_argument("--cuts", action="store_true",
                   help="pressure along --cut-y at --time (default 0.5), one curve per run, against Athena (Babbar & Ranocha 2026, Fig. 15)")
    p.add_argument("--cut-y", type=float, nargs="+", default=[0.3125, 0.4277],
                   help="--cuts: the lines y = const (default 0.3125 0.4277, the two with Athena data)")
    p.add_argument("--nop", type=int, nargs="+", default=[4], help="--cuts: polynomial order of the runs, one or one per run (default 4)")
    p.add_argument("--refine", type=int, default=4, help="--cuts: evaluation points per node interval along x (default 4)")
    p.add_argument("--pfield", default="p", help="--cuts: the pressure field (default p)")
    p.add_argument("--pscale", type=float, default=1.0,
                   help="--cuts: multiply the pressure by this (default 1; 0.2210 = gamma^2/(4 pi) for a run in the normalization of eq. 82)")
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
    if a.cuts:
        cuts(a, dirs, tags, labels)
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
