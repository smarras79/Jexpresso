#!/usr/bin/env python3
"""Figures of Leake, Linton & Antiochos (2010, ApJ 722:550) and of its Erratum (2011, ApJ 741:125)
from the output of problems/MHD/fluxEmergenceLeake2010.

RUN is an output directory (output/MHD/fluxEmergenceLeake2010/output-<date>, the one holding iter_N.pvtu);
a restarted run is DIR1,DIR2. Positions are plotted in Mm and times in s (code units: L0 = 150 km,
t0 = 23.1946 s), the density in g cm^-3 and B_z in G, as in the paper.

  python3 tools/plot_leake2010.py density RUN                 # Erratum Fig. 3 (Fig. 10 of the paper): log rho, field lines,
                                                              #   B_x arrows at t = 920, 1380, 1610, 1955, 2185, 2530 s
  python3 tools/plot_leake2010.py density RUN --times 920,1610 --label "Low Tension"
  python3 tools/plot_leake2010.py shear RUN                   # Erratum Figs. 1 and 8: B_z (G) and field lines
  python3 tools/plot_leake2010.py dome RUN [RUN2 ...]         # top of the dense dome on x = 0 vs time, against Erratum Fig. 3
  python3 tools/plot_leake2010.py axis RUN [RUN2 ...]         # Erratum Figs. 2(a), 7(a): height of the tube axis vs time
  python3 tools/plot_leake2010.py eshear RUN [RUN2 ...]       # Erratum Figs. 2(b,c), 5: E_shear above 1.2 and 5 Mm (Eq. 20)
  python3 tools/plot_leake2010.py profiles RUN                # Fig. 3(a,b): T, p, p_mag and beta at x = 0, t = 0

Every command takes --outdir (default figs_leake2010) and --format (png). Needs numpy and matplotlib;
it reads the VTU files with the reader of tools/plot_fluxemergence_vtu.py.
"""
import argparse, glob, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_fluxemergence_vtu import read_pvtu, pvtu_time, grid2d, flux_function_integrate  # noqa: E402

L0_MM = 0.15            # code length unit, Mm
T0_S = 23.1946          # code time unit, s (problems/MHD/fluxEmergenceLeake2010/user_flux.jl)
ERRATUM_FIG3_TIMES = [920, 1380, 1610, 1955, 2185, 2530]
ERRATUM_FIG1_TIMES = [920, 1610, 2070, 2254]

# Height of the axis of the low-tension tube (tube 2) emerging into a field-free corona, read off the
# vector data of the Erratum's Fig. 2(a) (green dashed line): t [s], height [Mm].
ERRATUM_FIG2A_TUBE2 = np.array([
    [0, -1.84], [200, -1.56], [408, -1.14], [593, -0.65], [801, -0.16], [917, -0.02], [1009, 0.13],
    [1194, 0.27], [1378, 0.41], [1494, 0.41], [1609, 0.48], [1956, 0.48], [2002, 0.41], [2187, 0.27],
    [2510, 0.13], [2533, 0.20], [2995, 0.20], [3503, 0.55], [3988, 0.83]])
# The same for the high-tension tube 3 (red dot-dashed line)
ERRATUM_FIG2A_TUBE3 = np.array([
    [0, -1.84], [200, -1.63], [408, -1.21], [593, -0.72], [801, -0.23], [917, -0.02], [1009, 0.06],
    [1194, 0.27], [1378, 0.34], [1609, 0.34], [1956, 0.27], [2187, 0.20], [2995, 0.20], [3503, 0.48],
    [3988, 0.97]])
# Initial atmosphere of the paper's Fig. 3(a), read off its vector data: y [Mm], log10 T [K], log10 p [dyn cm^-2]
PAPER_FIG3A = np.array([
    [-3.0, 4.732, 7.173], [-2.5, 4.665, 7.020], [-1.0, 4.340, 6.281], [-0.5, 4.137, 5.821], [0.0, 3.759, 4.938],
    [0.5, 3.767, 3.526], [1.0, 3.797, 2.167], [1.5, 3.895, 0.978], [2.0, 4.138, 0.159], [3.0, 5.035, -0.370],
    [5.0, 5.918, -0.458], [10.0, 5.933, -0.551]])
# Erratum Fig. 3 on x = 0: the highest point below 19 Mm where log10 rho > -13 (g cm^-3), the top of the
# dense emerging dome (inverted from the figure's colour scale, +-0.1 Mm): t [s], height [Mm]. At 2530 s the
# flux rope has moved to x = 4.5 Mm and the column is coronal above 3.65 Mm.
ERRATUM_FIG3_DOME_TOP = np.array([[920, 3.50], [1380, 7.50], [1610, 10.40], [1955, 13.05], [2185, 14.70], [2530, 3.65]])
# log10 beta at x = 0 for the medium quadrupole ("98 G", long dashes of Fig. 3(b))
PAPER_FIG3B_QUAD98 = np.array([
    [-3.0, 4.839], [-0.8, 4.394], [-0.5, 4.202], [0.0, 3.496], [0.5, 2.280], [1.0, 1.144], [1.5, 0.216],
    [2.0, -0.280], [2.5, -0.225], [3.0, 0.415], [3.5, 1.422], [4.0, 0.048], [5.0, -0.704], [6.0, -1.035],
    [8.0, -1.359], [10.0, -1.510]])


def runs_of(spec):
    """[(t [s], iter_N.pvtu)] of RUN (DIR or DIR1,DIR2,...), by time; a later directory wins a repeated time."""
    out = {}
    for d in spec.split(","):
        items = [(pvtu_time(f), f) for f in glob.glob(os.path.join(d, "iter_*.pvtu"))]
        if not items:
            raise SystemExit(f"{d}: no iter_N.pvtu (give the output-<date> directory of a run)")
        for t, f in items:
            if t is not None:
                out[round(t * T0_S, 6)] = f
    return sorted(out.items())


def nearest(items, times):
    """The snapshot nearest each requested time [s]; warns when it is more than 1 s away."""
    sel = []
    for tq in times:
        t, f = min(items, key=lambda it: abs(it[0] - tq))
        if abs(t - tq) > 1.0:
            print(f"   note: no snapshot at {tq:g} s, using t = {t:.1f} s ({os.path.basename(f)})")
        sel.append((t, f))
    return sel


def load(path, names):
    g = read_pvtu(path, set(names) | {"Bx", "By"})
    g["X"] = g["xc"] * L0_MM
    g["Y"] = g["yc"] * L0_MM
    return g


def _sem_integration_matrix(nop):
    """S[j, k] = int_{-1}^{xi_j} l_k(xi) dxi on the nop+1 LGL nodes (exact for the degree-nop interpolant)."""
    N = nop
    x = -np.cos(np.pi * np.arange(N + 1) / N)
    for _ in range(100):                                   # LGL nodes: roots of (1 - x^2) P_N'(x)
        P = np.polynomial.legendre.Legendre.basis(N)
        r = np.sort(np.real(P.deriv().roots()))
        x = np.concatenate([[-1.0], r, [1.0]])
        break
    V = np.vander(x, N + 1, increasing=True)               # l_k = sum_m C[m, k] x^m
    C = np.linalg.inv(V)
    S = np.zeros((N + 1, N + 1))
    for j in range(N + 1):
        for m in range(N + 1):
            S[j] += C[m] * (x[j] ** (m + 1) - (-1.0) ** (m + 1)) / (m + 1)
    return S


def flux_function_sem(g, Bx, By, nop=4):
    """A_z with B_x = dA/dy, B_y = -dA/dx on the tensor-product LGL lattice of LK10.geo (elements every nop
    nodes): A = -int B_y dx along the bottom row, then + int B_x dy up every column, each element's piece
    integrated exactly from its degree-nop interpolant (a trapezoidal rule would carry its error through the
    tube's axis up the whole column)."""
    bx, by = grid2d(g, Bx), grid2d(g, By)
    if bx is None or (len(g["xc"]) - 1) % nop or (len(g["yc"]) - 1) % nop:
        return None
    S = _sem_integration_matrix(nop)
    X, Y = g["xc"], g["yc"]
    ny, nx = bx.shape
    A = np.zeros((ny, nx))
    for e in range(0, nx - 1, nop):                        # bottom row, element by element
        h = 0.5 * (X[e + nop] - X[e])
        A[0, e:e + nop + 1] = A[0, e] - h * (S @ by[0, e:e + nop + 1])
    for e in range(0, ny - 1, nop):                        # up the columns
        h = 0.5 * (Y[e + nop] - Y[e])
        A[e:e + nop + 1, :] = A[e, :][None, :] + h * (S @ bx[e:e + nop + 1, :])
    return A


def a_grid(g):
    A = flux_function_sem(g, g["fields"]["Bx"], g["fields"]["By"])
    if A is None:
        A = grid2d(g, flux_function_integrate(g, g["fields"]["Bx"], g["fields"]["By"]))
    if A is None:
        raise SystemExit("the nodes do not form a tensor-product lattice (a mesh other than LK10.geo?)")
    return A * L0_MM  # flux per unit length in B0 Mm (only the spacing matters)


def field_lines(ax, g, levels, color="w", lw=0.7):
    A = a_grid(g)
    ax.contour(g["X"], g["Y"], A, levels=levels, colors=color, linewidths=lw, linestyles="solid")
    return A


def a_levels(g, xlim, ylim, n):
    """n equally spaced levels of A_z over the window of the first snapshot ("constant contour intervals")."""
    A = a_grid(g)
    mx = (g["X"] >= xlim[0]) & (g["X"] <= xlim[1])
    my = (g["Y"] >= ylim[0]) & (g["Y"] <= ylim[1])
    W = A[np.ix_(my, mx)]
    lo, hi = np.nanmin(W), np.nanmax(W)
    d = (hi - lo) / n
    return lo + d * (np.arange(n) + 0.5)


def bx_arrows(ax, g, heights, xlen=2.5):
    """White arrows at x = 0 pointing along B_x at the given heights [Mm] (the paper's arrows)."""
    ix0 = np.argmin(np.abs(g["X"]))
    Bx = grid2d(g, g["fields"]["Bx"])
    for h in heights:
        iy = np.argmin(np.abs(g["Y"] - h))
        s = np.sign(Bx[iy, ix0]) or 1.0
        ax.annotate("", xy=(s * xlen / 2, h), xytext=(-s * xlen / 2, h),
                    arrowprops=dict(arrowstyle="simple,head_length=0.7,head_width=0.7,tail_width=0.3",
                                    fc="w", ec="k", lw=0.6), zorder=4)


def _layout(fig, axs):
    """Margins in inches, so that one panel or six keep their labels and room for the colour bar."""
    w, h = fig.get_size_inches()
    fig.subplots_adjust(left=0.62 / w, right=1 - 0.12 / w, top=1 - 0.1 / h, bottom=1.05 / h,
                        wspace=0.2, hspace=0.14)


def _cbar_axes(fig, axs):
    """A horizontal colour bar under the panels, spanning them."""
    w, h = fig.get_size_inches()
    x0 = axs[-1, 0].get_position().x0
    x1 = axs[-1, -1].get_position().x1
    return fig.add_axes([x0 + 0.15 * (x1 - x0), 0.28 / h, 0.7 * (x1 - x0), 0.12 / h])


def panel_grid(n):
    ncol = 3 if n > 4 else n
    nrow = int(np.ceil(n / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(3.4 * ncol + 0.7, 3.0 * nrow + 1.2), squeeze=False)
    return fig, axs


def cmd_density(a):
    items = runs_of(a.run)
    times = [float(t) for t in a.times.split(",")] if a.times else ERRATUM_FIG3_TIMES
    sel = nearest(items, times)
    fig, axs = panel_grid(len(sel))
    levels = None
    for k, (t, f) in enumerate(sel):
        ax = axs.flat[k]
        g = load(f, ["log10rho_cgs"])
        if levels is None:
            levels = a_levels(load(items[0][1], []), a.xlim, a.ylim, a.nlines) if a.levels_from_t0 else \
                     a_levels(g, a.xlim, a.ylim, a.nlines)
        Z = grid2d(g, g["fields"]["log10rho_cgs"])
        im = ax.pcolormesh(g["X"], g["Y"], np.clip(Z, a.clim[0], a.clim[1]), cmap="jet",
                           vmin=a.clim[0], vmax=a.clim[1], shading="gouraud", rasterized=True)
        field_lines(ax, g, levels)
        if a.arrows:
            bx_arrows(ax, g, a.arrows)
        ax.set_xlim(a.xlim); ax.set_ylim(a.ylim); ax.set_aspect("equal")
        lab = f"({chr(97 + k)}) {a.label}" if (k == 0 and a.label) else f"({chr(97 + k)})"
        ax.text(0.03, 0.96, lab, color="w", transform=ax.transAxes, va="top", fontsize=9, weight="bold")
        ax.text(0.97, 0.96, f"t= {t:.0f}s", color="w", transform=ax.transAxes, va="top", ha="right",
                fontsize=9, weight="bold")
        ax.tick_params(labelsize=8)
    for ax in axs.flat[len(sel):]:
        ax.axis("off")
    for ax in axs[:, 0]:
        ax.set_ylabel("Y (Mm)", fontsize=8)
    for ax in axs[-1, :]:
        ax.set_xlabel("X (Mm)", fontsize=8)
    _layout(fig, axs)
    cb = fig.colorbar(im, cax=_cbar_axes(fig, axs), orientation="horizontal")
    cb.set_ticks(np.arange(2 * np.ceil(a.clim[0] / 2), a.clim[1] + 1e-9, 2.0))
    cb.set_label(r"log $\rho$ (g/cm$^3$)", fontsize=9)
    out = os.path.join(a.outdir, f"density.{a.format}")
    fig.savefig(out, dpi=a.dpi); plt.close(fig)
    print("wrote", out)


def cmd_shear(a):
    items = runs_of(a.run)
    times = [float(t) for t in a.times.split(",")] if a.times else ERRATUM_FIG1_TIMES
    sel = nearest(items, times)
    fig, axs = panel_grid(len(sel))
    levels = None
    for k, (t, f) in enumerate(sel):
        ax = axs.flat[k]
        g = load(f, ["Bz_G"])
        if levels is None:
            levels = a_levels(g, a.xlim, a.ylim, a.nlines)
        Z = grid2d(g, g["fields"]["Bz_G"])
        im = ax.pcolormesh(g["X"], g["Y"], np.clip(Z, a.bclim[0], a.bclim[1]), cmap="jet",
                           vmin=a.bclim[0], vmax=a.bclim[1], shading="gouraud", rasterized=True)
        field_lines(ax, g, levels)
        ax.set_xlim(a.xlim); ax.set_ylim(a.ylim); ax.set_aspect("equal")
        ax.text(0.03, 0.04, f"({chr(97 + k)}) {a.label if k == 0 else ''}", color="w", transform=ax.transAxes, fontsize=8)
        ax.text(0.97, 0.04, f"{t:.0f}s", color="w", transform=ax.transAxes, ha="right", fontsize=8)
        ax.tick_params(labelsize=8)
    for ax in axs.flat[len(sel):]:
        ax.axis("off")
    _layout(fig, axs)
    cb = fig.colorbar(im, cax=_cbar_axes(fig, axs), orientation="horizontal")
    cb.set_label(r"$B_z$ (G)", fontsize=9)
    out = os.path.join(a.outdir, f"shear.{a.format}")
    fig.savefig(out, dpi=a.dpi); plt.close(fig)
    print("wrote", out)


def axis_height(g):
    """Height [Mm] of the tube axis: where B_x on x = 0 changes sign from - to + going up (Sec. 3.2: "defined
    at the point where B_x|x=0 changes sign"). Of all such crossings below 10 Mm, the one nearest the peak of
    the tube's axial field B_z on x = 0 (far from the tube B_x is ~0 and noise crosses zero too); nan if none."""
    ix0 = np.argmin(np.abs(g["X"]))
    Bx = grid2d(g, g["fields"]["Bx"])[:, ix0]
    Bz = np.abs(grid2d(g, g["fields"]["Bz"])[:, ix0])
    Y = g["Y"]
    yz = Y[np.argmax(np.where(Y < 10.0, Bz, -1.0))]
    best = np.nan
    for k in range(len(Y) - 1):
        if Y[k] > 10.0:
            break
        if Bx[k] < 0.0 <= Bx[k + 1]:
            yc = Y[k] - Bx[k] * (Y[k + 1] - Y[k]) / (Bx[k + 1] - Bx[k])
            if np.isnan(best) or abs(yc - yz) < abs(best - yz):
                best = yc
    return best


def cmd_axis(a):
    fig, ax = plt.subplots(figsize=(4.5, 4.2))
    for spec in a.runs:
        lab, _, path = spec.rpartition("=")
        items = runs_of(path)
        T, H = [], []
        for t, f in items:
            T.append(t); H.append(axis_height(load(f, ["Bz"])))
        np.savetxt(os.path.join(a.outdir, f"axis_{(lab or os.path.basename(path.rstrip('/'))).replace(' ', '_')}.csv"),
                   np.c_[T, H], header="t [s], axis height [Mm]", delimiter=",")
        H = np.array(H)
        if a.shift:   # start from the erratum's initial height (the discrete axis of tube 2 sits higher on coarse meshes)
            H = H - H[0] + ERRATUM_FIG2A_TUBE2[0, 1]
        ax.plot(T, H, "-o", ms=3, label=(lab or os.path.basename(path.rstrip("/"))) + (" (shifted)" if a.shift else ""))
    ax.plot(*ERRATUM_FIG2A_TUBE2.T, "g--", lw=1.2, label="Erratum Fig. 2(a), low tension")
    ax.plot(*ERRATUM_FIG2A_TUBE3.T, "r-.", lw=1.0, label="Erratum Fig. 2(a), high tension")
    ax.set_xlabel("Time (s)"); ax.set_ylabel("Height of axis (Mm)")
    ax.set_xlim(0, a.tmax); ax.set_ylim(-3, 12 if a.tall else 2)
    ax.legend(fontsize=7); ax.grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(a.outdir, f"axis.{a.format}")
    fig.savefig(out, dpi=a.dpi); plt.close(fig)
    print("wrote", out)


def dome_top(g, thr=-13.0, ymax=19.0):
    """Highest point on x = 0 below ymax [Mm] where log10 rho > thr (g cm^-3), above 2 Mm; nan if none."""
    ix0 = np.argmin(np.abs(g["X"]))
    L = grid2d(g, g["fields"]["log10rho_cgs"])[:, ix0]
    Y = g["Y"]
    k = np.where((Y > 2.0) & (Y < ymax) & (L > thr))[0]
    return Y[k.max()] if len(k) else np.nan


def cmd_dome(a):
    fig, ax = plt.subplots(figsize=(4.5, 3.8))
    for spec in a.runs:
        lab, _, path = spec.rpartition("=")
        items = runs_of(path)
        T, H = [], []
        for t, f in items:
            T.append(t); H.append(dome_top(load(f, ["log10rho_cgs"]), a.thr))
        name = (lab or os.path.basename(path.rstrip("/"))).replace(" ", "_")
        np.savetxt(os.path.join(a.outdir, f"dome_{name}.csv"), np.c_[T, H],
                   header=f"t [s], top of log10 rho > {a.thr} on x = 0 [Mm]", delimiter=",")
        ax.plot(T, H, "-o", ms=3, label=lab or os.path.basename(path.rstrip("/")))
    ax.plot(*ERRATUM_FIG3_DOME_TOP.T, "ks", mfc="none", ms=7, label="Erratum Fig. 3")
    ax.set_xlabel("Time (s)"); ax.set_ylabel(r"top of log $\rho$ > %g on x = 0 (Mm)" % a.thr)
    ax.set_xlim(0, 2700); ax.set_ylim(0, 19); ax.grid(alpha=0.3); ax.legend(fontsize=7)
    fig.tight_layout()
    out = os.path.join(a.outdir, f"dome.{a.format}")
    fig.savefig(out, dpi=a.dpi); plt.close(fig)
    print("wrote", out)


def eshear(g, y1):
    """E_shear (Eq. 20): int B_z^2 / int |B|^2 over -45 < x < 45 Mm, y1 < y < 90 Mm (trapezoids on the lattice)."""
    Bx, By = grid2d(g, g["fields"]["Bx"]), grid2d(g, g["fields"]["By"])
    Bz = grid2d(g, g["fields"]["Bz"])
    m = g["Y"] >= y1
    Y = g["Y"][m]
    num = np.trapezoid(np.trapezoid(Bz[m] ** 2, g["X"], axis=1), Y)
    den = np.trapezoid(np.trapezoid((Bx ** 2 + By ** 2 + Bz ** 2)[m], g["X"], axis=1), Y)
    return num / den if den > 0 else 0.0


def cmd_eshear(a):
    fig, axs = plt.subplots(2, 1, figsize=(4.5, 5.5), sharex=True)
    for spec in a.runs:
        lab, _, path = spec.rpartition("=")
        items = runs_of(path)
        T = [t for t, _ in items]
        E = np.array([[eshear(load(f, ["Bz"]), y1) for y1 in (1.2, 5.0)] for _, f in items])
        for k, ax in enumerate(axs):
            ax.plot(T, E[:, k], "-o", ms=3, label=lab or os.path.basename(path.rstrip("/")))
        np.savetxt(os.path.join(a.outdir, f"eshear_{(lab or os.path.basename(path.rstrip('/'))).replace(' ', '_')}.csv"),
                   np.c_[T, E], header="t [s], E_shear above 1.2 Mm, above 5 Mm", delimiter=",")
    for ax, y1 in zip(axs, (1.2, 5.0)):
        ax.set_ylabel(r"$E_{shear}$ above %g Mm" % y1); ax.grid(alpha=0.3); ax.set_ylim(0, 0.1)
    axs[0].legend(fontsize=7); axs[1].set_xlabel("Time (s)")
    fig.tight_layout()
    out = os.path.join(a.outdir, f"eshear.{a.format}")
    fig.savefig(out, dpi=a.dpi); plt.close(fig)
    print("wrote", out)


def cmd_profiles(a):
    items = runs_of(a.run)
    t, f = items[0]
    g = load(f, ["p", "T", "β", "Bz"])
    ix0 = np.argmin(np.abs(g["X"]))
    Y = g["Y"]
    col = lambda n: grid2d(g, g["fields"][n])[:, ix0]
    # code -> cgs: p0 = 10^4.927 dyn cm^-2, T0 = 5700 K, B0 = (4 pi p0)^1/2
    p0 = 10 ** 4.927
    p = col("p") * p0
    T = col("T") * 5700.0
    pm = 0.5 * (col("Bx") ** 2 + col("By") ** 2 + col("Bz") ** 2) * p0
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8, 4.5), sharey=True)
    a1.semilogx(T, Y, "k-", lw=0.8, label="T (K)")
    a1b = a1.twiny()
    a1b.semilogx(p, Y, "k--", label="gas pressure")
    a1b.semilogx(np.maximum(pm, 1e-30), Y, "k-.", label="magnetic pressure")
    a1.semilogx(10 ** PAPER_FIG3A[:, 1], PAPER_FIG3A[:, 0], "o", mfc="none", color="tab:red", ms=5, label="paper T")
    a1b.semilogx(10 ** PAPER_FIG3A[:, 2], PAPER_FIG3A[:, 0], "s", mfc="none", color="tab:blue", ms=5, label="paper p")
    a1.set_xlim(1e3, 1e6); a1b.set_xlim(1e-2, 1e8); a1.set_ylim(-3, 20)
    a1.set_xlabel("Temperature (K)"); a1b.set_xlabel("Pressure (dyne/cm$^2$)"); a1.set_ylabel("Height (Mm)")
    a1.legend(loc="center right", fontsize=7); a1b.legend(loc="upper right", fontsize=7)
    a2.semilogx(col("β"), Y, "k-", label="this run")
    a2.semilogx(10 ** PAPER_FIG3B_QUAD98[:, 1], PAPER_FIG3B_QUAD98[:, 0], "o", mfc="none", color="tab:red", ms=5,
                label="paper, 98 G quadrupole")
    a2.axvline(0.2, color="k", ls=":", lw=0.8)
    a2.set_xlim(1e-3, 1e6); a2.set_xlabel(r"Plasma $\beta$"); a2.legend(fontsize=7)
    fig.tight_layout()
    out = os.path.join(a.outdir, f"profiles.{a.format}")
    fig.savefig(out, dpi=a.dpi); plt.close(fig)
    print("wrote", out)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default="figs_leake2010")
    p.add_argument("--format", default="png")
    p.add_argument("--dpi", type=int, default=150)
    sub = p.add_subparsers(dest="cmd", required=True)
    for name, fn in (("density", cmd_density), ("shear", cmd_shear)):
        s = sub.add_parser(name); s.set_defaults(fn=fn)
        s.add_argument("run")
        s.add_argument("--times", default="", help="comma-separated times in s")
        s.add_argument("--label", default="Low Tension")
        s.add_argument("--xlim", type=lambda t: tuple(map(float, t.split(","))),
                       default=(-12.75, 12.75) if name == "density" else (-5.5, 5.5))
        s.add_argument("--ylim", type=lambda t: tuple(map(float, t.split(","))),
                       default=(-3.0, 19.5) if name == "density" else (-2.2, 8.2))
        s.add_argument("--nlines", type=int, default=45)
        s.add_argument("--levels-from-t0", action="store_true", help="field-line levels from the first snapshot")
        s.add_argument("--clim", type=lambda t: tuple(map(float, t.split(","))), default=(-16.0, -5.3))
        s.add_argument("--bclim", type=lambda t: tuple(map(float, t.split(","))), default=(0.0, 500.0))
        s.add_argument("--arrows", type=lambda t: [float(v) for v in t.split(",")] if t else [],
                       default=[3.0, 11.0, 16.0], help="heights [Mm] of the B_x arrows on x = 0 ('' for none)")
    s = sub.add_parser("axis"); s.set_defaults(fn=cmd_axis)
    s.add_argument("runs", nargs="+", help="RUN or LABEL=RUN")
    s.add_argument("--tmax", type=float, default=4600.0)
    s.add_argument("--tall", action="store_true", help="y axis up to 12 Mm, as the paper's panel")
    s.add_argument("--shift", action="store_true", help="shift each run to start at the erratum's -1.84 Mm")
    s = sub.add_parser("dome"); s.set_defaults(fn=cmd_dome)
    s.add_argument("runs", nargs="+", help="RUN or LABEL=RUN")
    s.add_argument("--thr", type=float, default=-13.0, help="log10 rho threshold [g cm^-3]")
    s = sub.add_parser("eshear"); s.set_defaults(fn=cmd_eshear)
    s.add_argument("runs", nargs="+", help="RUN or LABEL=RUN")
    s = sub.add_parser("profiles"); s.set_defaults(fn=cmd_profiles)
    s.add_argument("run")
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    a.fn(a)


if __name__ == "__main__":
    main()
