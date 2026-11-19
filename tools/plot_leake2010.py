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
  python3 tools/plot_leake2010.py mass RUN                    # the plasma lifted above 3 Mm: its mass budget (with the flow
                                                              #   or not), temperature, pressure, draining speed, heating time
  python3 tools/plot_leake2010.py null RUN [RUN2 ...]         # on x = 0 vs time: the top of the emerged flux (the null
                                                              #   above it), the flux between the tube axis and it, the
                                                              #   mass lifted, and the DynSGS coefficients at the null

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
V0_KMS = L0_MM * 1e3 / T0_S   # code velocity unit, km/s
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


def _record(d):
    """lk10_case.txt of an output directory (initialize.jl): (the directory it was restarted from, by name,
    the time [s] of the snapshot it continued), or (None, None)."""
    f = os.path.join(d, "lk10_case.txt")
    if not os.path.isfile(f):
        return None, None
    for line in open(f, encoding="utf-8", errors="replace"):
        m = re.match(r"#\s*restarted from\s+(.+?),\s*iter_\d+(?:,\s*t\s*=\s*([-+0-9.eE]+)\s*s)?", line.strip())
        if m:
            return os.path.basename(os.path.normpath(m.group(1))), (float(m.group(2)) if m.group(2) else None)
    return None, None


def _chain(d):
    """[(directory, last time used [s])] of the run that ENDS in d, oldest first: d and the runs it was
    restarted from (JEXPRESSO_LK10_RESTART_DIR; each restart names the run and the snapshot it continues in
    its lk10_case.txt). An earlier run is used up to the snapshot its successor continued, so a run
    restarted from the middle of a finished one (to try other numerics) is not mixed with what that one
    did afterwards. Runs restarted FROM d are not followed: give the last directory of a chain."""
    root, name = os.path.split(os.path.normpath(d))
    root = root or "."
    if not os.path.isdir(os.path.join(root, name)):
        raise SystemExit(f"{d}: not a directory (give the output-<date> directory of a run)")
    rec = {s: _record(os.path.join(root, s)) for s in os.listdir(root) if os.path.isdir(os.path.join(root, s))}
    chain = [(name, np.inf)]
    while True:
        par, tr = rec.get(chain[0][0], (None, None))
        if par not in rec or par in [c[0] for c in chain]:
            break
        chain.insert(0, (par, np.inf if tr is None else tr))
    kids = sorted(s for s, (q, _) in rec.items() if q == name)
    if kids:
        print(f"   note: {', '.join(kids)} restarted from {name} and {'is' if len(kids) == 1 else 'are'} not included "
              "(give the last directory of a chain to plot the whole run)")
    return [(os.path.join(root, s), t) for s, t in chain]


FOLLOW_CHAIN = True    # --no-chain: only the directories given


def runs_of(spec):
    """[(t [s], iter_N.pvtu)] of RUN (DIR or DIR1,DIR2,...), by time. Unless --no-chain, each DIR stands for
    the run that ends in it: the runs it was restarted from are read too, each up to the snapshot its
    successor continued (_chain). A later directory wins a repeated time."""
    dirs = []
    for d in spec.split(","):
        for c, tmax in (_chain(d) if FOLLOW_CHAIN else [(os.path.normpath(d), np.inf)]):
            if c not in [x[0] for x in dirs]:
                dirs.append((c, tmax))
    out = {}
    for d, tmax in dirs:
        items = [(pvtu_time(f), f) for f in glob.glob(os.path.join(d, "iter_*.pvtu"))]
        items = [(t * T0_S, f) for t, f in items if t is not None]
        if not items:
            raise SystemExit(f"{d}: no iter_N.pvtu with a time (give the output-<date> directory of a run)")
        used = [it for it in items if it[0] <= tmax + 0.5]
        note = "" if len(used) == len(items) else f"; used to {tmax:.0f} s, where the next run continued it"
        print(f"   {d}: {len(items)} snapshots, t = {min(items)[0]:.0f} to {max(items)[0]:.0f} s{note}")
        for t, f in used:
            out[round(t, 6)] = f
    return sorted(out.items())


def nearest(items, times, tol=1.0):
    """The snapshot at each requested time [s], within tol. A time the run has no snapshot at is left out
    (the deck writes every figure time, so the run did not get there); with --nearest the nearest is used."""
    sel, miss = [], []
    for tq in times:
        t, f = min(items, key=lambda it: abs(it[0] - tq))
        if abs(t - tq) <= tol or NEAREST:
            if abs(t - tq) > tol:
                print(f"   note: no snapshot at {tq:g} s, using t = {t:.1f} s ({os.path.basename(f)})")
            sel.append((t, f))
        else:
            miss.append(tq)
    if miss:
        print(f"   no snapshot at {', '.join(f'{m:g}' for m in miss)} s: the snapshots end at {items[-1][0]:.1f} s. "
              "If the run stopped there, its log says why; a run continued by restarts is followed through "
              "lk10_case.txt (or give all its directories: RUN=dir1,dir2).")
    if not sel:
        raise SystemExit("   nothing to plot")
    return sel


NEAREST = False        # --nearest: the nearest snapshot for a time the run has none at


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


MU_NAME = "mu_dsgs_ρu_ρv_ρw"   # the DynSGS momentum coefficient in the output (:plot_dsgs)
ETA_NAMES = ("mu_dsgs_Bx_By_Bz", "mu_dsgs_ρ_ρE_Bx_By_Bz_ψ")   # the B coefficient: its own (JEXPRESSO_LK10_ETA=induction), or shared


def cmd_diag(a):
    """What grows before a run fails: per snapshot, the largest speed |u|, fast speed c_f = ((gamma p + B^2)/rho)^1/2,
    CFL rate (|u| + c_f)/dx (dx the local node spacing; times the run's dt, the CFL number), current |curl B|,
    and the smallest p/p(t0), beta and T = p/rho (T_ph units), with where they are (Mm), the number of
    nodes at the temperature floor (--tfloor, 0.1 by default), the smallest rho/rho(t0), and the largest DynSGS
    momentum coefficient (code units), where the operator works hardest.
    Code units, as the run (README, "Units")."""
    items = runs_of(a.run)
    names = ["ρ", "p", "u", "v", "w", "Bx", "By", "Bz", MU_NAME]
    g0 = load(items[0][1], names)          # the reference: the run's first snapshot, whatever --tmin
    if a.tmin is not None:
        items = [it for it in items if it[0] >= a.tmin - 1e-6]
    G = lambda g, n: grid2d(g, g["fields"][n])
    p0 = G(g0, "p")
    r0 = G(g0, "ρ")
    if p0 is None:
        raise SystemExit("diag: the nodes do not form a tensor-product grid")
    xc, yc = g0["xc"], g0["yc"]
    dx = np.gradient(xc)[None, :]; dy = np.gradient(yc)[:, None]
    dmin = np.minimum(np.broadcast_to(dx, p0.shape), np.broadcast_to(dy, p0.shape))
    win = np.ones(p0.shape, bool)
    if a.ymin is not None:
        win &= (yc[:, None] * L0_MM >= a.ymin)
    at = lambda Z, k: f"({xc[k[1]] * L0_MM:6.2f},{yc[k[0]] * L0_MM:6.2f})"
    def arg(Z, fn):
        Zw = np.where(win, Z, np.nan)
        k = np.unravel_index(fn(Zw), Z.shape)
        return Zw[k], k
    cfl_hdr = "CFL" if a.dt else "(|u|+cf)/dx"
    print(f"{'t [s]':>8} {'max|u|':>9} {'at (x,y) Mm':>16} {'max cf':>9} {cfl_hdr:>11} {'at':>16} "
          f"{'max|j|':>9} {'at':>16} {'min p/p0':>9} {'at':>16} {'min beta':>9} {'at':>16} {'min T':>8} {'at':>16} {'#T floor':>8} {'min rho/rho0':>12} {'at':>16} {'max mu':>9} {'at':>16}")
    for t, f in items:
        g = load(f, names)
        rho, p = G(g, "ρ"), G(g, "p")
        u, v, w = G(g, "u"), G(g, "v"), G(g, "w")
        Bx, By, Bz = G(g, "Bx"), G(g, "By"), G(g, "Bz")
        B2 = Bx**2 + By**2 + Bz**2
        sp = np.sqrt(u**2 + v**2 + w**2)
        cf = np.sqrt(np.maximum(5.0 / 3.0 * p + B2, 0.0) / rho)
        rate = (sp + cf) / dmin
        jz = np.gradient(By, xc, axis=1) - np.gradient(Bx, yc, axis=0)
        jx = np.gradient(Bz, yc, axis=0); jy = -np.gradient(Bz, xc, axis=1)
        j = np.sqrt(jx**2 + jy**2 + jz**2)
        beta = np.where(B2 > 0, 2.0 * p / np.maximum(B2, 1e-300), np.inf)
        su, ku = arg(sp, np.nanargmax); sc, _ = arg(cf, np.nanargmax); sr, kr = arg(rate, np.nanargmax)
        sj, kj = arg(j, np.nanargmax); spp, kp = arg(p / p0, np.nanargmin); sb, kb = arg(beta, np.nanargmin)
        sT, kT = arg(p / rho, np.nanargmin)
        nT = int(np.sum(np.where(win, p / rho, np.inf) <= 1.001 * a.tfloor))
        sR, kR = arg(rho / r0, np.nanargmin)
        mu = G(g, MU_NAME) if MU_NAME in g["fields"] else np.full(p.shape, np.nan)
        sM, kM = arg(mu, np.nanargmax) if np.isfinite(mu).any() else (np.nan, (0, 0))
        print(f"{t:8.1f} {su * V0_KMS:7.1f}km/s {at(sp, ku):>16} {sc:9.1f} {sr * (a.dt or 1.0):11.3g} {at(rate, kr):>16} "
              f"{sj:9.3g} {at(j, kj):>16} {spp:9.3g} {at(p, kp):>16} {sb:9.3g} {at(beta, kb):>16} {sT:8.3g} {at(p, kT):>16} {nT:8d} {sR:12.3g} {at(p, kR):>16} {sM:9.3g} {at(p, kM):>16}")


def cmd_atmos(a):
    """The atmosphere away from the emerging region against height, at several times: medians over
    XMIN <= |x| <= XMAX Mm of log10 rho (g cm^-3), log10 T (K), and the DynSGS momentum coefficient (code
    units; dashed: the coefficient of the B slots, the DynSGS resistivity). Shows whether the chromosphere
    heats and the transition region rises (heating), or the corona fills with dense plasma at the
    chromospheric temperature (transport), and where the DynSGS acts."""
    items = runs_of(a.run)
    times = [float(t) for t in a.times.split(",")] if a.times else [0, 464, 920, 1160, 1380, 1610, 1955, 2530]
    sel = nearest(items, times)
    fig, axs = plt.subplots(1, 3, figsize=(11.5, 4.6), sharey=True)
    cmap = plt.get_cmap("viridis")
    rows = []
    for i, (t, f) in enumerate(sel):
        g = load(f, ["log10rho_cgs", "T", MU_NAME, *ETA_NAMES])
        cols = (np.abs(g["X"]) >= a.xmin) & (np.abs(g["X"]) <= a.xmax)
        Y = g["Y"]
        R = np.median(grid2d(g, g["fields"]["log10rho_cgs"])[:, cols], axis=1)
        T = np.median(grid2d(g, g["fields"]["T"])[:, cols], axis=1) * 5700.0
        M = np.median(grid2d(g, g["fields"][MU_NAME])[:, cols], axis=1) if MU_NAME in g["fields"] else np.full(Y.shape, np.nan)
        en = next((n for n in ETA_NAMES if n in g["fields"]), None)
        E = np.median(grid2d(g, g["fields"][en])[:, cols], axis=1) if en else np.full(Y.shape, np.nan)
        c = cmap(i / max(len(sel) - 1, 1))
        axs[0].plot(R, Y, color=c, lw=1.2, label=f"{t:.0f} s")
        axs[1].plot(np.log10(np.clip(T, 1.0, None)), Y, color=c, lw=1.2)
        axs[2].plot(M, Y, color=c, lw=1.2)
        axs[2].plot(E, Y, color=c, lw=0.9, ls="--")
        # where log rho crosses -13 (the transition region), and the chromospheric temperature at 2 Mm
        k = np.where((R[:-1] > -13) & (R[1:] <= -13) & (Y[:-1] > 0))[0]
        ytr = Y[k[0]] + (R[k[0]] + 13) / (R[k[0]] - R[k[0] + 1]) * (Y[k[0] + 1] - Y[k[0]]) if len(k) else np.nan
        band = (Y > 0) & (Y < 10)
        rows.append((t, ytr, np.interp(1.5, Y, T), np.interp(6.0, Y, R), np.nanmax(np.where(band, M, np.nan)),
                     np.nanmax(np.where(band, E, np.nan)) if en else np.nan))
    for ax in axs:
        ax.set_ylim(a.ylim); ax.grid(alpha=0.3); ax.tick_params(labelsize=8)
    axs[0].set_xlabel(r"log $\rho$ (g cm$^{-3}$)"); axs[1].set_xlabel(r"log $T$ (K)"); axs[2].set_xlabel(r"DynSGS $\mu$ (solid), resistivity $\eta$ (dashed), code units")
    axs[0].set_ylabel("Y (Mm)"); axs[0].legend(fontsize=7, loc="upper right")
    fig.suptitle(f"median over {a.xmin:g} <= |x| <= {a.xmax:g} Mm", fontsize=9)
    fig.tight_layout()
    out = os.path.join(a.outdir, f"atmos.{a.format}")
    fig.savefig(out, dpi=a.dpi); plt.close(fig)
    print(f"{'t [s]':>7} {'TR (log rho=-13) [Mm]':>22} {'T at 1.5 Mm / T(t0)':>20} {'log rho at 6 Mm':>16} {'max mu (0-10 Mm)':>17} {'max eta (0-10 Mm)':>18}")
    for t, ytr, t2, r6, mm, ee in rows:
        print(f"{t:7.0f} {ytr:22.2f} {t2 / rows[0][2]:20.2f} {r6:16.2f} {mm:17.3g} {ee:18.3g}")
    print("wrote", out)


def _bx_crossings(Y, B, sign):
    """Heights where B (on x = 0) changes sign, going up: from - to + (sign = +1) or + to - (sign = -1)."""
    out = []
    for i in range(len(Y) - 1):
        if (sign > 0 and B[i] < 0 <= B[i + 1]) or (sign < 0 and B[i] > 0 >= B[i + 1]):
            out.append(Y[i] - B[i] * (Y[i + 1] - Y[i]) / (B[i + 1] - B[i]))
    return out


def null_row(g, thr=-13.0):
    """On x = 0: the tube axis (B_x from - to + below 3 Mm), the top of the emerged flux system (the first height
    above the axis where B_x turns from + to -: the quadrupole's null, or the O-point of a rope formed in its
    current sheet), the in-plane flux between the two (A(top) - A(axis), per unit length, in 10^10 Mx/cm), the
    overlying flux from the top to y = 40 Mm, the dome top (log rho > thr), the mass in 3 < y < 20 Mm, |x| < 13 Mm
    (g/cm), and the largest DynSGS mu and eta within 1 Mm of the top (code units)."""
    ix0 = np.argmin(np.abs(g["X"]))
    Y = g["Y"]
    Bx = grid2d(g, g["fields"]["Bx"])[:, ix0]
    A = a_grid(g)[:, ix0] / L0_MM                       # code units B0 L0
    axis = axis_height(g)
    # above the axis and above the low chromosphere (small reversals there are not the null)
    tops = [y for y in _bx_crossings(Y, Bx, -1) if y > max(axis if np.isfinite(axis) else 0.0, 1.0) and y < 40.0]
    top = tops[0] if tops else np.nan
    Ai = lambda y: np.interp(y, Y, A)
    # B0 = (4 pi p0)^1/2 G with p0 = 10^4.927 dyn cm^-2, L0 = 1.5e7 cm
    unit = np.sqrt(4.0 * np.pi * 10 ** 4.927) * 1.5e7 / 1e10
    flux = (Ai(top) - Ai(axis)) * unit if np.isfinite(top) and np.isfinite(axis) else np.nan
    over = (Ai(top) - Ai(40.0)) * unit if np.isfinite(top) else np.nan
    L = grid2d(g, g["fields"]["log10rho_cgs"])
    m = ((Y > 3.0) & (Y < 20.0))[:, None] & (np.abs(g["X"]) < 13.0)[None, :]
    rho = np.where(m, 10.0 ** L, 0.0)
    mass = np.trapezoid(np.trapezoid(rho, g["X"] * 1e8, axis=1), Y * 1e8)
    near = np.abs(Y - top) < 1.0 if np.isfinite(top) else np.zeros_like(Y, bool)
    def mx(name):
        if name not in g["fields"] or not near.any():
            return np.nan
        return float(np.max(grid2d(g, g["fields"][name])[near, ix0]))
    en = next((n for n in ETA_NAMES if n in g["fields"]), None)
    return dict(axis=axis, top=top, flux=flux, over=over, dome=dome_top(g, thr), mass=mass,
                mu=mx(MU_NAME), eta=mx(en) if en else np.nan, ntop=len(tops))


def cmd_null(a):
    """The emerged flux system against time on x = 0 (null_row). If the flux between the axis and the top
    falls while the top stalls, the emerging flux is being reconnected away at the null as fast as it rises;
    if it holds, the top is held down (weight, a stronger overlying field)."""
    fig, axs = plt.subplots(1, 3, figsize=(12.5, 3.8))
    for spec in a.runs:
        lab, _, path = spec.rpartition("=")
        lab = lab or os.path.basename(path.rstrip("/"))
        rows = []
        for t, f in runs_of(path):
            g = load(f, ["log10rho_cgs", "Bz", MU_NAME, *ETA_NAMES])
            rows.append((t, null_row(g, a.thr)))
        T = np.array([t for t, _ in rows])
        col = lambda k: np.array([r[k] for _, r in rows])
        name = lab.replace(" ", "_")
        np.savetxt(os.path.join(a.outdir, f"null_{name}.csv"),
                   np.c_[T, col("axis"), col("top"), col("dome"), col("flux"), col("over"), col("mass"), col("mu"), col("eta")],
                   header="t [s], axis [Mm], top of the emerged flux (B_x + to - on x = 0) [Mm], dome top [Mm], "
                          "flux axis-top [1e10 Mx/cm], flux top-40 Mm [1e10 Mx/cm], mass in 3 < y < 20 Mm, |x| < 13 Mm "
                          "[g/cm], max mu near the top, max eta near the top", delimiter=",")
        print(f"--- {lab}")
        print("heights in Mm on x = 0; in-plane flux in 1e10 Mx/cm; mass in 3 < y < 20 Mm, |x| < 13 Mm, against the first snapshot")
        print(f"{'t [s]':>7} {'axis':>6} {'top':>6} {'dome':>6} {'flux axis-top':>13} {'overlying':>9} {'mass':>8} {'mu@top':>8} {'eta@top':>8}")
        m0 = rows[0][1]["mass"]
        for t, r in rows:
            print(f"{t:7.0f} {r['axis']:6.2f} {r['top']:6.2f} {r['dome']:6.2f} {r['flux']:13.4g} {r['over']:9.4g} "
                  f"{r['mass'] / m0:8.2f} {r['mu']:8.3g} {r['eta']:8.3g}" + ("  (+ a second sign change: rope)" if r["ntop"] > 1 else ""))
        axs[0].plot(T, col("top"), "-o", ms=3, label=f"{lab}: top of the emerged flux")
        axs[0].plot(T, col("dome"), "--", lw=1, label=f"{lab}: dome top")
        axs[1].plot(T, col("flux"), "-o", ms=3, label=f"{lab}: axis to top")
        axs[1].plot(T, col("over"), "--", lw=1, label=f"{lab}: top to 40 Mm")
        axs[2].plot(T, col("mass") / col("mass")[0], "-o", ms=3, label=lab)
    axs[0].plot(*ERRATUM_FIG3_DOME_TOP.T, "ks", mfc="none", ms=7, label="Erratum Fig. 3, dome top")
    axs[0].set_ylabel("height on x = 0 (Mm)"); axs[0].set_ylim(0, 19)
    axs[1].set_ylabel(r"in-plane flux (10$^{10}$ Mx cm$^{-1}$)")
    axs[2].set_ylabel(r"mass in 3 < y < 20 Mm, |x| < 13 Mm, / first snapshot"); axs[2].set_yscale("log")
    for ax in axs:
        ax.set_xlabel("Time (s)"); ax.set_xlim(0, 2700); ax.grid(alpha=0.3); ax.legend(fontsize=6)
    fig.tight_layout()
    out = os.path.join(a.outdir, f"null.{a.format}")
    fig.savefig(out, dpi=a.dpi); plt.close(fig)
    print("wrote", out)


def _edge_flux(F, X, Y, box):
    """Inflow [per unit length] of the flux vector F = (Fx, Fy) (grids (ny, nx)) through the edges of the box
    (x0, x1, y0, y1): through y = y0 and x = x0 counted inward, y = y1 and x = x1 outward; nearest lattice lines."""
    x0, x1, y0, y1 = box
    iy0, iy1 = np.argmin(np.abs(Y - y0)), np.argmin(np.abs(Y - y1))
    ix0, ix1 = np.argmin(np.abs(X - x0)), np.argmin(np.abs(X - x1))
    cx, cy = slice(ix0, ix1 + 1), slice(iy0, iy1 + 1)
    return (np.trapezoid(F[1][iy0, cx], X[cx]) - np.trapezoid(F[1][iy1, cx], X[cx])
            + np.trapezoid(F[0][cy, ix0], Y[cy]) - np.trapezoid(F[0][cy, ix1], Y[cy]))


def mass_row(g, p0):
    """The plasma in the box 3 < y < 20 Mm, |x| < 13 Mm (cmd_mass). p0: the first snapshot's pressure grid."""
    X, Y = g["X"], g["Y"]
    G = lambda n: grid2d(g, g["fields"][n])
    box = (-13.0, 13.0, 3.0, 20.0)
    m = ((Y > box[2]) & (Y < box[3]))[:, None] & ((X > box[0]) & (X < box[1]))[None, :]
    L = G("log10rho_cgs"); rho = 10.0 ** L                        # g cm^-3
    u, v = G("u"), G("v")
    T = G("T") * 5700.0
    cm = 1e8                                                      # cm per Mm
    integ = lambda f: np.trapezoid(np.trapezoid(np.where(m, f, 0.0), X * cm, axis=1), Y * cm)
    M = integ(rho)
    vs = V0_KMS * 1e5                                             # cm/s per code velocity
    Fin = _edge_flux((rho * u * vs, rho * v * vs), X * cm, Y * cm, (box[0] * cm, box[1] * cm, box[2] * cm, box[3] * cm))
    w = np.where(m, rho, 0.0)
    fr = lambda sel: integ(np.where(sel, rho, 0.0)) / M
    P = G("p"); ratio = P / np.where(p0 > 0, p0, np.nan)
    o = np.argsort(ratio[m]); cw = np.cumsum(w[m][o]); pmed = ratio[m][o][np.searchsorted(cw, 0.5 * cw[-1])]
    vy = integ(rho * v) / M * V0_KMS
    # heating of the dense plasma (log rho > -13) by the DynSGS, code units: rho mu |grad v|^2 and eta |curl B|^2
    d = m & (L > -13.0)
    Xc, Yc = g["xc"], g["yc"]
    dx = lambda f: np.gradient(f, Xc, axis=1)
    dy = lambda f: np.gradient(f, Yc, axis=0)
    w3 = G("w")
    gv2 = dx(u) ** 2 + dy(u) ** 2 + dx(v) ** 2 + dy(v) ** 2 + dx(w3) ** 2 + dy(w3) ** 2
    Bx, By, Bz = G("Bx"), G("By"), G("Bz")
    j2 = (dx(By) - dy(Bx)) ** 2 + dy(Bz) ** 2 + dx(Bz) ** 2
    mu = G(MU_NAME) if MU_NAME in g["fields"] else np.zeros_like(u)
    en = next((n for n in ETA_NAMES if n in g["fields"]), None)
    eta = G(en) if en else np.zeros_like(u)
    eint = np.sum(np.where(d, P, 0.0)) / (5.0 / 3.0 - 1.0)
    qv = np.sum(np.where(d, G("ρ") * mu * gv2, 0.0))
    qo = np.sum(np.where(d, eta * j2, 0.0))
    tv = eint / qv * T0_S if qv > 0 else np.inf
    to = eint / qo * T0_S if qo > 0 else np.inf
    return dict(M=M, Fin=Fin, cold=fr(T < 2e4), warm=fr((T >= 2e4) & (T < 2e5)), hot=fr(T >= 2e5),
                pmed=pmed, vy=vy, tv=tv, to=to)


def cmd_mass(a):
    """Where the plasma lifted above the transition region goes, and how hot it is. Per snapshot, over the box
    3 < y < 20 Mm, |x| < 13 Mm: the mass against the first snapshot's; its change since the previous snapshot,
    split into what the flow carried across the box's edges (advective, trapezoid in time) and the rest
    (non-advective: the DynSGS mass diffusion, the density floor), in units of the first snapshot's mass; the
    mass fractions below 2e4 K, 2e4-2e5 K and above; the mass-weighted median of p over the first snapshot's p
    at the same point (the initial atmosphere: > 1 is plasma hotter than pressure balance with it allows); the
    mass-weighted mean v_y (km/s, < 0: draining); and, for the plasma denser than 10^-13 g cm^-3, the time in
    which the DynSGS viscous (rho mu |grad v|^2) and Ohmic (eta |curl B|^2) heating, estimated from the snapshot,
    would supply its internal energy. Plasma that is hot and in no hurry to drain is held up by its own pressure."""
    items = runs_of(a.run)
    rows = []
    p0 = None
    for t, f in items:
        g = load(f, ["log10rho_cgs", "ρ", "p", "T", "u", "v", "w", "Bz", MU_NAME, *ETA_NAMES])
        if p0 is None:
            p0 = grid2d(g, g["fields"]["p"])
        rows.append((t, mass_row(g, p0)))
    M0 = rows[0][1]["M"]
    print("box 3 < y < 20 Mm, |x| < 13 Mm; masses in units of the first snapshot's (" + f"{M0:.3g} g/cm)")
    print(f"{'t [s]':>7} {'mass':>8} {'change':>8} {'advect.':>8} {'other':>8} | {'<2e4 K':>6} {'2e4-2e5':>7} {'>2e5 K':>6} | "
          f"{'p/p_init':>8} {'v_y km/s':>8} | {'heat visc s':>11} {'heat Ohm s':>10}")
    prev = None
    for t, r in rows:
        if prev is None:
            ch = ad = ot = np.nan
        else:
            tp, rp = prev
            ch = (r["M"] - rp["M"]) / M0
            ad = 0.5 * (r["Fin"] + rp["Fin"]) * (t - tp) / M0
            ot = ch - ad
        print(f"{t:7.0f} {r['M'] / M0:8.2f} {ch:8.2f} {ad:8.2f} {ot:8.2f} | {r['cold']:6.2f} {r['warm']:7.2f} {r['hot']:6.2f} | "
              f"{r['pmed']:8.3g} {r['vy']:8.2f} | {r['tv']:11.3g} {r['to']:10.3g}")
        prev = (t, r)


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
    p.add_argument("--no-chain", action="store_true", help="only the directories given, not the restarts of the run")
    p.add_argument("--nearest", action="store_true", help="the nearest snapshot for a time the run has none at")
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
    s = sub.add_parser("mass"); s.set_defaults(fn=cmd_mass)
    s.add_argument("run")
    s = sub.add_parser("null"); s.set_defaults(fn=cmd_null)
    s.add_argument("runs", nargs="+", help="RUN or LABEL=RUN")
    s.add_argument("--thr", type=float, default=-13.0, help="log10 rho threshold of the dome top [g cm^-3]")
    s = sub.add_parser("eshear"); s.set_defaults(fn=cmd_eshear)
    s.add_argument("runs", nargs="+", help="RUN or LABEL=RUN")
    s = sub.add_parser("profiles"); s.set_defaults(fn=cmd_profiles)
    s.add_argument("run")
    s = sub.add_parser("atmos"); s.set_defaults(fn=cmd_atmos)
    s.add_argument("run")
    s.add_argument("--times", default="", help="comma-separated times in s (default 0,464,920,1160,1380,1610,1955,2530)")
    s.add_argument("--xmin", type=float, default=8.0, help="|x| range of the median, Mm")
    s.add_argument("--xmax", type=float, default=12.75)
    s.add_argument("--ylim", type=lambda t: tuple(map(float, t.split(","))), default=(-1.0, 20.0))
    s = sub.add_parser("diag"); s.set_defaults(fn=cmd_diag)
    s.add_argument("run")
    s.add_argument("--dt", type=float, default=None, help="the run's time step in t0 (its log): print the CFL number")
    s.add_argument("--ymin", type=float, default=None, help="only y >= YMIN Mm (e.g. 2 for the corona)")
    s.add_argument("--tmin", type=float, default=None, help="only snapshots from TMIN s")
    s.add_argument("--tfloor", type=float, default=0.1, help="the run's temperature floor (T_ph): count the nodes at it")
    a = p.parse_args()
    global FOLLOW_CHAIN, NEAREST
    FOLLOW_CHAIN, NEAREST = not a.no_chain, a.nearest
    os.makedirs(a.outdir, exist_ok=True)
    a.fn(a)


if __name__ == "__main__":
    main()
