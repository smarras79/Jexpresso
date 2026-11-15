#!/usr/bin/env python3
"""Figures of the 2.5D cross-section of Leake, Linton & Antiochos (2022, ApJ 934:10), from the output of
problems/MHD/fluxEmergenceLeake2022_2D (its README, "Output and figures").

RUN is an output directory (output/MHD/fluxEmergenceLeake2022_2D/output-<date>, the one holding iter_N.pvtu); a
restarted run is DIR1,DIR2,... (snapshots merged by time, a later directory wins a repeated time). LABEL=RUN names a
run in `evol`. x is horizontal and h (the code's y, the paper's z) up, in L0; t in t0; B in B0. The rope's axial
field is the code's Bz; its twist and the dipole are the in-plane (Bx, By), drawn as field lines.

  ic      RUN                                   the t = 0 state on x = 0 and at the side against the reference
                                                atmosphere: T, rho, p and beta (the paper's Fig. 1)
  snap    RUN --var log10_rho --times 0,20,40   a field with the in-plane field lines (and --vectors) at each
                                                time, one file each; --grid: all of them in one figure
  column  RUN --times 0,10,20 [--x 0]           rho/rho_e, T, p/p_e, v, Bz and Bx along a vertical line
  evol    LABEL=RUN ... [--csv FILE]            time series, runs overlaid: the rope's axis and top on x = 0, its
                                                axial flux above h = 0 and 10, the transition region on x = 0,
                                                kinetic and magnetic energy above h = 0, max T, min rho/rho_e,
                                                max |v|, nodes at the positivity floors

snap variables: log10_rho (default), rho, T, log10_T, p, log10_p, beta, log10_beta, Bx, By, Bz, B, vA, u, v, w,
speed, mu (the largest DynSGS coefficient), and rho_rel, p_rel, T_rel: log10 of rho, p, T over the reference
atmosphere (the case's closed form, without the rope; --mf and --cz-gradient if the run changed :lk2d_mf or
:lk2d_cz_gradient). --rho-floor and --p-floor are the run's LK2D_RHO_FLOOR and LK2D_P_FLOOR (default 0.2), for the
floor-node counts of `evol` and the floor lines of `column`.

Every iter_N.pvtu of a directory counts, at the time written in it, whether simulation.pvd lists it or not.
Needs numpy and matplotlib (it uses the reader and the flux function of tools/plot_fluxemergence_vtu.py).
Run any command with -h for its options.
"""
import argparse, csv, glob, os, sys
import xml.etree.ElementTree as ET
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_fluxemergence_vtu import read_pvtu, pvtu_time, grid2d, resolve_cmap   # noqa: E402
import matplotlib.pyplot as plt                                          # noqa: E402  (Agg: set by that module)
from matplotlib.colors import Normalize                                  # noqa: E402

GAMMA = 5.0 / 3.0
TCOR = 150.0                  # T_cor / T_ph
T_TR = np.sqrt(TCOR)          # 12.25: the middle of the transition region, at h = 15 in the reference
FIELDS = ("ρ", "p", "T", "u", "v", "w", "Bx", "By", "Bz")


# ---------------------------------------------------------------------------------------------------------
# Snapshots, reference atmosphere, one snapshot on the node grid
# ---------------------------------------------------------------------------------------------------------
def snapshots(run):
    """[(t, iter_N.pvtu)] of RUN (DIR or DIR1,DIR2,...), by time; a later directory wins a repeated time.
    Every iter_N.pvtu counts, at the time written in it (simulation.pvd only for a file without one), so files
    copied in from a restart, or written after the .pvd was copied, are not lost."""
    found = {}
    for d in [s for s in run.split(",") if s]:
        if not os.path.isdir(d):
            raise SystemExit(f"{d}: not a directory (give the output-<date> directory of a run)")
        listed = {}
        pvd = os.path.join(d, "simulation.pvd")
        if os.path.isfile(pvd):
            try:
                for ds in ET.parse(pvd).getroot().iter("DataSet"):
                    listed[os.path.basename(ds.get("file"))] = float(ds.get("timestep"))
            except (ET.ParseError, TypeError, ValueError):
                pass
        files = glob.glob(os.path.join(d, "iter_*.pvtu"))
        if not files:
            raise SystemExit(f"{d}: no iter_N.pvtu (give the output-<date> directory of a run)")
        for f in files:
            try:
                t = pvtu_time(f)
            except (ET.ParseError, OSError, ValueError):
                t = None                      # still being written
            t = listed.get(os.path.basename(f)) if t is None else t
            if t is None:
                print(f"   {f}: no time in it or in simulation.pvd, skipped")
                continue
            found[round(t, 9)] = (t, f)
    return [found[k] for k in sorted(found)]


def pick_times(snaps, spec):
    """Snapshots nearest to the requested times ('all', or comma-separated numbers)."""
    if not snaps:
        raise SystemExit("no snapshots")
    if spec in (None, "all"):
        return snaps
    ts = np.array([t for t, _ in snaps])
    out = []
    for tok in spec.split(","):
        want = float(tok)
        k = int(np.argmin(np.abs(ts - want)))
        if abs(ts[k] - want) > 1e-6 * max(1.0, abs(want)):
            print(f"   t = {want:g}: nearest snapshot is t = {ts[k]:g}")
        out.append(snaps[k])
    return out


def reference(h, mf=1.0, cz=1.0):
    """T, p, rho of the reference atmosphere at heights h (initialize.jl: hydrostatic with g = 1, rho(0) = 1,
    p = rho T / m_f; T = 1 - a (gamma-1)/gamma h below h = 0, 1 to 10, 150^((h-10)/10) to 20, 150 above)."""
    h = np.asarray(h, float)
    T, p = np.empty_like(h), np.empty_like(h)
    c, al, p0 = cz * (GAMMA - 1.0) / GAMMA, np.log(TCOR) / 10.0, 1.0 / mf
    p10 = p0 * np.exp(-mf * 10.0)
    p20 = p10 * np.exp(-mf * (1.0 - np.exp(-al * 10.0)) / al)
    m = h <= 0.0
    T[m] = 1.0 - c * h[m]
    p[m] = p0 * T[m] ** (mf / c) if c > 0 else p0 * np.exp(-mf * h[m])
    m = (h > 0.0) & (h < 10.0)
    T[m], p[m] = 1.0, p0 * np.exp(-mf * h[m])
    m = (h >= 10.0) & (h < 20.0)
    T[m], p[m] = TCOR ** ((h[m] - 10.0) / 10.0), p10 * np.exp(-mf * (1.0 - np.exp(-al * (h[m] - 10.0))) / al)
    m = h >= 20.0
    T[m], p[m] = TCOR, p20 * np.exp(-mf * (h[m] - 20.0) / TCOR)
    return T, p, mf * p / T


def dct_basis(n):
    """Orthonormal DCT-II matrix Q[k, i] = s_k cos(pi k (2i+1) / 2n) and the eigenvalues 2 - 2 cos(pi k / n) of the
    path-graph Laplacian it diagonalizes."""
    k, i = np.arange(n)[:, None], np.arange(n)[None, :]
    Q = np.cos(np.pi * k * (2 * i + 1) / (2 * n)) * np.sqrt(2.0 / n)
    Q[0] /= np.sqrt(2.0)
    return Q, 2.0 - 2.0 * np.cos(np.pi * np.arange(n) / n)


def flux_function(X, Z, bx, by):
    """Flux function A of the in-plane field on the (ny, nx) node grid, dA/dz = Bx and dA/dx = -By in the least-
    squares sense over the grid's edges (trapezoidal increments), A = 0 at the bottom-left corner. Integrating
    along paths instead carries the discrete div B of the rope along the path: up a column through the rope, it
    puts a step of several field lines' worth into the corona above it. The normal equations are the grid
    graph's Laplacian, solved exactly by cosine transforms."""
    gz = 0.5 * (bx[1:, :] + bx[:-1, :]) * np.diff(Z)[:, None]      # vertical edges (i, j) -> (i+1, j)
    gx = -0.5 * (by[:, 1:] + by[:, :-1]) * np.diff(X)[None, :]     # horizontal edges (i, j) -> (i, j+1)
    c = np.zeros_like(bx)
    c[1:, :] += gz; c[:-1, :] -= gz
    c[:, 1:] += gx; c[:, :-1] -= gx
    (Qz, lz), (Qx, lx) = dct_basis(len(Z)), dct_basis(len(X))
    den = lz[:, None] + lx[None, :]
    den[0, 0] = 1.0
    ah = (Qz @ c @ Qx.T) / den
    ah[0, 0] = 0.0
    A = Qz.T @ ah @ Qx
    return A - A[0, 0]


class Snap:
    """One snapshot on its tensor-product node grid: fields per node, the reference atmosphere, the flux
    function A of (Bx, By) (A = 0 at the domain's bottom-left corner) and the linear sub-cells."""

    def __init__(self, path, a, need_mu=False):
        g = read_pvtu(path, None if need_mu else set(FIELDS))
        missing = [n for n in FIELDS if n not in g["fields"]]
        if missing:
            raise SystemExit(f"{path}: no field(s) {', '.join(missing)} (is this a fluxEmergenceLeake2022_2D run?)")
        self.path, self.g, self.t = path, g, g["t"]
        self.x, self.h = g["x"], g["y"]
        self.f = g["fields"]
        self.Te, self.pe, self.rhoe = reference(self.h, a.mf, a.cz_gradient)
        if len(g["xc"]) * len(g["yc"]) != len(self.x):
            raise SystemExit(f"{path}: the nodes are not a tensor-product grid (not a LK22_2D.geo mesh?)")
        A = flux_function(g["xc"], g["yc"], grid2d(g, self.f["Bx"]), grid2d(g, self.f["By"]))
        self.A = A[g["iy"], g["ix"]]
        q = g["quads"]
        X, Y = self.x[q], self.h[q]
        self.q, self.hc = q, Y.mean(1)
        self.area = 0.5 * np.abs((X[:, 0] * Y[:, 1] - X[:, 1] * Y[:, 0]) + (X[:, 1] * Y[:, 2] - X[:, 2] * Y[:, 1])
                                 + (X[:, 2] * Y[:, 3] - X[:, 3] * Y[:, 2]) + (X[:, 3] * Y[:, 0] - X[:, 0] * Y[:, 3]))

    def var(self, name):
        f = self.f
        if name.startswith("mu"):
            mus = [v for k, v in f.items() if k.startswith("mu_dsgs")]
            if not mus:
                raise SystemExit(f"{self.path}: no DynSGS coefficient (mu_dsgs_*) in the output")
            return np.max(mus, axis=0)
        speed2 = f["u"] ** 2 + f["v"] ** 2 + f["w"] ** 2
        B2 = f["Bx"] ** 2 + f["By"] ** 2 + f["Bz"] ** 2
        lg = lambda v: np.log10(np.maximum(v, 1e-300))
        table = {"rho": lambda: f["ρ"], "log10_rho": lambda: lg(f["ρ"]), "T": lambda: f["T"],
                 "log10_T": lambda: lg(f["T"]), "p": lambda: f["p"], "log10_p": lambda: lg(f["p"]),
                 "beta": lambda: 2.0 * f["p"] / np.maximum(B2, 1e-300),
                 "log10_beta": lambda: lg(2.0 * f["p"] / np.maximum(B2, 1e-300)),
                 "Bx": lambda: f["Bx"], "By": lambda: f["By"], "Bz": lambda: f["Bz"], "B": lambda: np.sqrt(B2),
                 "vA": lambda: np.sqrt(B2 / np.maximum(f["ρ"], 1e-300)),
                 "u": lambda: f["u"], "v": lambda: f["v"], "w": lambda: f["w"], "speed": lambda: np.sqrt(speed2),
                 "rho_rel": lambda: lg(f["ρ"] / self.rhoe), "p_rel": lambda: lg(f["p"] / self.pe),
                 "T_rel": lambda: lg(f["T"] / self.Te)}
        if name not in table:
            raise SystemExit(f"unknown variable '{name}'; one of {', '.join(table)}, mu")
        return table[name]()

    def column(self, x0, *arrays):
        """h and the arrays on the node column nearest x = x0 (bottom to top)."""
        g = self.g
        ix = int(np.argmin(np.abs(g["xc"] - x0)))
        return (g["yc"],) + tuple(grid2d(g, v)[:, ix] for v in arrays)

    def integral(self, v, mask=None):
        c = v[self.q].mean(1) * self.area
        return c.sum() if mask is None else c[mask].sum()


VAR_LABEL = {"log10_rho": r"$\log_{10}\rho$", "rho": r"$\rho$", "T": "$T$", "log10_T": r"$\log_{10}T$", "p": "$p$",
             "log10_p": r"$\log_{10}p$", "beta": r"$\beta$", "log10_beta": r"$\log_{10}\beta$",
             "Bx": "$B_x$", "By": "$B_y$ (vertical)", "Bz": "$B_z$ (axial)", "B": r"$|\mathbf{B}|$",
             "vA": "$v_A$", "u": "$u$", "v": "$v$ (vertical)", "w": "$w$ (axial)", "speed": r"$|\mathbf{v}|$",
             "mu": "DynSGS coefficient", "rho_rel": r"$\log_{10}(\rho/\rho_e)$", "p_rel": r"$\log_{10}(p/p_e)$",
             "T_rel": r"$\log_{10}(T/T_e)$"}
DIVERGING = {"Bx", "By", "Bz", "u", "v", "w", "rho_rel", "p_rel", "T_rel"}
CMAP = {"T": "inferno", "log10_T": "inferno", "speed": "magma", "mu": "viridis", "beta": "viridis",
        "log10_beta": "viridis", "B": "viridis", "vA": "viridis"}


# ---------------------------------------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------------------------------------
def cmd_ic(a):
    snaps = snapshots(a.run)
    t, f = pick_times(snaps, "0")[0]
    s = Snap(f, a)
    hh = np.linspace(-30.0, 60.0, 1000)
    Te, pe, rhoe = reference(hh, a.mf, a.cz_gradient)
    fig, axs = plt.subplots(1, 4, figsize=(16, 4.2))
    for x0, ls, lab in ((0.0, "-", "x = 0 (the rope)"), (a.xside, "--", f"x = {a.xside:g}")):
        h, T, p, rho, Bx, By, Bz = s.column(x0, s.f["T"], s.f["p"], s.f["ρ"], s.f["Bx"], s.f["By"], s.f["Bz"])
        beta = 2.0 * p / np.maximum(Bx ** 2 + By ** 2 + Bz ** 2, 1e-300)
        for ax, v in zip(axs, (T, rho, p, beta)):
            ax.semilogy(h, v, ls, lw=1.4, label=lab)
        if x0 != 0.0:
            m = (h >= -30) & (h <= 60)
            Tr, pr, rr = reference(h[m], a.mf, a.cz_gradient)
            dev = [np.max(np.abs(v[m] / r - 1.0)) for v, r in ((T, Tr), (p, pr), (rho, rr))]
            print(f"   {f}: t = {t:g}; at x = {x0:g}, -30 <= h <= 60, |T/T_e - 1|, |p/p_e - 1|, |rho/rho_e - 1| "
                  f"<= {dev[0]:.1e}, {dev[1]:.1e}, {dev[2]:.1e}")
    for ax, v, name in zip(axs, (Te, rhoe, pe, None), ("$T$", r"$\rho$", "$p$", r"$\beta$")):
        if v is not None:
            ax.semilogy(hh, v, "k:", lw=1.0, label="reference")
        ax.set_xlabel("h / L0"); ax.set_title(name); ax.set_xlim(-30, 60); ax.grid(alpha=0.3)
    axs[0].legend(fontsize=8)
    fig.suptitle(rf"t = {t:g} $t_0$: the atmosphere on x = 0 and at the side", fontsize=11)
    fig.tight_layout()
    out = os.path.join(a.outdir, f"{a.name or 'ic'}.{a.format}")
    fig.savefig(out, dpi=a.dpi); plt.close(fig)
    print(f"   -> {out}")


def color_range(a, vals):
    if a.clim:
        return tuple(a.clim)
    lo, hi = min(np.nanmin(v) for v in vals), max(np.nanmax(v) for v in vals)
    if a.var in DIVERGING:
        m = max(abs(lo), abs(hi), 1e-300)
        return -m, m
    return lo, (hi if hi > lo else lo + 1.0)


def draw_snap(ax, s, a, V, levels, vmin, vmax, cmap):
    g = s.g
    Z = grid2d(g, V)
    ax.set_rasterization_zorder(0)
    ax.contourf(g["xc"], g["yc"], np.clip(Z, vmin, vmax), levels=np.linspace(vmin, vmax, 129), cmap=cmap,
                norm=Normalize(vmin, vmax), zorder=-1, antialiased=False)
    if levels is not None and len(levels):
        ax.contour(g["xc"], g["yc"], grid2d(g, s.A), levels=levels, colors=a.line_color, linewidths=a.linewidth,
                   linestyles="solid")
    if a.vectors:
        x, y = s.x, s.h
        (X0, X1), (Y0, Y1) = a.xlim, a.ylim
        xs = X0 + (X1 - X0) * np.arange(1, a.vec_n[0] + 1) / (a.vec_n[0] + 1)
        ys = Y0 + (Y1 - Y0) * np.arange(1, a.vec_n[1] + 1) / (a.vec_n[1] + 1)
        ip = np.array([np.argmin((x - px) ** 2 + (y - py) ** 2) for py in ys for px in xs])
        U, W = s.f["u"][ip], s.f["v"][ip]
        ref = a.vec_ref if a.vec_ref else max(np.max(np.hypot(U, W)), 1e-300)
        keep = np.hypot(U, W) >= 0.02 * ref
        kw = dict(angles="xy", scale_units="xy", scale=ref / (0.08 * (X1 - X0)), color=a.vec_color, width=0.002,
                  headwidth=4, headlength=4, headaxislength=3.5, zorder=3)
        ax.quiver(x[ip][keep], y[ip][keep], U[keep], W[keep], **kw)
        x0, y0 = X0 + 0.03 * (X1 - X0), Y0 + 0.04 * (Y1 - Y0)
        ax.quiver([x0], [y0], [ref], [0.0], **kw)
        ax.text(x0 + 0.1 * (X1 - X0), y0, f"{ref:.3g} $v_0$", color=a.vec_color, fontsize=8, va="center", zorder=3)
    ax.set_xlim(*a.xlim); ax.set_ylim(*a.ylim); ax.set_aspect("equal")


def cmd_snap(a):
    picks = pick_times(snapshots(a.run), a.times)
    snaps = [Snap(f, a, need_mu=a.var.startswith("mu")) for _, f in picks]
    (X0, X1), (Y0, Y1) = a.xlim, a.ylim
    wins = [(s.x >= X0) & (s.x <= X1) & (s.h >= Y0) & (s.h <= Y1) for s in snaps]
    vals = [s.var(a.var) for s in snaps]
    vmin, vmax = color_range(a, [v[w] for v, w in zip(vals, wins)])
    cmap = resolve_cmap(a.cmap or CMAP.get(a.var, "RdBu_r" if a.var in DIVERGING else "RainbowDesaturated"))
    # field lines: contours of A at fixed values, so that a line follows the same flux surface at every time (A is
    # gauged at the domain's bottom-left corner, where the field stays the dipole's). Two families set at the first
    # time: --nlines across the window's ambient field (|Bz| < 1e-3, outside the rope) and --nrope inside the rope,
    # from its axis (the largest |Bz|) to its edge (|Bz| = --edge); or with --dA one line per dA of flux everywhere
    # (the rope's twist flux, B_t w / 2 = 6.25, is a hundred times the dipole's across the window).
    levels = None
    if a.fieldlines:
        s0, w0 = snaps[0], wins[0]
        if a.dA:
            lo = min(s.A[w].min() for s, w in zip(snaps, wins))
            hi = max(s.A[w].max() for s, w in zip(snaps, wins))
            levels = (np.arange(np.floor(lo / a.dA), np.ceil(hi / a.dA)) + 0.5) * a.dA
        else:
            Bz0, levels = np.abs(s0.f["Bz"]), []
            amb = w0 & (Bz0 < 1e-3)
            lo = hi = None
            if amb.any() and a.nlines > 0:
                lo, hi = s0.A[amb].min(), s0.A[amb].max()
                levels += list(lo + (np.arange(a.nlines) + 0.5) * (hi - lo) / a.nlines)
            k = int(np.argmax(Bz0))
            ring = (Bz0 > 0.5 * a.edge) & (Bz0 < 1.5 * a.edge)
            if Bz0[k] > a.edge and ring.any() and a.nrope > 0:
                # from the axis to the rope's edge, or to the ambient field's last value if that comes first: a
                # level the ambient field also takes would draw a second line in it
                A_axis, A_edge = s0.A[k], np.median(s0.A[ring])
                if lo is not None and A_axis > hi:
                    A_edge = max(A_edge, hi)
                elif lo is not None and A_axis < lo:
                    A_edge = min(A_edge, lo)
                levels += list(A_axis + (A_edge - A_axis) * (np.arange(a.nrope) + 0.5) / a.nrope)
            levels = np.unique(levels)
    if a.vec_color is None:
        a.vec_color = "k" if (a.var in DIVERGING and not a.cmap) else "white"
    if a.line_color is None:
        a.line_color = "0.7" if getattr(cmap, "name", "") in ("inferno", "magma") else "k"
    label = VAR_LABEL.get(a.var, a.var)
    W, H = X1 - X0, Y1 - Y0
    if a.grid:
        nc = min(a.ncols, len(snaps)); nr = -(-len(snaps) // nc)
        fig, axs = plt.subplots(nr, nc, figsize=(a.width, a.width / nc * H / W * nr + 0.6), squeeze=False)
        for k, (s, V) in enumerate(zip(snaps, vals)):
            ax = axs[k // nc][k % nc]
            draw_snap(ax, s, a, V, levels, vmin, vmax, cmap)
            ax.set_title(rf"$t = {s.t:g}\,t_0$", fontsize=10, loc="left")
            if k % nc == 0:
                ax.set_ylabel("h / L0")
            if k // nc == nr - 1:
                ax.set_xlabel("x / L0")
        for k in range(len(snaps), nr * nc):
            axs[k // nc][k % nc].axis("off")
        fig.colorbar(plt.cm.ScalarMappable(Normalize(vmin, vmax), cmap), ax=axs, shrink=0.85, pad=0.02, label=label)
        out = os.path.join(a.outdir, f"{a.name or 'snap_' + a.var}.{a.format}")
        fig.savefig(out, dpi=a.dpi, bbox_inches="tight"); plt.close(fig)
        print(f"   {len(snaps)} snapshots -> {out}")
        return
    for s, V in zip(snaps, vals):
        fig, ax = plt.subplots(figsize=(a.width, a.width * H / W + 0.8))
        draw_snap(ax, s, a, V, levels, vmin, vmax, cmap)
        ax.set_xlabel("x / L0"); ax.set_ylabel("h / L0")
        ax.set_title(rf"{label}   $t = {s.t:g}\,t_0$", fontsize=11)
        fig.colorbar(plt.cm.ScalarMappable(Normalize(vmin, vmax), cmap), ax=ax, shrink=0.85, pad=0.02)
        out = os.path.join(a.outdir, f"{a.name or 'snap_' + a.var}_t{s.t:g}.{a.format}")
        fig.savefig(out, dpi=a.dpi, bbox_inches="tight"); plt.close(fig)
        print(f"   {s.path}: t = {s.t:g} -> {out}")


def cmd_column(a):
    picks = pick_times(snapshots(a.run), a.times)
    fig, axs = plt.subplots(2, 3, figsize=(15, 7.5), sharex=True)
    colors = plt.cm.viridis(np.linspace(0.0, 0.9, len(picks)))
    for (t, f), c in zip(picks, colors):
        s = Snap(f, a)
        h, rho, T, p, v, Bz, Bx, rhoe, pe = s.column(a.x, s.f["ρ"], s.f["T"], s.f["p"], s.f["v"], s.f["Bz"],
                                                     s.f["Bx"], s.rhoe, s.pe)
        m = (h >= a.hlim[0]) & (h <= a.hlim[1])
        for ax, val in zip(axs.flat, (rho / rhoe, T, p / pe, v, Bz, Bx)):
            ax.plot(h[m], val[m], color=c, lw=1.3, label=rf"$t = {t:g}$")
        print(f"   {f}: t = {t:g}")
    hh = np.linspace(*a.hlim, 600)
    axs[0, 1].plot(hh, reference(hh, a.mf, a.cz_gradient)[0], "k:", lw=1.0, label="reference")
    for ax, name, logy in zip(axs.flat, (r"$\rho/\rho_e$", "$T$", "$p/p_e$", "$v$ (vertical)", "$B_z$ (axial)", "$B_x$"),
                              (True, True, True, False, False, False)):
        ax.set_title(name, fontsize=11); ax.grid(alpha=0.3)
        if logy:
            ax.set_yscale("log")
    axs[0, 0].axhline(a.rho_floor, color="r", lw=0.8, ls="--")
    axs[0, 2].axhline(a.p_floor, color="r", lw=0.8, ls="--")
    for ax in axs[1]:
        ax.set_xlabel("h / L0")
    axs[0, 1].legend(fontsize=8)
    fig.suptitle(f"x = {a.x:g}", fontsize=11)
    fig.tight_layout()
    out = os.path.join(a.outdir, f"{a.name or 'column_x' + format(a.x, 'g')}.{a.format}")
    fig.savefig(out, dpi=a.dpi); plt.close(fig)
    print(f"   -> {out}")


EVOL = [("axis", "rope axis h (O-point on x = 0)"), ("top", "rope top h (|Bz| = edge on x = 0)"),
        ("flux0", "axial flux above h = 0 (fraction)"), ("flux10", "axial flux above h = 10 (fraction)"),
        ("h_tr", "T = 12.25 on x = 0: h"), ("Ek", "kinetic energy, h > 0"), ("dEm", "magnetic energy, h > 0, change"),
        ("Tmax", "max T"), ("rho_min", "min rho/rho_e"), ("vmax", "max |v|"), ("floors", "nodes at the floors"),
        ("Bz_max", "max |Bz|")]


def metrics(s, a):
    f = s.f
    h, Bz, T, Bx = s.column(0.0, np.abs(f["Bz"]), f["T"], f["Bx"])
    # the axis: the O-point of the in-plane field on x = 0 (By = 0 there by symmetry), i.e. the sign change of Bx
    # inside the rope (|Bz| >= edge) nearest the |Bz|-weighted centre of the column. The largest |Bz| is no
    # good: the axial field flattens, and takes two maxima, as the rope presses against the photosphere.
    rope = Bz >= a.edge
    axis = np.nan
    if rope.any():
        hc = np.sum(h[rope] * Bz[rope]) / np.sum(Bz[rope])
        j = np.nonzero(rope[:-1] & rope[1:] & (np.sign(Bx[:-1]) * np.sign(Bx[1:]) < 0))[0]
        if len(j):
            z = h[j] + Bx[j] / (Bx[j] - Bx[j + 1]) * (h[j + 1] - h[j])
            axis = z[np.argmin(np.abs(z - hc))]
        else:
            axis = hc
    up = np.nonzero(Bz >= a.edge)[0]
    top = np.nan
    if len(up):
        j = up[-1]
        top = h[j] if j == len(h) - 1 else h[j] + (Bz[j] - a.edge) / (Bz[j] - Bz[j + 1]) * (h[j + 1] - h[j])
    m = np.nonzero((h > 2.0) & (T >= T_TR))[0]
    h_tr = np.nan
    if len(m) and m[0] > 0:
        j = m[0]
        h_tr = h[j - 1] + (T_TR - T[j - 1]) / (T[j] - T[j - 1]) * (h[j] - h[j - 1])
    phi = s.integral(f["Bz"])
    sgn = 1.0 if phi >= 0 else -1.0
    tot = abs(phi) if phi != 0 else np.nan
    flux0 = sgn * s.integral(f["Bz"], s.hc > 0.0) / tot
    flux10 = sgn * s.integral(f["Bz"], s.hc > 10.0) / tot
    above = s.hc > 0.0
    Ek = s.integral(0.5 * f["ρ"] * (f["u"] ** 2 + f["v"] ** 2 + f["w"] ** 2), above)
    Em = s.integral(0.5 * (f["Bx"] ** 2 + f["By"] ** 2 + f["Bz"] ** 2), above)
    rr, pr = f["ρ"] / s.rhoe, f["p"] / s.pe
    floors = int(np.count_nonzero((rr <= a.rho_floor * 1.001) | (pr <= a.p_floor * 1.001)))
    return dict(t=s.t, axis=axis, top=top, flux0=flux0, flux10=flux10, h_tr=h_tr, Ek=Ek, Em=Em,
                Tmax=float(np.max(f["T"])), rho_min=float(np.min(rr)),
                vmax=float(np.sqrt(np.max(f["u"] ** 2 + f["v"] ** 2 + f["w"] ** 2))), floors=floors,
                Bz_max=float(np.max(np.abs(f["Bz"]))))


def cmd_evol(a):
    runs = []
    for spec in a.runs:
        lab, _, run = spec.rpartition("=")
        lab = lab or os.path.basename(run.split(",")[-1].rstrip("/"))
        rows = []
        for t, f in pick_times(snapshots(run), a.times):
            try:
                rows.append(metrics(Snap(f, a), a))
            except (ValueError, OSError) as e:
                print(f"   {f}: skipped ({e})")
        if not rows:
            continue
        for r in rows:
            r["dEm"] = r["Em"] - rows[0]["Em"]
        runs.append((lab, rows))
        print(f"== {lab}   ({run})")
        print("      t    axis     top  flux>0  flux>10   h_TR       Ek      dEm     Tmax  rho/rho_e    |v|  floors  |Bz|max"
              "      (axis: O-point, top: |Bz| = edge, both on x = 0)")
        for r in rows:
            print(f"  {r['t']:6.1f} {r['axis']:7.2f} {r['top']:7.2f} {r['flux0']:7.3f} {r['flux10']:8.3f} {r['h_tr']:6.2f} "
                  f"{r['Ek']:8.2e} {r['dEm']:8.2e} {r['Tmax']:8.1f} {r['rho_min']:10.3f} {r['vmax']:6.2f} {r['floors']:7d} "
                  f"{r['Bz_max']:8.3f}")
    if not runs:
        raise SystemExit("no snapshot could be read")
    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            keys = ["t"] + [k for k, _ in EVOL]
            w.writerow(["run"] + keys)
            for lab, rows in runs:
                for r in rows:
                    w.writerow([lab] + [r[k] for k in keys])
        print(f"   table -> {a.csv}")
    panels = [("axis", "top"), ("flux0", "flux10"), ("h_tr",), ("Ek",), ("dEm",), ("Tmax",), ("rho_min",), ("vmax",),
              ("floors",)]
    titles = ["rope on x = 0: axis (O-point, solid), top (dashed)", "axial flux above h = 0 (solid), 10 (dashed)",
              "T = 12.25 on x = 0 (reference: h = 15)", "kinetic energy above h = 0", "magnetic energy above h = 0, change",
              "max T", r"min $\rho/\rho_e$", r"max $|\mathbf{v}|$", "nodes at the positivity floors"]
    fig, axs = plt.subplots(3, 3, figsize=(16, 11), sharex=True)
    for (lab, rows), c in zip(runs, plt.rcParams["axes.prop_cycle"].by_key()["color"] * 4):
        t = [r["t"] for r in rows]
        for ax, keys in zip(axs.flat, panels):
            for k, ls in zip(keys, ("-", "--")):
                y = [r[k] if (k != "Ek" or r[k] > 0) else np.nan for r in rows]
                ax.plot(t, y, ls, color=c, marker=".", ms=4, lw=1.3, label=lab if ls == "-" else None)
    for ax, title in zip(axs.flat, titles):
        ax.set_title(title, fontsize=10); ax.grid(alpha=0.3)
    axs[0, 0].axhline(0.0, color="k", lw=0.6)
    axs[1, 0].set_yscale("log")
    axs[2, 0].axhline(a.rho_floor, color="r", lw=0.8, ls="--")
    for ax in axs[2]:
        ax.set_xlabel("t / t0")
    axs[0, 0].legend(fontsize=8)
    fig.tight_layout()
    out = os.path.join(a.outdir, f"{a.name or 'evol'}.{a.format}")
    fig.savefig(out, dpi=a.dpi); plt.close(fig)
    print(f"   -> {out}")


# ---------------------------------------------------------------------------------------------------------
def floats(text):
    return [float(v) for v in text.split(",")]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(s):
        s.add_argument("--outdir", default="figs_leake2022_2d", help="output directory (default figs_leake2022_2d)")
        s.add_argument("--name", default=None, help="file name stem (default per command)")
        s.add_argument("--format", default="png", help="png, pdf, svg, ...")
        s.add_argument("--dpi", type=int, default=150)
        s.add_argument("--mf", type=float, default=1.0, help="the run's :lk2d_mf (default 1)")
        s.add_argument("--cz-gradient", type=float, default=1.0, help="the run's :lk2d_cz_gradient (default 1)")
        s.add_argument("--rho-floor", type=float, default=0.2, help="the run's LK2D_RHO_FLOOR (default 0.2)")
        s.add_argument("--p-floor", type=float, default=0.2, help="the run's LK2D_P_FLOOR (default 0.2)")

    s = sub.add_parser("ic", help="the t = 0 atmosphere against the reference")
    s.add_argument("run"); s.add_argument("--xside", type=float, default=-200.0, help="the side column (default -200)")
    common(s); s.set_defaults(func=cmd_ic)

    s = sub.add_parser("snap", help="a field and the in-plane field lines at each time")
    s.add_argument("run")
    s.add_argument("--var", default="log10_rho", help="log10_rho (default), rho_rel, T, Bz, beta, speed, ... (see -h of the tool)")
    s.add_argument("--times", default="all", help="comma-separated times (nearest snapshots), or all (default)")
    s.add_argument("--xlim", type=floats, default=[-60.0, 60.0], help="x range, e.g. -60,60 (default)")
    s.add_argument("--ylim", type=floats, default=[-30.0, 90.0], help="h range, e.g. -30,90 (default)")
    s.add_argument("--clim", type=floats, default=None, help="color range, e.g. -1,1 (default: all times' range)")
    s.add_argument("--cmap", default=None, help="matplotlib colormap or RainbowDesaturated (default per variable)")
    s.add_argument("--fieldlines", action=argparse.BooleanOptionalAction, default=True)
    s.add_argument("--nlines", type=int, default=20,
                   help="field lines of the ambient (dipole) field across the window at the first time (default 20)")
    s.add_argument("--nrope", type=int, default=8, help="field lines inside the rope at the first time (default 8)")
    s.add_argument("--edge", type=float, default=0.1, help="|Bz| of the rope's edge (default 0.1 B0)")
    s.add_argument("--dA", type=float, default=None, help="instead: one field line per dA of in-plane flux (B0 L0)")
    s.add_argument("--linewidth", type=float, default=0.6)
    s.add_argument("--line-color", default=None, help="field lines' color (default black, light grey on inferno/magma)")
    s.add_argument("--vectors", action="store_true", help="the in-plane velocity")
    s.add_argument("--vec-n", type=int, nargs=2, default=(24, 20), help="arrows in x and h (default 24 20)")
    s.add_argument("--vec-ref", type=float, default=1.0, help="reference speed (default 1; 0 = the largest)")
    s.add_argument("--vec-color", default=None, help="default black on the diverging maps, else white")
    s.add_argument("--grid", action="store_true", help="all times in one figure")
    s.add_argument("--ncols", type=int, default=3, help="columns of --grid (default 3)")
    s.add_argument("--width", type=float, default=9.0, help="figure width in inches (default 9)")
    common(s); s.set_defaults(func=cmd_snap)

    s = sub.add_parser("column", help="profiles along a vertical line")
    s.add_argument("run")
    s.add_argument("--times", default="all")
    s.add_argument("--x", type=float, default=0.0, help="the line's x (nearest node column; default 0)")
    s.add_argument("--hlim", type=floats, default=[-30.0, 60.0], help="h range (default -30,60)")
    common(s); s.set_defaults(func=cmd_column)

    s = sub.add_parser("evol", help="time series of the emergence, runs overlaid")
    s.add_argument("runs", nargs="+", help="LABEL=RUN (or RUN); RUN may be DIR1,DIR2,...")
    s.add_argument("--times", default="all")
    s.add_argument("--edge", type=float, default=0.1, help="|Bz| of the rope's edge (default 0.1 B0; Bt = 5)")
    s.add_argument("--csv", default=None, help="also write the table to this file")
    common(s); s.set_defaults(func=cmd_evol)

    argv = sys.argv[1:]
    # let comma lists start with a minus sign: --xlim -60,60 -> --xlim=-60,60
    for k in range(len(argv) - 1):
        if argv[k] in ("--xlim", "--ylim", "--clim", "--hlim") and argv[k + 1].startswith("-"):
            argv[k], argv[k + 1] = argv[k] + "=" + argv[k + 1], ""
    a = p.parse_args([v for v in argv if v != ""])
    os.makedirs(a.outdir, exist_ok=True)
    a.func(a)


if __name__ == "__main__":
    main()
