#!/usr/bin/env python3
"""Figures of Leake, Linton & Antiochos (2014, ApJ 787:46) from the output of problems/MHD/fluxEmergenceLeake2014.

RUN is an output directory (output/MHD/fluxEmergenceLeake2014/output-<date>, the one holding iter_N.pvtu);
a restarted run can be given as DIR1,DIR2 (snapshots merged by time). LABEL=RUN names a curve; the labels
SD, MD, WD, ND get the papers' colors and line styles.

  ic      [LABEL=RUN ...]                          Fig. 1(a): T, p and beta(0,0,z) at t = 0 (analytic; runs overlaid)
  slice   RUN --plane y=0 --var j --times ...      Fig. 8(a-c); --plane z=0 --var Bz --vectors: Fig. 5(d-f);
                                                   --fieldlines --nulls: the (Bx, Bz) lines and O/X-points of Figs. 6-7
  line    RUN --var w --times ...                  Fig. 8(d-f): a field along the line x = y = 0
  height  LABEL=RUN ...                            Fig. 11 (O-point height in y = 0), Fig. 12 (its speed), Fig. 13 (--mono)
  axflux  LABEL=RUN ...                            Fig. 14(b): axial flux above z = 0 and 10 over the total, in y = 0
  render  RUN --preset NAME --times ...            3D field lines, Figs. 1(b), 2-7, 10 (needs pyvista)
  fig2    RUN --times 5,40,60                      Fig. 2: side view (a-c) and top view (d-f) at each time, one file
                                                   each, Bz at z = 0 on the paper's rainbow and ticks; --grid: all in
                                                   one figure (needs pyvista)

Units are the papers' (L0, t0, B0, v0, j0 = B0/(mu0 L0)): the case writes its output in them.
Needs numpy and matplotlib; render and fig2 also pyvista (pip install pyvista). Run any command with -h for its options.
"""
import argparse, glob, os, sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True          # no tools/__pycache__ in the work tree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_fluxemergence_vtu import read_vtu, _cluster, pvtu_time   # noqa: E402

import numpy as np                      # noqa: E402
import matplotlib                       # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt         # noqa: E402

trapz = getattr(np, "trapezoid", None) or np.trapz      # numpy >= 2 / < 2

# Leake, Linton & Antiochos (2014), Figs. 11-12: SD red solid, MD green dashed, WD blue dot-dashed, ND black dotted
PAPER_STYLE = {"SD": ("tab:red", "-"), "MD": ("tab:green", "--"), "WD": ("tab:blue", "-."), "ND": ("k", ":")}
MONO_LS = ["-", "--", "-.", ":", (0, (5, 1, 1, 1, 1, 1))]

# the case's constants (initialize.jl)
GAMMA, TCOR = 5.0 / 3.0, 150.0
ZT, RT, BT, ZD, LAM = -12.0, 2.5, 5.0, -100.0, 10.0
RATIO = {"SD": 0.13, "MD": 0.1, "WD": 0.067, "ND": 0.0}
ZTOP = 209.84530500030144


# ---------------------------------------------------------------------------------------------------------
# Output: snapshots, lattice reader, planes, derivatives
# ---------------------------------------------------------------------------------------------------------
def snapshots(run):
    """[(t, iter_N.pvtu)] of RUN (DIR or DIR1,DIR2,...), by time; a later directory wins a repeated time."""
    found = {}
    for d in run.split(","):
        pvd, items = os.path.join(d, "simulation.pvd"), []
        if os.path.isfile(pvd):
            for ds in ET.parse(pvd).getroot().iter("DataSet"):
                f = os.path.join(d, ds.get("file"))
                if os.path.isfile(f):
                    items.append((float(ds.get("timestep")), f))
        if not items:
            items = [(pvtu_time(f), f) for f in glob.glob(os.path.join(d, "iter_*.pvtu"))]
        if not items:
            raise SystemExit(f"{d}: no iter_N.pvtu (give the output-<date> directory of a run)")
        for t, f in items:
            if t is not None:
                found[round(t, 9)] = (t, f)
    return [found[k] for k in sorted(found)]


def pick_times(snaps, spec):
    """Snapshots nearest to the requested times ('all', or comma-separated numbers)."""
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


def read_lattice(path, keep, box=None):
    """Point fields `keep` of iter_N.pvtu on the tensor-product node lattice of LK.geo, restricted to
    box = ((x0, x1), (y0, y1), (z0, z1)) (None: all). Returns (x, y, z), {name: array[nz, ny, nx]}."""
    root = ET.parse(path).getroot()
    pieces = [os.path.join(os.path.dirname(path), p.get("Source")) for p in root.iter("Piece")]
    Ps, Vs = [], {k: [] for k in keep}
    for f in pieces:
        d = read_vtu(f)
        P = next(v for (s, _), v in d.items() if s == "Points").reshape(-1, 3)
        m = np.ones(len(P), bool)
        if box is not None:                     # padded: a lattice plane on the box edge must not split in roundoff
            pad = 1e-6 * max(np.abs(P).max(), 1.0)
            for a, (lo, hi) in enumerate(box):
                m &= (P[:, a] >= lo - pad) & (P[:, a] <= hi + pad)
        if not m.any():
            continue
        Ps.append(np.array(P[m]))
        for k in keep:
            if ("PointData", k) not in d:
                names = sorted(n for (s, n) in d if s == "PointData")
                raise SystemExit(f"{f}: no point field '{k}'. Available: {', '.join(names)}")
            Vs[k].append(np.array(d[("PointData", k)][m], float))
    if not Ps:
        raise SystemExit(f"{path}: no node inside the box {box}")
    P = np.concatenate(Ps)
    tol = 1e-7 * max(np.ptp(P, axis=0).max(), 1.0)
    (ix, xc), (iy, yc), (iz, zc) = (_cluster(P[:, a], tol) for a in range(3))
    shape = (len(zc), len(yc), len(xc))
    flat = (iz * shape[1] + iy) * shape[2] + ix
    hit = np.zeros(np.prod(shape), bool)
    hit[flat] = True
    if not hit.all():
        raise SystemExit(f"{path}: the nodes do not form a tensor-product lattice (is the mesh from LK.geo?)")
    F = {}
    for k in keep:
        a = np.empty(hit.size)
        a[flat] = np.concatenate(Vs[k])
        F[k] = a.reshape(shape)
    return (xc, yc, zc), F


AX = {"x": 2, "y": 1, "z": 0}           # array axis of each coordinate (arrays are [z, y, x])


def at_plane(coords, arr, axis, value):
    """arr on the plane axis = value (linear between the two bracketing lattice planes).
    Returns (h, v, A[nv, nh]): y-plane (x, z), z-plane (x, y), x-plane (y, z)."""
    a = AX[axis]
    c = coords["xyz".index(axis)]
    if not c[0] - 1e-9 <= value <= c[-1] + 1e-9:
        raise SystemExit(f"{axis} = {value} is outside the data [{c[0]:g}, {c[-1]:g}]")
    k = int(np.clip(np.searchsorted(c, value), 1, len(c) - 1)) if len(c) > 1 else 0
    if len(c) == 1 or abs(c[k] - value) < 1e-9 * max(1.0, abs(value)):
        A = np.take(arr, k, axis=a)
    elif abs(c[k - 1] - value) < 1e-9 * max(1.0, abs(value)):
        A = np.take(arr, k - 1, axis=a)
    else:
        w = (value - c[k - 1]) / (c[k] - c[k - 1])
        A = (1 - w) * np.take(arr, k - 1, axis=a) + w * np.take(arr, k, axis=a)
    x, y, z = coords
    return {"y": (x, z, A), "z": (x, y, A), "x": (y, z, A)}[axis]


def ddx(f, c, axis):
    return np.gradient(f, c, axis=axis) if f.shape[axis] > 1 else np.zeros_like(f)


def current(coords, F):
    """(jx, jy, jz) = curl B on the lattice (second order, nonuniform spacing), in j0 = B0/(mu0 L0)."""
    x, y, z = coords
    Bx, By, Bz = F["Bx"], F["By"], F["Bz"]
    return (ddx(Bz, y, 1) - ddx(By, z, 0), ddx(Bx, z, 0) - ddx(Bz, x, 2), ddx(By, x, 2) - ddx(Bx, y, 1))


DERIVED = {"j": ("Bx", "By", "Bz"), "jx": ("Bx", "By", "Bz"), "jy": ("Bx", "By", "Bz"), "jz": ("Bx", "By", "Bz"),
           "vh": ("u", "v"), "vmag": ("u", "v", "w")}
LABEL = {"j": r"$|\mathbf{j}|/j_0$", "jx": r"$j_x/j_0$", "jy": r"$j_y/j_0$", "jz": r"$j_z/j_0$",
         "Bx": r"$B_x/B_0$", "By": r"$B_y/B_0$", "Bz": r"$B_z/B_0$", "Bmag": r"$|\mathbf{B}|/B_0$",
         "u": r"$v_x/v_0$", "v": r"$v_y/v_0$", "w": r"$v_z/v_0$", "vh": r"$|v_h|/v_0$", "vmag": r"$|\mathbf{v}|/v_0$",
         "ρ": r"$\rho/\rho_0$", "rho": r"$\rho/\rho_0$", "log10ρ": r"$\log_{10}\rho/\rho_0$",
         "p": r"$p/P_0$", "T": r"$T/T_0$", "β": r"$\beta$", "beta": r"$\beta$", "ψ": r"$\psi$"}
ALIAS = {"rho": "ρ", "beta": "β", "psi": "ψ", "log10rho": "log10ρ"}


def needed(var):
    var = ALIAS.get(var, var)
    return DERIVED.get(var, (var,))


def evaluate(var, coords, F):
    var = ALIAS.get(var, var)
    if var in ("j", "jx", "jy", "jz"):
        jx, jy, jz = current(coords, F)
        return {"j": np.sqrt(jx**2 + jy**2 + jz**2), "jx": jx, "jy": jy, "jz": jz}[var]
    if var == "vh":
        return np.hypot(F["u"], F["v"])
    if var == "vmag":
        return np.sqrt(F["u"]**2 + F["v"]**2 + F["w"]**2)
    return F[var]


def inplane_lines(ax, h, v, Bh, Bv, density, color, lw=0.6):
    """In-plane field lines of (Bh, Bv)[nv, nh] (matplotlib streamplot on a uniform resampling of the lattice):
    traced, so they show the topology (O- and X-points) of a 3D cut, where no flux function exists."""
    nh = 500
    hu = np.linspace(h[0], h[-1], nh)
    vu = np.linspace(v[0], v[-1], max(2, int(nh * np.ptp(v) / np.ptp(h))))
    def resample(F):
        G = np.array([np.interp(hu, h, row) for row in F])
        return np.array([np.interp(vu, v, col) for col in G.T]).T
    ax.streamplot(hu, vu, resample(Bh), resample(Bv), density=density, color=color, linewidth=lw,
                  arrowstyle="-", integration_direction="both")


def nulls(h, v, Bh, Bv):
    """Nulls of the in-plane field (Bh, Bv)[nv, nh] on the lattice (h, v), bilinear in each cell:
    [(h0, v0, det J)], O-point where det J > 0, X-point where det J < 0."""
    def straddles(f):
        c = np.stack([f[:-1, :-1], f[:-1, 1:], f[1:, :-1], f[1:, 1:]])
        return (c.min(0) <= 0) & (c.max(0) >= 0)
    out = []
    for k, i in np.argwhere(straddles(Bh) & straddles(Bv)):
        a = Bh[k, i], Bh[k, i + 1], Bh[k + 1, i], Bh[k + 1, i + 1]
        b = Bv[k, i], Bv[k, i + 1], Bv[k + 1, i], Bv[k + 1, i + 1]
        bil = lambda c, s, t: c[0] * (1 - s) * (1 - t) + c[1] * s * (1 - t) + c[2] * (1 - s) * t + c[3] * s * t
        ds_ = lambda c, t: (c[1] - c[0]) * (1 - t) + (c[3] - c[2]) * t
        dt_ = lambda c, s: (c[2] - c[0]) * (1 - s) + (c[3] - c[1]) * s
        s = t = 0.5
        for _ in range(40):
            fa, fb = bil(a, s, t), bil(b, s, t)
            J = ds_(a, t) * dt_(b, s) - dt_(a, s) * ds_(b, t)
            if J == 0:
                break
            us, ut = (fa * dt_(b, s) - fb * dt_(a, s)) / J, (ds_(a, t) * fb - ds_(b, t) * fa) / J
            s, t = s - us, t - ut
            if abs(us) + abs(ut) < 1e-12:
                break
        scale = max(np.abs(a).max(), np.abs(b).max(), 1e-300)
        if not (-1e-6 <= s <= 1 + 1e-6 and -1e-6 <= t <= 1 + 1e-6) or \
           max(abs(bil(a, s, t)), abs(bil(b, s, t))) > 1e-8 * scale:
            continue
        dh, dv = h[i + 1] - h[i], v[k + 1] - v[k]
        detJ = (ds_(a, t) * dt_(b, s) - dt_(a, s) * ds_(b, t)) / (dh * dv)
        p = (h[i] + s * dh, v[k] + t * dv, detJ)
        if not any(abs(p[0] - q[0]) < 1e-6 * dh and abs(p[1] - q[1]) < 1e-6 * dv for q in out):
            out.append(p)
    return out


def unpair(found, dmin):
    """Drop O-X pairs closer than dmin (nearest first): the null pairs of grid-scale noise."""
    O = [n for n in found if n[2] > 0]
    X = [n for n in found if n[2] <= 0]
    pairs = sorted((np.hypot(o[0] - x[0], o[1] - x[1]), i, j) for i, o in enumerate(O) for j, x in enumerate(X))
    dO, dX = set(), set()
    for d, i, j in pairs:
        if d >= dmin:
            break
        if i not in dO and j not in dX:
            dO.add(i); dX.add(j)
    return [o for i, o in enumerate(O) if i not in dO] + [x for j, x in enumerate(X) if j not in dX]


def opoint(h, v, Bh, Bv, hwin, vlim, dmin):
    """The highest O-point with |h| <= hwin and vlim[0] <= v <= vlim[1] (the papers' coronal flux-rope axis),
    after unpair(dmin). Returns (h0, v0, all retained nulls in the window)."""
    win = [n for n in unpair(nulls(h, v, Bh, Bv), dmin) if abs(n[0]) <= hwin and vlim[0] <= n[1] <= vlim[1]]
    O = [n for n in win if n[2] > 0]
    h0, v0 = max(O, key=lambda n: n[1])[:2] if O else (np.nan, np.nan)
    return h0, v0, win


def runs_arg(items):
    """['SD=dir', 'dir2'] -> [('SD', 'dir'), ('dir2', 'dir2')]"""
    out = []
    for it in items:
        lab, _, d = it.partition("=") if "=" in it else (it, "", it)
        out.append((lab, d))
    return out


def style(label, i, mono):
    if mono:
        return "k", MONO_LS[i % len(MONO_LS)]
    if label in PAPER_STYLE:
        return PAPER_STYLE[label]
    return f"C{i % 10}", "-"


def savefig(fig, out, dpi):
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f" -> {out}")


# ---------------------------------------------------------------------------------------------------------
# ic: Fig. 1(a)
# ---------------------------------------------------------------------------------------------------------
def atmosphere(z, a=1.0, MF=1.0):
    """(T, p) of the cases' hydrostatic atmosphere (p = rho T/MF, rho(0) = T(0) = 1)."""
    al = np.log(TCOR) / 10.0
    p0, p10 = 1.0 / MF, np.exp(-10.0 * MF) / MF
    p20 = p10 * np.exp(-MF * (1.0 - np.exp(-10.0 * al)) / al)
    T = np.select([z <= 0, z < 10, z < 20], [1.0 - a * (GAMMA - 1) / GAMMA * z, np.ones_like(z),
                                              TCOR ** ((z - 10.0) / 10.0)], TCOR)
    with np.errstate(over="ignore", invalid="ignore"):
        p = np.select([z <= 0, z < 10, z < 20],
                      [p0 * np.maximum(T, 1e-300) ** (MF * GAMMA / (a * (GAMMA - 1))), p0 * np.exp(-MF * z),
                       p10 * np.exp(-MF * (1.0 - np.exp(-al * (z - 10.0))) / al)],
                      p20 * np.exp(-MF * (z - 20.0) / TCOR))
    return T, p


def tube_axis(z):
    """(B_x, B_y, p1) of the tube on x = 0."""
    zeta = z - ZT
    By = BT * np.exp(-zeta**2 / RT**2)
    p1 = 0.5 * BT**2 * np.exp(-2 * zeta**2 / RT**2) * ((0.5 * RT**2 - zeta**2) / RT**2 - 1.0)
    return zeta / RT * By, By, p1


def dipole_axis(z, Bd):
    return 2.0 * Bd / (z - ZD) ** 3            # B_x of the dipole on x = 0 (B_z = 0 there)


def flux_ratio(Bdabs, ztop, n=20001, defn="whole"):
    """Phi_dip/Phi_tube on x = y = 0 (initialize.jl): Phi_tube of the whole tube (B_t R) or its top half."""
    z = np.linspace(ZT, ztop, n)
    Bx = tube_axis(z)[0] + dipole_axis(z, -Bdabs)
    i0 = int(np.argmax(Bx > 0))
    neg = np.nonzero(Bx[i0:] < 0)[0]
    if len(neg) == 0:
        return 0.0
    i = i0 + neg[0]
    Pt = trapz(Bx[:i], z[:i]) if defn == "top" else BT * RT
    return abs(trapz(Bx[i:], z[i:])) / Pt


def solve_dipole(ratio, ztop, defn="whole"):
    if ratio <= 0:
        return 0.0
    lo, hi = 0.0, 1.0e5
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if flux_ratio(mid, ztop, defn=defn) < ratio else (lo, mid)
    return -0.5 * (lo + hi)


def cmd_ic(a):
    z = np.linspace(-30.0, a.ztop, 4001)
    T, p0 = atmosphere(z, a.cz_gradient, a.mf)
    Bxt, Byt, p1 = tube_axis(z)
    fig, ax = plt.subplots(figsize=(a.width, 0.75 * a.width))
    ax.semilogy(z, T, "k-", lw=1.5, label=r"$T/T_0$")
    ax.semilogy(z, p0, "k--", lw=1.5, label=r"$p/P_0$")
    for lab in ("SD", "MD", "WD"):
        Bd = solve_dipole(RATIO[lab], a.ztop, a.flux_def)
        beta = 2.0 * (p0 + p1) / np.maximum((Bxt + dipole_axis(z, Bd)) ** 2 + Byt**2, 1e-300)
        ax.semilogy(z, beta, color=PAPER_STYLE[lab][0], lw=1.2, label=rf"$\beta$ {lab} ($B_d = {Bd:.0f}$)")
        print(f"   {lab}: B_d = {Bd:.1f}, min beta on the axis above z = 0: {beta[z > 0].min():.3g}")
    for i, (lab, run) in enumerate(runs_arg(a.runs)):
        t, f = snapshots(run)[0]
        (x, y, zc), F = read_lattice(f, ("T", "p", "β"), box=((-1e-3, 1e-3), (-1e-3, 1e-3), (-1e9, 1e9)))
        col = PAPER_STYLE.get(lab, (f"C{i}", ""))[0]
        sl = (slice(None), 0, 0)
        ax.semilogy(zc, F["β"][sl], "o", ms=2.5, color=col, label=rf"{lab} nodes ($t={t:g}$): $\beta$, $T$, $p$")
        ax.semilogy(zc, F["T"][sl], "o", ms=2.5, color=col)
        ax.semilogy(zc, F["p"][sl], "o", ms=2.5, color=col)
    ax.axvline(0, color="0.7", lw=0.6)
    ax.set_xlim(z[0], z[-1])
    ax.set_xlabel(r"$z/L_0$")
    ax.legend(fontsize=8, ncol=2)
    ax.set_title(r"Initial atmosphere and $\beta(x=0, y=0, z)$")
    savefig(fig, os.path.join(a.outdir, f"leake_fig1a.{a.format}"), a.dpi)


# ---------------------------------------------------------------------------------------------------------
# slice: Figs. 2(d-f), 5(d-f), 6-7 (in-plane lines), 8(a-c)
# ---------------------------------------------------------------------------------------------------------
def parse_plane(text):
    ax, _, val = text.partition("=")
    if ax not in AX or not val:
        raise SystemExit(f"--plane {text}: expected x=VALUE, y=VALUE or z=VALUE")
    return ax, float(val)


def lim(text):
    return None if not text else tuple(float(s) for s in text.split(","))


def read_box(axis, value, margin, xlim, ylim, zlim):
    big = (-1e9, 1e9)
    box = [xlim or big, ylim or big, zlim or big]
    box = [(lo - margin, hi + margin) for lo, hi in box]
    box["xyz".index(axis)] = (value - margin, value + margin)
    return box


def cmd_slice(a):
    axis, value = parse_plane(a.plane)
    var = ALIAS.get(a.var, a.var)
    inplane = {"y": ("Bx", "Bz", "u", "w"), "z": ("Bx", "By", "u", "v"), "x": ("By", "Bz", "v", "w")}[axis]
    keep = set(needed(var))
    if a.fieldlines or a.nulls:
        keep |= set(inplane[:2])
    if a.vectors:
        keep |= set(inplane[2:])
    snaps = pick_times(snapshots(a.run), a.times)
    xlim, ylim, zlim = lim(a.xlim), lim(a.ylim), lim(a.zlim)
    hlim, vlim = {"y": (xlim, zlim), "z": (xlim, ylim), "x": (ylim, zlim)}[axis]
    hname, vname = {"y": ("x", "z"), "z": ("x", "y"), "x": ("y", "z")}[axis]
    panels = []
    for t, f in snaps:
        coords, F = read_lattice(f, sorted(keep), read_box(axis, value, a.margin, xlim, ylim, zlim))
        h, v, A = at_plane(coords, evaluate(var, coords, F), axis, value)
        extra = {k: at_plane(coords, F[k], axis, value)[2] for k in keep & set(inplane)}
        panels.append((t, h, v, A, extra))
        print(f"   {f}: t = {t:g}, {var} in [{np.nanmin(A):.3g}, {np.nanmax(A):.3g}]")

    def window(h, v, A):
        mh = np.ones(len(h), bool) if hlim is None else (h >= hlim[0]) & (h <= hlim[1])
        mv = np.ones(len(v), bool) if vlim is None else (v >= vlim[0]) & (v <= vlim[1])
        return mh, mv

    if a.clim:
        vmin, vmax = lim(a.clim)
    else:
        vals = np.concatenate([A[np.ix_(*window(h, v, A)[::-1])].ravel() for _, h, v, A, _ in panels])
        vmin, vmax = np.nanmin(vals), np.nanmax(vals)
        if a.log:
            vmin = max(vmin, vmax * 1e-4)
    if a.symmetric:
        vmax = max(abs(vmin), abs(vmax)); vmin = -vmax
    cmap = a.cmap or ("RdBu_r" if a.symmetric or var in ("Bx", "By", "Bz", "u", "v", "w", "jx", "jy", "jz")
                      else "inferno")
    if cmap == "RdBu_r" and not a.symmetric and not a.clim:
        vmax = max(abs(vmin), abs(vmax)); vmin = -vmax
    norm = matplotlib.colors.LogNorm(vmin, vmax) if a.log else matplotlib.colors.Normalize(vmin, vmax)

    n = len(panels)
    hspan = np.ptp(hlim) if hlim else np.ptp(panels[0][1])
    vspan = np.ptp(vlim) if vlim else np.ptp(panels[0][2])
    pw = a.width / n
    fig, axs = plt.subplots(1, n, figsize=(a.width, pw * vspan / hspan + 0.9), squeeze=False, sharey=True)
    for ax, (t, h, v, A, ex) in zip(axs[0], panels):
        mh, mv = window(h, v, A)
        hh, vv = h[mh], v[mv]
        S = A[np.ix_(mv, mh)]
        ax.set_rasterization_zorder(0)
        ax.pcolormesh(hh, vv, np.maximum(S, vmin) if a.log else S, shading="gouraud", cmap=cmap, norm=norm, zorder=-1)
        if a.fieldlines:
            inplane_lines(ax, hh, vv, ex[inplane[0]][np.ix_(mv, mh)], ex[inplane[1]][np.ix_(mv, mh)],
                          a.fieldlines, a.line_color)
        if a.nulls:
            for p, q, d in unpair(nulls(hh, vv, ex[inplane[0]][np.ix_(mv, mh)], ex[inplane[1]][np.ix_(mv, mh)]),
                                  a.pair_dist):
                ax.plot(p, q, "o" if d > 0 else "x", ms=7, mfc="none", mec="lime" if d > 0 else "cyan", mew=1.8)
        if a.vectors and not np.any(np.hypot(ex[inplane[2]], ex[inplane[3]]) > 0):
            print(f"   t = {t:g}: the in-plane velocity is zero, no arrows")
        elif a.vectors:
            Uh, Uv = ex[inplane[2]][np.ix_(mv, mh)], ex[inplane[3]][np.ix_(mv, mh)]
            ih = np.unique(np.searchsorted(hh, np.linspace(hh[0], hh[-1], a.nvec)).clip(0, len(hh) - 1))
            iv = np.unique(np.searchsorted(vv, np.linspace(vv[0], vv[-1], a.nvec)).clip(0, len(vv) - 1))
            ref = a.vec_ref or float(np.nanmax(np.hypot(Uh, Uv)))
            scale = a.vec_scale or ref / (0.8 * np.ptp(hh) / a.nvec)     # ref spans 0.8 of the arrow spacing
            Q = ax.quiver(hh[ih], vv[iv], Uh[np.ix_(iv, ih)], Uv[np.ix_(iv, ih)], color=a.vec_color,
                          angles="xy", scale_units="xy", scale=scale, width=0.003)
            ax.quiverkey(Q, 0.92, 1.025, ref, f"{ref:.2g} $v_0$", labelpos="W", coordinates="axes",
                         fontproperties={"size": 8})
        ax.set_aspect("equal")
        ax.set_xlim(hh[0], hh[-1]); ax.set_ylim(vv[0], vv[-1])
        ax.set_xlabel(rf"${hname}/L_0$")
        ax.set_title(rf"$t = {t:g}\,t_0$", fontsize=10, loc="left")
    axs[0][0].set_ylabel(rf"${vname}/L_0$")
    cb = fig.colorbar(plt.cm.ScalarMappable(norm, cmap), ax=axs[0].tolist(), shrink=0.85, pad=0.02)
    cb.set_label(LABEL.get(var, var) + f"  at ${axis} = {value:g}\\,L_0$")
    tag = a.name or f"slice_{var}_{axis}{value:g}"
    savefig(fig, os.path.join(a.outdir, f"{tag}.{a.format}"), a.dpi)


# ---------------------------------------------------------------------------------------------------------
# line: Fig. 8(d-f)
# ---------------------------------------------------------------------------------------------------------
def cmd_line(a):
    var = ALIAS.get(a.var, a.var)
    snaps = pick_times(snapshots(a.run), a.times)
    zlim = lim(a.zlim)
    m = a.margin
    box = ((a.x - m, a.x + m), (a.y - m, a.y + m), zlim if zlim else (-1e9, 1e9))
    n = len(snaps)
    fig, axs = plt.subplots(1, n, figsize=(a.width, a.height), squeeze=False, sharey=True)
    for ax, (t, f) in zip(axs[0], snaps):
        coords, F = read_lattice(f, needed(var), box)
        x, y, z = coords
        Q = evaluate(var, coords, F)
        prof = at_plane(coords, Q, "x", a.x)[2]                   # [nz, ny] on x = a.x; now y = a.y
        k = int(np.clip(np.searchsorted(y, a.y), 1, len(y) - 1)) if len(y) > 1 else 0
        if len(y) == 1 or abs(y[k] - a.y) < 1e-9:
            col = prof[:, k]
        else:
            w = (a.y - y[k - 1]) / (y[k] - y[k - 1])
            col = (1 - w) * prof[:, k - 1] + w * prof[:, k]
        ax.plot(z, col, "k-", lw=1.2)
        ax.axhline(0, color="0.7", lw=0.6)
        ax.set_xlabel(r"$z/L_0$")
        ax.set_title(rf"$t = {t:g}\,t_0$", fontsize=10)
        if zlim:
            ax.set_xlim(*zlim)
        print(f"   {f}: t = {t:g}, {var} on (x, y) = ({a.x:g}, {a.y:g}) in [{col.min():.3g}, {col.max():.3g}]")
    axs[0][0].set_ylabel(LABEL.get(var, var) + rf" on $x={a.x:g}, y={a.y:g}$")
    savefig(fig, os.path.join(a.outdir, f"{a.name or 'line_' + var}.{a.format}"), a.dpi)


# ---------------------------------------------------------------------------------------------------------
# height: Figs. 11-13
# ---------------------------------------------------------------------------------------------------------
def cmd_height(a):
    runs = runs_arg(a.runs)
    os.makedirs(a.outdir, exist_ok=True)
    series = []
    for lab, run in runs:
        rows = []
        for t, f in pick_times(snapshots(run), a.times):
            if a.trange and not (a.trange[0] <= t <= a.trange[1]):
                continue
            coords, F = read_lattice(f, ("Bx", "Bz"), read_box("y", 0.0, a.margin, (-a.xwin, a.xwin), None, None))
            x, z, Bx = at_plane(coords, F["Bx"], "y", 0.0)
            Bz = at_plane(coords, F["Bz"], "y", 0.0)[2]
            xo, zo, win = opoint(x, z, Bx, Bz, a.xwin, (a.zmin, a.zmax), a.pair_dist)
            nO, nX = sum(1 for n in win if n[2] > 0), sum(1 for n in win if n[2] <= 0)
            rows.append((t, zo, xo, nO, nX))
            print(f"   {lab}: t = {t:8.3f}  O-point z = {zo:8.3f}  x = {xo:7.3f}   ({nO} O, {nX} X in the window)")
        rows = np.array(rows)
        csv = os.path.join(a.outdir, f"leake_height_{lab}.csv")
        np.savetxt(csv, rows, header="t z_O x_O n_O n_X   (highest O-point of (Bx, Bz) in y = 0, |x| <= "
                   f"{a.xwin:g}, {a.zmin:g} <= z <= {a.zmax:g})", fmt="%.6g")
        print(f" -> {csv}")
        series.append((lab, rows))

    fig, ax = plt.subplots(figsize=(a.width, 0.8 * a.width))
    for i, (lab, r) in enumerate(series):
        c, ls = style(lab, i, a.mono)
        ax.plot(r[:, 0], r[:, 1], color=c, ls=ls, lw=1.5, label=lab)
    ax.axhline(a.damp_z, color="0.6", lw=0.8, ls=(0, (2, 2)))
    ax.set_xlabel(r"$t/t_0$"); ax.set_ylabel(r"O-point height $z/L_0$")
    ax.legend(fontsize=9)
    savefig(fig, os.path.join(a.outdir, f"leake_fig11_height.{a.format}"), a.dpi)

    fig, ax = plt.subplots(figsize=(a.width, 0.8 * a.width))
    for i, (lab, r) in enumerate(series):
        m = np.isfinite(r[:, 1])
        if m.sum() < 3:
            continue
        t, z = r[m, 0], r[m, 1]
        if a.smooth > 1:
            z = np.convolve(np.pad(z, a.smooth // 2, mode="edge"), np.ones(a.smooth) / a.smooth, "valid")[:len(t)]
        w = np.gradient(z, t)
        s = (t >= a.speed_range[0]) & (t <= a.speed_range[1])
        c, ls = style(lab, i, a.mono)
        ax.plot(t[s], w[s], color=c, ls=ls, lw=1.5, label=lab)
    ax.set_xlabel(r"$t/t_0$"); ax.set_ylabel(r"$\partial(z/L_0)/\partial(t/t_0)$")
    ax.legend(fontsize=9)
    savefig(fig, os.path.join(a.outdir, f"leake_fig12_speed.{a.format}"), a.dpi)


# ---------------------------------------------------------------------------------------------------------
# axflux: Fig. 14(b)
# ---------------------------------------------------------------------------------------------------------
def cmd_axflux(a):
    runs = runs_arg(a.runs)
    heights = [float(s) for s in a.heights.split(",")]
    os.makedirs(a.outdir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(a.width, 0.7 * a.width))
    for i, (lab, run) in enumerate(runs):
        rows = []
        for t, f in pick_times(snapshots(run), a.times):
            coords, F = read_lattice(f, ("By",), read_box("y", 0.0, a.margin, None, None, None))
            x, z, By = at_plane(coords, F["By"], "y", 0.0)
            if a.positive:
                By = np.maximum(By, 0.0)
            Iz = trapz(By, x, axis=1)                                   # int B_y dx, per height
            C = np.concatenate([[0.0], np.cumsum(0.5 * (Iz[1:] + Iz[:-1]) * np.diff(z))])
            tot = C[-1]
            rows.append([t] + [(tot - np.interp(h, z, C)) / tot for h in heights])
            print(f"   {lab}: t = {t:8.3f}  " + "  ".join(f"z>{h:g}: {r:.4f}" for h, r in zip(heights, rows[-1][1:])))
        rows = np.array(rows)
        np.savetxt(os.path.join(a.outdir, f"leake_axflux_{lab}.csv"), rows, fmt="%.6g",
                   header="t " + " ".join(f"ratio_z>{h:g}" for h in heights) + "   (int B_y dx dz in y = 0)")
        c, _ = style(lab, i, a.mono)
        for j, h in enumerate(heights):
            ax.plot(rows[:, 0], rows[:, 1 + j], color=c, ls=MONO_LS[j % len(MONO_LS)], lw=1.5,
                    label=f"{lab}, $z > {h:g}$")
    ax.set_xlabel(r"$t/t_0$"); ax.set_ylabel(r"$\Phi_y(z > z_h)/\Phi_y$")
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8, ncol=2)
    savefig(fig, os.path.join(a.outdir, f"leake_fig14b_axflux.{a.format}"), a.dpi)


# ---------------------------------------------------------------------------------------------------------
# render, fig2: 3D field lines (pyvista)
# ---------------------------------------------------------------------------------------------------------
# The papers' color scale (Tecplot's rainbow): blue - cyan - green (B_z = 0) - yellow - red. In Fig. 2 the
# colorbar ticks -0.5, -0.2, 0.0, 0.2, 0.5 are equally spaced: the map is piecewise linear in B_z between them.
RAINBOW = ("#0000ff", "#00ffff", "#00ff00", "#ffff00", "#ff0000")
FIG2_TICKS = (-0.5, -0.2, 0.0, 0.2, 0.5)


def colormap(name):
    if name in (None, "paper", "rainbow"):
        return matplotlib.colors.LinearSegmentedColormap.from_list("leake_rainbow", RAINBOW)
    return plt.get_cmap(name)


def bz_map(v, ticks):
    """B_z on the colorbar coordinate [-1, 1]: ticks equally spaced, linear in between, clamped outside."""
    return np.interp(v, ticks, np.linspace(-1.0, 1.0, len(ticks)))


def tick_text(ticks):
    tenths = all(abs(v * 10 - round(v * 10)) < 1e-9 for v in ticks)
    return [f"{v:.1f}" if tenths else f"{v:.3g}" for v in ticks]


# Camera ("render"): "yz" looks along -x with the tube (along y) left to right (Fig. 2); "xz" looks along +y
# across the tube (Figs. 1b, 10). dseeds = nx,ny,x0,x1,ymax: bottom-boundary seeds at x in [x0, x1], |y| <= ymax
# (a dipole line from (x, z = -30) peaks at z = -5 for x = 50, 30 for 80, 109 for 127). bz = color range of Bz
# at z = 0; bzticks = colorbar ticks, equally spaced on the bar (piecewise-linear map).
# fig2 only: sidebox / topbox = display boxes of the side (a-c) and top (d-f) views, sideplane = extent of the Bz
# plane in the side views; side = (azimuth, elevation, visible width in L0 at x = 0, z of the view center).
_EMERGE = dict(base="yz", azimuth=30.0, elevation=8.0, zoom=1.4, box=(-45, 45, -45, 45, -25, 55),
               dseeds=(3, 5, 70, 95, 40), bz=0.5)
PRESETS = {
    # Fig. 1(b): green dipole lines from the bottom, multicolored tube lines from the side boundaries,
    # transparent Bz plane (the paper's +-7.5 G)
    "ic":        dict(dipole="green", tube="multi", axis=False, surface=0.6, iso=None, base="xz", azimuth=35.0,
                      elevation=12.0, zoom=1.3, box=(-100, 100, -128, 128, -30, 200), dseeds=(5, 3, 30, 120, 60),
                      bz=0.006),
    # Fig. 2(a-c) (--view top, or the fig2 command: d-f): one gray dipole line from the bottom boundary at each
    # y = -20, -10, 0, 10, 20, multicolored tube lines and the black tube-axis line from the side boundaries,
    # Bz at z = 0 on the paper's rainbow with ticks -0.5, -0.2, 0, 0.2, 0.5
    "emergence": dict(dipole="#707070", tube="multi", axis=True, surface=1.0, iso=None,
                      **{**_EMERGE, "dseeds": (1, 5, 70, 70, 20), "bzticks": FIG2_TICKS,
                         "sidebox": (-100, 100, -50, 50, -30, 150), "sideplane": (-40, 40, -50, 50, -1, 1),
                         "topbox": (-40, 40, -40, 40, -30, 200), "side": (15.0, 6.0, 90.0, 0.0)}),
    # Fig. 3: gray dipole lines, multicolored tube lines, red isosurface |j| = 0.004 above z = 50
    "arcade":    dict(dipole="gray", tube="multi", axis=False, surface=None, iso=(0.004, 50.0),
                      **{**_EMERGE, "box": (-100, 100, -100, 100, -25, 150), "dseeds": (4, 5, 60, 120, 60),
                         "elevation": 15.0}),
    # Fig. 4: gray dipole lines, side-boundary circles of radius 0.5, 1, 2 (orange, blue, black)
    "dipoles":   dict(dipole="gray", tube="rings", axis=False, surface=None, iso=None,
                      **{**_EMERGE, "box": (-80, 80, -128, 128, -25, 100), "dseeds": (4, 5, 60, 120, 60),
                         "elevation": 15.0}),
    # Fig. 5(a-c): Bz at z = 0, black (purple) line from the tube axis on y = max (min)
    "rotation":  dict(dipole=None, tube=None, axis="pair", surface=1.0, iso=None, **{**_EMERGE, "bz": 0.2}),
    # Fig. 6: as rotation plus the orange line through the O-point in y = 0
    "fluxrope":  dict(dipole=None, tube=None, axis="pair", surface=1.0, iso=None, opoint="line",
                      **{**_EMERGE, "box": (-60, 60, -60, 60, -25, 80), "elevation": 15.0, "bz": 0.2}),
    # Fig. 7: as rotation plus the red isosurface |j| = 0.2 above z = 10
    "current":   dict(dipole=None, tube=None, axis="pair", surface=1.0, iso=(0.2, 10.0),
                      **{**_EMERGE, "box": (-60, 60, -60, 60, -25, 80), "elevation": 15.0, "bz": 0.2}),
    # Fig. 10: gray dipole lines, orange line through the O-point, blue lines within 10 L0 of it, Bz at z = 0 (+-0.2)
    "eruption":  dict(dipole="gray", tube=None, axis=False, surface=1.0, iso=None, opoint="disk", base="xz",
                      azimuth=35.0, elevation=12.0, zoom=1.3, box=(-100, 100, -100, 100, -30, 200),
                      dseeds=(4, 3, 40, 120, 60), bz=0.2),
}


def need_pyvista():
    try:
        import pyvista as pv
    except ImportError:
        raise SystemExit("render and fig2 need pyvista: pip install pyvista")
    pv.OFF_SCREEN = True
    return pv


def preset_params(a):
    """The preset of a, with the command-line overrides."""
    P = dict(PRESETS[a.preset])
    for k in ("dipole", "tube", "axis", "surface"):
        v = getattr(a, k)
        if v is not None:
            P[k] = {"none": None, "off": None, "true": True, "false": False}.get(v.lower(), v)
    if a.iso:
        lv, zz = lim(a.iso)
        P["iso"] = (lv, zz)
    if a.no_iso:
        P["iso"] = None
    for k, d in (("azimuth", "azimuth"), ("elevation", "elevation"), ("zoom", "zoom"), ("bz_clim", "bz")):
        if getattr(a, k, None) is not None:
            P[d] = getattr(a, k)
    return P


def bz_ticks(P, a):
    """Colorbar ticks of Bz: --bz-ticks, else the preset's (unless --bz-clim), else None (+-bz, linear)."""
    if a.bz_ticks:
        return lim(a.bz_ticks)
    if a.bz_clim is None and P.get("bzticks"):
        return tuple(P["bzticks"])
    return None


def scene(pv, f, P, a, dseeds):
    """The grid (B, Bz and |j| if needed) of snapshot f and its field lines over the whole domain:
    (grid, [(lines, color or 'multi', name, number of seeds)]). The same seeds at every time."""
    (x, y, z), F = read_lattice(f, ["Bx", "By", "Bz"])
    grid = pv.RectilinearGrid(x, y, z)
    grid.point_data["B"] = np.stack([F["Bx"].ravel(), F["By"].ravel(), F["Bz"].ravel()], axis=1)
    grid.point_data["Bz"] = F["Bz"].ravel()
    if P.get("iso"):
        jx, jy, jz = current((x, y, z), F)
        J = np.sqrt(jx**2 + jy**2 + jz**2)
        J[z < P["iso"][1], :, :] = 0.0
        grid.point_data["J"] = J.ravel()
        del jx, jy, jz, J
    eps = 1e-3 * (z[-1] - z[0])
    ymax, ymin = y[-1] - eps, y[0] + eps
    seeds = []                                                 # (seeds, color or 'multi', name)
    if P.get("dipole"):
        nx_, ny_, x0_, x1_, ys = dseeds
        S = np.array([(xx, yy, z[0] + eps) for yy in np.linspace(-ys, ys, int(ny_))
                      for xx in np.linspace(x0_, x1_, int(nx_))])
        seeds.append((S, P["dipole"], "dipole"))
    sides = [ymax, ymin] if a.sides == "both" else ([ymax] if a.sides == "max" else [ymin])
    if P.get("tube") == "multi":                               # colored by the angle of the seed around the axis
        S = np.array([(r * np.cos(th), ys, ZT + r * np.sin(th)) for ys in sides for r in lim(a.tube_radii)
                      for th in np.linspace(0, 2 * np.pi, a.tube_n, endpoint=False)])
        seeds.append((S, "multi", "tube"))
    elif P.get("tube") == "rings":
        for r, col in zip((0.5, 1.0, 2.0), ("orange", "royalblue", "k")):
            S = np.array([(r * np.cos(th), ys, ZT + r * np.sin(th)) for ys in sides
                          for th in np.linspace(0, 2 * np.pi, a.tube_n, endpoint=False)])
            seeds.append((S, col, f"ring{r:g}"))
    if P.get("axis") is True:
        seeds.append((np.array([(0.0, s, ZT) for s in sides]), "k", "axis"))
    elif P.get("axis") == "pair":
        seeds.append((np.array([(0.0, ymax, ZT)]), "k", "axis+"))
        seeds.append((np.array([(0.0, ymin, ZT)]), "purple", "axis-"))
    if P.get("opoint"):
        iy0 = int(np.argmin(np.abs(y)))
        xo, zo, _ = opoint(x, z, F["Bx"][:, iy0, :], F["Bz"][:, iy0, :], a.xwin, (a.zmin, a.zmax), a.pair_dist)
        if np.isfinite(zo):
            print(f"   O-point in y = {y[iy0]:g}: (x, z) = ({xo:.3f}, {zo:.3f})")
            seeds.append((np.array([(xo, y[iy0], zo)]), "orange", "opoint"))
            if P["opoint"] == "disk":
                S = np.array([(xo + r * np.cos(th), y[iy0], zo + r * np.sin(th)) for r in (2.5, 5.0, 7.5, 10.0)
                              for th in np.linspace(0, 2 * np.pi, 8, endpoint=False)])
                seeds.append((S, "royalblue", "disk"))
        else:
            print("   no O-point in the window (--xwin, --zmin, --zmax)")
    traced = []
    for S, col, name in seeds:
        sl = grid.streamlines_from_source(pv.PolyData(S), vectors="B", integration_direction="both",
                                          initial_step_length=0.5, max_step_length=2.0, max_steps=a.max_steps,
                                          max_length=a.max_length, terminal_speed=1e-12)
        if sl.n_points == 0:
            print(f"   {name}: no line from {len(S)} seeds")
            continue
        if col == "multi":
            for data in (sl.point_data, sl.cell_data):
                if "SeedIds" in data:
                    data["hue"] = (np.asarray(data["SeedIds"]) % a.tube_n) / max(a.tube_n - 1, 1)
        traced.append((sl, col, name, len(S)))
    return grid, traced


def draw(pv, pl, grid, traced, P, a, box, cmap, ticks, scalar_bar=True, plane_box=None):
    """Add the field lines, the Bz plane (clipped to plane_box, default box) and the |j| isosurface, clipped to box,
    to the plotter pl. Returns the Bz colorbar ticks used."""
    def clip(m, b=None):
        b = b or box
        if not b:
            return m
        m = m.clip_box(list(b), invert=False)
        if m.n_points and not isinstance(m, pv.PolyData):      # older pyvista: clip_box gives a grid
            m = m.extract_surface()
        return m
    fg = "k" if getattr(a, "background", "white") == "white" else "w"
    for sl, col, name, ns in traced:
        s = clip(sl)
        if s.n_points == 0:
            continue
        geo = s.tube(radius=a.tube_width) if a.tube_width > 0 else s
        if col == "multi":
            pl.add_mesh(geo, scalars="hue", cmap=cmap, clim=(0.0, 1.0), show_scalar_bar=False, smooth_shading=True)
        else:
            pl.add_mesh(geo, color=col, smooth_shading=True)
        print(f"   {name}: {ns} seeds, {s.n_cells} line pieces in the box")
    if P.get("surface"):
        sf = clip(grid.slice(normal="z", origin=(0.0, 0.0, a.surface_z)), plane_box)
        if ticks is None:
            c = P["bz"] or float(np.abs(sf["Bz"]).max()) or 1.0
            ticks = (-c, 0.0, c)
        sf["Bz_bar"] = bz_map(np.asarray(sf["Bz"]), ticks)
        pos = np.linspace(-1.0, 1.0, len(ticks))
        pl.add_mesh(sf, scalars="Bz_bar", cmap=cmap, clim=(-1.0, 1.0), opacity=float(P["surface"]), lighting=False,
                    show_scalar_bar=scalar_bar,
                    annotations={float(m): s for m, s in zip(pos, tick_text(ticks))},
                    scalar_bar_args=dict(title=f"Bz/B0 at z = {a.surface_z:g}", vertical=True, n_labels=0, color=fg))
    if P.get("iso"):
        iso = grid.contour([P["iso"][0]], scalars="J")
        if box and iso.n_points:
            iso = clip(iso)
        if iso.n_points:
            pl.add_mesh(iso, color="red", opacity=0.6, smooth_shading=True)
            print(f"   |j| = {P['iso'][0]:g} above z = {P['iso'][1]:g}: {iso.n_cells} triangles")
        else:
            print(f"   |j| never reaches {P['iso'][0]:g} above z = {P['iso'][1]:g}")
    if a.outline and box:
        pl.add_mesh(pv.Box(box).outline(), color="gray")
    return ticks


def cmd_render(a):
    pv = need_pyvista()
    P = preset_params(a)
    box = lim(a.box) if a.box else P["box"]
    dseeds = lim(a.dipole_seeds) if a.dipole_seeds else P["dseeds"]
    cmap, ticks = colormap(a.bz_cmap), bz_ticks(P, a)
    os.makedirs(a.outdir, exist_ok=True)
    for t, f in pick_times(snapshots(a.run), a.times):
        print(f"   {f}: t = {t:g}")
        grid, traced = scene(pv, f, P, a, dseeds)
        pl = pv.Plotter(off_screen=True, window_size=(a.size[0], a.size[1]))
        pl.set_background(a.background)
        draw(pv, pl, grid, traced, P, a, box, cmap, ticks)
        if a.view == "top":
            pl.view_xy()
            P["zoom"] = a.zoom if a.zoom is not None else 1.0
        else:
            pl.view_yz(negative=True) if P["base"] == "yz" else pl.view_xz()
            pl.camera.azimuth = P["azimuth"]
            pl.camera.elevation = P["elevation"]
        pl.camera.zoom(P["zoom"])
        if a.title:
            pl.add_text(f"t = {t:g} t0", font_size=12, color="k" if a.background == "white" else "w")
        out = os.path.join(a.outdir, f"{a.name or 'render_' + a.preset}_t{t:g}_{a.view}.png")
        pl.screenshot(out)
        pl.close()
        print(f" -> {out}")


def side_camera(pl, azimuth, elevation, width, zc, aspect, flip=False, view_angle=12.0):
    """Perspective camera on the -x side (+x with flip) looking at (0, 0, zc), the tube (along y) across the
    image, turned by azimuth about z and raised by elevation (degrees); width L0 visible across at x = 0."""
    az, el = np.radians(azimuth), np.radians(elevation)
    d = width / (2.0 * np.tan(np.radians(view_angle) / 2.0) * aspect)
    u = np.array([(1.0 if flip else -1.0) * np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])
    fp = np.array([0.0, 0.0, zc])
    pl.camera_position = [tuple(fp + d * u), tuple(fp), (0.0, 0.0, 1.0)]
    pl.camera.view_angle = view_angle
    pl.reset_camera_clipping_range()


def top_camera(pl, box, zabove):
    """Orthographic view from above framing exactly [x0, x1] x [y0, y1] (the window must have their aspect)."""
    x0, x1, y0, y1 = box[:4]
    c = (0.5 * (x0 + x1), 0.5 * (y0 + y1))
    pl.camera_position = [(c[0], c[1], zabove), (c[0], c[1], 0.0), (0.0, 1.0, 0.0)]
    pl.enable_parallel_projection()
    pl.camera.parallel_scale = 0.5 * (y1 - y0)
    pl.reset_camera_clipping_range()


def cmd_fig2(a):
    """Fig. 2 of the paper: for each time, the field lines seen from the side (a-c) and from above (d-f) with Bz
    at z = 0, each view its own figure; the same seeds at every time. --grid also puts them in one 2 x N figure."""
    pv = need_pyvista()
    P = preset_params(a)
    dseeds = lim(a.dipole_seeds) if a.dipole_seeds else P["dseeds"]
    sbox = lim(a.box) if a.box else P.get("sidebox", P["box"])
    splane = lim(a.side_plane) if a.side_plane else P.get("sideplane", sbox)
    tbox = lim(a.top_box) if a.top_box else P.get("topbox", P["box"])
    saz, sel, swid, szc = P.get("side", (P["azimuth"], P["elevation"], float(np.ptp(sbox[2:4])), 0.0))
    saz = saz if a.side_azimuth is None else a.side_azimuth
    sel = sel if a.side_elevation is None else a.side_elevation
    swid = swid if a.side_width is None else a.side_width
    szc = szc if a.side_zc is None else a.side_zc
    cmap, ticks0 = colormap(a.bz_cmap), bz_ticks(P, a)
    W, H = a.size
    x0, x1, y0, y1 = tbox[:4]
    Wt = int(round(H * (x1 - x0) / (y1 - y0)))
    views = ("side", "top") if a.view == "both" else (a.view,)
    os.makedirs(a.outdir, exist_ok=True)
    tag = a.name or "leake_fig2"
    letters = "abcdefghijklmnopqrstuvwxyz"

    def side_panel(ax, t, img, ticks, letter):
        ax.imshow(img)
        ax.set_xticks([]); ax.set_yticks([])
        ax.text(0.03, 0.975, (f"({letter})  " if letter else "") + rf"$t = {t:g}\,t_0$", transform=ax.transAxes,
                va="top", ha="left", fontsize=10)
        if P.get("surface"):
            cax = ax.inset_axes([0.25, 0.80, 0.55, 0.035])
            cb = ax.figure.colorbar(matplotlib.cm.ScalarMappable(matplotlib.colors.Normalize(-1.0, 1.0), cmap),
                                    cax=cax, orientation="horizontal", ticks=np.linspace(-1.0, 1.0, len(ticks)))
            cb.ax.set_xticklabels(tick_text(ticks))
            cb.ax.xaxis.set_ticks_position("top")
            cb.ax.tick_params(labelsize=7, length=2, pad=1)

    def top_panel(ax, t, img, letter, ylabel=True):
        ax.imshow(img, extent=(x0, x1, y0, y1), origin="upper")
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_locator(matplotlib.ticker.MultipleLocator(10))
            axis.set_minor_locator(matplotlib.ticker.MultipleLocator(2))
        ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=7)
        ax.set_xlabel(r"$x/L_0$", fontsize=9)
        if ylabel:
            ax.set_ylabel(r"$y/L_0$", fontsize=9)
        else:
            ax.tick_params(labelleft=False)
        ax.text(0.03, 0.975, (f"({letter})  " if letter else "") + rf"$t = {t:g}\,t_0$", transform=ax.transAxes,
                va="top", ha="left", fontsize=10)

    panels = []
    for t, f in pick_times(snapshots(a.run), a.times):
        print(f"   {f}: t = {t:g}")
        grid, traced = scene(pv, f, P, a, dseeds)
        side = top = None
        ticks = ticks0
        if "side" in views:
            pl = pv.Plotter(off_screen=True, window_size=(W, H))
            pl.set_background("white")
            ticks = draw(pv, pl, grid, traced, P, a, sbox, cmap, ticks0, scalar_bar=False, plane_box=splane)
            side_camera(pl, saz, sel, swid, szc, W / H, a.side_flip, a.view_angle)
            side = pl.screenshot(return_img=True)
            pl.close()
            fig, ax = plt.subplots(figsize=(a.width, a.width * H / W))
            side_panel(ax, t, side, ticks, None)
            savefig(fig, os.path.join(a.outdir, f"{tag}_side_t{t:g}.{a.format}"), a.dpi)
        if "top" in views:
            pl = pv.Plotter(off_screen=True, window_size=(Wt, H))
            pl.set_background("white")
            ticks = draw(pv, pl, grid, traced, P, a, tbox, cmap, ticks0, scalar_bar=False)
            top_camera(pl, tbox, float(grid.bounds[5]) + 50.0)
            top = pl.screenshot(return_img=True)
            pl.close()
            fig, ax = plt.subplots(figsize=(a.width, a.width * (y1 - y0) / (x1 - x0)))
            top_panel(ax, t, top, None)
            if P.get("surface"):                                 # colorbar under the axes
                cb = fig.colorbar(matplotlib.cm.ScalarMappable(matplotlib.colors.Normalize(-1.0, 1.0), cmap),
                                  ax=ax, orientation="vertical", shrink=0.8, pad=0.02,
                                  ticks=np.linspace(-1.0, 1.0, len(ticks)))
                cb.ax.set_yticklabels(tick_text(ticks))
                cb.ax.tick_params(labelsize=7)
                cb.set_label(rf"$B_z/B_0$ at $z = {a.surface_z:g}$", fontsize=9)
            savefig(fig, os.path.join(a.outdir, f"{tag}_top_t{t:g}.{a.format}"), a.dpi)
        panels.append((t, side, top, ticks))

    if not a.grid:
        return
    n = len(panels)
    rows = [v for v in ("side", "top") if v in views]
    ratio = {"side": H / W, "top": (y1 - y0) / (x1 - x0)}
    cw = a.width                                                   # each column as wide as a single view
    fig = plt.figure(figsize=(cw * n, cw * sum(ratio[r] for r in rows) * 1.06 + 0.4))
    gs = fig.add_gridspec(len(rows), n, height_ratios=[ratio[r] for r in rows], wspace=0.06, hspace=0.03)
    for i, (t, side, top, ticks) in enumerate(panels):
        for k, r in enumerate(rows):
            ax = fig.add_subplot(gs[k, i])
            letter = letters[(k * n + i) % 26]
            side_panel(ax, t, side, ticks, letter) if r == "side" else top_panel(ax, t, top, letter, i == 0)
    savefig(fig, os.path.join(a.outdir, f"{tag}.{a.format}"), a.dpi)


# ---------------------------------------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(s, fmt=True):
        s.add_argument("--outdir", default="figs_leake2014")
        s.add_argument("--dpi", type=int, default=200)
        s.add_argument("--width", type=float, default=6.0, help="figure width in inches")
        if fmt:
            s.add_argument("--format", default="png", help="png, pdf, svg, ...")

    s = sub.add_parser("ic", help="Fig. 1(a)")
    s.add_argument("runs", nargs="*", help="LABEL=RUN: overlay the nodal T, p, beta on x = y = 0 at t = 0")
    s.add_argument("--ztop", type=float, default=ZTOP)
    s.add_argument("--cz-gradient", type=float, default=1.0, help=":lk_cz_gradient of the deck")
    s.add_argument("--mf", type=float, default=1.0, help=":lk_mf of the deck (p = rho T/m_f)")
    s.add_argument("--flux-def", choices=("whole", "top"), default="whole", help=":lk_flux_def of the deck")
    common(s)
    s.set_defaults(func=cmd_ic)

    s = sub.add_parser("slice", help="a field on a plane (Figs. 2d-f, 5d-f, 6-7, 8a-c)")
    s.add_argument("run")
    s.add_argument("--plane", default="y=0", help="y=0 (x-z cut), z=0 (surface), x=VALUE")
    s.add_argument("--var", default="j", help="any output field (rho, u, v, w, p, Bx, By, Bz, T, beta, Bmag, ...), "
                   "or j, jx, jy, jz (curl B), vh (horizontal speed), vmag")
    s.add_argument("--times", default=None, help="comma-separated times (nearest snapshots); default all")
    s.add_argument("--xlim"); s.add_argument("--ylim"); s.add_argument("--zlim")
    s.add_argument("--clim", help="vmin,vmax")
    s.add_argument("--log", action="store_true", help="log color scale")
    s.add_argument("--symmetric", action="store_true")
    s.add_argument("--cmap")
    s.add_argument("--fieldlines", type=float, default=0.0,
                   help="draw the in-plane field lines, with this streamplot density (e.g. 1.5)")
    s.add_argument("--line-color", default="w")
    s.add_argument("--nulls", action="store_true", help="mark O-points (circles) and X-points (crosses)")
    s.add_argument("--pair-dist", type=float, default=3.0, help="drop O-X pairs closer than this (grid noise)")
    s.add_argument("--vectors", action="store_true", help="in-plane velocity arrows")
    s.add_argument("--nvec", type=int, default=20)
    s.add_argument("--vec-scale", type=float, default=None, help="v0 per L0 of arrow length (default: --vec-ref "
                   "spans 0.8 of the arrow spacing)")

    s.add_argument("--vec-ref", type=float, default=None)
    s.add_argument("--vec-color", default="k")
    s.add_argument("--margin", type=float, default=4.0, help="read this far either side of the plane (for derivatives)")
    s.add_argument("--name", help="output file name (no extension)")
    common(s)
    s.set_defaults(func=cmd_slice, width=12.0)

    s = sub.add_parser("line", help="a field along a vertical line (Fig. 8d-f)")
    s.add_argument("run")
    s.add_argument("--var", default="w")
    s.add_argument("--x", type=float, default=0.0); s.add_argument("--y", type=float, default=0.0)
    s.add_argument("--times", default=None)
    s.add_argument("--zlim")
    s.add_argument("--height", type=float, default=3.0)
    s.add_argument("--margin", type=float, default=4.0)
    s.add_argument("--name")
    common(s)
    s.set_defaults(func=cmd_line, width=10.0)

    s = sub.add_parser("height", help="O-point height in y = 0 and its speed (Figs. 11-13)")
    s.add_argument("runs", nargs="+", help="LABEL=RUN (SD, MD, WD, ND get the paper's styles)")
    s.add_argument("--times", default=None)
    s.add_argument("--trange", type=lambda t: tuple(float(v) for v in t.split(",")))
    s.add_argument("--xwin", type=float, default=30.0, help="O-points with |x| <= xwin only")
    s.add_argument("--zmin", type=float, default=-30.0); s.add_argument("--zmax", type=float, default=180.0)
    s.add_argument("--damp-z", type=float, default=180.0, help="marks the top damping layer")
    s.add_argument("--speed-range", type=lambda t: tuple(float(v) for v in t.split(",")), default=(120.0, 165.0))
    s.add_argument("--smooth", type=int, default=1, help="running mean over this many snapshots before d/dt")
    s.add_argument("--pair-dist", type=float, default=3.0, help="drop O-X pairs closer than this (grid noise)")
    s.add_argument("--mono", action="store_true", help="black lines, styles in order (Fig. 13)")
    s.add_argument("--margin", type=float, default=4.0)
    common(s)
    s.set_defaults(func=cmd_height)

    s = sub.add_parser("axflux", help="axial-flux ratio above given heights in y = 0 (Fig. 14b)")
    s.add_argument("runs", nargs="+")
    s.add_argument("--heights", default="0,10")
    s.add_argument("--times", default=None)
    s.add_argument("--positive", action="store_true", help="count B_y > 0 only")
    s.add_argument("--mono", action="store_true")
    s.add_argument("--margin", type=float, default=4.0)
    common(s)
    s.set_defaults(func=cmd_axflux)

    def scene_args(s):
        """Options shared by render and fig2: what is drawn (seeds, plane, isosurface) and how."""
        s.add_argument("run")
        s.add_argument("--preset", choices=sorted(PRESETS), default="emergence")
        s.add_argument("--dipole", help="override: color of the bottom-boundary lines, or none")
        s.add_argument("--dipole-seeds", help="nx,ny,x0,x1,ymax: bottom seeds at x in [x0, x1], |y| <= ymax "
                       "(default: the preset's; emergence: 1,5,70,70,20 = one line at each y = -20..20 step 10)")
        s.add_argument("--tube", help="override: multi, rings or none")
        s.add_argument("--tube-radii", default="1,2", help="seed circles around the tube axis on the side boundaries")
        s.add_argument("--tube-n", type=int, default=8, help="seeds per circle (multi: colored by their angle)")
        s.add_argument("--sides", choices=("both", "max", "min"), default="both")
        s.add_argument("--axis", help="override: true, pair or none")
        s.add_argument("--surface", help="override: opacity of the Bz plane, or none")
        s.add_argument("--surface-z", type=float, default=0.0)
        s.add_argument("--bz-clim", type=float, help="color range +-c of Bz at z = 0, linear (default: the preset's)")
        s.add_argument("--bz-ticks", help="colorbar ticks of Bz, equally spaced on the bar, the map piecewise "
                       "linear between them (emergence: -0.5,-0.2,0,0.2,0.5 as in Fig. 2)")
        s.add_argument("--bz-cmap", default="paper",
                       help="paper: the papers' rainbow blue-cyan-green(0)-yellow-red (also the tube lines); "
                            "or any matplotlib name (jet, RdBu_r, ...)")
        s.add_argument("--iso", help="override: LEVEL,ZMIN of the |j| isosurface")
        s.add_argument("--no-iso", action="store_true")
        s.add_argument("--xwin", type=float, default=30.0); s.add_argument("--zmin", type=float, default=-30.0)
        s.add_argument("--zmax", type=float, default=180.0)
        s.add_argument("--pair-dist", type=float, default=3.0)
        s.add_argument("--tube-width", type=float, default=0.3, help="field-line tube radius in L0 (0: lines)")
        s.add_argument("--max-steps", type=int, default=20000); s.add_argument("--max-length", type=float, default=3000.0)
        s.add_argument("--outline", action="store_true", help="outline the display box")
        s.add_argument("--name")

    s = sub.add_parser("render", help="3D field lines (pyvista): Figs. 1b, 2-7, 10, one image per time")
    scene_args(s)
    s.add_argument("--times", default=None)
    s.add_argument("--view", choices=("side", "top"), default="side")
    s.add_argument("--azimuth", type=float); s.add_argument("--elevation", type=float)
    s.add_argument("--zoom", type=float)
    s.add_argument("--box", help="display box x0,x1,y0,y1,z0,z1 (default: the preset's)")
    s.add_argument("--size", type=lambda t: tuple(int(v) for v in t.split(",")), default=(1600, 1100))
    s.add_argument("--background", default="white")
    s.add_argument("--title", action="store_true")
    s.add_argument("--outdir", default="figs_leake2014")
    s.set_defaults(func=cmd_render)

    s = sub.add_parser("fig2", help="Fig. 2: side view (a-c) and top view (d-f) at each time, one file each")
    scene_args(s)
    s.add_argument("--times", default="5,40,60", help="comma-separated times (nearest snapshots), or all")
    s.add_argument("--view", choices=("both", "side", "top"), default="both")
    s.add_argument("--grid", action="store_true", help="also all views in one figure, a column per time (Fig. 2)")
    s.add_argument("--box", help="display box of the side views x0,x1,y0,y1,z0,z1 (default: the preset's sidebox)")
    s.add_argument("--top-box", help="display box of the top views; its x, y range is the axes (default: "
                   "-40,40,-40,40,-30,200)")
    s.add_argument("--side-azimuth", type=float, help="turn of the side camera about z, degrees (default 15)")
    s.add_argument("--side-elevation", type=float, help="height of the side camera above z = zc, degrees (default 6)")
    s.add_argument("--side-width", type=float, help="L0 visible across the side view at x = 0 (default 90)")
    s.add_argument("--side-zc", type=float, help="z of the side-view center (default 0)")
    s.add_argument("--side-flip", action="store_true", help="side camera on the +x side (mirrors y)")
    s.add_argument("--side-plane", help="extent x0,x1,y0,y1,z0,z1 of the Bz plane in the side views "
                   "(default: -40,40,-50,50,-1,1)")
    s.add_argument("--view-angle", type=float, default=12.0, help="side camera view angle, degrees (small: "
                   "nearly orthographic, as in the paper)")
    s.add_argument("--size", type=lambda t: tuple(int(v) for v in t.split(",")), default=(1000, 1000),
                   help="pixels of each side view W,H; the top views are H high")
    common(s)
    s.set_defaults(func=cmd_fig2, width=5.0)

    # let comma lists start with a minus sign: --xlim -30,30 -> --xlim=-30,30
    argv = sys.argv[1:]
    for i in range(len(argv) - 1, 0, -1):
        if argv[i - 1].startswith("--") and "=" not in argv[i - 1] and argv[i][:1] == "-" and \
           len(argv[i]) > 1 and (argv[i][1].isdigit() or argv[i][1] == "."):
            argv[i - 1:i + 1] = [argv[i - 1] + "=" + argv[i]]
    a = p.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
