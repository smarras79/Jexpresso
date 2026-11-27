#!/usr/bin/env python3
"""Figures of Leake, Linton & Antiochos (2010, ApJ 722:550) and of its Erratum (2011, ApJ 741:125)
from the output of problems/MHD/fluxEmergenceLeake2010.

RUN is an output directory (output/MHD/fluxEmergenceLeake2010/output-<date>, the one holding iter_N.pvtu);
a restarted run is DIR1,DIR2 (DIR1@920,DIR2: DIR1 only up to 920 s). A limit that starts with a minus sign is
written with "=": --ylim=-3,19.5. Positions are plotted in Mm and times in s (code units: L0 = 150 km,
t0 = 23.1946 s), the density in g cm^-3 and B_z in G, as in the paper.

  python3 tools/plot_leake2010.py density RUN                 # Erratum Fig. 3 (Fig. 10 of the paper): log rho, field lines,
                                                              #   B_x arrows at t = 920, 1380, 1610, 1955, 2185, 2530 s
  python3 tools/plot_leake2010.py density RUN --times 920,1610 --label "Low Tension"
  python3 tools/plot_leake2010.py shear RUN                   # Erratum Figs. 1 and 8: B_z (G) and field lines
  python3 tools/plot_leake2010.py dsgs RUN                    # the DynSGS coefficient of the momentum equations at those
                                                              #   times (--which induction | energy: of those equations)
  python3 tools/plot_leake2010.py movie RUN                   # one image per snapshot and an animated GIF of them
                                                              #   (--field density | momentum | induction | energy)
  python3 tools/plot_leake2010.py dome RUN [RUN2 ...]         # top of the dense dome on x = 0 vs time, against Erratum Fig. 3
  python3 tools/plot_leake2010.py axis RUN [RUN2 ...]         # Erratum Figs. 2(a), 7(a): height of the tube axis vs time
  python3 tools/plot_leake2010.py eshear RUN [RUN2 ...]       # Erratum Figs. 2(b,c), 5: E_shear above 1.2 and 5 Mm (Eq. 20)
  python3 tools/plot_leake2010.py profiles RUN                # Fig. 3(a,b): T, p, p_mag and beta at x = 0, t = 0
  python3 tools/plot_leake2010.py forces RUN                  # Erratum Fig. 4: vertical forces on x = 0 at 2185 s (gravity,
                                                              #   magnetic tension and pressure, gas pressure, their sum);
                                                              #   --rope: on the line through the flux rope, wherever it is
  python3 tools/plot_leake2010.py mass RUN                    # the plasma lifted above 3 Mm: its mass budget (with the flow
                                                              #   or not), temperature, pressure, draining speed, heating time
  python3 tools/plot_leake2010.py drain RUN                   # why the lifted plasma drains or does not: how cold it is,
                                                              #   the field it sits on (flat, weak), its speed along the
                                                              #   field and how much of its weight its pressure carries
  python3 tools/plot_leake2010.py fastest RUN --times 2070,2087   # where the fastest waves are (|v| + c_f): position,
                                                              #   speeds, density, field, temperature (a step that collapses)
  python3 tools/plot_leake2010.py column RUN --times 1380,1610 --x 0   # the profiles on a vertical line as a table: density,
                                                              #   temperature, pressure, field, beta, velocity, DynSGS
  python3 tools/plot_leake2010.py null RUN [RUN2 ...]         # on x = 0 vs time: the top of the emerged flux (the null
                                                              #   above it), the flux between the tube axis and it, the
                                                              #   mass lifted, and the DynSGS coefficients at the null
  python3 tools/plot_leake2010.py noise RUN --times 2185,2530  # stripes and grain: the node-to-node part of log rho in each
                                                              #   element (along y, along x) and the nodes held at the floors
  python3 tools/plot_leake2010.py core RUN [RUN2 ...]         # the plasma at the tube axis vs time: height, density against
                                                              #   the atmosphere's (buoyancy), B_z/rho against the initial
                                                              #   tube's (kept in ideal MHD), the in-plane flux around it
  python3 tools/plot_leake2010.py lines RUN [RUN2 ...]        # Erratum Fig. 3 line by line: its field lines are contours of
                                                              #   A_z at known levels, so each is a line of the run; their
                                                              #   heights on x = 0 in the run and in the figure, the tube
                                                              #   axis, the flux reconnected at the null
  python3 tools/plot_leake2010.py --paper-lines density RUN   # any map with exactly the erratum's field lines (not N
                                                              #   equally spaced ones): the panels then compare line for line

Every command takes --outdir (default figs_leake2010), --format (png), --dpi, --lw (the width of the field
lines) and --prefix, before or after the command name. A figure is <prefix><variable>.<format>, by default
flux_emergence_Leake2010_density.png, ..._mu_dsgs_momentum.png, ..._Bz.png, ..._dome.png, and so on. Needs numpy and matplotlib;
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
# The field lines of the erratum's Fig. 3 are contours of A_z at a constant interval. Interval and levels were
# found from panel (a) (920 s) against this case's own quadrupole: its analytic A_z (A_z -> 0 at infinity) at the
# 83 places where the figure's lines cross x = +-8, +-9, ..., +-12 Mm falls on PAPER_A0 - k PAPER_DA with an rms
# misfit of 2 G Mm (the figure's vertical scale, checked on its density far from the tube, is good to 0.07 Mm,
# which is 2.4 G Mm on the levels). The quadrupole's null has A_z = 820.49 G Mm: line 0 hugs the separatrix, and
# lines 1, 2, ... are crossed twice on x = 0, in the emerged flux below the null and in the overlying field.
PAPER_DA = 61.35    # G Mm (6.135e9 Mx/cm)
PAPER_A0 = 820.0    # G Mm
# Heights [Mm] on x = 0 of the lines of the emerged flux (below the null), by line number, read off the figure
# (+-0.03 Mm). Lines 1 (1955 s), 2 (2185 s), 7 (1380 s) and 9 (1610-2185 s) are behind a B_x arrow there: taken
# from x = +-2.3 Mm and the offsets of their neighbours (+-0.1 Mm). At 1380 s lines 14-16 merge (+-0.06 Mm). The
# lines missing at the top have reconnected at the null: line 1 between 1955 and 2185 s, 2 and 3 by 2530 s.
ERRATUM_FIG3_LINES = {
    1380: {1: 6.09, 2: 5.37, 3: 4.80, 4: 4.28, 5: 3.82, 6: 3.39, 7: 3.02, 8: 2.70, 9: 2.39, 10: 2.10, 11: 1.87,
           12: 1.64, 13: 1.38, 14: 1.16, 15: 0.98, 16: 0.80, 17: 0.60},
    1610: {1: 8.76, 2: 7.67, 3: 6.70, 4: 5.89, 5: 5.17, 6: 4.54, 7: 3.99, 8: 3.48, 9: 3.05, 10: 2.70, 11: 2.33,
           12: 2.01, 13: 1.75, 14: 1.52, 15: 1.26, 16: 0.98, 17: 0.72},
    1955: {1: 11.45, 2: 9.40, 3: 7.90, 4: 6.75, 5: 5.72, 6: 4.91, 7: 4.22, 8: 3.65, 9: 3.12, 10: 2.73, 11: 2.33,
           12: 1.98, 13: 1.70, 14: 1.44, 15: 1.18, 16: 0.92, 17: 0.66},
    2185: {2: 11.1, 3: 8.88, 4: 7.27, 5: 6.03, 6: 5.03, 7: 4.22, 8: 3.59, 9: 3.01, 10: 2.53, 11: 2.16, 12: 1.81,
           13: 1.49, 14: 1.24, 15: 1.01, 16: 0.78, 17: 0.49},
    2530: {4: 9.17, 5: 7.30, 6: 5.95, 7: 4.86, 8: 3.99, 9: 3.30, 10: 2.73, 11: 2.27, 12: 1.93, 13: 1.55, 14: 1.26,
           15: 1.01, 16: 0.75, 17: 0.40}}
# The same lines in the overlying field (above the null) on x = 0, and the O-point of the tube (its axis)
ERRATUM_FIG3_OVERLYING = {920: {1: 9.74, 2: 13.68}, 1380: {1: 9.57, 2: 13.07, 3: 16.98},
                          1610: {1: 11.64, 2: 14.63, 3: 18.19}, 1955: {1: 13.51, 2: 16.0, 3: 19.31}, 2185: {2: 16.41}}
ERRATUM_FIG3_AXIS = {1380: 0.46, 1610: 0.52, 1955: 0.46, 2185: 0.29, 2530: 0.23}

# The original paper's Fig. 9(a) (ApJ 722:550, 2010; before the erratum removed the stray resistivity at the
# origin): height of the centre of the three tubes in a field-free corona, read exactly from the vector paths
# of its EPS (arXiv 1007.5484v2, axis_shear_run2a4a5a.eps), every 100 s. To ~1400 s it agrees with the
# erratum's Fig. 2(a) (the stray resistivity acts only once the centre's crescent sheet forms): all three
# centres settle at 0.27-0.41 Mm by 1200-1400 s.
PAPER_FIG9A_TUBE1 = np.array([[0, -1.84], [100, -1.77], [200, -1.63], [300, -1.48], [400, -1.24], [500, -1.00], [600, -0.70], [700, -0.46], [800, -0.22], [900, -0.06], [1000, 0.06], [1100, 0.20], [1200, 0.23], [1300, 0.27], [1400, 0.34], [1500, 0.34], [1600, 0.34], [1700, 0.34], [1800, 0.34], [1900, 0.27], [2000, 0.27], [2100, 0.20], [2200, 0.20], [2300, 0.20], [2400, 0.27], [2500, 0.27], [2600, 0.27], [2700, 0.27], [2800, 0.34], [2900, 0.34], [3000, 0.34], [3100, 0.34], [3200, 0.27], [3300, 0.27], [3400, 0.27], [3500, 0.27], [3600, 0.27], [3700, 0.27], [3800, 0.27], [3900, 0.34], [4000, 0.34]])
PAPER_FIG9A_TUBE2 = np.array([[0, -1.84], [100, -1.74], [200, -1.58], [300, -1.41], [400, -1.17], [500, -0.93], [600, -0.63], [700, -0.32], [800, -0.15], [900, 0.01], [1000, 0.13], [1100, 0.20], [1200, 0.27], [1300, 0.34], [1400, 0.43], [1500, 0.55], [1600, 0.76], [1700, 1.09], [1800, 1.57], [1900, 2.18], [2000, 2.73], [2100, 3.14], [2200, 3.82], [2300, 4.25], [2400, 4.63], [2500, 6.49], [2600, 7.00], [2700, 7.31], [2800, 7.37], [2900, 7.57], [3000, 7.97], [3100, 8.40], [3200, 8.57], [3300, 8.57], [3400, 8.45], [3500, 8.36], [3600, 8.39], [3700, 8.57], [3800, 8.85], [3900, 8.99], [4000, 9.06]])
PAPER_FIG9A_TUBE3 = np.array([[0, -1.84], [100, -1.77], [200, -1.63], [300, -1.48], [400, -1.24], [500, -1.00], [600, -0.70], [700, -0.39], [800, -0.22], [900, -0.06], [1000, 0.11], [1100, 0.20], [1200, 0.27], [1300, 0.27], [1400, 0.34], [1500, 0.34], [1600, 0.34], [1700, 0.34], [1800, 0.34], [1900, 0.27], [2000, 0.27], [2100, 0.20], [2200, 0.20], [2300, 0.27], [2400, 0.27], [2500, 0.27], [2600, 0.27], [2700, 0.34], [2800, 0.34], [2900, 0.34], [3000, 0.34], [3100, 0.34], [3200, 0.34], [3300, 0.34], [3400, 0.39], [3500, 0.55], [3600, 0.62], [3700, 0.68], [3800, 0.71], [3900, 0.83], [4000, 0.83]])
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
    successor continued (_chain). A later directory wins a repeated time, except with its FIRST snapshot: that
    is the state it was restarted from, written with a DynSGS coefficient evaluated before the sensor had a
    history (0.05 of the true one on REF 1), so the snapshot of the run it continues is kept there. DIR@T uses
    DIR only up to T seconds: BASE@920,NEW is BASE to 920 s and then NEW, for a NEW whose lk10_case.txt does not
    lead back to BASE."""
    dirs = []
    for d in spec.split(","):
        d, _, cap = d.partition("@")                       # DIR@T: this directory only up to T seconds
        chain = _chain(d) if FOLLOW_CHAIN else [(os.path.normpath(d), np.inf)]
        if cap:
            chain[-1] = (chain[-1][0], float(cap))
        for c, tmax in chain:
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
        t0 = min(items)[0]
        for t, f in used:
            if t == t0 and round(t, 6) in out:
                continue                                   # the first snapshot of a restart: keep the parent's
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
    g["dir"] = os.path.dirname(os.path.abspath(path))
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


_GL600 = None


def quadrupole_az(x, y, qd):
    """A_z [G Mm] of the case's quadrupole at (x, y) [Mm] (initialize.jl, lk10_quadrupole): harmonic above the
    source line y = -7.5 Mm, on which A_z = q_d (x/L0)^2 [1 - (x/x_a)^2] for |x| <= x_a = 21 Mm (q_d in G cm), and
    zero at infinity; 600-point Gauss-Legendre quadrature of its Poisson integral."""
    global _GL600
    if _GL600 is None:
        _GL600 = np.polynomial.legendre.leggauss(600)
    xa, ys = 21.0, -7.5
    xi, w = _GL600[0] * xa, _GL600[1] * xa
    x, y = np.asarray(x, float), np.asarray(y, float)
    h = (y - ys)[..., None]
    As = w * (qd * 1e-8) * (xi / L0_MM) ** 2 * (1.0 - (xi / xa) ** 2)
    return (As * h / ((x[..., None] - xi) ** 2 + h * h)).sum(-1) / np.pi


def a_gauged(g):
    """A_z [G Mm] on the lattice in the gauge of the initial state (A_z -> 0 at infinity), in which the erratum's
    lines are PAPER_A0 - k PAPER_DA; None for a run without the quadrupole. The flux function of a_grid (zero at
    the bottom left corner) is shifted so that the top of the column x = 0 has the quadrupole's value: that node
    is on the closed top wall, where v = 0 keeps A_z at its initial value, 90 Mm from the tube."""
    d = g.get("dir", "")
    if _case_value(d, "field", "quadrupole") != "quadrupole":
        return None
    A = a_grid(g) * B0_GAUSS
    ix0 = np.argmin(np.abs(g["X"]))
    return A + (quadrupole_az(g["X"][ix0], g["Y"][-1], float(_case_value(d, "qd", "5.85e7"))) - A[-1, ix0])


LINE_WIDTH = 0.9       # field lines, points (--lw)
PAPER_LINES = False    # --paper-lines: the field lines of the maps are the erratum's (PAPER_A0 - k PAPER_DA)
PREFIX = "flux_emergence_Leake2010_"   # every figure is <PREFIX><variable>.<format> (--prefix)


def fig_path(a, variable, ext=None):
    return os.path.join(a.outdir, f"{PREFIX}{variable}.{ext or a.format}")


def _window(g, xlim, ylim):
    """Row and column slices of the lattice covering the window xlim x ylim [Mm] and one node beyond each side. A map
    is drawn there only: a REF 3 snapshot (1.8 M nodes) drawn whole, smooth-shaded, is a few million coloured
    triangles per panel, and four panels of it exhaust a login node's memory."""
    X, Y = g["X"], g["Y"]
    i0 = max(int(np.searchsorted(X, xlim[0])) - 1, 0); i1 = min(int(np.searchsorted(X, xlim[1])) + 1, len(X) - 1)
    j0 = max(int(np.searchsorted(Y, ylim[0])) - 1, 0); j1 = min(int(np.searchsorted(Y, ylim[1])) + 1, len(Y) - 1)
    return slice(j0, j1 + 1), slice(i0, i1 + 1)


def field_lines(ax, g, levels, color="w", lw=None, alpha=1.0, win=None):
    if isinstance(levels, str):            # "paper": the erratum's own lines, PAPER_A0 - k PAPER_DA (--paper-lines)
        A = a_gauged(g)
        if A is None:
            raise SystemExit("--paper-lines: the erratum's levels (its Fig. 3) are those of a quadrupole run")
        k = np.arange(np.floor((PAPER_A0 - np.nanmax(A)) / PAPER_DA), np.ceil((PAPER_A0 - np.nanmin(A)) / PAPER_DA) + 1)
        levels = np.sort(PAPER_A0 - k * PAPER_DA)
    else:
        A = a_grid(g)
    rs, cs = win if win is not None else (slice(None), slice(None))   # A over the whole lattice, drawn in the window
    ax.contour(g["X"][cs], g["Y"][rs], A[rs, cs], levels=levels, colors=color, linewidths=LINE_WIDTH if lw is None else lw,
               linestyles="solid", alpha=alpha)
    return A


# The colour scale of the DynSGS coefficient in the group's figures: dark navy where there is none, through
# cyan and cream to dark red where it is largest (ParaView's "Cool to Warm (Extended)" control points, here
# interpolated in RGB). The density and B_z maps keep Leake et al.'s rainbow scale, to be laid beside theirs.
_C2W = [(0.0, 0.0, 0.0, 0.34902), (0.03125, 0.039216, 0.062745, 0.380392), (0.0625, 0.062745, 0.117647, 0.411765),
        (0.09375, 0.090196, 0.184314, 0.45098), (0.125, 0.12549, 0.262745, 0.501961), (0.15625, 0.160784, 0.337255, 0.541176),
        (0.1875, 0.2, 0.396078, 0.568627), (0.21875, 0.239216, 0.454902, 0.6), (0.25, 0.286275, 0.521569, 0.65098),
        (0.28125, 0.337255, 0.592157, 0.701961), (0.3125, 0.388235, 0.654902, 0.74902), (0.34375, 0.466667, 0.737255, 0.819608),
        (0.375, 0.572549, 0.819608, 0.878431), (0.40625, 0.654902, 0.866667, 0.909804), (0.4375, 0.752941, 0.917647, 0.941176),
        (0.46875, 0.823529, 0.956863, 0.968627), (0.5, 0.941176, 0.984314, 0.988235), (0.5, 0.988235, 0.960784, 0.901961),
        (0.52, 0.988235, 0.945098, 0.85098), (0.54, 0.980392, 0.898039, 0.784314), (0.5625, 0.968627, 0.835294, 0.698039),
        (0.59375, 0.94902, 0.733333, 0.588235), (0.625, 0.929412, 0.65098, 0.509804), (0.65625, 0.909804, 0.564706, 0.435294),
        (0.6875, 0.878431, 0.458824, 0.352941), (0.71875, 0.839216, 0.388235, 0.286275), (0.75, 0.760784, 0.294118, 0.211765),
        (0.78125, 0.701961, 0.211765, 0.168627), (0.8125, 0.65098, 0.156863, 0.129412), (0.84375, 0.6, 0.094118, 0.094118),
        (0.875, 0.54902, 0.066667, 0.098039), (0.90625, 0.501961, 0.05098, 0.12549), (0.9375, 0.45098, 0.054902, 0.172549),
        (0.96875, 0.4, 0.054902, 0.192157), (1.0, 0.34902, 0.070588, 0.211765)]
DSGS_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list("dsgs_cool_to_warm", [(x, (r, g, b)) for x, r, g, b in _C2W], N=512)
INK = "#1f2933"        # titles

# What a map shows: the output field(s) it reads (the first one present), whether its log10 is taken, the
# colour map and its default limits, the colour-bar label, the colour of the text and field lines drawn on it.
MAPS = {
    "density": dict(names=("log10rho_cgs",), log=False, cmap="jet", clim=(-16.0, -5.3), ink="w", line="w", alpha=1.0,
                    label=r"log$_{10}$ $\rho$ (g cm$^{-3}$)", tick=2.0,
                    title=r"Density: log$_{10}$ $\rho$ (g cm$^{-3}$), with the field lines and the direction of $B_x$"),
    # The DynSGS coefficients. The solver has ONE coefficient mu per element, from the largest normalized residual
    # of the eight equations (rho, rho u, rho v, E, rho w, Bx, By, Bz), and writes one output field per group of
    # equations whose applied coefficient differs, named after them (src/io/write_output.jl):
    #   mu_dsgs_ρu_ρv_ρw   momentum:            mu + the paper's explicit viscosity (0.013 above the photosphere)
    #   mu_dsgs_Bx_By_Bz   induction:           its own coefficient eta from the induction residuals alone
    #                                           (JEXPRESSO_LK10_ETA=induction) + the anomalous resistivity where active
    #   mu_dsgs_ρ_ρE_ψ     mass, energy, psi:   mu itself
    # (JEXPRESSO_LK10_ETA=all: the induction equations share mu, field mu_dsgs_ρ_ρE_Bx_By_Bz_ψ.)
    "momentum":  dict(names=("mu_dsgs_ρu_ρv_ρw", "mu_dsgs"), log=True, cmap=DSGS_CMAP, clim=(-2.0, 2.0), ink="w", line="w",
                      alpha=0.45, plate=True, var="mu_dsgs_momentum", tick=1.0,
                      label=r"log$_{10}$ $\mu$ ($\mu$ in $L_0v_0$ = 970 km$^2$ s$^{-1}$)",
                      title=r"DynSGS $\mu$ of the momentum equations ($\rho u$, $\rho v$, $\rho w$)",
                      lin=dict(clim=(0.0, 5.0), label=r"$\mu$ (in $L_0v_0$ = 970 km$^2$ s$^{-1}$)")),
    "induction": dict(names=("mu_dsgs_Bx_By_Bz", "mu_dsgs_ρ_ρE_Bx_By_Bz_ψ", "mu_dsgs"), log=True, cmap=DSGS_CMAP, clim=(-4.0, 1.0),
                      ink="w", line="w", alpha=0.45, plate=True, var="mu_dsgs_induction", tick=1.0,
                      label=r"log$_{10}$ $\eta$ ($\eta$ in $L_0v_0$ = 970 km$^2$ s$^{-1}$)",
                      title=r"DynSGS $\eta$ of the induction equations ($B_x$, $B_y$, $B_z$)",
                      lin=dict(clim=(0.0, 0.5), label=r"$\eta$ (in $L_0v_0$ = 970 km$^2$ s$^{-1}$)")),
    "energy":    dict(names=("mu_dsgs_ρ_ρE_ψ", "mu_dsgs_ρ_ρE_Bx_By_Bz_ψ", "mu_dsgs_ρE_ψ", "mu_dsgs"), log=True, cmap=DSGS_CMAP,
                      clim=(-2.0, 2.0), ink="w", line="w", alpha=0.45, plate=True, var="mu_dsgs_mass_energy", tick=1.0,
                      label=r"log$_{10}$ $\mu$ ($\mu$ in $L_0v_0$ = 970 km$^2$ s$^{-1}$)",
                      title=r"DynSGS $\mu$ of the mass, energy and $\psi$ equations ($\rho$, $E$, $\psi$)",
                      lin=dict(clim=(0.0, 5.0), label=r"$\mu$ (in $L_0v_0$ = 970 km$^2$ s$^{-1}$)")),
}
MAP_ALIAS = {"mu": "momentum", "eta": "induction", "rho": "energy"}


def map_spec(name, linear=False):
    """The map `name` of MAPS; with linear = True a DynSGS coefficient itself instead of its log10."""
    spec = dict(MAPS[MAP_ALIAS.get(name, name)])
    spec.setdefault("var", name)
    if linear and "lin" in spec:
        spec.update(spec["lin"], log=False, tick=None)
    return spec


def cbar_ticks(cb, clim, spec):
    if spec.get("tick"):
        cb.set_ticks(np.arange(spec["tick"] * np.ceil(clim[0] / spec["tick"]), clim[1] + 1e-9, spec["tick"]))
    cb.set_label(spec["label"], fontsize=9); cb.ax.tick_params(labelsize=8)


def map_field(g, spec):
    """The grid of a map's quantity (MAPS): the first of its output fields the snapshot has, log10 if asked."""
    name = next((n for n in spec["names"] if n in g["fields"]), None)
    if name is None:
        raise SystemExit(f"the snapshots have none of the fields {', '.join(spec['names'])} "
                         "(the DynSGS coefficients are written with :plot_dsgs => true)")
    Z = grid2d(g, g["fields"][name])
    return np.log10(np.maximum(Z, 1e-300)) if spec["log"] else Z


def corner_text(ax, x, y, text, spec, size=9, **kw):
    """A label in a corner of a panel, in the map's ink; on the DynSGS scale it gets a dark plate, since the
    tone under it can be anything from navy to cream."""
    plate = dict(boxstyle="round,pad=0.25", fc="#0b1030", ec="none", alpha=0.7) if spec.get("plate") else None
    ax.text(x, y, text, color=spec["ink"], transform=ax.transAxes, fontsize=size, weight="bold", bbox=plate, **kw)


def draw_map(ax, g, spec, clim, levels, xlim, ylim, arrows=()):
    """One panel: the map's quantity in its colours, the field lines, the B_x arrows; returns the image."""
    rs, cs = _window(g, xlim, ylim)
    Z = map_field(g, spec)[rs, cs]
    im = ax.pcolormesh(g["X"][cs], g["Y"][rs], np.clip(Z, clim[0], clim[1]), cmap=spec["cmap"], vmin=clim[0], vmax=clim[1],
                       shading="gouraud", rasterized=True)
    if levels is not None:
        field_lines(ax, g, levels, color=spec["line"], alpha=spec["alpha"], win=(rs, cs))
    if arrows:
        bx_arrows(ax, g, arrows)
    ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_aspect("equal")
    ax.tick_params(labelsize=8)
    return im


def a_levels(g, xlim, ylim, n):
    """n equally spaced levels of A_z over the window of the first snapshot ("constant contour intervals"); with
    --paper-lines the erratum's own levels instead (field_lines draws them in the gauge of a_gauged)."""
    if PAPER_LINES:
        return "paper"
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


def _layout(fig, axs, title=""):
    """Margins in inches, so that one panel or six keep their labels, the title naming what the colours show,
    and room for the colour bar with its label."""
    w, h = fig.get_size_inches()
    fig.subplots_adjust(left=0.62 / w, right=1 - 0.12 / w, top=1 - (0.42 if title else 0.1) / h, bottom=1.30 / h,
                        wspace=0.2, hspace=0.14)
    if title:
        fig.suptitle(title, y=1 - 0.13 / h, fontsize=11, weight="bold", color="#1f2933")


def _cbar_axes(fig, axs):
    """A horizontal colour bar under the panels, spanning them, high enough for its ticks and label."""
    w, h = fig.get_size_inches()
    x0 = axs[-1, 0].get_position().x0
    x1 = axs[-1, -1].get_position().x1
    return fig.add_axes([x0 + 0.15 * (x1 - x0), 0.60 / h, 0.7 * (x1 - x0), 0.13 / h])


def panel_grid(n):
    ncol = 3 if n > 4 else n
    nrow = int(np.ceil(n / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(3.4 * ncol + 0.7, 3.0 * nrow + 1.75), squeeze=False)
    return fig, axs


def cmd_density(a):
    spec = MAPS["density"]
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
        im = draw_map(ax, g, spec, a.clim, levels, a.xlim, a.ylim, a.arrows)
        lab = f"({chr(97 + k)}) {a.label}" if (k == 0 and a.label) else f"({chr(97 + k)})"
        ax.text(0.03, 0.96, lab, color="w", transform=ax.transAxes, va="top", fontsize=9, weight="bold")
        ax.text(0.97, 0.96, f"t= {t:.0f}s", color="w", transform=ax.transAxes, va="top", ha="right",
                fontsize=9, weight="bold")
    for ax in axs.flat[len(sel):]:
        ax.axis("off")
    for ax in axs[:, 0]:
        ax.set_ylabel("Y (Mm)", fontsize=8)
    for ax in axs[-1, :]:
        ax.set_xlabel("X (Mm)", fontsize=8)
    _layout(fig, axs, spec["title"])
    cb = fig.colorbar(im, cax=_cbar_axes(fig, axs), orientation="horizontal")
    cb.set_ticks(np.arange(2 * np.ceil(a.clim[0] / 2), a.clim[1] + 1e-9, 2.0))
    cb.set_label(spec["label"], fontsize=9); cb.ax.tick_params(labelsize=8)
    out = fig_path(a, "density")
    fig.savefig(out, dpi=a.dpi); plt.close(fig)
    print("wrote", out)


def cmd_dsgs(a):
    """A DynSGS coefficient at the times of the density figure, in the same window, with the field lines, as
    log10 (or --linear) in code units L0 v0 = 970 km^2/s, on the colour scale of the group's DynSGS figures.
    --which says whose (MAPS): the momentum equations' (default; the element's mu plus the paper's explicit
    viscosity, 0.013 above the photosphere, so navy is "no DynSGS"), the induction equations' (eta), or the
    mass, energy and psi equations' (mu itself)."""
    spec = map_spec(a.which, a.linear)
    clim = a.clim or spec["clim"]
    items = runs_of(a.run)
    times = [float(t) for t in a.times.split(",")] if a.times else ERRATUM_FIG3_TIMES
    sel = nearest(items, times)
    fig, axs = panel_grid(len(sel))
    for k, (t, f) in enumerate(sel):
        ax = axs.flat[k]
        g = load(f, list(spec["names"]))
        levels = None if a.nlines <= 0 else a_levels(g, a.xlim, a.ylim, a.nlines)
        im = draw_map(ax, g, spec, clim, levels, a.xlim, a.ylim)
        corner_text(ax, 0.03, 0.96, f"({chr(97 + k)})", spec, va="top")
        corner_text(ax, 0.97, 0.96, f"t= {t:.0f}s", spec, va="top", ha="right")
    for ax in axs.flat[len(sel):]:
        ax.axis("off")
    for ax in axs[:, 0]:
        ax.set_ylabel("Y (Mm)", fontsize=8)
    for ax in axs[-1, :]:
        ax.set_xlabel("X (Mm)", fontsize=8)
    _layout(fig, axs, spec["title"])
    cb = fig.colorbar(im, cax=_cbar_axes(fig, axs), orientation="horizontal")
    cbar_ticks(cb, clim, spec)
    out = fig_path(a, spec["var"])
    fig.savefig(out, dpi=a.dpi); plt.close(fig)
    print("wrote", out)


def cmd_movie(a):
    """One image per snapshot of the run (every output time, or --tmin/--tmax/--every), all with the same
    window, colour limits and field-line levels (those of the first snapshot of the run, so that a line is the
    same flux surface in every frame), and an animated GIF of them in which the simulated time runs evenly
    (a frame stays in proportion to the time to the next one: the snapshots are not equally spaced)."""
    spec = map_spec(a.field, a.linear)
    clim = a.clim or spec["clim"]
    items = runs_of(a.run)
    levels = None if a.nlines <= 0 else a_levels(load(items[0][1], []), a.xlim, a.ylim, a.nlines)
    sel = [(t, f) for t, f in items if a.tmin - 0.5 <= t <= a.tmax + 0.5][::max(a.every, 1)]
    if not sel:
        raise SystemExit("   no snapshot in that range of times")
    var = spec["var"]
    fdir = os.path.join(a.outdir, f"frames_{var}")
    os.makedirs(fdir, exist_ok=True)
    arrows = a.arrows if spec["var"] == "density" else []
    files = []
    for k, (t, f) in enumerate(sel):
        g = load(f, list(spec["names"]))
        fig = plt.figure(figsize=(5.4, 6.0))
        ax = fig.add_axes([0.13, 0.215, 0.82, 0.715])
        fig.suptitle(spec["title"], y=0.975, fontsize=9.5, weight="bold", color="#1f2933")
        im = draw_map(ax, g, spec, clim, levels, a.xlim, a.ylim, arrows)
        if a.label:
            corner_text(ax, 0.03, 0.97, a.label, spec, size=10, va="top")
        corner_text(ax, 0.97, 0.97, f"t = {t:.0f} s", spec, size=10, va="top", ha="right")
        ax.set_xlabel("X (Mm)", fontsize=9); ax.set_ylabel("Y (Mm)", fontsize=9)
        cb = fig.colorbar(im, cax=fig.add_axes([0.20, 0.095, 0.68, 0.022]), orientation="horizontal")
        cbar_ticks(cb, clim, spec)
        out = os.path.join(fdir, f"{PREFIX}{var}_{k + 1:04d}_t{t:05.0f}s.png")
        fig.savefig(out, dpi=a.dpi); plt.close(fig)
        files.append(out)
        print(f"   {os.path.basename(out)}", flush=True)
    print(f"wrote {len(files)} images in {fdir}")
    if a.no_gif or len(files) < 2:
        return
    try:
        from PIL import Image
    except ImportError:
        print("   no GIF: the Python package Pillow is missing (pip install pillow), or join the images yourself, e.g.\n"
              f"   ffmpeg -framerate {a.fps:g} -pattern_type glob -i '{fdir}/*.png' {PREFIX}{var}.mp4")
        return
    T = np.array([t for t, _ in sel])
    gaps = np.append(np.diff(T), np.median(np.diff(T)))
    base = 1000.0 / a.fps                                  # ms for a snapshot interval of the median length
    dur = [int(max(20, round(base * d / np.median(np.diff(T))))) for d in gaps]
    dur[-1] = int(max(dur[-1], 3 * base))                 # hold the last frame
    # a GIF has 256 colours a frame: each frame gets its own best 256 (one palette for all bands the gradients)
    frames = [Image.open(f).convert("RGB").quantize(colors=256, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE)
              for f in files]
    out = fig_path(a, var, "gif")
    frames[0].save(out, save_all=True, append_images=frames[1:], duration=dur, loop=0, optimize=False)
    print(f"wrote {out} ({len(frames)} frames, {os.path.getsize(out) / 1e6:.1f} MB)")


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
        rs, cs = _window(g, a.xlim, a.ylim)
        Z = grid2d(g, g["fields"]["Bz_G"])[rs, cs]
        im = ax.pcolormesh(g["X"][cs], g["Y"][rs], np.clip(Z, a.bclim[0], a.bclim[1]), cmap="jet",
                           vmin=a.bclim[0], vmax=a.bclim[1], shading="gouraud", rasterized=True)
        field_lines(ax, g, levels, win=(rs, cs))
        ax.set_xlim(a.xlim); ax.set_ylim(a.ylim); ax.set_aspect("equal")
        ax.text(0.03, 0.04, f"({chr(97 + k)}) {a.label if k == 0 else ''}", color="w", transform=ax.transAxes, fontsize=8)
        ax.text(0.97, 0.04, f"{t:.0f}s", color="w", transform=ax.transAxes, ha="right", fontsize=8)
        ax.tick_params(labelsize=8)
    for ax in axs.flat[len(sel):]:
        ax.axis("off")
    _layout(fig, axs, r"Shear field $B_z$ (G), with the field lines")
    cb = fig.colorbar(im, cax=_cbar_axes(fig, axs), orientation="horizontal")
    cb.set_label(r"$B_z$ (G)", fontsize=9); cb.ax.tick_params(labelsize=8)
    out = fig_path(a, "Bz")
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
    out = fig_path(a, "axis")
    fig.savefig(out, dpi=a.dpi); plt.close(fig)
    print("wrote", out)


def tube_axis_initial(tube, mcz=1.3):
    """B_z/rho [G per photospheric density] at the axis of the initial tube (initialize.jl, lk10_tube_field, the
    tube in thermal equilibrium at y_t = -1.8 Mm): what the plasma at the axis keeps in ideal 2.5D MHD."""
    g = G_EXACT
    B0, a, c1, c2, yt = 7800.0 / B0_GAUSS, 2.0, 1.0, 10.0, -12.0
    T = 1.0 - g * yt / (mcz + 1.0)
    k0, c = 2.0 / a ** 2, 1.0 / a
    I = {1: 1.0 / k0, 2: 1.0 / k0 - 2.0 / (k0 + 1.0 / c1 ** 2) + 1.0 / (k0 + 2.0 / c1 ** 2),
         3: 1.0 / (k0 + 2.0 / c2 ** 2)}.get(tube)
    if I is None:
        return np.nan
    p = T ** (mcz + 1.0) - 0.5 * B0 ** 2 + 0.5 * c * c * B0 * B0 * I
    return 7800.0 / (p / T)


def core_row(g, mcz=1.3):
    """The plasma at the tube axis (the O-point of the in-plane field on x = 0): height [Mm], log10 rho, its
    density over the reference atmosphere's at that height, T [K], B_z [G], B_z/rho [G per photospheric density],
    the in-plane flux around the axis (A_z far above minus A_z at the axis, G Mm), and over the tube (B_z above
    2 % of its peak) int B_z^2 dA [G^2 Mm^2] and the flux-weighted mean of B_z/rho. In ideal 2.5D MHD the plasma
    at the O-point keeps its B_z/rho (no in-plane field there to shear it), the in-plane flux around the axis is
    kept, and so are the axial flux and the mass inside every closed field line; the density deficit against
    the atmosphere is the buoyancy."""
    ix0 = np.argmin(np.abs(g["X"])); Y = g["Y"]
    ya = axis_height(g)
    L = grid2d(g, g["fields"]["log10rho_cgs"])
    BZ = grid2d(g, g["fields"]["Bz"]) * B0_GAUSS
    A = a_grid(g)[:, ix0] * B0_GAUSS
    if not np.isfinite(ya):
        return None
    at = lambda F: float(np.interp(ya, Y, F[:, ix0]))
    rho = 10.0 ** at(L) / 2.02e-7
    rhoe = _ref_atmos(np.array([ya / L0_MM]), mcz)[0][0]
    R = 10.0 ** L / 2.02e-7
    tube = BZ > 0.02 * BZ.max()
    w = np.outer(_trapw(Y), _trapw(g["X"]))
    T = at(grid2d(g, g["fields"]["T"])) * 5700.0 if "T" in g["fields"] else np.nan
    # the DynSGS coefficients at the axis (L0 v0), when the snapshot has them: momentum (with the case's nu),
    # field (resistivity), mass
    cof = lambda n: at(grid2d(g, g["fields"][n])) if n in g["fields"] else np.nan
    en = next((n for n in ETA_NAMES if n in g["fields"]), None)
    return dict(axis=ya, logrho=at(L), dens=rho / rhoe, T=T, bz=at(BZ), q=at(BZ) / rho,
                flux=A[-1] - float(np.interp(ya, Y, A)), ebz=float((w * BZ ** 2)[tube].sum()),
                qmean=float((w * BZ * BZ / R)[tube].sum() / (w * BZ)[tube].sum()),
                mu=cof(MU_NAME), eta=cof(en) if en else np.nan, mur=cof("mu_dsgs_ρ"))


def cmd_core(a):
    """What the tube keeps on its way up (core_row), snapshot by snapshot, against the initial tube and against
    the erratum's axis (its Fig. 2(a), field-free corona, and the O-point of its Fig. 3)."""
    fig, axs = plt.subplots(1, 3, figsize=(12.5, 3.8))
    for spec in a.runs:
        lab, _, path = spec.rpartition("=")
        lab = lab or os.path.basename(path.rstrip("/"))
        items = runs_of(path)
        d0 = os.path.dirname(os.path.abspath(items[0][1]))
        mcz = float(_case_value(d0, "mcz", "1.3"))
        tube = int(float(_case_value(d0, "tube", "2")))
        q0 = tube_axis_initial(tube, mcz) if _case_value(d0, "tube_eq", "thermal") == "thermal" else np.nan
        rows = [(t, core_row(load(f, ["Bz", "log10rho_cgs", "T", MU_NAME, *ETA_NAMES, "mu_dsgs_ρ"]), mcz)) for t, f in items]
        rows = [(t, r) for t, r in rows if r is not None]
        if not rows:
            print(f"--- {lab}: no tube axis found"); continue
        T = np.array([t for t, _ in rows]); col = lambda k: np.array([r[k] for _, r in rows])
        first = rows[0][1]
        np.savetxt(os.path.join(a.outdir, f"core_{lab.replace(' ', '_')}.csv"),
                   np.c_[T, col("axis"), col("logrho"), col("dens"), col("T"), col("bz"), col("q"), col("flux"), col("ebz"), col("qmean")],
                   header="t [s], axis [Mm], log10 rho, rho/rho_e, T [K], B_z [G], B_z/rho [G per rho_ph], in-plane flux "
                          "around the axis [G Mm], int B_z^2 dA over the tube [G^2 Mm^2], flux-weighted mean B_z/rho", delimiter=",")
        print(f"--- {lab}: the plasma at the tube axis; B_z/rho of the initial tube {q0:.0f} G per photospheric density "
              f"(tube {tube}), kept at the axis in ideal MHD")
        print("    (read between the nodes, the axis of the initial state has 0.96 of it on REF 1 and 0.95 on REF 2 and 4)")
        print("    and the DynSGS coefficients at the axis [L0 v0 = 970 km^2/s]: momentum (with the case's nu), field, mass")
        print(f"{'t [s]':>7} {'axis Mm':>8} {'log rho':>8} {'rho/rho_e':>9} {'T K':>8} {'Bz G':>7} {'Bz/rho':>8} {'/initial':>8} "
              f"{'in-plane flux G Mm':>19} {'/first':>6} {'int Bz^2':>9} {'mean Bz/rho':>11} | {'mu':>7} {'eta':>7} {'mu_rho':>7}")
        for t, r in rows:
            print(f"{t:7.0f} {r['axis']:8.2f} {r['logrho']:8.2f} {r['dens']:9.3f} {r['T']:8.0f} {r['bz']:7.0f} {r['q']:8.0f} "
                  f"{r['q'] / q0:8.2f} {r['flux']:19.1f} {r['flux'] / first['flux']:6.3f} {r['ebz'] / first['ebz']:9.3f} "
                  f"{r['qmean'] / first['qmean']:11.3f} | {r['mu']:7.3g} {r['eta']:7.3g} {r['mur']:7.3g}")
        axs[0].plot(T, col("axis"), "-o", ms=3, label=lab)
        axs[1].plot(T, col("q") / q0, "-o", ms=3, label=lab)
        axs[2].plot(T, col("dens"), "-o", ms=3, label=lab)
    axs[0].plot(*ERRATUM_FIG2A_TUBE2.T, "g--", lw=1.2, label="Erratum Fig. 2(a), low tension")
    axs[0].plot(list(ERRATUM_FIG3_AXIS), list(ERRATUM_FIG3_AXIS.values()), "ks", mfc="none", ms=6, label="Erratum Fig. 3, O-point")
    axs[0].plot(*PAPER_FIG9A_TUBE3.T, ":", color="0.45", lw=1.0, label="Paper Fig. 9(a), high tension (tube 3)")
    axs[0].set_ylabel("height of the axis (Mm)"); axs[0].set_ylim(-2.0, 1.0)
    axs[1].axhline(1.0, color="k", lw=0.6)
    axs[1].set_ylabel(r"$B_z/\rho$ at the axis / initial tube (1 in ideal MHD)"); axs[1].set_ylim(0, 1.2)
    axs[2].axhline(1.0, color="k", lw=0.6)
    axs[2].set_ylabel(r"$\rho$ at the axis / atmosphere at that height"); axs[2].set_yscale("log")
    for ax in axs:
        ax.set_xlabel("Time (s)"); ax.set_xlim(0, a.tmax); ax.grid(alpha=0.3); ax.legend(fontsize=6)
    fig.tight_layout()
    out = fig_path(a, "core")
    fig.savefig(out, dpi=a.dpi); plt.close(fig)
    print("wrote", out)


def cmd_lare(a):
    """Our DynSGS momentum coefficient against the viscosity of the paper's code for the same flow. LareXd's tensor
    viscosity (lagran.F90 of 2.6-2.9, built without Q_MONO) gives each component s_ij of the traceless strain rate
    the stress L rho (0.1 c_f + 0.5 L |s_ij|) s_ij in every cell: a kinematic viscosity 0.05 c_f L + 0.25 L^2 |s_ij|
    (here with the node's largest |s_ij|, from finite differences on the lattice: an upper bound where the flow is
    resolved), L its cell, --lare-dx Mm (0.0375, the paper's cell near the surface). At each requested time, over
    the regions around the tube (its axis as in 'axis'): in the tube (|B| > 300 G); above it, the field-free
    atmosphere it pushes into (|x| < 3 Mm, from 0.5 Mm above the axis to 3 Mm, |B| <= 30 G); beside it in the
    convection zone (3 < |x| < 8 Mm, -1.5 < y < 0 Mm): the median of each, the percentiles of their ratio node by
    node, and the flow there (speed, Mach number at the 90th percentile). The coefficient is the momentum slot's
    (mu_dsgs_rhou_rhov_rhow, with the case's factors); with JEXPRESSO_LK10_VISC=velocity that is still the whole
    viscosity, mu_rho on the momentum plus (mu - mu_rho)_+ on the velocity."""
    items = runs_of(a.run)
    Lc = a.lare_dx / L0_MM
    for tt in a.times:
        t, f = min(items, key=lambda it: abs(it[0] - tt))
        g = load(f, ["Bz", "p", "ρ", "u", "v", MU_NAME])
        if MU_NAME not in g["fields"]:
            raise SystemExit(f"lare: {f} has no {MU_NAME} (a run without :plot_dsgs)")
        X, Y = g["xc"], g["yc"]
        G = lambda n: grid2d(g, g["fields"][n])
        R, p, u, v, mu = G("ρ"), G("p"), G("u"), G("v"), G(MU_NAME)
        B2 = G("Bx") ** 2 + G("By") ** 2 + G("Bz") ** 2
        cs = np.sqrt(np.maximum(GAMMA * p / R, 0.0)); cf = np.sqrt(cs ** 2 + B2 / R)
        ux, uy = np.gradient(u, X, axis=1), np.gradient(u, Y, axis=0)
        vx, vy = np.gradient(v, X, axis=1), np.gradient(v, Y, axis=0)
        smax = np.maximum(np.maximum(np.abs(2 / 3 * ux - 1 / 3 * vy), np.abs(2 / 3 * vy - 1 / 3 * ux)), np.abs(0.5 * (uy + vx)))
        nuL = 0.05 * cf * Lc + 0.25 * Lc * Lc * smax
        Bt = np.sqrt(B2) * B0_GAUSS
        ya = axis_height(g)
        XX, YY = np.meshgrid(g["X"], g["Y"])
        spd = np.sqrt(u ** 2 + v ** 2)
        print(f"--- {os.path.basename(os.path.dirname(f))}, t = {t:.0f} s, tube axis at {ya:.2f} Mm; viscosities in L0 v0 "
              f"(970 km^2/s); Lare's cell {a.lare_dx * 1000:.1f} km")
        print(f"{'region':34s} {'nodes':>6} | {'ours':>7} {'Lare':>7} | {'ours/Lare p10':>13} {'p50':>6} {'p90':>6} | "
              f"{'|v| km/s':>8} {'Mach p90':>8}")
        regs = {"in the tube, |B| > 300 G": Bt > 300.0,
                "above it, |x|<3, axis+0.5 to 3 Mm": (np.abs(XX) < 3) & (YY > ya + 0.5) & (YY < 3) & (Bt <= 30.0),
                "beside it, 3<|x|<8, -1.5<y<0 Mm": (np.abs(XX) > 3) & (np.abs(XX) < 8) & (YY > -1.5) & (YY < 0)}
        for k, m in regs.items():
            if m.sum() == 0:
                print(f"{k:34s} {0:6d}"); continue
            r = np.percentile((mu / nuL)[m], [10, 50, 90])
            print(f"{k:34s} {m.sum():6d} | {np.median(mu[m]):7.4f} {np.median(nuL[m]):7.4f} | {r[0]:13.1f} {r[1]:6.1f} {r[2]:6.1f} | "
                  f"{np.median(spd[m]) * V0_KMS:8.2f} {np.percentile((spd / np.maximum(cs, 1e-12))[m], 90):8.2f}")


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
    out = fig_path(a, "dome")
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
    out = fig_path(a, "eshear")
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
    out = fig_path(a, "atmos")
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
    out = fig_path(a, "null")
    fig.savefig(out, dpi=a.dpi); plt.close(fig)
    print("wrote", out)


def lines_row(g):
    """On x = 0, in the gauge of a_gauged: the tube axis and the top of the emerged flux (as null_row), A_z there,
    and the heights of the erratum's field lines PAPER_A0 - k PAPER_DA: in the emerged flux (between the axis and
    the top, where A_z rises with height) and in the overlying field (from the top to 40 Mm, where it falls)."""
    A2 = a_gauged(g)
    if A2 is None:
        raise SystemExit("lines: the erratum's levels (its Fig. 3) are those of a quadrupole run")
    ix0 = np.argmin(np.abs(g["X"]))
    Y, A = g["Y"], A2[:, ix0]
    Bx = grid2d(g, g["fields"]["Bx"])[:, ix0]
    axis = axis_height(g)
    tops = [y for y in _bx_crossings(Y, Bx, -1) if y > max(axis if np.isfinite(axis) else 0.0, 1.0) and y < 40.0]
    top = tops[0] if tops else np.nan
    out = dict(axis=axis, top=top, Atop=np.nan, Aaxis=np.nan, dome={}, over={}, ntop=len(tops))
    if not (np.isfinite(top) and np.isfinite(axis)):
        return out
    Ai = lambda y: float(np.interp(y, Y, A))
    out["Atop"], out["Aaxis"] = Ai(top), Ai(axis)

    def cross(level, lo, hi, rising):
        m = (Y > lo) & (Y < hi)
        y, v = np.r_[lo, Y[m], hi], np.r_[Ai(lo), A[m], Ai(hi)]
        for i in range(len(y) - 1):
            if (v[i] < level <= v[i + 1]) if rising else (v[i] > level >= v[i + 1]):
                return y[i] + (level - v[i]) * (y[i + 1] - y[i]) / (v[i + 1] - v[i])
        return np.nan
    for k in range(0, 30):
        lev = PAPER_A0 - k * PAPER_DA
        if out["Aaxis"] < lev < out["Atop"]:
            out["dome"][k] = cross(lev, axis, top, True)
        if lev < out["Atop"] and k <= 6:
            out["over"][k] = cross(lev, top, 40.0, False)
    return out


def cmd_lines(a):
    """The erratum's Fig. 3 line by line. Its field lines are contours of A_z at a constant interval, at levels
    known in this case's gauge (PAPER_DA, PAPER_A0), so each line of the figure is a line of the run: this prints
    the height on x = 0 of every one of them in the run beside the height read off the figure, at the figure's
    times, with the height of the tube axis, the top of the emerged flux and the flux reconnected there (the fall
    of A_z at the top from the null's 820.49 G Mm); writes them for every snapshot (CSV) and draws both."""
    cmap = plt.get_cmap("viridis")
    col = lambda k: cmap((k - 1) / 16.0)
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.6))
    ptimes = sorted(ERRATUM_FIG3_LINES)
    want = [float(t) for t in a.times.split(",")] if a.times else ptimes
    for ir, spec in enumerate(a.runs):
        lab, _, path = spec.rpartition("=")
        lab = lab or os.path.basename(path.rstrip("/"))
        ls = ("-", "--", ":", "-.")[ir % 4]
        rows = [(t, lines_row(load(f, ["Bz"]))) for t, f in runs_of(path)]
        T = np.array([t for t, _ in rows])
        H = np.array([[r["dome"].get(k, np.nan) for k in range(1, 18)] for _, r in rows])
        O = np.array([[r["over"].get(k, np.nan) for k in range(1, 4)] for _, r in rows])
        name = lab.replace(" ", "_")
        np.savetxt(os.path.join(a.outdir, f"lines_{name}.csv"),
                   np.c_[T, [r["axis"] for _, r in rows], [r["top"] for _, r in rows], [r["Atop"] for _, r in rows], H, O],
                   header="t [s], axis [Mm], top of the emerged flux [Mm], A_z there [G Mm], heights on x = 0 [Mm] of the "
                          "erratum's lines 1..17 in the emerged flux, and of lines 1..3 in the overlying field", delimiter=",")
        print(f"--- {lab}: the erratum's field lines (A_z = {PAPER_A0:g} - {PAPER_DA:g} k G Mm) on x = 0, heights in Mm")
        for tq in want:
            i = int(np.argmin(np.abs(T - tq)))
            if abs(T[i] - tq) > 1.0 and not NEAREST:
                print(f"   no snapshot at {tq:g} s (the snapshots run from {T[0]:.0f} to {T[-1]:.0f} s)")
                continue
            r = rows[i][1]
            tp = min(ptimes, key=lambda t: abs(t - tq))
            pl = ERRATUM_FIG3_LINES[tp] if abs(tp - tq) <= 1.0 else {}
            po = ERRATUM_FIG3_OVERLYING.get(tp, {}) if abs(tp - tq) <= 1.0 else {}
            pa = ERRATUM_FIG3_AXIS.get(tp, np.nan) if abs(tp - tq) <= 1.0 else np.nan
            print(f"   t = {T[i]:.0f} s: axis {r['axis']:.2f} (erratum {pa:.2f}); top of the emerged flux {r['top']:.2f}"
                  f"{' (+ a rope in the sheet)' if r['ntop'] > 1 else ''}; reconnected there {820.49 - r['Atop']:.0f} G Mm "
                  f"({(820.49 - r['Atop']) / PAPER_DA:.1f} lines); between axis and top {r['Atop'] - r['Aaxis']:.0f} G Mm")
            print("     line    A_z   emerged flux: run  erratum   diff  |  overlying: run  erratum   diff")
            for k in range(1, 18):
                hr, hp = r["dome"].get(k, np.nan), pl.get(k, np.nan)
                orr, op = r["over"].get(k, np.nan), po.get(k, np.nan)
                if not (np.isfinite(hr) or np.isfinite(hp) or np.isfinite(orr) or np.isfinite(op)):
                    continue
                f = lambda v, w=7: f"{v:{w}.2f}" if np.isfinite(v) else " " * (w - 1) + "-"
                print(f"     {k:4d}  {PAPER_A0 - k * PAPER_DA:6.0f}             {f(hr)}  {f(hp)} {f(hr - hp, 6)}  |"
                      f"         {f(orr) if k <= 6 else '':>7}  {f(op) if k <= 6 else '':>7} {f(orr - op, 6) if k <= 6 else ''}")
            # the run minus the erratum, line by line
            ks = [k for k in range(1, 18) if np.isfinite(r["dome"].get(k, np.nan)) and k in pl]
            if ks:
                axs[1].plot(ks, [r["dome"][k] - pl[k] for k in ks], ls, marker="o", ms=3,
                            color=plt.get_cmap("plasma")(ptimes.index(tp) / len(ptimes)), label=f"{lab}, {tp} s")
        for k in range(1, 18):
            axs[0].plot(T, H[:, k - 1], ls, color=col(k), lw=1.0, label=lab if k == 1 else None)
        axs[0].plot(T, [r["axis"] for _, r in rows], ls, color="k", lw=1.4)
    for tp in ptimes:
        for k, h in ERRATUM_FIG3_LINES[tp].items():
            axs[0].plot(tp, h, "s", color=col(k), mfc="none", ms=5)
        axs[0].plot(tp, ERRATUM_FIG3_AXIS[tp], "ks", ms=5)
    axs[0].plot([], [], "s", color="0.4", mfc="none", ms=5, label="Erratum Fig. 3 (same colours)")
    axs[0].plot([], [], "ks", ms=5, label="tube axis (black)")
    axs[0].set_xlabel("Time (s)"); axs[0].set_ylabel("height on x = 0 of field lines 1 (top) to 17 (Mm)")
    axs[0].set_xlim(900, 2600); axs[0].set_ylim(0, 12); axs[0].grid(alpha=0.3); axs[0].legend(fontsize=7, loc="upper left")
    axs[1].axhline(0.0, color="k", lw=0.6)
    axs[1].set_xlabel("field line (1: next to the null; 17: next to the axis)"); axs[1].set_ylabel("run - erratum (Mm)")
    axs[1].grid(alpha=0.3)
    if axs[1].get_legend_handles_labels()[0]:
        axs[1].legend(fontsize=7)
    fig.tight_layout()
    out = fig_path(a, "lines")
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
    which the DynSGS viscous (rho mu |grad v|^2) and Ohmic (eta |curl B|^2) dissipation, estimated from the
    snapshot, would supply its internal energy if all of it were heat (JEXPRESSO_LK10_HEAT=1; with 'shock' most of
    it is not: `heat` gives what is). Plasma that is hot and in no hurry to drain is held up by its own pressure."""
    items = runs_of(a.run)
    rows = []
    p0 = M0 = None
    if a.ref:
        # the initial atmosphere and the unit of mass from the first snapshot of another run (the one this
        # run was restarted from, when its directory is not found through lk10_case.txt)
        t0, f0 = runs_of(a.ref)[0]
        g0 = load(f0, ["log10rho_cgs", "ρ", "p", "T", "u", "v", "w", "Bz", MU_NAME, *ETA_NAMES])
        p0 = grid2d(g0, g0["fields"]["p"])
        M0 = mass_row(g0, p0)["M"]
    for t, f in items:
        g = load(f, ["log10rho_cgs", "ρ", "p", "T", "u", "v", "w", "Bz", MU_NAME, *ETA_NAMES])
        if p0 is None:
            p0 = grid2d(g, g["fields"]["p"])
        rows.append((t, mass_row(g, p0)))
    whose = "the reference run's first snapshot" if M0 is not None else "the first snapshot"
    M0 = rows[0][1]["M"] if M0 is None else M0
    print(f"box 3 < y < 20 Mm, |x| < 13 Mm; masses in units of {whose}'s ({M0:.3g} g/cm)")
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


G_CODE = 0.98273        # gravity in code units, g_sun L0 / v0^2 (README, Units)


def _wmedian(x, w):
    """Weighted median of x (1D) with weights w; nan if there is no weight."""
    if x.size == 0 or not np.any(w > 0):
        return np.nan
    o = np.argsort(x); cw = np.cumsum(w[o])
    return float(x[o][min(np.searchsorted(cw, 0.5 * cw[-1]), x.size - 1)])


def drain_row(g, thr=-13.0, ymin=3.0, tfloor=0.1, rhofac=1.0, viscfac=1.0, vform="momentum"):
    """The dense plasma (log10 rho > thr) of the box ymin < y < 20 Mm, |x| < 13 Mm (cmd_drain)."""
    X, Y = g["X"], g["Y"]
    G = lambda n: grid2d(g, g["fields"][n])
    L = G("log10rho_cgs"); rho = 10.0 ** L
    d = ((Y > ymin) & (Y < 20.0))[:, None] & ((X > -13.0) & (X < 13.0))[None, :] & (L > thr)
    cm = 1e8
    integ = lambda f: np.trapezoid(np.trapezoid(f, X * cm, axis=1), Y * cm)
    M = integ(np.where(d, rho, 0.0))
    if not np.any(d) or M <= 0:
        return dict(M=0.0, floor=np.nan, cold=np.nan, warm=np.nan, beta=np.nan, flat=np.nan, vpar=np.nan, sup=np.nan,
                    vsup=np.nan, nsup=np.nan, ssup=np.nan, tsup=np.nan, mud=np.nan, mufl=np.nan, muoff=np.nan, T=np.nan)
    fr = lambda sel: integ(np.where(d & sel, rho, 0.0)) / M
    T = G("T")                                                    # in T_ph = 5700 K
    P, r = G("p"), G("ρ")                                         # code units
    Bx, By, Bz = G("Bx"), G("By"), G("Bz")
    Bp = np.sqrt(Bx * Bx + By * By)
    beta = P / np.maximum(0.5 * (Bp * Bp + Bz * Bz), 1e-300)
    sg = np.where(By < 0, -1.0, 1.0)                              # the in-plane unit vector, taken upward
    bx = sg * Bx / np.maximum(Bp, 1e-300); by = sg * By / np.maximum(Bp, 1e-300)
    incl = by > np.sin(np.radians(15.0))                          # field lines steeper than 15 degrees
    Xc, Yc = g["xc"], g["yc"]
    dpds = bx * _sem_d(P, Xc, 1) + by * _sem_d(P, Yc, 0)         # along the field, upward (element derivatives)
    sup = -dpds / np.maximum(r * G_CODE * by, 1e-300)             # 1: hydrostatic along the field; 0: nothing holds it
    u, v = G("u"), G("v")
    vdown = -(u * bx + v * by) * V0_KMS                           # km/s down the field line
    # the force of the diffusion terms on the VELOCITY along the field, rho dv/dt = R_rhov - v R_rho, as the same share:
    # the momentum slots diffuse rho v/rho_e with rho_e (mu + nu) (:dsgs_ref_weight; the coefficient written, after
    # JEXPRESSO_LK10_VISCFAC), the mass slot rho/rho_e with rho_e rho_factor mu. Split into the case's nu alone, and
    # the shear part div(rho (mu + nu) grad v), which a viscosity acting on the velocity would give: the rest is
    # momentum carried by the diffusion out of dense plasma moving as a block, with no mass to match it.
    if MU_NAME in g["fields"]:
        mu = G(MU_NAME)
        mud = np.maximum(mu - NU_CODE, 0.0)
        re = _ref_atmos(Yc, 1.3)[0][:, None] * np.ones_like(u)
        L2 = lambda c, q: _sem_dcd(c, q, Xc, 1) + _sem_dcd(c, q, Yc, 0)
        mur = (rhofac / viscfac) * mud                          # the mass slot's coefficient
        Rr = L2(re * mur, r / re)
        along = lambda fx, fy: (bx * fx + by * fy) / np.maximum(r * G_CODE * by, 1e-300)
        if vform.startswith("velocity") and _sem_ok(Xc) and _sem_ok(Yc):
            # JEXPRESSO_LK10_VISC=velocity (SGS.jl, _dsgs_velocity_viscosity_2d!): the momentum is diffused with the
            # mass slot's coefficient, the rest of the viscosity is a stress rho_e s_min (mu - mu_rho)+ grad v
            sr = r / re
            smin = np.array([[sr[j:j + 5, i:i + 5].min() for i in range(0, sr.shape[1] - 1, 4)]
                             for j in range(0, sr.shape[0] - 1, 4)])            # per element (nely, nelx)
            fx, fy = _el_factor(smin, Xc, Yc)
            S2 = lambda c, q: _sem_dcd(c, q, Xc, 1, cel=fx) + _sem_dcd(c, q, Yc, 0, cel=fy)
            ex = np.maximum(mu - mur, 0.0)
            vsup = along(L2(re * mur, r * u / re) - u * Rr + S2(re * ex, u), L2(re * mur, r * v / re) - v * Rr + S2(re * ex, v))
            nsup = along(S2(re * np.minimum(NU_CODE, ex), u), S2(re * np.minimum(NU_CODE, ex), v))
            ssup = along(L2(r * mur, u) + S2(re * ex, u), L2(r * mur, v) + S2(re * ex, v))
        else:
            vsup = along(L2(re * mu, r * u / re) - u * Rr, L2(re * mu, r * v / re) - v * Rr)
            nsup = along(L2(re * NU_CODE, r * u / re), L2(re * NU_CODE, r * v / re))
            ssup = along(L2(r * mu, u), L2(r * mu, v))
    else:
        vsup = nsup = ssup = mud = np.full_like(u, np.nan)
    # cell weights of the mass (for the medians)
    wx = np.gradient(X); wy = np.gradient(Y)
    w = rho * wy[:, None] * wx[None, :]
    k = d & incl
    return dict(M=M, floor=fr(T < 1.2 * tfloor), cold=fr(T * 5700.0 < 2e4), warm=fr(T * 5700.0 >= 2e4),
                beta=fr(beta > 1.0), flat=fr(~incl),
                vpar=float(np.sum(w[k] * vdown[k]) / np.sum(w[k])) if np.any(k) else np.nan,
                sup=_wmedian(sup[k], w[k]), vsup=_wmedian(vsup[k], w[k]), nsup=_wmedian(nsup[k], w[k]),
                ssup=_wmedian(ssup[k], w[k]), tsup=_wmedian((sup + vsup)[k], w[k]), mud=_wmedian(mud[d], w[d]),
                mufl=_wmedian(mud[d & (T < 1.2 * tfloor)], w[d & (T < 1.2 * tfloor)]),
                muoff=_wmedian(mud[d & (T >= 1.2 * tfloor)], w[d & (T >= 1.2 * tfloor)]),
                T=_wmedian(T[d] * 5700.0, w[d]))


def cmd_drain(a):
    """Why the plasma lifted into the corona drains, or does not. In the erratum's Fig. 3 the dome loses most of
    its plasma between 1610 and 1955 s, down its arched field lines, and only then rises to 12-13 Mm; a dome that
    keeps its plasma keeps its weight and stays where the field can carry it. Per snapshot, for the plasma denser
    than 10^thr g cm^-3 in the box ymin < y < 20 Mm, |x| < 13 Mm (fractions are of its mass):
      mass      in units of the first snapshot's whole box (as `mass`; --ref RUN0 for a restarted run);
      T med     its median temperature [K];  floor: within 20 % of the temperature floor;  <2e4, >=2e4 K;
      beta>1    in field too weak to guide it (the layer between the emerged field and the field above it);
      flat      on field lines within 15 degrees of the horizontal, where gravity has no component to drain it;
      v down    its mean speed down the field lines [km/s], on the steeper field lines (free fall from 10 Mm: 74);
      support   the median share of its weight along the field that its pressure gradient carries, on the same
                field lines: 1 is hydrostatic (it is held by its own pressure, it is too warm to drain), 0 is
                nothing in the way (it falls, or something that is not pressure holds it: viscosity);
      visc      the same share carried by the diffusion terms' force on the velocity (rho dv/dt = R_rhov - v R_rho,
                the DynSGS and the case's nu, with the code's operators; > 0: a drag on plasma that falls), and of it:
                nu, the case's explicit viscosity alone; shear, div(rho mu grad v), what a viscosity acting on the
                velocity would give (the rest is momentum the diffusion carries out of dense plasma moving as a
                block, with no mass diffusion to match it);  p+visc: pressure and visc, whose shortfall from 1 is
                what accelerates it down the field;
      mu        the median DynSGS momentum coefficient applied in this plasma [L0 v0 = 970 km^2/s] (after
                JEXPRESSO_LK10_VISCFAC; the case's nu, 0.0133, not included);  @floor, off: the same in its part at
                the temperature floor and in the rest (the floor resets the energy each step, which a sensor built
                on the residual of the equations reads as one).
    Nothing along a field line but pressure, gravity and viscosity acts on the plasma, so a dome that does not
    drain shows up here as one of: warm and support ~ 1; flat; beta > 1; or cold, steep, support ~ 0 and slow."""
    items = runs_of(a.run)
    M0 = None
    forms = {}
    if a.ref:
        g0 = load(runs_of(a.ref)[0][1], ["log10rho_cgs", "ρ", "p", "T", "u", "v", "Bz"])
        M0 = drain_row(g0, -99.0, a.ymin, a.tfloor)["M"]
    rows = []
    for t, f in items:
        g = load(f, ["log10rho_cgs", "ρ", "p", "T", "u", "v", "Bz", MU_NAME])
        if M0 is None:
            M0 = drain_row(g, -99.0, a.ymin, a.tfloor)["M"]
        d = os.path.dirname(f)
        rows.append((t, drain_row(g, a.thr, a.ymin, a.tfloor, float(_case_value(d, "rho_factor", "1.0")),
                                  float(_case_value(d, "visc_factor", "1.0")), _case_value(d, "visc_form", "momentum"))))
        forms[os.path.basename(os.path.normpath(d))] = _case_value(d, "visc_form", "momentum")
    print(f"plasma with log10 rho > {a.thr:g} in {a.ymin:g} < y < 20 Mm, |x| < 13 Mm; mass in units of the "
          f"{'reference run' if a.ref else 'first snapshot'}'s whole box ({M0:.3g} g/cm); fractions of its mass")
    if any(v.startswith("velocity") for v in forms.values()):
        print("   viscosity on the velocity (JEXPRESSO_LK10_VISC=velocity or nodal; evaluated with the s_min weight, a lower bound "
              "of the stress for nodal) in: " + ", ".join(k for k, v in forms.items() if v.startswith("velocity"))
              + "; visc, nu and shear are those of that operator there")
    print(f"{'t [s]':>7} {'mass':>8} | {'T med K':>8} {'floor':>6} {'<2e4 K':>6} {'>=2e4':>6} | {'beta>1':>6} {'flat':>6} | "
          f"{'v down':>7} {'support':>7} {'visc':>6} {'nu':>6} {'shear':>6} {'p+visc':>6} | {'mu':>7} {'@floor':>7} {'off':>7}")
    for t, r in rows:
        print(f"{t:7.0f} {r['M'] / M0:8.2f} | {r['T']:8.3g} {r['floor']:6.2f} {r['cold']:6.2f} {r['warm']:6.2f} | "
              f"{r['beta']:6.2f} {r['flat']:6.2f} | {r['vpar']:7.2f} {r['sup']:7.2f} {r['vsup']:6.2f} {r['nsup']:6.2f} "
              f"{r['ssup']:6.2f} {r['tsup']:6.2f} | {r['mud']:7.3g} {r['mufl']:7.3g} {r['muoff']:7.3g}")
    np.savetxt(os.path.join(a.outdir, "drain.csv"),
               np.array([[t, r["M"] / M0, r["T"], r["floor"], r["cold"], r["warm"], r["beta"], r["flat"], r["vpar"], r["sup"],
                          r["vsup"], r["nsup"], r["ssup"], r["tsup"], r["mud"], r["mufl"], r["muoff"]] for t, r in rows]),
               header="t [s], mass, median T [K], at the floor, < 2e4 K, >= 2e4 K, beta > 1, flat field, v down the field [km/s], "
                      "pressure support, viscous support, of it the case's nu, the shear part, pressure + viscous, "
                      "median DynSGS momentum coefficient [L0 v0], the same at the temperature floor, off it", delimiter=",")
    print("wrote", os.path.join(a.outdir, "drain.csv"))


GAMMA = 5.0 / 3.0
NU_CODE = 0.01326                       # the case's explicit kinematic viscosity above the photosphere (user_flux.jl, LK10_NU)
ETA0_CODE, JCRIT_CODE = 1.18e-4, 126.3  # its anomalous resistivity eta0 max(0, |curl B|/jcrit - 1) (user_flux.jl), code units
HEAT_TERMS = ["p dV", "nu", "DSGS heat", "anom eta", "E diff", "rho diff", "relax", "GLM psi"]
CH_DEFAULT = 274.655                    # the GLM speed of the tube + quadrupole case (its log: "GLM: c_h = ...")
TAU0_CODE, TAU_EXP, RELAX_RHOMIN = 23.0 / 23.1946, 1.43, 1.0e-5   # the energy relaxation (user_source.jl)
G_EXACT = 2.74e4 * 1.5e7 / (1.380649e-16 * 5700.0 / (1.125 * 1.67262192e-24))   # gravity, code (user_flux.jl, LK10_G)
LOST_TERMS = ["DSGS visc", "DSGS eta"]
LARE_VISC1, LARE_VISC2 = 0.1, 0.5       # LareXd's shock-viscosity coefficients (control.f90 of 2.6-2.9, the defaults)
LARE_DX_MM = 0.048                      # its cell in the figure window: the paper's grid, 0.0375 Mm at the surface (REF 4's mean node gap)


def _case_value(d, key, default):
    """The value of `key = value` in lk10_case.txt of the output directory d, or default."""
    f = os.path.join(d, "lk10_case.txt")
    if os.path.isfile(f):
        for line in open(f, encoding="utf-8", errors="replace"):
            k, s, v = line.partition("=")
            if s and k.strip() == key:
                return v.strip()
    return default


def _trapw(x):
    """Trapezoid weights of the nodes x (the lattice is not uniform)."""
    dx = np.diff(x)
    w = np.zeros_like(x, dtype=float)
    w[:-1] += 0.5 * dx
    w[1:] += 0.5 * dx
    return w


def _lgl(nop):
    """LGL nodes, weights and differentiation matrix D[i, j] = l_j'(x_i) of degree nop on [-1, 1]."""
    P = np.polynomial.legendre.Legendre.basis(nop)
    x = np.concatenate([[-1.0], np.sort(np.real(P.deriv().roots())), [1.0]])
    w = 2.0 / (nop * (nop + 1) * P(x) ** 2)
    c = np.array([np.prod([x[i] - x[j] for j in range(nop + 1) if j != i]) for i in range(nop + 1)])
    D = np.zeros((nop + 1, nop + 1))
    for i in range(nop + 1):
        for j in range(nop + 1):
            if i != j:
                D[i, j] = c[i] / (c[j] * (x[i] - x[j]))
        D[i, i] = -np.sum(D[i, :])
    return x, w, D


def _sem_ok(coords, nop=4):
    return (len(coords) - 1) % nop == 0 and len(coords) > nop


def _sem_d(f, coords, axis, nop=4):
    """d f/d(coord) along `axis` on the mesh's tensor-product LGL lattice (elements every nop nodes, affine), element
    by element with the LGL differentiation matrix and averaged at the nodes two elements share, as the code's
    element derivatives assembled (DSS); np.gradient on a lattice that is not one."""
    if not _sem_ok(coords, nop):
        return np.gradient(f, coords, axis=axis)
    _, _, D = _lgl(nop)
    fa = np.moveaxis(f, axis, -1)
    out = np.zeros_like(fa); cnt = np.zeros(fa.shape[-1])
    for e in range(0, fa.shape[-1] - 1, nop):
        sl = slice(e, e + nop + 1)
        out[..., sl] += fa[..., sl] @ D.T / (0.5 * (coords[e + nop] - coords[e]))
        cnt[sl] += 1.0
    return np.moveaxis(out / cnt, -1, axis)


def _sem_dcd(c, q, coords, axis, nop=4, cel=None):
    """d/dx (c dq/dx) along `axis` as the code's CG weak form: -sum_e D^T W (c D q) over the elements, assembled
    and divided by the lumped LGL mass. On this lattice of rectangles the 2D operator div(c grad q) is exactly the
    sum of this along x on every row and along y on every column. cel: a factor of c that is constant in each
    element, one value per line and per element along `axis` (_el_factor)."""
    if not _sem_ok(coords, nop):
        return np.gradient(c * np.gradient(q, coords, axis=axis), coords, axis=axis)
    _, w, D = _lgl(nop)
    qa = np.moveaxis(q, axis, -1); ca = np.moveaxis(c, axis, -1)
    ce = None if cel is None else np.moveaxis(cel, axis, -1)
    A = np.zeros_like(qa); m = np.zeros(qa.shape[-1])
    for e in range(0, qa.shape[-1] - 1, nop):
        sl = slice(e, e + nop + 1)
        h = 0.5 * (coords[e + nop] - coords[e])
        F = ca[..., sl] * (qa[..., sl] @ D.T / h) * w      # w_k c_k (dq/dx)_k
        if ce is not None:
            F = F * ce[..., e // nop][..., None]
        A[..., sl] -= F @ D                                 # -sum_k w_k c_k (dq/dx)_k D[k, i]  (dphi_i/dx = D[k,i]/h, J = h)
        m[sl] += w * h
    return np.moveaxis(A / m, -1, axis)


def _el_factor(f, xc, yc, nop=4):
    """For a factor that is constant in each element (f: (nely, nelx)), what the assembled operator sees on each
    lattice line: (per row and x-element, per y-element and column). A line shared by two elements gets their
    mean weighted by the element sizes across it, as the sum of the two elements' contributions over the lumped mass."""
    hy = np.diff(yc[::nop]); hx = np.diff(xc[::nop])
    def lines(g, h, n):                      # g: (nel_across, ...) -> (n_lines_across, ...)
        out = np.empty((n,) + g.shape[1:])
        for j in range(n):
            e = min(j // nop, g.shape[0] - 1)
            if j % nop == 0 and 0 < j < n - 1:
                out[j] = (h[e - 1] * g[e - 1] + h[e] * g[e]) / (h[e - 1] + h[e])
            else:
                out[j] = g[e]
        return out
    fx = lines(f, hy, len(yc))                               # (ny, nelx): the x-pass, row by row
    fy = lines(f.T, hx, len(xc)).T                           # (nely, nx): the y-pass, column by column
    return fx, fy


def _ref_atmos(yc, mcz):
    """The case's reference atmosphere at the heights yc [L0], as initialize.jl builds it (lk10_temperature,
    lk10_pressure): (rho_e, p_e) in code units. Above the photosphere T = 1 + (T_cor - 1)/2 [tanh((y - y_tr)/w_tr) + 1]
    in hydrostatic balance from p(0) = 1; below it the polytrope T = 1 - g y/(m + 1), p = T^(m + 1)."""
    g = G_EXACT
    tcor, ytr, wtr = 8.6e5 / 5700.0, 25.0, 5.0
    Tup = lambda y: 1.0 + 0.5 * (tcor - 1.0) * (np.tanh((y - ytr) / wtr) + 1.0)
    ymax = max(float(np.max(yc)), 0.002)
    ys = np.linspace(0.0, ymax, int(np.ceil(ymax / 0.002)) + 1)
    f = 1.0 / Tup(np.maximum(ys, 1e-300))
    lnp = np.concatenate([[0.0], -g * np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(ys))])
    up = yc > 0.0
    T = np.where(up, Tup(np.maximum(yc, 0.0)), 1.0 - g * np.minimum(yc, 0.0) / (mcz + 1.0))
    p = np.where(up, np.exp(np.interp(np.maximum(yc, 0.0), ys, lnp)), np.maximum(T, 1e-300) ** (mcz + 1.0))
    return p / T, p


def heat_fields(g, thr, ymin, ymax, rhofac, mcz, ch, relax=True, relax_ymin=-np.inf, mode="shock", lare_dx=LARE_DX_MM):
    """Every term of the temperature equation of the plasma of one snapshot, as the code applies it with
    JEXPRESSO_LK10_HEAT=mode (shock, viscous or a fraction), on the mesh's own lattice (cmd_heat). The DynSGS
    diffuses each slot q as div(rho_e mu grad((q - q_e)/rho_e)) (:dsgs_ref_weight, the case's user_primitives!), q_e
    the reference state (rho_e, p_e of the atmosphere, at rest), and its heat is SGS.jl's
    _dsgs_energy_consistent_2d!. Also what LareXd's default tensor shock viscosity would heat with this flow
    (lagran.F90 of 2.6-2.9, Q_MONO off): sigma_ij = L rho (visc1 c_f + visc2 L |s_ij|) s_ij, s the traceless
    strain rate, in every cell, heating by sigma_ij dv_i/dx_j; L its cell, lare_dx [Mm]."""
    X, Y = g["X"], g["Y"]
    G = lambda n: grid2d(g, g["fields"][n])
    L = G("log10rho_cgs")
    box = ((Y > ymin) & (Y < ymax))[:, None] & ((X > -13.0) & (X < 13.0))[None, :]
    dense = box & (L > thr)
    xc, yc = g["xc"], g["yc"]
    dx = lambda f: _sem_d(f, xc, 1)
    dy = lambda f: _sem_d(f, yc, 0)
    div2 = lambda c, f: _sem_dcd(c, f, xc, 1) + _sem_dcd(c, f, yc, 0)      # div(c grad f), the code's weak form
    r, p, T = G("ρ"), G("p"), G("T")
    u, v, w = G("u"), G("v"), G("w")
    Bx, By, Bz = G("Bx"), G("By"), G("Bz")
    if MU_NAME not in g["fields"]:
        raise SystemExit(f"heat: the snapshots have no {MU_NAME} (the DynSGS coefficient; the case writes it by default)")
    mom = G(MU_NAME)                                     # the momentum coefficient written: the DynSGS mu + the case's nu
    mud = np.maximum(mom - NU_CODE, 0.0)                 # the DynSGS mu, also that of rho (x rho_factor) and of e
    en = next((n for n in ETA_NAMES if n in g["fields"]), None)
    etaB = G(en) if en else np.zeros_like(u)            # the field's: the DynSGS eta + the anomalous one
    re1, pe1 = _ref_atmos(yc, mcz)
    re, pe = re1[:, None] * np.ones_like(u), pe1[:, None] * np.ones_like(u)
    e, ee = p / (GAMMA - 1.0), pe / (GAMMA - 1.0)
    div = dx(u) + dy(v)
    # grad u : grad(rho u / rho_e), as the code's dissipation (Dm, Dp of _dsgs_energy_consistent_2d!)
    gm = sum(dx(c) * dx(r * c / re) + dy(c) * dy(r * c / re) for c in (u, v, w))
    Dp = np.maximum(re * NU_CODE * gm, 0.0)              # the case's viscosity: heats in full
    Dm = np.maximum(re * (mud + NU_CODE) * gm, 0.0)      # everything the momentum diffusion dissipates
    shock = np.minimum(r * mud * np.minimum(div, 0.0) ** 2, np.maximum(Dm - Dp, 0.0))
    j2 = (dx(By) - dy(Bx)) ** 2 + dx(Bz) ** 2 + dy(Bz) ** 2
    eta_an = ETA0_CODE * np.maximum(np.sqrt(j2) / JCRIT_CODE - 1.0, 0.0)
    visc_all = np.maximum(Dm - Dp, 0.0)                  # all the DynSGS viscosity dissipates
    eta_dsgs = np.maximum(etaB - eta_an, 0.0) * j2       # and its resistivity
    if mode == "shock":
        applied, lost_v, lost_e = shock, visc_all - shock, eta_dsgs
    elif mode == "viscous":
        applied, lost_v, lost_e = visc_all, 0.0 * visc_all, eta_dsgs
    else:                                                # a fraction f of all of it (ψ diffusion not counted)
        fh = float(mode)
        applied, lost_v, lost_e = fh * (visc_all + eta_dsgs), (1.0 - fh) * visc_all, (1.0 - fh) * eta_dsgs
    # LareXd's default tensor shock viscosity on this flow (2.5D: no z derivatives; its qxx, qyy, qxy, qxz, qyz)
    Lc = lare_dx / L0_MM
    cf = np.sqrt(np.maximum(GAMMA * p + Bx ** 2 + By ** 2 + Bz ** 2, 0.0) / np.maximum(r, 1e-300))
    ux, uy, vx, vy, wx, wy = dx(u), dy(u), dx(v), dy(v), dx(w), dy(w)
    strains = ((0.5 * (uy + vx), uy + vx), (0.5 * wx, wx), (0.5 * wy, wy),
               ((2.0 * ux - vy) / 3.0, ux), ((2.0 * vy - ux) / 3.0, vy))
    lare = sum(Lc * r * (LARE_VISC1 * cf + LARE_VISC2 * Lc * np.abs(sij)) * sij * dij for sij, dij in strains)
    lare = np.maximum(lare, 0.0)
    mrn = next((n for n in ("mu_dsgs_ρ",) if n in g["fields"]), None)
    mur = G(mrn) if mrn else rhofac * mud                # the mass slot's coefficient (its own field where it has one:
    Rr = div2(re * mur, r / re)                          # JEXPRESSO_LK10_RHOFAC < 1, _RHOSENS), else rho_factor x mu
    # the energy relaxation to the internal energy the reference atmosphere has at the height of the same density
    yt = np.linspace(min(float(np.min(yc)), -200.0), float(np.max(yc)), 20001)
    rt, pt = _ref_atmos(yt, mcz)
    T0r = np.interp(np.log(np.maximum(r, 1e-300)), np.log(rt[::-1]), (pt / rt)[::-1])
    on = relax & (r > RELAX_RHOMIN) & (yc[:, None] > relax_ymin)
    qrel = np.where(on, -(np.maximum(r, 0.0) ** TAU_EXP / TAU0_CODE) * (p - r * T0r) / (GAMMA - 1.0), 0.0)
    # GLM: the damping of psi (user_source.jl) takes its energy 1/2 psi^2 out of psi but not out of E, so that
    # the pressure p = (gamma - 1)(E - ... - 1/2 psi^2) gains it
    dh = min(float(np.min(np.diff(xc))), float(np.min(np.diff(yc))))
    psi = G("ψ") if "ψ" in g["fields"] else np.zeros_like(u)
    q = {"p dV": -p * div,                               # compression heats, expansion cools: ideal
         "nu": Dp,
         "DSGS heat": applied,                           # the DynSGS heat as the run's mode applies it
         "anom eta": eta_an * j2,                        # the case's anomalous resistivity, in full
         "E diff": div2(re * mud, (e - ee) / re),        # the DynSGS diffusion of the internal energy departure
         "rho diff": -(e / r) * Rr,                      # it leaves e in place: T = p/rho rises where mass is taken away
         "relax": qrel,
         "GLM psi": (ch / dh) * psi * psi,
         "DSGS visc": lost_v,                            # dissipated, not heated
         "DSGS eta": lost_e,
         "Lare visc": lare}                              # what Lare's default tensor shock viscosity would heat
    W = _trapw(yc)[:, None] * _trapw(xc)[None, :]
    return dict(X=X, dense=dense, W=W, r=r, T=T, u=u, v=v, q=q, dx=dx, dy=dy, dh=dh,
                Mbox=np.sum(np.where(box, W * r, 0.0)))     # all the plasma of the box: the unit of mass, as drain's


KPER100 = (GAMMA - 1.0) * 5700.0 / T0_S * 100.0         # code dT/dt -> K per 100 s


def heat_means(F, sel):
    """Mass, mass-weighted mean T [K] and each term as a mass-weighted mean dT/dt [K per 100 s] over the nodes sel."""
    m = np.sum(np.where(sel, F["W"] * F["r"], 0.0))
    if m <= 0:
        return None
    out = {"M": m, "T": np.sum(np.where(sel, F["W"] * F["r"] * F["T"], 0.0)) / m * 5700.0}
    for n, f in F["q"].items():
        out[n] = np.sum(np.where(sel, F["W"] * f, 0.0)) / m * KPER100
    return out


def heat_lagrangian(F1, F2, t1, t2, sel1, sel2, tcut=2.0e5, jump=0.3):
    """The observed DT/Dt between two close snapshots, dT/dt + v.grad T at each node (fields averaged over the pair),
    against the budget (the mean of the two snapshots' terms), as mass-weighted means [K per 100 s] over the nodes
    that are in the set in both, below tcut [K] in both and whose T changed by less than `jump` of itself: a front
    that sweeps over a node (the hot layer above the dome) is not an advection a finite difference can follow.
    Returns (observed, budget, share of the set's mass compared)."""
    T1, T2 = F1["T"], F2["T"]
    Tm = 0.5 * (T1 + T2)
    c = sel1 & sel2 & (np.maximum(T1, T2) * 5700.0 < tcut) & (np.abs(T2 - T1) < jump * Tm)
    rm = 0.5 * (F1["r"] + F2["r"])
    W = F2["W"]
    m = np.sum(np.where(c, W * rm, 0.0))
    mall = np.sum(np.where(sel2, W * F2["r"], 0.0))
    if m <= 0 or mall <= 0:
        return np.nan, np.nan, 0.0
    dTdt = (T2 - T1) / ((t2 - t1) / T0_S)
    adv = 0.5 * (F1["u"] + F2["u"]) * F2["dx"](Tm) + 0.5 * (F1["v"] + F2["v"]) * F2["dy"](Tm)
    obs = np.sum(np.where(c, W * rm * (dTdt + adv), 0.0)) / m * 5700.0 / T0_S * 100.0
    qs = sum(0.5 * (F1["q"][n] + F2["q"][n]) for n in HEAT_TERMS)
    bud = np.sum(np.where(c, W * qs, 0.0)) / m * KPER100
    return obs, bud, m / mall


def cmd_heat(a):
    """What heats the dense plasma, from the snapshots. For the plasma denser than 10^thr g cm^-3 in the box
    ymin < y < 20 Mm, |x| < 13 Mm, split into the center (|x| < --xsplit Mm, where the dome's field lies flat) and
    the flanks, each term of its temperature equation as the code applies it with the run's JEXPRESSO_LK10_HEAT
    (lk10_case.txt), at the instant of each snapshot, as the mass-weighted mean dT/dt [K per 100 s]:
      p dV      compression (heats) and expansion (cools), -p div v: ideal, the paper's too;
      nu        the case's explicit viscosity (0.01326 L0 v0 above the photosphere), in full;
      DSGS heat the DynSGS heat the run applies: with heating = shock its viscosity where the flow is compressed,
                rho mu (div v)^2 for div v < 0; with viscous all its viscous dissipation; with a fraction f, f of
                all it dissipates;
      anom eta  the case's anomalous resistivity, eta |curl B|^2, in full (zero below |curl B| = j_crit);
      E diff    the DynSGS diffusion of the internal energy's departure from the reference atmosphere;
      rho diff  the DynSGS mass diffusion (rho_factor x mu, lk10_case.txt): it leaves the internal energy where it
                is, so T = p/rho rises where it takes mass away and falls where it brings it;
      relax     the energy relaxation (user_source.jl) to the reference's temperature at the same density;
      GLM psi   the GLM damping of psi, whose energy 1/2 psi^2 the pressure gains (alpha_p c_h/dh psi^2; --ch);
      sum       of these: DT/Dt of the plasma, but for the temperature floor (JEXPRESSO_LK10_TFLOOR), which a
                snapshot does not show.
    Then what the DynSGS dissipates and the run does not heat (heating = 1 would): its viscosity (outside
    compression, with shock), and its resistivity. Last, 'Lare visc': what LareXd's default tensor shock viscosity
    (the paper's code; lagran.F90 of 2.6-2.9, built without Q_MONO) would heat with the same flow: it acts in every
    cell, expanding ones too, sigma_ij = L rho (0.1 c_f + 0.5 L |s_ij|) s_ij on the traceless strain rate s, and all
    of what it dissipates is heat; L is its cell, --lare-dx Mm (0.048: the paper's grid in the figure window, a
    factor ~1.3 either way). Where two snapshots are close (the run's off-grid outputs, 12-17 s apart),
    the observed DT/Dt = dT/dt + v.grad T, node by node, is set against the mean of the two budgets on the nodes
    that stay in the set, below 2e5 K, and change by less than 30 % (fronts sweeping over a node excluded): what
    the budget misses (the floor, and its own errors) is their difference. The gradients are finite differences on
    the mesh's own LGL lattice, so each term is an estimate to a factor ~2: enough to tell which one heats."""
    if a.ymin < 0.0:
        raise SystemExit("heat: --ymin must be >= 0: below the photosphere the case's viscosity is dynamic (nu/rho_e)")
    items = runs_of(a.run)
    names = ["log10rho_cgs", "ρ", "p", "T", "u", "v", "w", "Bz", "ψ", MU_NAME, "mu_dsgs_ρ", *ETA_NAMES]
    seen, rows, pairs = {}, [], []
    prev = None
    for t, f in items:
        d = os.path.dirname(f)
        if d not in seen:
            seen[d] = (float(_case_value(d, "rho_factor", "1.0")) if a.rhofac is None else a.rhofac,
                       _case_value(d, "heating", "not recorded"), float(_case_value(d, "mcz", "1.3")),
                       _case_value(d, "relax", "true") == "true", float(_case_value(d, "relax_ymin", "-Inf")))
        rf, hm, mcz, rl, rly = seen[d]
        hm = "shock" if hm == "not recorded" else hm
        F = heat_fields(load(f, names), a.thr, a.ymin, a.ymax, rf, mcz, a.ch, rl, rly, hm, a.lare_dx)   # relax_ymin: in L0
        regs = {"center": F["dense"] & (np.abs(F["X"]) < a.xsplit)[None, :],
                "flanks": F["dense"] & (np.abs(F["X"]) >= a.xsplit)[None, :]}
        rows.append((t, {k: heat_means(F, s) for k, s in regs.items()}))
        M0 = F["Mbox"] if not rows[:-1] else M0
        if prev is not None and 0.0 < t - prev[0] <= a.pair and prev[1]["T"].shape == F["T"].shape:
            pairs.append((prev[0], t, {k: heat_lagrangian(prev[1], F, prev[0], t, prev[2][k], s) for k, s in regs.items()}))
        prev = ({k: F[k] for k in ("T", "u", "v", "r", "q", "dx", "dy", "W")}, regs)
        prev = (t, prev[0], prev[1])
    if not rows:
        raise SystemExit("heat: no snapshots")
    print(f"plasma with log10 rho > {a.thr:g} in {a.ymin:g} < y < {a.ymax:g} Mm, |x| < 13 Mm: each term of DT/Dt as applied "
          f"with the run's heating, mass-weighted mean [K per 100 s]; GLM: c_h = {a.ch:g}, dh = {F['dh']:.4g} (lattice); "
          f"Lare visc with its cell {a.lare_dx:g} Mm")
    for d in seen:
        if _case_value(d, "visc_form", "momentum").startswith("velocity"):
            print(f"   {os.path.basename(os.path.normpath(d))}: viscosity on the velocity (JEXPRESSO_LK10_VISC=velocity); the nu and DSGS "
                  f"heat columns below are those of the momentum form, an upper bound there (rho for rho_e s_min)")
    for d, (rf, h, *_) in seen.items():
        print(f"   {os.path.basename(os.path.normpath(d))}: heating = {h}"
              + (" (read as shock)" if h == "not recorded" else "") + f", rho_factor = {rf:g}")
    hdr = (f"{'t [s]':>7} {'mass':>7} {'T K':>8} | " + " ".join(f"{n:>8}" for n in HEAT_TERMS)
           + f" {'sum':>8} | " + " ".join(f"{n:>9}" for n in LOST_TERMS) + f" | {'Lare visc':>9}")
    csv = []
    for ir, (reg, label) in enumerate((("center", f"|x| < {a.xsplit:g} Mm"), ("flanks", f"{a.xsplit:g} <= |x| < 13 Mm"))):
        print(f"--- {reg}: {label}; mass in units of the first snapshot's whole box (as drain's: center + flanks = its mass)")
        print(hdr)
        for t, r in rows:
            x = r[reg]
            if x is None:
                continue
            s = sum(x[n] for n in HEAT_TERMS)
            print(f"{t:7.0f} {x['M'] / M0:7.3g} {x['T']:8.3g} | " + " ".join(f"{x[n]:8.3g}" for n in HEAT_TERMS)
                  + f" {s:8.3g} | " + " ".join(f"{x[n]:9.3g}" for n in LOST_TERMS) + f" | {x['Lare visc']:9.3g}")
            csv.append([ir, t, x["M"], x["T"], *(x[n] for n in HEAT_TERMS), s, *(x[n] for n in LOST_TERMS), x["Lare visc"]])
        for t1, t2, ob in pairs:
            o, b, fr = ob[reg]
            if not fr > 0:
                continue
            print(f"   {t1:.0f}-{t2:.0f} s, on the {100 * fr:.0f} % of the mass that changed smoothly: observed DT/Dt {o:9.3g}, "
                  f"budget {b:9.3g}, not in the budget {o - b:9.3g}")
    out = os.path.join(a.outdir, "heat.csv")
    np.savetxt(out, np.array(csv), delimiter=",",
               header="region (0 center, 1 flanks), t [s], mass [code], mean T [K], "
                      + ", ".join(HEAT_TERMS) + ", sum, " + ", ".join(LOST_TERMS) + ", Lare visc (K per 100 s)")
    print("wrote", out)


B0_GAUSS = 1030.6       # code field unit, G
P_PH_CGS = 10 ** 4.927  # code pressure unit, dyn cm^-2


def cmd_column(a):
    """The state on a vertical line, as a table: at each requested time (a snapshot the run has; --nearest for the
    nearest one) and every --dy Mm of height between --ymin and --ymax, at the mesh column nearest x = --x:
    log10 rho [g cm^-3], T [K], p [dyn cm^-2], the in-plane field and B_z [G], beta, v_x and v_y [km/s], the
    angle of the in-plane field to the horizontal [deg], and the DynSGS coefficients of the snapshot for the
    momentum and for the field [L0 v0 = 970 km^2/s]. Where the dense plasma is, how warm it is and what field it
    sits in, in numbers that can be pasted."""
    items = runs_of(a.run)
    sel = nearest(items, a.times)
    for t, f in sel:
        g = load(f, ["log10rho_cgs", "p", "T", "u", "v", "Bz", MU_NAME, *ETA_NAMES])
        X, Y = g["X"], g["Y"]
        G = lambda n: grid2d(g, g["fields"][n])
        i = int(np.argmin(np.abs(X - a.x)))
        L, T, P = G("log10rho_cgs")[:, i], G("T")[:, i] * 5700.0, G("p")[:, i]
        Bx, By, Bz = G("Bx")[:, i], G("By")[:, i], G("Bz")[:, i]
        u, v = G("u")[:, i] * V0_KMS, G("v")[:, i] * V0_KMS
        mu = G(MU_NAME)[:, i] if MU_NAME in g["fields"] else np.full(len(Y), np.nan)
        en = next((n for n in ETA_NAMES if n in g["fields"]), None)
        eta = G(en)[:, i] if en else np.full(len(Y), np.nan)
        Bp = np.hypot(Bx, By)
        beta = P / np.maximum(0.5 * (Bp ** 2 + Bz ** 2), 1e-300)
        ang = np.degrees(np.arctan2(np.abs(By), np.abs(Bx)))
        print(f"--- t = {t:.0f} s, x = {X[i]:.2f} Mm ({os.path.basename(f)})")
        print(f"{'y Mm':>6} {'log rho':>8} {'T K':>9} {'p cgs':>9} {'|Bp| G':>8} {'Bx G':>8} {'Bz G':>8} {'beta':>9} "
              f"{'vx km/s':>8} {'vy km/s':>8} {'angle':>6} {'mu':>8} {'eta':>8}")
        for yq in np.arange(a.ymin, a.ymax + 1e-9, a.dy):
            j = int(np.argmin(np.abs(Y - yq)))
            print(f"{Y[j]:6.2f} {L[j]:8.2f} {T[j]:9.3g} {P[j] * P_PH_CGS:9.3g} {Bp[j] * B0_GAUSS:8.2f} {Bx[j] * B0_GAUSS:8.2f} "
                  f"{Bz[j] * B0_GAUSS:8.2f} {beta[j]:9.3g} {u[j]:8.2f} {v[j]:8.2f} {ang[j]:6.1f} {mu[j]:8.3g} {eta[j]:8.3g}")


def cmd_fastest(a):
    """Where the fastest wave is: at each requested time (default: the last three snapshots of the run), the --n
    places with the largest |v| + c_f, at least --sep Mm apart, with c_f = sqrt((gamma p + B^2)/rho) as the solver
    computes it. Printed for each: position [Mm], |v| and c_f [km/s], (|v| + c_f)/v0, log10 rho [g cm^-3], |B| [G],
    T [K], beta. A run whose step collapses (JEXPRESSO_LK10_CFL) or that blows up names its region here."""
    items = runs_of(a.run)
    sel = nearest(items, a.times) if a.times else items[-3:]
    for t, f in sel:
        g = load(f, ["ρ", "p", "u", "v", "w", "Bz", "T", "log10rho_cgs"])
        F = g["fields"]
        rho = np.maximum(F["ρ"], 1e-300)
        B2 = F["Bx"] ** 2 + F["By"] ** 2 + F["Bz"] ** 2
        vel = np.sqrt(F["u"] ** 2 + F["v"] ** 2 + (F["w"] ** 2 if "w" in F else 0.0))
        cf = np.sqrt((5.0 / 3.0 * np.maximum(F["p"], 0.0) + B2) / rho)
        spd = vel + cf
        x, y = g["x"] * L0_MM, g["y"] * L0_MM
        picked = []
        for i in np.argsort(-spd):
            if all(np.hypot(x[i] - x[j], y[i] - y[j]) >= a.sep for j in picked):
                picked.append(i)
            if len(picked) == a.n:
                break
        print(f"--- t = {t:.1f} s ({os.path.basename(f)}): the {len(picked)} fastest places, {a.sep:g} Mm apart or more")
        print(f"{'x Mm':>8} {'y Mm':>8} {'|v| km/s':>10} {'c_f km/s':>10} {'(|v|+c_f)/v0':>13} {'log rho':>8} "
              f"{'|B| G':>9} {'T K':>9} {'beta':>9}")
        for i in picked:
            beta = F["p"][i] / max(0.5 * B2[i], 1e-300)
            print(f"{x[i]:8.2f} {y[i]:8.2f} {vel[i] * V0_KMS:10.4g} {cf[i] * V0_KMS:10.4g} {spd[i]:13.4g} "
                  f"{F['log10rho_cgs'][i]:8.2f} {np.sqrt(B2[i]) * B0_GAUSS:9.4g} {F['T'][i] * 5700.0:9.3g} {beta:9.3g}")


def cmd_massbudget(a):
    """The mass budget node by node on a vertical line, for the stripes: at the snapshot nearest --time and on the mesh
    column nearest each --x, for every node between --ylim (an element face is marked |):
      s, s next   rho/rho_e at that snapshot and at the next one;
      -d(ru)/dx, -d(rv)/dy   the advection of mass along x and along y, with the code's element derivatives
                  averaged at the faces two elements share;
      diffusion   the DynSGS mass diffusion, div(rho_e mu_rho grad(rho/rho_e)), with the mass slot's coefficient
                  of the snapshot;
      sum, observed   their sum, and the change of ln rho to the next snapshot;
    all as rates of ln rho per 100 s (-100: the node would lose a factor e of its density in 100 s), then the
    node-to-node part of log10(rho/rho_e) of the node's element along y [dex] (what JEXPRESSO_LK10_RHOSENS reads) and
    the mass slot's coefficient. A stripe in the making is one row of nodes, at an element face, whose mass leaves
    through -d(rv)/dy alone while the rows 50 km above and below it keep theirs, and to which the diffusion returns
    a fraction (README, "The density sensor"). The next snapshot should be close (the off-grid outputs, 12-18 s
    after 1380, 1610, 1955, 2070, 2185 s): over a longer interval the budget of one instant is not the change."""
    items = runs_of(a.run)
    k = int(np.argmin([abs(t - a.time) for t, _ in items]))
    if k + 1 >= len(items):
        raise SystemExit(f"massbudget: the snapshot nearest {a.time:g} s is the last one; give an earlier time")
    (t1, f1), (t2, f2) = items[k], items[k + 1]
    d0 = os.path.dirname(os.path.abspath(f1))
    mcz = float(_case_value(d0, "mcz", "1.3"))
    g0 = read_pvtu(f1, None)
    mus = sorted(n for n in g0["fields"] if n.startswith("mu_dsgs"))
    mu_m = next((n for n in mus if "ρu" in n), mus[0] if mus else None)
    mu_r = next((n for n in mus if n == "mu_dsgs" or n.startswith("mu_dsgs_ρ_") or n == "mu_dsgs_ρ"), mu_m)
    g = load(f1, ["ρ", "u", "v"] + ([mu_r] if mu_r else []))
    g2 = load(f2, ["ρ"])
    G = lambda n: grid2d(g, g["fields"][n])
    r, u, v = G("ρ"), G("u"), G("v")
    mur = G(mu_r) if mu_r else np.zeros_like(r)
    r2 = grid2d(g2, g2["fields"]["ρ"])
    xc, yc, X, Y = g["xc"], g["yc"], g["X"], g["Y"]
    if r2.shape != r.shape or not (_sem_ok(xc) and _sem_ok(yc)):
        raise SystemExit("massbudget: needs two snapshots on the same tensor-product LGL lattice")
    re = _ref_atmos(yc, mcz)[0][:, None] * np.ones_like(r)
    per100 = 100.0 / T0_S
    ax = -_sem_d(r * u, xc, 1) / r * per100
    ay = -_sem_d(r * v, yc, 0) / r * per100
    df = (_sem_dcd(re * mur, r / re, xc, 1) + _sem_dcd(re * mur, r / re, yc, 0)) / r * per100
    obs = np.log(np.maximum(r2, 1e-300) / r) / (t2 - t1) * 100.0
    hy = _highmode(np.log10(r / re), axis=0)
    print(f"--- t = {t1:.0f} s ({os.path.basename(f1)}) against the change to {t2:.0f} s ({os.path.basename(f2)}, {t2 - t1:.0f} s later); "
          f"rates of ln rho per 100 s; mass coefficient: {mu_r}")
    for xq in a.x:
        i = int(np.argmin(np.abs(X - xq)))
        print(f"   x = {X[i]:.2f} Mm (node column {i}: {'an element face' if i % 4 == 0 else 'inside an element'})")
        print(f"   {'y Mm':>6} {'s':>9} {'s next':>9} | {'-d(ru)/dx':>9} {'-d(rv)/dy':>9} {'diffusion':>9} {'sum':>8} | {'observed':>8} | "
              f"{'n2n dex':>7} {'mu_rho':>7}")
        for j in range(len(Y)):
            if a.ylim[0] <= Y[j] <= a.ylim[1]:
                print(f" {'|' if j % 4 == 0 else ' '} {Y[j]:6.3f} {r[j, i] / re[j, i]:9.3f} {r2[j, i] / re[j, i]:9.3f} | {ax[j, i]:9.1f} {ay[j, i]:9.1f} "
                      f"{df[j, i]:9.1f} {ax[j, i] + ay[j, i] + df[j, i]:8.1f} | {obs[j, i]:8.1f} | {hy[j, i]:7.3f} {mur[j, i]:7.4f}")


# The erratum's Fig. 3 (Leake et al. 2011, ApJ 741:125) by height: its colours turned back into log10 rho
# [g cm^-3] through its colour bar (white field lines and arrows left out), median over the centre (|x| <= 3 Mm)
# and over the flanks (3 < |x| <= 8 Mm) at each height of ERRATUM_Y [Mm]; nan where too few pixels are not
# white. Read off the published figure (a factor ~0.1 dex in a smooth region); for cmd_profile.
ERRATUM_Y = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0]
ERRATUM_PROFILE = {
    920: ([-10.78, -12.44, -13.41, -14.32, -14.32, -14.32, -14.32, -14.32, -14.31, -14.32, -14.32, -14.32, -14.32, -14.32],
          [-11.62, -13.02, -13.91, -14.20, -14.21, -14.32, -14.32, -14.32, -14.33, -14.33, -14.43, -14.32, -14.32, -14.32]),
    1380: ([np.nan, -10.52, -10.62, -10.52, -10.19, -10.95, -14.42, -14.32, -14.32, -14.32, -14.32, np.nan, -14.32, -14.32],
           [-9.95, -10.51, -10.41, -10.40, -11.62, -14.32, -14.32, -14.32, -14.32, -14.32, -14.32, -14.32, -14.32, -14.32]),
    1610: ([np.nan, -11.20, -11.38, -11.41, -11.43, -11.36, -11.14, -10.96, -11.40, -14.43, -14.32, -14.32, -14.32, -14.32],
           [-10.32, -10.62, -10.65, -11.20, -11.30, -11.13, -10.99, -11.23, -14.43, -14.43, -14.43, -14.35, -14.32, -14.32]),
    1955: ([-10.74, -11.38, -12.27, -12.59, -12.62, -12.54, -12.44, -12.21, -11.99, -12.11, -12.37, -14.42, -14.36, -14.32],
           [-10.42, -10.72, -11.08, -11.38, -11.60, -11.59, -12.07, -12.21, -12.21, -13.04, -14.43, -14.32, -14.32, -14.32]),
    2185: ([-10.95, -11.59, -12.54, -12.70, -12.75, -13.28, -13.30, -13.19, -12.98, -12.72, -12.21, -12.46, -13.02, -14.04],
           [-10.65, -10.94, -11.28, -11.50, -11.79, -12.02, -12.19, -12.69, -13.00, -13.17, -14.43, -14.46, -14.43, -14.43]),
    2530: ([-10.94, -12.15, -13.11, -13.62, -13.90, -14.18, -14.42, -14.44, -14.48, -14.78, -14.74, -14.75, -14.62, -14.57],
           [-10.74, -11.19, -11.51, -11.85, -12.26, -12.93, -13.90, -14.22, -14.59, -14.81, -14.56, -14.60, -14.49, -14.43]),
}


def cmd_profile(a):
    """The density of the dome by height, against the erratum's Fig. 3. At each requested time, for bands 0.5 Mm
    thick around y = 2, 3, ..., 15 Mm: the run's median log10 rho [g cm^-3] over the centre (|x| <= 3 Mm) and over
    the flanks (3 < |x| <= 8 Mm), the erratum's at that height and time (ERRATUM_PROFILE: its colours turned back
    into log10 rho, medians over the same x ranges), their difference (run - erratum, dex; negative: the run is
    thinner), the run's 10th percentile (how far its thinnest nodes go), and the share of its nodes held at the
    density floor (0.2 p_e(y)/T_cor) and at the temperature floor (lk10_case.txt's tfloor). Whether, and where, the
    run's dome drains more than the paper's, with the floors that a thinner dome meets kept apart."""
    items = runs_of(a.run)
    sel = nearest(items, a.times)
    d0 = os.path.dirname(os.path.abspath(sel[0][1]))
    mcz = float(_case_value(d0, "mcz", "1.3"))
    tfl = float(_case_value(d0, "tfloor", "0.1"))
    TCOR = 8.6e5 / 5700.0
    for t, f in sel:
        g = load(f, ["log10rho_cgs", "ρ", "T"])
        X, Y = g["X"], g["Y"]
        G = lambda n: grid2d(g, g["fields"][n])
        L, R, T = G("log10rho_cgs"), G("ρ"), G("T")
        _, pe = _ref_atmos(Y / L0_MM, mcz)
        atR = R <= (0.2 * pe / TCOR)[:, None] * (1.0 + 1e-6)
        atT = T <= tfl * (1.0 + 1e-6)
        te = min(ERRATUM_PROFILE, key=lambda k: abs(k - t))
        ec, ef = ERRATUM_PROFILE[te] if abs(te - t) <= 1.0 else ([np.nan] * len(ERRATUM_Y),) * 2
        print(f"--- t = {t:.0f} s ({os.path.basename(f)}): median log10 rho [g cm^-3] by height"
              + (f"; erratum's Fig. 3 at {te} s" if abs(te - t) <= 1.0 else "; no panel of the erratum at this time"))
        print(f"{'y Mm':>5} | {'centre |x|<=3 Mm: run':>21} {'erratum':>8} {'diff':>6} {'p10':>7} {'rho fl':>7} {'T fl':>6}"
              f" | {'flanks 3<|x|<=8: run':>20} {'erratum':>8} {'diff':>6} {'p10':>7} {'rho fl':>7} {'T fl':>6}")
        for k, yq in enumerate(ERRATUM_Y):
            band = np.abs(Y - yq) <= 0.25
            cells = []
            for msk, ev in (((np.abs(X) <= 3.0), ec[k]), ((np.abs(X) > 3.0) & (np.abs(X) <= 8.0), ef[k])):
                Wm = np.ix_(band, msk)
                v = L[Wm]
                if v.size == 0:
                    cells.append(f"{'':>21} {'':>8} {'':>6} {'':>7} {'':>7} {'':>6}")
                    continue
                m = float(np.median(v))
                cells.append(f"{m:21.2f} {ev:8.2f} {m - ev:6.2f} {float(np.percentile(v, 10)):7.2f} "
                             f"{atR[Wm].mean():7.1%} {atT[Wm].mean():6.1%}")
            print(f"{yq:5.1f} | {cells[0]} | {cells[1]}")


def _highmode(F, nop=4, axis=0):
    """Amplitude of the highest Legendre mode (degree nop) of F within each element along `axis`, on the tensor
    lattice of LK10.geo (elements every nop nodes): the node-to-node part of a field, which the element cannot
    tell from noise. Returns an array of F's shape with each element's value at its nodes (the larger of the
    two estimates at a shared face)."""
    x, _, _ = _lgl(nop)
    V = np.polynomial.legendre.legvander(x, nop)         # V[i, k] = P_k(x_i)
    Vi = np.linalg.inv(V)
    G = np.moveaxis(np.asarray(F, float), axis, 0)
    n = G.shape[0]
    out = np.zeros_like(G)
    for e in range(0, n - 1, nop):
        c = np.tensordot(Vi[nop], G[e:e + nop + 1], axes=(0, 0))    # coefficient of P_nop, per line
        out[e:e + nop + 1] = np.maximum(out[e:e + nop + 1], np.abs(c)[None, ...])
    return np.moveaxis(out, 0, axis)


def cmd_noise(a):
    """Stripes and grain in a map: at each time, the node-to-node part of log10 rho (the amplitude of the highest
    Legendre mode in each element, along y and along x, in dex), in the window (--xlim, --ylim), and the nodes held
    at the positivity floors (rho at 0.2 p_e/T_cor, T at 0.1 T_ph). Prints, per row of elements, the largest and
    the mean amplitude and the share of nodes at a floor, at either and at each; a figure with the density, the two
    amplitudes and the floor nodes (grey: the temperature floor, black: the density floor). Horizontal stripes that
    are a property of the discretization show as a y amplitude that is large where the stripes are and small along
    x; plasma at the density floor as floor nodes; real layers of plasma as neither. Rows of nodes at the
    temperature floor one element apart in height, over a smooth density, are the pressure of a drained low-beta
    region set by the aliasing of the total-energy equation (README, Positivity: The internal energy where it is a
    small share)."""
    items = runs_of(a.run)
    sel = nearest(items, a.times)
    d0 = os.path.dirname(os.path.abspath(sel[0][1]))
    mcz = float(_case_value(d0, "mcz", "1.3"))
    tfl = float(_case_value(d0, "tfloor", "0.1"))
    TCOR = 8.6e5 / 5700.0
    fig, axs = plt.subplots(len(sel), 4, figsize=(16, 3.6 * len(sel)), squeeze=False)
    for r, (t, f) in enumerate(sel):
        g = load(f, ["log10rho_cgs", "ρ", "T"])
        X, Y = g["X"], g["Y"]
        G = lambda n: grid2d(g, g["fields"][n])
        L, R, T = G("log10rho_cgs"), G("ρ"), G("T")
        ny_, nx_ = L.shape
        hy = _highmode(L, axis=0); hx = _highmode(L, axis=1)
        _, pe = _ref_atmos(Y / L0_MM, mcz)
        rfl = 0.2 * pe / TCOR                                # the density floor of the run, per height (code units)
        atR = R <= rfl[:, None] * (1.0 + 1e-6)              # held at the density floor
        atT = T <= tfl * (1.0 + 1e-6)                        # held at the temperature floor
        atfl = atR | atT
        mx = (X >= a.xlim[0]) & (X <= a.xlim[1]); my = (Y >= a.ylim[0]) & (Y <= a.ylim[1])
        W = np.ix_(my, mx)
        print(f"--- t = {t:.0f} s ({os.path.basename(f)}), x in {a.xlim}, y in {a.ylim} Mm: node-to-node part of log10 rho "
              f"[dex], along y and along x, and the share of nodes at a floor")
        print(f"   whole window: along y max {hy[W].max():.3f} mean {hy[W].mean():.4f}; along x max {hx[W].max():.3f} "
              f"mean {hx[W].mean():.4f}; at a floor {atfl[W].mean():.3%} ({atfl[W].sum()} nodes: "
              f"{atT[W].sum()} at the temperature floor, {atR[W].sum()} at the density floor)")
        print(f"{'y Mm':>12} {'y: max':>8} {'mean':>8} {'x: max':>8} {'mean':>8} {'at floor':>9} {'T floor':>8} {'rho floor':>9}")
        js = np.where(my)[0]
        for e0 in range(js[0] - js[0] % 4, js[-1], 4):
            rows = slice(e0, min(e0 + 5, ny_))
            sub = (rows, mx)
            print(f"{Y[e0]:5.2f}-{Y[min(e0 + 4, ny_ - 1)]:5.2f} {hy[sub].max():8.3f} {hy[sub].mean():8.4f} "
                  f"{hx[sub].max():8.3f} {hx[sub].mean():8.4f} {atfl[sub].mean():9.2%} {atT[sub].mean():8.2%} {atR[sub].mean():9.2%}")
        ext = dict(xlim=a.xlim, ylim=a.ylim)
        for k, (Z, ttl, cm, lim) in enumerate(((L, r"log$_{10}\rho$", "jet", (-16, -5.3)),
                                                (hy, "node-to-node part along y (dex)", "magma", (0, a.vmax)),
                                                (hx, "node-to-node part along x (dex)", "magma", (0, a.vmax)),
                                                (np.where(atR, 1.0, np.where(atT, 0.5, 0.0)),
                                                 "nodes at a floor: grey T, black rho", "Greys", (0, 1)))):
            ax = axs[r, k]
            rs, cs = _window(g, a.xlim, a.ylim)
            im = ax.pcolormesh(X[cs], Y[rs], Z[rs, cs], cmap=cm, vmin=lim[0], vmax=lim[1], shading="nearest", rasterized=True)
            ax.set_xlim(a.xlim); ax.set_ylim(a.ylim); ax.set_aspect("equal")
            ax.set_title(f"{ttl}, t = {t:.0f} s", fontsize=9); fig.colorbar(im, ax=ax, shrink=0.8)
    fig.tight_layout()
    out = fig_path(a, "noise")
    fig.savefig(out, dpi=a.dpi); plt.close(fig)
    print("wrote", out)


def cmd_elements(a):
    """Which quantity carries a pattern that repeats element by element (the stripes `noise` finds): along the
    vertical line of nodes nearest each --x, for every element it crosses within --ylim, the element's smallest
    rho/rho_e (s_min, which sets the viscosity of the velocity form, JEXPRESSO_LK10_VISC=velocity: the stress
    rho_e s_min (mu - mu_rho)), its DynSGS coefficients (as the run wrote them, one value per element), that
    stress over rho_e, s_min (mu - mu_rho), and its ratio to the element's below (a jump between neighbouring
    elements, which the velocity form's stress has where the plasma drains), and on the line the node-to-node part (the highest
    Legendre mode in the element) of log10 rho, log10 T and the two velocity components, and the nodes of the
    element at the density and temperature floors."""
    items = runs_of(a.run)
    sel = nearest(items, a.times)
    d0 = os.path.dirname(os.path.abspath(sel[0][1]))
    mcz = float(_case_value(d0, "mcz", "1.3"))
    tfl = float(_case_value(d0, "tfloor", "0.1"))
    TCOR = 8.6e5 / 5700.0
    V0 = L0_MM * 1.0e3 / 23.1946                                  # code velocity unit, km/s
    for t, f in sel:
        g0 = read_pvtu(f, None)
        mus = sorted(k for k in g0["fields"] if k.startswith("mu_dsgs"))
        mu_m = next((k for k in mus if "ρu" in k), mus[0] if mus else None)          # the momentum slots'
        mu_r = next((k for k in mus if k == "mu_dsgs" or k.startswith("mu_dsgs_ρ_") or k == "mu_dsgs_ρ"), mu_m)
        g = load(f, ["ρ", "T", "u", "v"] + ([mu_m] if mu_m else []) + ([mu_r] if mu_r and mu_r != mu_m else []))
        X, Y = g["X"], g["Y"]
        G = lambda n: grid2d(g, g["fields"][n])
        R, T, U, W = G("ρ"), G("T"), G("u") * V0, G("v") * V0
        MM = G(mu_m) if mu_m else np.zeros_like(R)
        MR = G(mu_r) if mu_r else np.zeros_like(R)
        rho_e, pe = _ref_atmos(g["yc"], mcz)
        S = R / rho_e[:, None]
        atR = R <= 0.2 * pe[:, None] / TCOR * (1.0 + 1e-6)
        atT = T <= tfl * (1.0 + 1e-6)
        LR, LT = np.log10(np.maximum(R, 1e-300)), np.log10(np.maximum(T, 1e-300))
        h = np.linalg.inv(np.polynomial.legendre.legvander(_lgl(4)[0], 4))[4]       # P_4 coefficient from 5 nodes
        print(f"--- t = {t:.0f} s ({os.path.basename(f)}); DynSGS fields: viscosity {mu_m}, mass {mu_r}; "
              f"velocities km/s, node-to-node parts on the line, logs in dex")
        for x0 in a.x:
            i0 = int(np.argmin(np.abs(X - x0)))
            ic = min(i0 - i0 % 4, len(X) - 5)
            print(f"   x = {X[i0]:.2f} Mm (element columns {X[ic]:.2f} to {X[ic + 4]:.2f} Mm)")
            print(f"{'y Mm':>13} {'s_min':>9} {'mu':>8} {'mu_rho':>8} {'s(mu-mu_r)':>10} {'/below':>8} {'log rho':>8} "
                  f"{'log T':>7} {'u':>7} {'v':>7} {'rho fl':>6} {'T fl':>5}")
            js = np.where((Y >= a.ylim[0]) & (Y <= a.ylim[1]))[0]
            cprev = None
            for e in range(js[0] - js[0] % 4, js[-1] - 3, 4):
                blk = (slice(e, e + 5), slice(ic, ic + 5))
                smin = S[blk].min()
                jm, im = e + 2, ic + 2                               # an interior node: the element's own values
                mu, mur = MM[jm, im], MR[jm, im]
                c = smin * max(mu - mur, 0.0)                       # the stress over rho_e: no stratification in it
                rat = c / cprev if (cprev is not None and cprev > 0) else np.nan
                cprev = c
                col = slice(e, e + 5)
                hm = lambda F: abs(float(h @ F[col, i0]))
                print(f"{Y[e]:6.2f}-{Y[e + 4]:5.2f} {smin:9.3g} {mu:8.3g} {mur:8.3g} {c:10.3g} {rat:8.3g} {hm(LR):8.3f} "
                      f"{hm(LT):7.3f} {hm(U):7.2f} {hm(W):7.2f} {int(atR[blk].sum()):6d} {int(atT[blk].sum()):5d}")


def find_opoints(A, X, Y, xlim=(-12.75, 12.75), ylim=(3.0, 19.5), k=6):
    """O-points of the in-plane field in a window: nodes where the flux function A (grid (ny, nx)) is the largest or
    the smallest of its (2k+1)^2 neighbours (an X-point is a saddle of A and is not one). Returns
    [(x, y, depth)], deepest first; depth = |A - mean of A on the square 2k nodes away|, the flux it encloses
    at that scale, which ranks a rope above a ripple."""
    ny, nx = A.shape
    hi = np.full_like(A, -np.inf); lo = np.full_like(A, np.inf)
    P = np.pad(A, k, mode="edge")
    for dj in range(2 * k + 1):
        for di in range(2 * k + 1):
            W = P[dj:dj + ny, di:di + nx]
            np.maximum(hi, W, out=hi); np.minimum(lo, W, out=lo)
    ext = (A >= hi) | (A <= lo)
    inwin = ((Y > ylim[0]) & (Y < ylim[1]))[:, None] & ((X > xlim[0]) & (X < xlim[1]))[None, :]
    out = []
    for j, i in zip(*np.where(ext & inwin)):
        j0, j1, i0, i1 = max(j - 2 * k, 0), min(j + 2 * k, ny - 1), max(i - 2 * k, 0), min(i + 2 * k, nx - 1)
        ring = np.concatenate([A[j0, i0:i1 + 1], A[j1, i0:i1 + 1], A[j0:j1 + 1, i0], A[j0:j1 + 1, i1]])
        out.append((float(X[i]), float(Y[j]), float(abs(A[j, i] - ring.mean()))))
    out.sort(key=lambda o: -o[2])
    # one entry per rope: drop those within 0.5 Mm of a deeper one
    keep = []
    for o in out:
        if all(np.hypot(o[0] - q[0], o[1] - q[1]) > 0.5 for q in keep):
            keep.append(o)
    return keep


FORCE_UNIT = 10 ** 4.927 / 1.5e7 / 1e-5      # code force density p_ph/L0, in 1e-5 dyn cm^-3
G_SUN = 2.74e4                               # cm s^-2


def vertical_forces(g, x):
    """The four vertical force densities of the erratum's Fig. 4 on the lattice column nearest x [Mm], in
    1e-5 dyn cm^-3 (heights Y [Mm]): gravity -rho g, magnetic tension (B.grad)B_y/4pi, magnetic pressure
    -d(B^2/8pi)/dy, gas pressure -dP/dy (2.5D: no z derivatives). Derivatives by differences on the LGL nodes."""
    G = lambda n: grid2d(g, g["fields"][n])
    Xc, Yc = g["xc"], g["yc"]
    Bx, By, Bz, P = G("Bx"), G("By"), G("Bz"), G("p")
    dx = lambda f: np.gradient(f, Xc, axis=1)
    dy = lambda f: np.gradient(f, Yc, axis=0)
    i = int(np.argmin(np.abs(g["X"] - x)))
    grav = -(10.0 ** G("log10rho_cgs")[:, i]) * G_SUN / 1e-5
    tens = (Bx * dx(By) + By * dy(By))[:, i] * FORCE_UNIT
    magp = -dy(0.5 * (Bx ** 2 + By ** 2 + Bz ** 2))[:, i] * FORCE_UNIT
    gasp = -dy(P)[:, i] * FORCE_UNIT
    return g["X"][i], g["Y"], grav, tens, magp, gasp


def cmd_forces(a):
    """Erratum Fig. 4 (Fig. 11 of the paper): the vertical forces along a vertical line at one time, against height,
    and their sum. The paper's is x = 0 at 2185 s, through its flux rope at 13.7 Mm, over 10-17 Mm (the defaults).
    --rope puts the line through the flux rope of this run at that time, wherever it is (the deepest O-point of the
    field in the corona, find_opoints), with 7 Mm of height around it; --x, --time, --ylim choose by hand."""
    items = runs_of(a.run)
    (t, f), = nearest(items, [a.time])
    g = load(f, ["log10rho_cgs", "p", "Bz"])
    x, ylim = a.x, a.ylim
    if a.rope:
        ropes = find_opoints(a_grid(g), g["X"], g["Y"])
        if not ropes:
            raise SystemExit(f"   no O-point of the field between 3 and 19.5 Mm at t = {t:.0f} s: no flux rope then; "
                             "give --x and --ylim, or another --time")
        print("   O-points of the field above 3 Mm (x, y [Mm], enclosed flux [B0 Mm]), deepest first:")
        for xo, yo, d in ropes[:5]:
            print(f"      ({xo:6.2f}, {yo:5.2f})   {d:.3g}")
        x = ropes[0][0]
        ylim = (ropes[0][1] - 3.5, ropes[0][1] + 3.5)
    xi, Y, grav, tens, magp, gasp = vertical_forces(g, x)
    m = (Y >= ylim[0]) & (Y <= ylim[1])
    tot = grav + tens + magp + gasp
    np.savetxt(fig_path(a, "forces", "csv"), np.c_[Y[m], grav[m], tens[m], magp[m], gasp[m], tot[m]], delimiter=",",
               header=f"x = {xi:.3f} Mm, t = {t:.1f} s; height [Mm], -rho g, (B.grad)B_y/4pi, -d(B^2/8pi)/dy, -dP/dy, total "
                      "[1e-5 dyn cm^-3]")
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.4, 6.2), sharex=True, gridspec_kw=dict(height_ratios=[3.2, 1], hspace=0.30))
    ink = "#0b0b0b"
    for v, c, ls, lab in ((grav, "#2a78d6", (0, (5, 2)), r"$-\rho g$"),
                          (tens, "#eb6834", (0, (5, 1.5, 1, 1.5)), r"$+\mathbf{B}\cdot\nabla B_y/4\pi$"),
                          (magp, "#1baf7a", (0, (5, 1.5, 1, 1.5, 1, 1.5)), r"$-\partial_y(B^2/8\pi)$"),
                          (gasp, ink, "solid", r"$-\partial_y P$")):
        a1.plot(Y[m], v[m], color=c, ls=ls, lw=1.8, label=lab, solid_capstyle="round")
    a2.plot(Y[m], tot[m], color=ink, lw=1.8)
    top = 1.15 * max(np.nanmax(np.abs(np.c_[grav[m], tens[m], magp[m], gasp[m]])), 1e-12)
    for ax, lim in ((a1, top), (a2, 1.15 * max(np.nanmax(np.abs(tot[m])), 1e-12))):
        ax.axhline(0.0, color="#c9c8c2", lw=0.8, zorder=0)
        ax.set_ylim(-lim, lim); ax.set_xlim(ylim); ax.tick_params(labelsize=8)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    a1.set_ylabel(r"Force (10$^{-5}$ dyn cm$^{-3}$)", fontsize=9)
    a1.legend(fontsize=9, frameon=False, loc="upper right", handlelength=3.2)
    a1.set_title("(a) the four forces", loc="left", fontsize=9, color=INK)
    a2.set_title("(b) their sum", loc="left", fontsize=9, color=INK)
    a2.set_ylabel(r"(10$^{-5}$ dyn cm$^{-3}$)", fontsize=9); a2.set_xlabel("Height (Mm)", fontsize=9)
    fig.suptitle(f"Vertical forces on the line x = {xi:.1f} Mm at t = {t:.0f} s", y=0.975, fontsize=11, weight="bold", color=INK)
    fig.subplots_adjust(left=0.13, right=0.97, top=0.90, bottom=0.09)
    out = fig_path(a, "forces")
    fig.savefig(out, dpi=a.dpi); plt.close(fig)
    print("wrote", out)


def cmd_balance(a):
    """What holds the dome: the vertical force balance on a vertical line, as a table. At each requested time and
    every --dy Mm between --ymin and --ymax on the mesh column nearest x = --x: log10 rho, B_x [G], v_y [km/s], and
    the four vertical force densities of the erratum's Fig. 4 [1e-5 dyn cm^-3] (gravity -rho g, the magnetic tension
    (B.grad)B_y/4pi, the magnetic pressure force -d(B^2/8pi)/dy, the gas pressure force -dP/dy) with their sum, the
    derivatives being those of the elements. Below it, each force integrated over three layers [dyn cm^-2]: from
    --ybase (3 Mm: the lifted plasma) to the top of the dense plasma (log10 rho > --thr), from there to where B_x
    changes sign (the top of the emerged field on this line), and from there to ymax (the field above). A dome at rest has a sum near zero in
    each; the terms say what carries the weight of the lifted plasma and what presses on it from above."""
    items = runs_of(a.run)
    for t, f in nearest(items, a.times):
        g = load(f, ["log10rho_cgs", "p", "v", "Bz"])
        G = lambda n: grid2d(g, g["fields"][n])
        Xc, Yc, X, Y = g["xc"], g["yc"], g["X"], g["Y"]
        i = int(np.argmin(np.abs(X - a.x)))
        Bx, By, Bz, P, L = G("Bx"), G("By"), G("Bz"), G("p"), G("log10rho_cgs")
        dx = lambda q: _sem_d(q, Xc, 1)
        dy = lambda q: _sem_d(q, Yc, 0)
        grav = -(10.0 ** L[:, i]) * G_SUN / 1e-5
        tens = (Bx * dx(By) + By * dy(By))[:, i] * FORCE_UNIT
        magp = -dy(0.5 * (Bx ** 2 + By ** 2 + Bz ** 2))[:, i] * FORCE_UNIT
        gasp = -dy(P)[:, i] * FORCE_UNIT
        tot = grav + tens + magp + gasp
        vy = G("v")[:, i] * V0_KMS
        bx = Bx[:, i] * B0_GAUSS
        col = (Y >= a.ymin) & (Y <= a.ymax)
        dense = col & (L[:, i] > a.thr) & (Y > 1.0)
        ydome = float(Y[dense].max()) if np.any(dense) else np.nan
        k = np.where(col & (Y > 1.0) & (bx > 0))[0]            # the emerged field: B_x > 0 above the photosphere
        yflux = float(Y[k].max()) if k.size else np.nan
        print(f"--- t = {t:.0f} s, x = {X[i]:.2f} Mm ({os.path.basename(f)}): top of log10 rho > {a.thr:g} at {ydome:.2f} Mm, "
              f"B_x changes sign at {yflux:.2f} Mm")
        print(f"{'y Mm':>6} {'log rho':>8} {'Bx G':>8} {'vy km/s':>8} | {'gravity':>9} {'tension':>9} {'mag pres':>9} {'gas pres':>9} {'sum':>9}   [1e-5 dyn cm^-3]")
        for y in np.arange(a.ymin, a.ymax + 1e-9, a.dy):
            j = int(np.argmin(np.abs(Y - y)))
            print(f"{Y[j]:6.2f} {L[j, i]:8.2f} {bx[j]:8.2f} {vy[j]:8.2f} | {grav[j]:9.3g} {tens[j]:9.3g} {magp[j]:9.3g} {gasp[j]:9.3g} {tot[j]:9.3g}")
        w = _trapw(Y) * 1e3                                    # (1e-5 dyn cm^-3) x Mm -> dyn cm^-2
        layers = [("ybase to the top of the dense plasma", a.ybase, ydome), ("from there to the sign change of B_x", ydome, yflux),
                  ("from there to ymax", yflux, a.ymax)]
        print(f"   integrals [dyn cm^-2]:{'':24} {'gravity':>9} {'tension':>9} {'mag pres':>9} {'gas pres':>9} {'sum':>9}")
        for name, y0, y1 in layers:
            if not (np.isfinite(y0) and np.isfinite(y1)) or y1 <= y0:
                print(f"   {name:<45} (none)")
                continue
            m = (Y >= y0) & (Y <= y1)
            I = [float(np.sum(w[m] * q[m])) for q in (grav, tens, magp, gasp, tot)]
            print(f"   {name + f' ({y0:.2f}-{y1:.2f} Mm)':<45} " + " ".join(f"{v:9.3g}" for v in I))


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
    out = fig_path(a, "profiles")
    fig.savefig(out, dpi=a.dpi); plt.close(fig)
    print("wrote", out)


def main():
    global FOLLOW_CHAIN, NEAREST, LINE_WIDTH, PREFIX, PAPER_LINES
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default="figs_leake2010")
    p.add_argument("--format", default="png")
    p.add_argument("--dpi", type=int, default=150)
    p.add_argument("--no-chain", action="store_true", help="only the directories given, not the restarts of the run")
    p.add_argument("--nearest", action="store_true", help="the nearest snapshot for a time the run has none at")
    p.add_argument("--lw", type=float, default=LINE_WIDTH, help="width of the field lines, in points")
    p.add_argument("--prefix", default=PREFIX, help="start of every figure's file name (default flux_emergence_Leake2010_)")
    p.add_argument("--paper-lines", action="store_true",
                   help="the field lines of the maps are the erratum's own (A_z = 820 - 61.35 k G Mm; quadrupole runs)")
    # the same options after the command name (they keep the value given before it, if any)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--outdir", default=argparse.SUPPRESS)
    common.add_argument("--format", default=argparse.SUPPRESS)
    common.add_argument("--dpi", type=int, default=argparse.SUPPRESS)
    common.add_argument("--no-chain", action="store_true", default=argparse.SUPPRESS)
    common.add_argument("--nearest", action="store_true", default=argparse.SUPPRESS)
    common.add_argument("--lw", type=float, default=argparse.SUPPRESS)
    common.add_argument("--prefix", default=argparse.SUPPRESS)
    common.add_argument("--paper-lines", action="store_true", default=argparse.SUPPRESS)
    sub = p.add_subparsers(dest="cmd", required=True)
    for name, fn in (("density", cmd_density), ("shear", cmd_shear)):
        s = sub.add_parser(name, parents=[common]); s.set_defaults(fn=fn)
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
    pair = lambda t: tuple(map(float, t.split(",")))
    s = sub.add_parser("dsgs", parents=[common]); s.set_defaults(fn=cmd_dsgs)
    s.add_argument("run")
    s.add_argument("--which", choices=("momentum", "induction", "energy", "mu", "eta", "rho"), default="momentum",
                   help="whose coefficient: the momentum equations' (mu), the induction equations' (eta), "
                        "or the mass, energy and psi equations' (rho)")
    s.add_argument("--times", default="", help="comma-separated times in s")
    s.add_argument("--xlim", type=pair, default=(-12.75, 12.75))
    s.add_argument("--ylim", type=pair, default=(-3.0, 19.5))
    s.add_argument("--clim", type=pair, default=None, help="limits of log10 of the coefficient (default -2,2; eta: -4,1)")
    s.add_argument("--linear", action="store_true", help="the coefficient itself, not its log10 (default limits 0,5; eta: 0,0.5)")
    s.add_argument("--nlines", type=int, default=45, help="field lines (0: none)")
    s = sub.add_parser("movie", parents=[common]); s.set_defaults(fn=cmd_movie)
    s.add_argument("run")
    s.add_argument("--field", choices=tuple(MAPS) + tuple(MAP_ALIAS), default="density",
                   help="density, or a DynSGS coefficient: momentum (mu), induction (eta), energy (rho)")
    s.add_argument("--fps", type=float, default=4.0, help="snapshots per second of the GIF (at the usual snapshot interval)")
    s.add_argument("--tmin", type=float, default=-np.inf, help="first time [s]")
    s.add_argument("--tmax", type=float, default=np.inf, help="last time [s]")
    s.add_argument("--every", type=int, default=1, help="every N-th snapshot")
    s.add_argument("--xlim", type=pair, default=(-12.75, 12.75))
    s.add_argument("--ylim", type=pair, default=(-3.0, 19.5))
    s.add_argument("--clim", type=pair, default=None, help="colour limits (default: the field's)")
    s.add_argument("--linear", action="store_true", help="a DynSGS coefficient itself, not its log10")
    s.add_argument("--nlines", type=int, default=45, help="field lines (0: none)")
    s.add_argument("--arrows", type=lambda t: [float(v) for v in t.split(",")] if t else [],
                   default=[3.0, 11.0, 16.0], help="heights [Mm] of the B_x arrows of the density frames ('' for none)")
    s.add_argument("--label", default="", help="text in the top left corner of each frame")
    s.add_argument("--no-gif", action="store_true", help="the images only")
    s = sub.add_parser("axis", parents=[common]); s.set_defaults(fn=cmd_axis)
    s.add_argument("runs", nargs="+", help="RUN or LABEL=RUN")
    s.add_argument("--tmax", type=float, default=4600.0)
    s.add_argument("--tall", action="store_true", help="y axis up to 12 Mm, as the paper's panel")
    s.add_argument("--shift", action="store_true", help="shift each run to start at the erratum's -1.84 Mm")
    s = sub.add_parser("dome", parents=[common]); s.set_defaults(fn=cmd_dome)
    s.add_argument("runs", nargs="+", help="RUN or LABEL=RUN")
    s.add_argument("--thr", type=float, default=-13.0, help="log10 rho threshold [g cm^-3]")
    s = sub.add_parser("forces", parents=[common]); s.set_defaults(fn=cmd_forces)
    s.add_argument("run")
    s.add_argument("--time", type=float, default=2185.0, help="time [s] (default 2185, the paper's)")
    s.add_argument("--x", type=float, default=0.0, help="the vertical line, x [Mm] (default 0)")
    s.add_argument("--ylim", type=lambda t: tuple(map(float, t.split(","))), default=(10.0, 17.0), help="heights [Mm] (default 10,17)")
    s.add_argument("--rope", action="store_true", help="the line through this run's flux rope at that time, 7 Mm around it")
    s = sub.add_parser("balance", parents=[common]); s.set_defaults(fn=cmd_balance)
    s.add_argument("run")
    s.add_argument("--times", type=lambda t: [float(x) for x in t.split(",")], default=[1380.0, 1610.0, 1955.0, 2185.0],
                   help="times [s], comma-separated (default 1380,1610,1955,2185)")
    s.add_argument("--x", type=float, default=0.0, help="the vertical line, x [Mm] (default 0)")
    s.add_argument("--ymin", type=float, default=0.0); s.add_argument("--ymax", type=float, default=14.0)
    s.add_argument("--dy", type=float, default=0.5, help="height step of the table [Mm]")
    s.add_argument("--thr", type=float, default=-13.0, help="log10 rho of the dense plasma [g cm^-3]")
    s.add_argument("--ybase", type=float, default=3.0, help="bottom of the first layer of the integrals [Mm]")
    s = sub.add_parser("mass", parents=[common]); s.set_defaults(fn=cmd_mass)
    s.add_argument("run")
    s.add_argument("--ref", default="", help="run whose first snapshot (t = 0) gives the initial pressure and the unit "
                                             "of mass, for a restarted run listed from its restart time on")
    s = sub.add_parser("column", parents=[common]); s.set_defaults(fn=cmd_column)
    s.add_argument("run")
    s.add_argument("--times", type=lambda t: [float(x) for x in t.split(",")], default=[1380.0, 1610.0, 1955.0, 2185.0],
                   help="times [s], comma-separated (default 1380,1610,1955,2185)")
    s.add_argument("--x", type=float, default=0.0, help="the vertical line, x [Mm] (default 0)")
    s.add_argument("--ymin", type=float, default=0.0); s.add_argument("--ymax", type=float, default=16.0)
    s.add_argument("--dy", type=float, default=0.5, help="height step of the table [Mm]")
    s = sub.add_parser("fastest", parents=[common]); s.set_defaults(fn=cmd_fastest)
    s.add_argument("run")
    s.add_argument("--times", type=lambda t: [float(x) for x in t.split(",")], default=None,
                   help="times [s], comma-separated (default: the last three snapshots)")
    s.add_argument("--n", type=int, default=5, help="how many places (default 5)")
    s.add_argument("--sep", type=float, default=0.5, help="their least distance apart [Mm] (default 0.5)")
    s = sub.add_parser("massbudget", parents=[common]); s.set_defaults(fn=cmd_massbudget)
    s.add_argument("run")
    s.add_argument("--time", type=float, default=2185.0, help="the snapshot nearest this time [s] (default 2185); the next one gives the change")
    s.add_argument("--x", type=lambda t: [float(v) for v in t.split(",")], default=[-2.0, 2.0],
                   help="the vertical lines [Mm], comma-separated (default -2,2)")
    s.add_argument("--ylim", type=lambda t: tuple(map(float, t.split(","))), default=(3.0, 5.7))
    s = sub.add_parser("profile", parents=[common]); s.set_defaults(fn=cmd_profile)
    s.add_argument("run")
    s.add_argument("--times", type=lambda t: [float(x) for x in t.split(",")], default=[1955.0, 2185.0, 2530.0],
                   help="times [s], comma-separated (default 1955,2185,2530; the erratum has 920, 1380, 1610 too)")
    s = sub.add_parser("noise", parents=[common]); s.set_defaults(fn=cmd_noise)
    s.add_argument("run")
    s.add_argument("--times", type=lambda t: [float(x) for x in t.split(",")], default=[1955.0, 2185.0, 2530.0],
                   help="times [s], comma-separated (default 1955,2185,2530)")
    s.add_argument("--xlim", type=lambda t: tuple(map(float, t.split(","))), default=(-6.0, 6.0))
    s.add_argument("--ylim", type=lambda t: tuple(map(float, t.split(","))), default=(3.0, 14.0))
    s.add_argument("--vmax", type=float, default=0.3, help="top of the colour scale of the node-to-node maps [dex]")
    s = sub.add_parser("elements", parents=[common]); s.set_defaults(fn=cmd_elements)
    s.add_argument("run")
    s.add_argument("--times", type=lambda t: [float(x) for x in t.split(",")], default=[2185.0, 2530.0],
                   help="times [s], comma-separated (default 2185,2530)")
    s.add_argument("--x", type=lambda t: [float(v) for v in t.split(",")], default=[-2.0, 0.0, 2.0],
                   help="the vertical lines [Mm], comma-separated (default -2,0,2)")
    s.add_argument("--ylim", type=lambda t: tuple(map(float, t.split(","))), default=(3.0, 9.5))
    s = sub.add_parser("drain", parents=[common]); s.set_defaults(fn=cmd_drain)
    s.add_argument("run")
    s.add_argument("--ref", default="", help="run whose first snapshot (t = 0) gives the unit of mass, for a restarted run")
    s.add_argument("--thr", type=float, default=-13.0, help="log10 rho above which the plasma counts as lifted [g cm^-3]")
    s.add_argument("--ymin", type=float, default=3.0, help="bottom of the box [Mm]")
    s.add_argument("--tfloor", type=float, default=0.1, help="the run's temperature floor in T_ph (JEXPRESSO_LK10_TFLOOR)")
    s = sub.add_parser("heat", parents=[common]); s.set_defaults(fn=cmd_heat)
    s.add_argument("run")
    s.add_argument("--thr", type=float, default=-13.0, help="log10 rho above which the plasma counts as lifted [g cm^-3]")
    s.add_argument("--ymin", type=float, default=3.0, help="bottom of the box [Mm] (above the photosphere)")
    s.add_argument("--ymax", type=float, default=20.0, help="top of the box [Mm]")
    s.add_argument("--ch", type=float, default=CH_DEFAULT, help="the run's GLM speed c_h, code units (its log: 'GLM: c_h = ...')")
    s.add_argument("--xsplit", type=float, default=2.5, help="|x| [Mm] between the center and the flanks (default 2.5)")
    s.add_argument("--rhofac", type=float, default=None,
                   help="the factor on the DynSGS mass diffusion (default: rho_factor of each directory's lk10_case.txt)")
    s.add_argument("--pair", type=float, default=20.0, help="snapshots at most this far apart [s] are compared with the budget")
    s.add_argument("--lare-dx", type=float, default=LARE_DX_MM,
                   help="the cell of the paper's Lare grid for the 'Lare visc' column [Mm] (default 0.048)")
    s = sub.add_parser("null", parents=[common]); s.set_defaults(fn=cmd_null)
    s.add_argument("runs", nargs="+", help="RUN or LABEL=RUN")
    s.add_argument("--thr", type=float, default=-13.0, help="log10 rho threshold of the dome top [g cm^-3]")
    s = sub.add_parser("core", parents=[common]); s.set_defaults(fn=cmd_core)
    s.add_argument("runs", nargs="+", help="RUN or LABEL=RUN")
    s.add_argument("--tmax", type=float, default=2600.0, help="limit of the time axis [s]")
    s = sub.add_parser("lare", parents=[common]); s.set_defaults(fn=cmd_lare)
    s.add_argument("run")
    s.add_argument("--times", type=lambda t: [float(x) for x in t.split(",")], default=[464.0, 696.0, 920.0, 1044.0, 1380.0],
                   help="times [s], comma-separated (default 464,696,920,1044,1380: the rise, the arrival at the surface, the emergence)")
    s.add_argument("--lare-dx", type=float, default=0.0375, help="the paper's Lare cell [Mm] (default 0.0375, its cell at the surface)")
    s = sub.add_parser("lines", parents=[common]); s.set_defaults(fn=cmd_lines)
    s.add_argument("runs", nargs="+", help="RUN or LABEL=RUN")
    s.add_argument("--times", default="", help="comma-separated times in s of the tables (default: those of the erratum's Fig. 3)")
    s = sub.add_parser("eshear", parents=[common]); s.set_defaults(fn=cmd_eshear)
    s.add_argument("runs", nargs="+", help="RUN or LABEL=RUN")
    s = sub.add_parser("profiles", parents=[common]); s.set_defaults(fn=cmd_profiles)
    s.add_argument("run")
    s = sub.add_parser("atmos", parents=[common]); s.set_defaults(fn=cmd_atmos)
    s.add_argument("run")
    s.add_argument("--times", default="", help="comma-separated times in s (default 0,464,920,1160,1380,1610,1955,2530)")
    s.add_argument("--xmin", type=float, default=8.0, help="|x| range of the median, Mm")
    s.add_argument("--xmax", type=float, default=12.75)
    s.add_argument("--ylim", type=lambda t: tuple(map(float, t.split(","))), default=(-1.0, 20.0))
    s = sub.add_parser("diag", parents=[common]); s.set_defaults(fn=cmd_diag)
    s.add_argument("run")
    s.add_argument("--dt", type=float, default=None, help="the run's time step in t0 (its log): print the CFL number")
    s.add_argument("--ymin", type=float, default=None, help="only y >= YMIN Mm (e.g. 2 for the corona)")
    s.add_argument("--tmin", type=float, default=None, help="only snapshots from TMIN s")
    s.add_argument("--tfloor", type=float, default=0.1, help="the run's temperature floor (T_ph): count the nodes at it")
    a = p.parse_args()
    FOLLOW_CHAIN, NEAREST, LINE_WIDTH, PREFIX, PAPER_LINES = not a.no_chain, a.nearest, a.lw, a.prefix, a.paper_lines
    os.makedirs(a.outdir, exist_ok=True)
    a.fn(a)


if __name__ == "__main__":
    main()
