# Figure of the one-page DynSGS-cutoff note: state L2 error of the smooth Alfven
# wave at t = 2 against #DOFs, N = 5, 6, 7: Galerkin, DynSGS without the cutoff,
# DynSGS with the cutoff 1e-3 (the deck's default). Reads the case's error store.
import glob, math, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

E, OUT = sys.argv[1], sys.argv[2]
R = {}
for f in glob.glob(os.path.join(E, "*.dat")):
    meta, l2 = {}, {}
    for ln in open(f, encoding="utf-8"):
        t = ln.split()
        if not t: continue
        if t[0] == "#": meta.update(dict(s.split("=", 1) for s in t if "=" in s))
        else: l2[t[0]] = float(t[1])
    R[(int(meta["nop"]), int(meta["nelx"]), meta["visc"])] = math.sqrt(sum(v*v for v in l2.values()))

# the case figure's colours (a categorical palette checked for colour-vision
# deficiency on neighbouring curves): N = 5 aqua, 6 yellow, 7 magenta
COL = {5: "#1baf7a", 6: "#eda100", 7: "#e87ba4"}
MK  = {5: "^", 6: "*", 7: "h"}
INK, MUTED = "#0b0b0b", "#52514e"

plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6})
fig, ax = plt.subplots(figsize=(3.35, 2.75))
for N in (5, 6, 7):
    ms = sorted({k[1] for k in R if k[0] == N})
    x  = [(m*N)**2 for m in ms]
    g  = [R.get((N, m, "galerkin"))   for m in ms]
    nc = [R.get((N, m, "dsgs_nocut")) for m in ms]
    ct = [R.get((N, m, "dsgs"))       for m in ms]
    c = COL[N]
    ax.plot(x, nc, "-",  color=c, lw=1.5, marker=MK[N], ms=4.5, mfc=c, mec=c, zorder=3)
    ax.plot(x, g,  "--", color=c, lw=1.2, marker=MK[N], ms=8,   mfc="white", mec=c, mew=1.0, zorder=2)
    ax.plot(x, ct, "none", marker=MK[N], ms=3.6, mfc=c, mec=c, zorder=4)
    # direct label at the fine end of the Galerkin curve (text in ink; the
    # marker beside it carries the colour)
    ax.annotate(f"N = {N}", (x[-1], g[-1]), xytext=(7, 0), textcoords="offset points",
                ha="left", va="center", fontsize=7.5, color=INK)
# slope guides beside the curves they describe: 4 above DynSGS without the
# cutoff (N = 6), 8 below Galerkin N = 7
for p, x0, y0, lab in ((4, 9216.0, 2.2e-7, "4"), (8, 3136.0, 2.0e-10, "8")):
    xs = [x0, x0*4.0]; ys = [y0, y0*4.0**(-p/2)]
    ax.plot(xs, ys, ":", color=MUTED, lw=1.0)
    ax.annotate(lab, (xs[1], ys[1]), xytext=(3, 0), textcoords="offset points",
                fontsize=7.5, color=MUTED, va="center")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(3.0e2, 1.3e5)
ax.set_xlabel(r"#DOFs $=(n_{el}N)^2$")
ax.set_ylabel(r"$L^2$ error of the conserved state, $t=2$")
ax.grid(True, which="major", color="#d9d8d4", lw=0.5)
ax.set_axisbelow(True)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
# legend of the three treatments, drawn in neutral ink
from matplotlib.lines import Line2D
h = [Line2D([], [], color=MUTED, lw=1.5, marker="o", ms=4, mfc=MUTED, mec=MUTED, label="DynSGS, no cutoff"),
     Line2D([], [], color=MUTED, lw=0, marker="o", ms=3.4, mfc=MUTED, mec=MUTED, label=r"DynSGS, cutoff $10^{-3}$"),
     Line2D([], [], color=MUTED, lw=1.2, ls="--", marker="o", ms=7, mfc="white", mec=MUTED, label="Galerkin")]
ax.legend(handles=h, loc="upper right", fontsize=7, frameon=False, handlelength=2.4)
fig.tight_layout(pad=0.3)
fig.savefig(OUT, dpi=220)
print("wrote", OUT)
