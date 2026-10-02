#!/usr/bin/env python3
"""Merge and plot the 3D periodic Poisson solver comparison.

    python3 tools/poisson3d_benchmark/plot3d.py OUTDIR [--direct]

Reads every OUTDIR/parts/*/results.csv (one per solver) and, if present,
OUTDIR/results.csv; a configuration (solver, d, ne, N) seen twice keeps its
last occurrence. Writes OUTDIR/results.csv, OUTDIR/results.md and, per SEM
order N, into OUTDIR/assets (light and dark SVG):

  ppb3d_cost_N<N>      setup + solve vs unknowns n, one curve per solver, with
                       the reference slopes n (optimal iterative) and n^2
                       (3D sparse direct with nested dissection)
  ppb3d_total_N<N>     time-to-solution (assembly + rhs + setup + solve) vs n
  ppb3d_memory_N<N>    peak resident memory of the run vs n (MPI: summed over ranks)
  ppb3d_iters_N<N>     CG iterations vs n (AMG and Jacobi preconditioners)
  ppb3d_error_N<N>     L-inf error vs n (the SEM solvers share one curve)

and the figures of the 2D benchmark, per mesh level (ne^d elements, the
order N varying; drawn when a level has two or more orders):

  ppb3d_error_vs_solve_time_ne<ne>   L-inf error vs time of the solve step
  ppb3d_error_vs_cost_ne<ne>         L-inf error vs setup + solve
  ppb3d_error_vs_total_time_ne<ne>   L-inf error vs time-to-solution
  ppb3d_error_vs_order_ne<ne>        L-inf error vs SEM order N
  ppb3d_error_vs_Ng                  L-inf error vs unknowns per direction ne*N,
                                     all levels (one SEM curve per mesh; Fourier
                                     solvers; the predicted r^(ne*N/2))

The direct solvers (sem, sc_direct, mumps) stay in the tables but are left
out of the figures unless --direct is given.

No dependencies (uses the SVG chart of ../periodic_poisson_benchmark/plot.py).
"""
import csv, glob, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "periodic_poisson_benchmark"))
import plot as P  # noqa: E402

SOLVERS = [  # key, label, (colour index, marker)
    ("sem",        "SEM direct",  (0, "circle")),
    ("sem_amg",    "SEM AMG",     (1, "square")),
    ("pmg_amg",    "p-MG + AMG",  (0, "circle")),
    ("pmg_gmg",    "p-MG + GMG",  (2, "triangle")),
    ("sem_jacobi", "SEM Jacobi",  (6, "square")),
    ("sc_direct",  "SC direct",   (2, "triangle")),
    ("sc_amg",     "SC AMG",      (3, "diamond")),
    ("ps",         "pseudo-spectral", (4, "tridown")),
    ("fft",        "FFT",         (5, "ring")),
    # MPI solvers (mpi/bench3d_mpi.jl), same discretisation as the SEM solvers
    ("mumps",      "MUMPS (MPI)",        (0, "circle")),
    ("boomeramg",  "BoomerAMG-CG (MPI)", (1, "square")),
    ("jacobi",     "Jacobi-CG (MPI)",    (6, "square")),
]
SEM = {"sem", "sem_amg", "pmg_amg", "pmg_gmg", "sem_jacobi", "sc_direct", "sc_amg", "mumps", "boomeramg", "jacobi"}
# direct solvers: kept in results.csv / results.md, left out of the figures
# (never used for large problems) unless --direct is given; they share the
# colours of the p-multigrid solvers
DIRECT = {"sem", "sc_direct", "mumps"}
LABEL = {k: l for k, l, _ in SOLVERS}


def merge(outdir):
    files = sorted(glob.glob(os.path.join(outdir, "parts", "*", "results.csv")))
    top = os.path.join(outdir, "results.csv")
    if os.path.isfile(top):
        files = [top] + files
    byconf, hdr = {}, None
    for f in files:
        with open(f) as fh:
            rd = csv.DictReader(fh)
            hdr = hdr or rd.fieldnames
            for r in rd:
                byconf[(r["solver"], r["d"], r["ne"], r["nop"])] = r
            hdr = hdr + [k for k in rd.fieldnames if k not in hdr]
    rows = sorted(byconf.values(), key=lambda r: (int(r["d"]), int(r["nop"]), int(r["n"] or 0),
                                                  [k for k, _, _ in SOLVERS].index(r["solver"])))
    if not rows:
        sys.exit(f"plot3d.py: no results in {outdir}")
    with open(top, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=hdr, restval="")
        w.writeheader()
        w.writerows(rows)
    print(f"merged {len(rows)} configurations from {len(files)} file(s) -> {top}")
    return rows


def fnum(x, default=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def t_str(x):
    x = fnum(x)
    if x <= 0:
        return "—"
    return f"{x*1e3:.3g} ms" if x < 1 else (f"{x:.3g} s" if x < 600 else f"{x/60:.3g} min")


def write_md(path, rows):
    with open(path, "w") as fh:
        fh.write("| solver | d | elements | N | unknowns n | CG its | ‖e‖∞ | assembly | setup | solve | time-to-solution | factor nnz | peak memory | threads (julia/BLAS) | MPI ranks | status |\n")
        fh.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|\n")
        for r in rows:
            ok = r["status"] == "ok"
            fh.write("| {} | {} | {}^{} | {} | {:,} | {} | {} | {} | {} | {} | **{}** | {} | {} GB | {}/{} | {} | {} |\n".format(
                LABEL.get(r["solver"], r["solver"]), r["d"], r["ne"], r["d"], r["nop"], int(fnum(r["n"])),
                r["iters"] if ok and fnum(r["iters"]) > 0 else "—",
                f"{fnum(r['linf']):.2e}" if ok else "—",
                t_str(r["assembly"]) if ok else "—", t_str(r["setup"]) if ok else "—",
                t_str(r["solve"]) if ok else "—", t_str(r["total"]) if ok else "—",
                f"{int(fnum(r['factor_nnz'])):,}" if fnum(r["factor_nnz"]) > 0 else "—",
                r.get("maxrss_gb", ""), r.get("julia_threads", ""), r.get("blas_threads", ""), r.get("nranks") or "1", r["status"]).replace(",", " "))
    print("wrote", path)


def decades(vals):
    vals = [v for v in vals if v > 0]
    return (math.floor(math.log10(min(vals))), math.ceil(math.log10(max(vals))))


def ref_line(x0, y0, slope, xs, ydec):
    """y = y0 (x/x0)^slope over xs, clipped to the y range"""
    pts = [(x, y0 * (x / x0) ** slope) for x in xs]
    return [(x, y) for x, y in pts if 10.0 ** ydec[0] <= y <= 10.0 ** ydec[1]]


def figure(assets, name, title, sub, desc, series, ylabel, ydec, refs=()):
    xs = sorted({p[0] for _, pts, _ in series for p in pts})
    if len(xs) < 2:
        return
    xgrid = [xs[0] * (xs[-1] / xs[0]) ** (i / 60) for i in range(61)]
    allser = list(series)
    for (x0, y0, slope, lab) in refs:
        rl = ref_line(x0, y0, slope, xgrid, ydec)
        if len(rl) > 1:
            allser.append((lab, rl, (None, None, True)))
    yticks_time = "seconds" in ylabel
    P.write(assets, name, lambda th: P.chart(
        th, title, sub, desc, allser, (xs[0], xs[-1]), True, P.log_ticks(xs),
        "unknowns n (log scale)", ylabel, ydec, end_labels=False, legend_cols=5))


def plots(rows, assets):
    os.makedirs(assets, exist_ok=True)
    ok = [r for r in rows if r["status"] == "ok"]
    for d in sorted({r["d"] for r in ok}):
        for N in sorted({int(r["nop"]) for r in ok if r["d"] == d}):
            rs = [r for r in ok if r["d"] == d and int(r["nop"]) == N]
            sup = {"2": "²", "3": "³"}[d]
            sub = f"−∇²u = f on [0,2π]{sup}, fully periodic · SEM order N = {N} (Fourier grids (ne·N){sup}) · second-run timings"

            def ser(key, solvers=None, transform=lambda r, v: v):
                out = []
                for k, lab, (c, m) in SOLVERS:
                    if solvers is not None and k not in solvers:
                        continue
                    pts = sorted((int(fnum(r["n"])), transform(r, fnum(r[key]))) for r in rs if r["solver"] == k)
                    pts = [p for p in pts if p[1] > 0]
                    if pts:
                        out.append((lab, pts, (c, m, False)))
                return out

            sfx = f"_d{d}_N{N}" if d != "3" else f"_N{N}"
            cost = ser("setup", transform=lambda r, v: v + fnum(r["solve"]))
            if cost:
                yd = decades([p[1] for _, pts, _ in cost for p in pts])
                refs = []
                for keys, slope, lab in ((("sem_amg", "boomeramg"), 1.0, "∝ n"), (("sem", "mumps"), 2.0, "∝ n²")):
                    first = [pts[0] for l, pts, _ in cost if l in [LABEL[k] for k in keys]]
                    if first:
                        refs.append((first[0][0], first[0][1], slope, lab))
                figure(assets, "ppb3d_cost" + sfx, f"Solver cost vs unknowns, N = {N}", sub,
                       f"Setup plus solve wall-clock versus number of unknowns for every solver at SEM order {N}, "
                       "with dashed reference slopes n and n squared.", cost, "setup + solve, seconds (log scale)", yd, refs)
            tot = ser("total")
            if tot:
                figure(assets, "ppb3d_total" + sfx, f"Time-to-solution vs unknowns, N = {N}", sub,
                       f"Assembly plus right-hand side plus setup plus solve versus unknowns at SEM order {N}.",
                       tot, "time-to-solution, seconds (log scale)", decades([p[1] for _, pts, _ in tot for p in pts]))
            mem = ser("maxrss_gb")
            if mem:
                figure(assets, "ppb3d_memory" + sfx, f"Peak memory vs unknowns, N = {N}", sub,
                       f"Peak resident memory of the run versus unknowns at SEM order {N}.",
                       mem, "peak memory, GB (log scale)", decades([p[1] for _, pts, _ in mem for p in pts]))
            its = ser("iters", solvers={"sem_amg", "pmg_amg", "pmg_gmg", "sem_jacobi", "sc_amg", "boomeramg", "jacobi"})
            if its:
                figure(assets, "ppb3d_iters" + sfx, f"CG iterations vs unknowns, N = {N}", sub,
                       f"Conjugate-gradient iterations to a relative preconditioned residual of 1e-12 versus unknowns at SEM order {N}.",
                       its, "CG iterations (log scale)", decades([p[1] for _, pts, _ in its for p in pts]))
            err = []
            for k, lab, (c, m) in SOLVERS:
                if k not in ("sem", "ps", "fft"):
                    continue
                pts = sorted((int(fnum(r["n"])), fnum(r["linf"])) for r in rs if r["solver"] == k or
                             (k == "sem" and r["solver"] in SEM))
                byn = {}
                for n, e in pts:
                    byn.setdefault(n, e)
                pts = [(n, e) for n, e in sorted(byn.items()) if e > 0]
                if pts:
                    err.append(("SEM (all solvers)" if k == "sem" else lab, pts, (c, m, False)))
            if err:
                figure(assets, "ppb3d_error" + sfx, f"Error vs unknowns, N = {N}", sub,
                       f"L-infinity error against the exact solution versus unknowns at SEM order {N}; the SEM solvers share one curve.",
                       err, "L∞ error (log scale)", decades([p[1] for _, pts, _ in err for p in pts]))


# ---------------------------------------------------------------------------
# The figures of the 2D benchmark (periodic_poisson_benchmark/plot.py), for
# a sweep over mesh levels (elements per side ne) and orders N:
#   per level:  error vs time of the solve step, error vs setup + solve,
#               error vs time-to-solution, error vs order
#   all levels: error vs unknowns per direction, ∛n = ne·N
# ---------------------------------------------------------------------------
SUPD = {"2": "²", "3": "³"}


def _sem_label(present):
    k = [s for s in present if s in SEM]
    return f"SEM (all {len(k)} solvers)" if len(k) > 1 else LABEL[k[0]]


def level_figures(rs, assets, d, ne):
    """The figures of one mesh level (ne^d elements), the order N varying."""
    by = {k: sorted([r for r in rs if r["solver"] == k], key=lambda r: int(r["nop"])) for k, _, _ in SOLVERS}
    present = [k for k, _, _ in SOLVERS if by[k]]
    if not present or len({r["nop"] for r in rs}) < 2:
        return
    sup = SUPD[d]
    sub = (f"−∇²u = f on [0,2π]{sup}, fully periodic · {ne}{sup} elements · second-run timings")
    ylab = "L∞ error vs exact solution"
    ydec = decades([fnum(r["linf"]) for r in rs])
    nops = sorted({int(r["nop"]) for r in rs})
    sfx = f"_ne{ne}" + ("" if d == "3" else f"_d{d}")

    def series(val):
        out = []
        for k, lab, (c, m) in SOLVERS:
            pts = [(val(r), fnum(r["linf"])) for r in by[k] if val(r) > 0 and fnum(r["linf"]) > 0]
            if pts:
                out.append((lab, pts, (c, m, False)))
        return out

    # the SEM solvers solve the same discrete system: one error curve for them
    err, sem_done = [], False
    for k, lab, (c, m) in SOLVERS:
        if not by[k]:
            continue
        if k in SEM:
            if sem_done:
                continue
            sem_done, lab = True, _sem_label(present)
        err.append((lab, [(int(r["nop"]), fnum(r["linf"])) for r in by[k]], (c, m, False)))
    P.write(assets, "ppb3d_error_vs_order" + sfx, lambda th: P.chart(
        th, "Error vs order", sub,
        f"L-infinity error versus SEM polynomial order N on {ne}^{d} elements (one SEM curve: its solvers give the "
        f"same discrete solution); the Fourier solvers on the ({ne}N)^{d} grid, the same unknowns.",
        err, (nops[0], nops[-1]), False, [(n, str(n)) for n in nops],
        f"SEM order N   (Fourier solvers: ({ne}N){sup} grid, same unknowns)", ylab, ydec))

    for key, val, name, title, xl in (
            ("solve", lambda r: fnum(r["solve"]), "ppb3d_error_vs_solve_time", "Error vs time of the solve step",
             "solve-step wall-clock (log scale)"),
            ("cost", lambda r: fnum(r["setup"]) + fnum(r["solve"]), "ppb3d_error_vs_cost",
             "Error vs setup + solve", "setup + solve wall-clock (log scale)"),
            ("total", lambda r: fnum(r["total"]), "ppb3d_error_vs_total_time", "Error vs time-to-solution",
             "time-to-solution incl. assembly (log scale)")):
        ser = series(val)
        ts = [p[0] for _, pts, _ in ser for p in pts]
        if not ts:
            continue
        e0, e1 = math.floor(math.log10(min(ts))), math.ceil(math.log10(max(ts)))
        P.write(assets, name + sfx, lambda th, ser=ser, title=title, xl=xl, e0=e0, e1=e1: P.chart(
            th, title, sub,
            f"L-infinity error versus {xl} on {ne}^{d} elements, one curve per solver over the orders N.",
            ser, (10.0 ** e0, 10.0 ** e1), True,
            [(10.0 ** e, P.time_label(e)) for e in range(e0, e1 + 1)], xl, ylab, ydec))


def ng_figure(rs, assets, d, r):
    """Error versus unknowns per direction, ne·N, for all mesh levels: one SEM
    curve per mesh, the Fourier solvers of all levels on one curve each (they
    depend on ne·N only), and the predicted Fourier error r^(ne·N/2)."""
    sup = SUPD[d]
    root = {"2": "√n", "3": "∛n"}[d]
    ng = lambda q: int(q["ne"]) * int(q["nop"])
    levels = sorted({int(q["ne"]) for q in rs})
    marks = ["circle", "square", "triangle", "diamond", "tridown", "ring"]
    sem_rows = [q for q in rs if q["solver"] in SEM]
    series = []
    for i, ne in enumerate(levels):
        byng = {}
        for q in sem_rows:                          # the SEM solvers agree: one point per (ne, N)
            if int(q["ne"]) == ne:
                byng.setdefault(ng(q), fnum(q["linf"]))
        pts = sorted((x, y) for x, y in byng.items() if y > 0)
        if len(pts) >= 2:
            series.append((f"SEM, {ne}{sup} el.", pts, (0, marks[i % len(marks)], False)))
    for key, lab, (c, m) in [x for x in SOLVERS if x[0] in ("ps", "fft")]:
        byng = {ng(q): fnum(q["linf"]) for q in rs if q["solver"] == key and fnum(q["linf"]) > 0}
        if byng:
            series.append((lab, sorted(byng.items()), (c, m, False)))
    allpts = [p for _, pts, _ in series for p in pts]
    if len(series) < 2 and not any(len(pts) > 1 for _, pts, _ in series):
        return
    x0, x1 = min(p[0] for p in allpts), max(p[0] for p in allpts)
    ydec = decades([p[1] for p in allpts])
    pred = [(n, r ** (n / 2)) for n in range(x0, x1 + 1, max(1, (x1 - x0) // 200))
            if r ** (n / 2) >= 10.0 ** ydec[0]]
    if len(pred) > 1:
        series.append((f"{r}^({root}/2)", pred, (None, None, True)))
    step = next(s for s in (16, 32, 64, 128, 256) if (x1 - x0) / s <= 8)
    xticks = [(n, str(n)) for n in range(0, x1 + 1, step) if n >= x0]
    name = "ppb3d_error_vs_Ng" + ("" if d == "3" else f"_d{d}")
    P.write(assets, name, lambda th: P.chart(
        th, "Error vs unknowns per direction, all mesh levels",
        f"−∇²u = f on [0,2π]{sup}, fully periodic · {root} = (elements per side) × N, the same for every solver",
        f"L-infinity error versus the number of unknowns per direction, (elements per side) times N, for every "
        f"mesh level: one SEM curve per mesh, the Fourier solvers of all levels, and the predicted Fourier error "
        f"{r}^({root}/2) (dashed).",
        series, (x0, x1), False, xticks,
        f"unknowns per direction, {root} = (elements per side) × N (linear scale)",
        "L∞ error vs exact solution", ydec, end_labels=False, legend_cols=4))


def level_plots(rows, assets):
    ok = [r for r in rows if r["status"] == "ok"]
    for d in sorted({r["d"] for r in ok}):
        rd = [r for r in ok if r["d"] == d]
        for ne in sorted({int(r["ne"]) for r in rd}):
            level_figures([r for r in rd if int(r["ne"]) == ne], assets, d, ne)
        rr = sorted({fnum(r["r"]) for r in rd})
        ng_figure(rd, assets, d, rr[0] if len(rr) == 1 else 0.5)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    outdir = sys.argv[1]
    rows = merge(outdir)
    write_md(os.path.join(outdir, "results.md"), rows)
    if "--direct" not in sys.argv[2:]:
        rows = [r for r in rows if r["solver"] not in DIRECT]
    plots(rows, os.path.join(outdir, "assets"))
    level_plots(rows, os.path.join(outdir, "assets"))
