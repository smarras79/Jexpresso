#!/usr/bin/env python3
"""Merge and plot the 3D periodic Poisson solver comparison.

    python3 tools/poisson3d_benchmark/plot3d.py OUTDIR

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

No dependencies (uses the SVG chart of ../periodic_poisson_benchmark/plot.py).
"""
import csv, glob, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "periodic_poisson_benchmark"))
import plot as P  # noqa: E402

SOLVERS = [  # key, label, (colour index, marker)
    ("sem",        "SEM direct",  (0, "circle")),
    ("sem_amg",    "SEM AMG",     (1, "square")),
    ("sem_jacobi", "SEM Jacobi",  (6, "square")),
    ("sc_direct",  "SC direct",   (2, "triangle")),
    ("sc_amg",     "SC AMG",      (3, "diamond")),
    ("ps",         "pseudo-sp.",  (4, "tridown")),
    ("fft",        "FFT",         (5, "ring")),
    # MPI solvers (mpi/bench3d_mpi.jl), same discretisation as the SEM solvers
    ("mumps",      "MUMPS (MPI)",        (0, "circle")),
    ("boomeramg",  "BoomerAMG-CG (MPI)", (1, "square")),
    ("jacobi",     "Jacobi-CG (MPI)",    (6, "square")),
]
SEM = {"sem", "sem_amg", "sem_jacobi", "sc_direct", "sc_amg", "mumps", "boomeramg", "jacobi"}
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
            its = ser("iters", solvers={"sem_amg", "sem_jacobi", "sc_amg", "boomeramg", "jacobi"})
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


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    outdir = sys.argv[1]
    rows = merge(outdir)
    write_md(os.path.join(outdir, "results.md"), rows)
    plots(rows, os.path.join(outdir, "assets"))
