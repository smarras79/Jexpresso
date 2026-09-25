#!/usr/bin/env python3
"""Draw the README figures of the periodic Poisson solver benchmark.

    python3 tools/periodic_poisson_benchmark/plot.py [results.csv] [assets_dir]

Reads results.csv (written by pipeline.jl, second-run timings) and writes, each
in a light and a dark variant (<assets_dir>/<name>.svg and <name>-dark.svg;
assets_dir defaults to the repository's assets/, the README figures):

  ppb_error_vs_dofs          L-inf error vs number of unknowns
  ppb_error_vs_order         L-inf error vs SEM order N (Fourier grids 16N)
  ppb_error_vs_solve_time    L-inf error vs time of the solve step
  ppb_error_vs_total_time    L-inf error vs time-to-solution incl. infrastructure

If the CSV holds several refinement levels (pipeline.jl levels = ...), those
four figures are drawn per level (suffix _L<level>), plus, at every order N
run on two or more levels, time versus unknowns under h-refinement:

  ppb_hrefine_solve_N<N>     solve step
  ppb_hrefine_cost_N<N>      setup + solve (the solver's cost)
  ppb_hrefine_total_N<N>     time-to-solution incl. infrastructure
  ppb_error_vs_Ng            error vs unknowns per direction, sqrt(n) = n_e*N, for
                             all levels (one SEM curve per mesh; pseudo-spectral,
                             FFT and the predicted 0.8^(sqrt(n)/2))

No dependencies.
"""
import csv, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

THEMES = {
    "light": dict(surface="#fcfcfb", text1="#0b0b0b", text2="#52514e", grid="#e4e3dd",
                  axis="#b5b3aa", s=["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]),
    "dark":  dict(surface="#1a1a19", text1="#ffffff", text2="#c3c2b7", grid="#34342f",
                  axis="#5d5c55", s=["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300"]),
}
MARKERS = ["circle", "square", "triangle", "diamond", "tridown", "ring"]
# solver key → legend label, in the fixed categorical order
NAMES = [("sem", "SEM direct"), ("sem_amg", "SEM AMG"), ("sc_direct", "SC direct"),
         ("sc_amg", "SC AMG"), ("ps", "pseudo-spectral"), ("fft", "FFT")]
SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
W, H = 680, 420
L, R, T, B = 78, 150, 108, 62


def read_rows(path):
    with open(path) as f:
        rows = [dict(solver=r["solver"], nop=int(r["nop"]), dofs=int(r["dofs"]),
                     level=int(r.get("level", 0) or 0), linf=float(r["linf"]),
                     solve=float(r["solve"]), cost=float(r["setup"]) + float(r["solve"]),
                     total=float(r["total"]))
                for r in csv.DictReader(f)]
    for r in rows:                                   # a single-level CSV has no nel column
        r["nel"] = round(math.sqrt(r["dofs"]) / r["nop"])
    return rows


def marker(kind, x, y, col, surface):
    ring = f'stroke="{surface}" stroke-width="2"'
    if kind == "circle":
        return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4.5" fill="{col}" {ring}/>'
    if kind == "square":
        return f'<rect x="{x-4.5:.1f}" y="{y-4.5:.1f}" width="9" height="9" rx="1.5" fill="{col}" {ring}/>'
    if kind == "triangle":
        return (f'<path d="M{x:.1f},{y-5.5:.1f} L{x+5.5:.1f},{y+4.5:.1f} L{x-5.5:.1f},{y+4.5:.1f} Z" '
                f'fill="{col}" {ring} stroke-linejoin="round"/>')
    if kind == "tridown":
        return (f'<path d="M{x:.1f},{y+5.5:.1f} L{x+5.5:.1f},{y-4.5:.1f} L{x-5.5:.1f},{y-4.5:.1f} Z" '
                f'fill="{col}" {ring} stroke-linejoin="round"/>')
    if kind == "diamond":
        return (f'<path d="M{x:.1f},{y-6:.1f} L{x+6:.1f},{y:.1f} L{x:.1f},{y+6:.1f} L{x-6:.1f},{y:.1f} Z" '
                f'fill="{col}" {ring} stroke-linejoin="round"/>')
    return (f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{surface}" stroke="{col}" stroke-width="2.5"/>')


def spread(ys, gap=15.0, lo=None, hi=None):
    """Label y-positions near their targets `ys`, at least `gap` apart.

    Overlapping labels are merged into groups, each group centred on the mean
    of its targets (repeated until no two groups overlap), so every label stays
    as close to its own line end as the spacing allows. Optional bounds keep
    the stack inside the plot.
    """
    order = sorted(range(len(ys)), key=lambda i: ys[i])
    groups = [[i] for i in order]                     # each: list of indices, top to bottom

    def place(g):
        c = sum(ys[i] for i in g) / len(g)
        top = c - (len(g) - 1) * gap / 2
        if lo is not None: top = max(top, lo)
        if hi is not None: top = min(top, hi - (len(g) - 1) * gap)
        return top

    merged = True
    while merged:
        merged = False
        for k in range(len(groups) - 1):
            a, b = groups[k], groups[k + 1]
            if place(a) + (len(a) - 1) * gap + gap > place(b):
                groups[k:k + 2] = [a + b]
                merged = True
                break
    out = [0.0] * len(ys)
    for g in groups:
        top = place(g)
        for j, i in enumerate(g):
            out[i] = top + j * gap
    return out


def style(th, slot):
    """slot: an index into the palette/markers, or (color index | None for grey,
    marker kind | None, dashed) -> (colour, marker kind or None, dashed)"""
    if isinstance(slot, tuple):
        c, kind, dash = slot
        return (th["text2"] if c is None else th["s"][c]), kind, dash
    return th["s"][slot], MARKERS[slot], False


def chart(th, title, subtitle, desc, series, xs, xlog, xticks, xlabel, ylabel, ydec,
          end_labels=True, legend_cols=3):
    """series: list of (label, [(x, y), ...], slot); y on a log axis spanning decades ydec."""
    ymin, ymax = ydec
    if xlog:
        lx0, lx1 = math.log10(xs[0]), math.log10(xs[1])
        px = lambda v: L + (math.log10(v) - lx0) / (lx1 - lx0) * (W - L - R)
    else:
        px = lambda v: L + (v - xs[0]) / (xs[1] - xs[0]) * (W - L - R)
    py = lambda v: T + (ymax - math.log10(v)) / (ymax - ymin) * (H - T - B)
    o = []
    a = o.append
    a(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
      f'role="img" aria-labelledby="t d" font-family="-apple-system, Segoe UI, Helvetica, Arial, sans-serif">')
    a(f'<title id="t">{title}</title><desc id="d">{desc}</desc>')
    a(f'<rect width="{W}" height="{H}" fill="{th["surface"]}"/>')
    a(f'<text x="{L}" y="28" font-size="16" font-weight="600" fill="{th["text1"]}">{title}</text>')
    a(f'<text x="{L}" y="48" font-size="12.5" fill="{th["text2"]}">{subtitle}</text>')
    step = 2 if ymax - ymin > 8 else 1
    for e in range(ymin, ymax + 1):
        y = py(10.0 ** e)
        a(f'<line x1="{L}" x2="{W-R}" y1="{y:.1f}" y2="{y:.1f}" stroke="{th["grid"]}" stroke-width="1"/>')
        if (e - ymin) % step == 0:
            a(f'<text x="{L-10}" y="{y+4:.1f}" font-size="12" text-anchor="end" fill="{th["text2"]}">'
              f'10{str(e).translate(SUP)}</text>')
    yb = H - B
    a(f'<line x1="{L}" x2="{W-R}" y1="{yb}" y2="{yb}" stroke="{th["axis"]}" stroke-width="1"/>')
    for v, lab in xticks:
        x = px(v)
        a(f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{yb}" y2="{yb+5}" stroke="{th["axis"]}" stroke-width="1"/>')
        a(f'<text x="{x:.1f}" y="{yb+20}" font-size="12" text-anchor="middle" fill="{th["text2"]}">{lab}</text>')
    a(f'<text x="{(L + W - R)/2:.1f}" y="{H-14}" font-size="12.5" text-anchor="middle" fill="{th["text2"]}">{xlabel}</text>')
    a(f'<text x="18" y="{(T + H - B)/2:.1f}" font-size="12.5" text-anchor="middle" fill="{th["text2"]}" '
      f'transform="rotate(-90 18 {(T + H - B)/2:.1f})">{ylabel}</text>')
    ends = []
    for i, (label, pts, slot) in enumerate(series):
        col, kind, dash = style(th, slot)
        da = ' stroke-dasharray="6 4"' if dash else ""
        a(f'<polyline points="{" ".join(f"{px(x):.1f},{py(y):.1f}" for x, y in pts)}" fill="none" '
          f'stroke="{col}" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"{da}/>')
        if kind:
            for x, y in pts:
                a(marker(kind, px(x), py(y), col, th["surface"]))
        ends.append((label, px(pts[-1][0]), py(pts[-1][1])))
    if end_labels:
        for (label, x, _), y in zip(ends, spread([e[2] for e in ends], lo=T - 4, hi=H - B + 4)):
            a(f'<text x="{x+12:.1f}" y="{y+4:.1f}" font-size="12.5" fill="{th["text1"]}">{label}</text>')
    if len(series) > 1:                        # legend: rows under the subtitle
        colw = 570 // legend_cols
        for i, (label, _, slot) in enumerate(series):
            col, kind, dash = style(th, slot)
            lx, ly = L + (i % legend_cols) * colw, 68 + (i // legend_cols) * 20
            da = ' stroke-dasharray="6 4"' if dash else ""
            a(f'<line x1="{lx}" x2="{lx+22}" y1="{ly}" y2="{ly}" stroke="{col}" stroke-width="2"{da}/>')
            if kind:
                a(marker(kind, lx + 11, ly, col, th["surface"]))
            a(f'<text x="{lx+30}" y="{ly+4}" font-size="12.5" fill="{th["text1"]}">{label}</text>')
    a('</svg>')
    return "\n".join(o) + "\n"


def write(assets, name, make):
    for mode, th in THEMES.items():
        path = os.path.join(assets, name + ("" if mode == "light" else "-dark") + ".svg")
        with open(path, "w") as f:
            f.write(make(th))
        print("wrote", path)


def time_label(e):
    return {-6: "1 µs", -5: "10 µs", -4: "100 µs", -3: "1 ms", -2: "10 ms", -1: "100 ms",
            0: "1 s", 1: "10 s", 2: "100 s"}.get(e, f"10{str(e).translate(SUP)} s")


def level_figures(rows, assets, suffix=""):
    """The four error figures of one mesh (one refinement level)."""
    by = {k: sorted([r for r in rows if r["solver"] == k], key=lambda r: r["nop"]) for k, _ in NAMES}
    names = [(k, lab, slot) for slot, (k, lab) in enumerate(NAMES) if by[k]]
    allr = [r for k, _, _ in names for r in by[k]]
    nel = allr[0]["nel"]
    sub = f"−∇²u = f on [0,2π]², doubly periodic · {nel}×{nel} elements · second-run timings"
    ydec = (math.floor(math.log10(min(r["linf"] for r in allr))),
            math.ceil(math.log10(max(r["linf"] for r in allr))))
    dofs = sorted({r["dofs"] for r in allr})
    nops = sorted({r["nop"] for r in allr})
    dticks = log_ticks(dofs)
    ylab = "L∞ error vs exact solution"

    def series(key):
        return [(lab, [(r[key], r["linf"]) for r in by[k]], slot) for k, lab, slot in names]

    # The four SEM solvers compute the same discrete solution, so their error
    # curves coincide (to 3-4 digits) and would hide one another: plot the
    # error against dofs/order once for the SEM, labelled as such.
    SEM4 = ("sem", "sem_amg", "sc_direct", "sc_amg")

    def error_series(key):
        out, sem_done = [], False
        for k, lab, slot in names:
            if k in SEM4:
                if sem_done:
                    continue
                sem_done, lab = True, "SEM (all 4 solvers)"
            out.append((lab, [(r[key], r["linf"]) for r in by[k]], slot))
        return out

    write(assets, "ppb_error_vs_dofs" + suffix, lambda th: chart(
        th, "Error vs unknowns", sub,
        "L-infinity error versus number of unknowns for the SEM (one curve: its four solvers give the same "
        "discrete solution), pseudo-spectral and FFT solvers.",
        error_series("dofs"), (dofs[0], dofs[-1]), True, dticks, "unknowns (log scale)", ylab, ydec))

    write(assets, "ppb_error_vs_order" + suffix, lambda th: chart(
        th, "Error vs order", sub,
        "L-infinity error versus SEM polynomial order N (one SEM curve for its four solvers, which give the "
        f"same discrete solution); the Fourier solvers are plotted at the order whose SEM grid has the same "
        f"number of unknowns (a {nel}N by {nel}N grid).",
        error_series("nop"), (nops[0], nops[-1]), False, [(n, str(n)) for n in nops],
        f"SEM order N   (Fourier solvers: {nel}N × {nel}N grid, same unknowns)", ylab, ydec))

    for key, name, title, xl in (
            ("solve", "ppb_error_vs_solve_time", "Error vs time of the solve step",
             "solve-step wall-clock (log scale)"),
            ("total", "ppb_error_vs_total_time", "Error vs time-to-solution (with infrastructure)",
             "time-to-solution incl. infrastructure (log scale)")):
        ts = [r[key] for r in allr]
        e0, e1 = math.floor(math.log10(min(ts))), math.ceil(math.log10(max(ts)))
        write(assets, name + suffix, lambda th, key=key, title=title, xl=xl, e0=e0, e1=e1: chart(
            th, title, sub,
            f"L-infinity error versus {xl} for the six solvers (four SEM solves, pseudo-spectral, FFT); "
            "each curve runs over increasing resolution.",
            series(key), (10.0 ** e0, 10.0 ** e1), True,
            [(10.0 ** e, time_label(e)) for e in range(e0, e1 + 1)], xl, ylab, ydec))


def hrefine_figures(rows, assets):
    """Time versus unknowns under h-refinement (one curve per solver) at each
    SEM order N that was run on at least two levels: where the curves of AMG
    and of the direct solves cross, if they do."""
    for n in sorted({r["nop"] for r in rows}):
        rn = [r for r in rows if r["nop"] == n]
        if len({r["level"] for r in rn}) < 2:
            continue
        by = {k: sorted([r for r in rn if r["solver"] == k], key=lambda r: r["level"]) for k, _ in NAMES}
        names = [(k, lab, slot) for slot, (k, lab) in enumerate(NAMES) if len(by[k]) >= 1]
        dofs = sorted({r["dofs"] for r in rn})
        for key, what in (("solve", "solve step"), ("cost", "setup + solve (no SEM infrastructure)"),
                          ("total", "time-to-solution incl. infrastructure")):
            ts = [r[key] for r in rn if r[key] > 0]
            if not ts:
                continue
            ydec = (math.floor(math.log10(min(ts))), math.ceil(math.log10(max(ts))))
            series = [(lab, [(r["dofs"], r[key]) for r in by[k] if r[key] > 0], slot) for k, lab, slot in names]
            series = [x for x in series if x[1]]
            write(assets, f"ppb_hrefine_{key}_N{n}", lambda th, key=key, what=what, series=series, ydec=ydec, dofs=dofs: chart(
                th, f"h-refinement at order N = {n}: {what}",
                "−∇²u = f on [0,2π]², doubly periodic · meshes 16²·4^L elements · second-run timings",
                f"Wall-clock of the {what} versus number of unknowns under uniform mesh refinement at SEM order {n}, "
                "for the six solvers.",
                series, (dofs[0], dofs[-1]), True, log_ticks(dofs), "unknowns (log scale)", "seconds (log scale)", ydec))


def ng_figure(rows, assets, r=0.8):
    """Error versus the number of unknowns per direction, sqrt(n) = n_e*N (the
    same for every solver at a given configuration: SEM nodes per direction =
    Fourier grid points per direction), for every level of a multi-level run:
    one SEM curve per mesh, the pseudo-spectral and FFT results of all levels
    merged (they depend on sqrt(n) only), and the predicted Fourier error
    r^(sqrt(n)/2). Each error is measured against the exact solution at that
    solver's own points. On a linear axis geometric decay is a straight line,
    until the round-off floor."""
    ng = lambda q: q["nel"] * q["nop"]
    levels = sorted({q["level"] for q in rows})
    series = []
    sem_marks = ["circle", "square", "triangle", "diamond", "tridown", "ring"]
    for i, L in enumerate(levels):
        pts = sorted((ng(q), q["linf"]) for q in rows if q["solver"] == "sem" and q["level"] == L)
        if pts:
            n = pts[0][0] // min(q["nop"] for q in rows if q["solver"] == "sem" and q["level"] == L)
            series.append((f"SEM, {n}² el.", pts, (0, sem_marks[i % len(sem_marks)], False)))
    for key, lab, slot in (("ps", "pseudo-spectral", 4), ("fft", "FFT", 5)):
        byng = {}
        for q in rows:
            if q["solver"] == key:
                byng[ng(q)] = q["linf"]              # the same N_g on two levels: same solve
        if byng:
            series.append((lab, sorted(byng.items()), slot))
    allpts = [p for _, pts, _ in series for p in pts]
    if not allpts:
        return
    x0, x1 = min(p[0] for p in allpts), max(p[0] for p in allpts)
    y0 = min(p[1] for p in allpts)
    ydec = (math.floor(math.log10(y0)), math.ceil(math.log10(max(p[1] for p in allpts))))
    pred = [(n, r ** (n / 2)) for n in range(x0, x1 + 1, max(1, (x1 - x0) // 200))
            if r ** (n / 2) >= 10.0 ** ydec[0]]
    if len(pred) > 1:
        series.append((f"{r}^(√n/2)", pred, (None, None, True)))
    step = 64 if x1 - x0 <= 640 else 128
    xticks = [(n, str(n)) for n in range(0, x1 + 1, step) if n >= x0]
    write(assets, "ppb_error_vs_Ng", lambda th: chart(
        th, "Error vs unknowns per direction, all mesh levels",
        "−∇²u = f on [0,2π]², doubly periodic · √n = (elements per side) × N, the same for every solver",
        "L-infinity error versus the number of unknowns per direction, the square root of n = (elements per "
        "side) times N, the same for every solver at a given configuration, for every mesh level: one SEM curve "
        "per mesh, the pseudo-spectral and FFT results of all levels, and the predicted Fourier error "
        f"{r}^(√n/2) (dashed).",
        series, (x0, x1), False, xticks, "unknowns per direction, √n = (elements per side) × N (linear scale)",
        "L∞ error vs exact solution", ydec, end_labels=False, legend_cols=4))


def log_ticks(dofs):
    t = [(d, f"{d:,}".replace(",", " ")) for d in dofs if d in (1024, 2304, 4096, 9216, 16384)]
    if len(t) >= 2 and max(dofs) <= 16384:
        return t
    lo, hi = math.floor(math.log10(dofs[0])), math.ceil(math.log10(dofs[-1]))
    t = [(10.0 ** e, f"10{str(e).translate(SUP)}") for e in range(lo, hi + 1) if dofs[0] <= 10.0 ** e <= dofs[-1]]
    return t or [(d, str(d)) for d in (dofs[0], dofs[-1])]


if __name__ == "__main__":
    import sys
    rows = read_rows(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "results.csv"))
    assets = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "assets")
    os.makedirs(assets, exist_ok=True)
    levels = sorted({r["level"] for r in rows})
    if len(levels) == 1:
        level_figures(rows, assets)                  # the README figures
    else:
        for l in levels:
            level_figures([r for r in rows if r["level"] == l], assets, f"_L{l}")
        hrefine_figures(rows, assets)
        ng_figure(rows, assets)
