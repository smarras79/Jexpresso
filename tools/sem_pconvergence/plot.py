#!/usr/bin/env python3
"""Draw the README figures for the doubly periodic Poisson comparison.

    python3 tools/sem_pconvergence/plot.py

Reads periodic_poisson_solvers.csv (written by sweep.jl) and writes, each in a
light and a dark variant (assets/<name>.svg and assets/<name>-dark.svg):

  sem_pconvergence_periodic_poisson      SEM error vs polynomial order
  periodic_poisson_error_vs_unknowns     SEM / pseudo-spectral / FFT error
  periodic_poisson_time_vs_unknowns      SEM / pseudo-spectral / FFT solve time

No dependencies.
"""
import csv, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

THEMES = {
    "light": dict(surface="#fcfcfb", text1="#0b0b0b", text2="#52514e", grid="#e4e3dd",
                  axis="#b5b3aa", s=["#2a78d6", "#eb6834", "#1baf7a"]),
    "dark":  dict(surface="#1a1a19", text1="#ffffff", text2="#c3c2b7", grid="#34342f",
                  axis="#5d5c55", s=["#3987e5", "#d95926", "#199e70"]),
}
MARKERS = ["circle", "square", "triangle"]
SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
W, H = 680, 420
L, R, T, B = 78, 150, 70, 62


def read_rows():
    with open(os.path.join(HERE, "periodic_poisson_solvers.csv")) as f:
        rows = []
        for r in csv.DictReader(f):
            rows.append(dict(solver=r["solver"], nop=int(r["nop"]), unknowns=int(r["unknowns"]),
                             linf=float(r["linf"]), l2rel=float(r["l2rel"]),
                             solve_s=float(r["solve_s"])))
        return rows


def marker(kind, x, y, col, surface):
    ring = f'stroke="{surface}" stroke-width="2"'
    if kind == "circle":
        return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4.5" fill="{col}" {ring}/>'
    if kind == "square":
        return f'<rect x="{x-4.5:.1f}" y="{y-4.5:.1f}" width="9" height="9" rx="1.5" fill="{col}" {ring}/>'
    return (f'<path d="M{x:.1f},{y-5.5:.1f} L{x+5.5:.1f},{y+4.5:.1f} L{x-5.5:.1f},{y+4.5:.1f} Z" '
            f'fill="{col}" {ring} stroke-linejoin="round"/>')


def chart(th, title, subtitle, desc, series, xs, xlog, xticks, xlabel, ylabel, ydec, legend="right"):
    """series: list of (label, [(x, y), ...]); y on a log axis spanning decades ydec."""
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
    for i, (label, pts) in enumerate(series):
        col, kind = th["s"][i], MARKERS[i]
        a(f'<polyline points="{" ".join(f"{px(x):.1f},{py(y):.1f}" for x, y in pts)}" fill="none" '
          f'stroke="{col}" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>')
        for x, y in pts:
            a(marker(kind, px(x), py(y), col, th["surface"]))
        a(f'<text x="{px(pts[-1][0])+12:.1f}" y="{py(pts[-1][1])+4:.1f}" font-size="12.5" '
          f'fill="{th["text1"]}">{label}</text>')
    if len(series) > 1:                        # legend, top corner inside the plot
        lx, ly = (W - R - 175 if legend == "right" else L + 16), T + 14
        for i, (label, _) in enumerate(series):
            col, yy = th["s"][i], ly + i * 20
            a(f'<line x1="{lx}" x2="{lx+22}" y1="{yy}" y2="{yy}" stroke="{col}" stroke-width="2"/>')
            a(marker(MARKERS[i], lx + 11, yy, col, th["surface"]))
            a(f'<text x="{lx+30}" y="{yy+4}" font-size="12.5" fill="{th["text1"]}">{label}</text>')
    a('</svg>')
    return "\n".join(o) + "\n"


def write(name, make):
    for mode, th in THEMES.items():
        path = os.path.join(ROOT, "assets", name + ("" if mode == "light" else "-dark") + ".svg")
        with open(path, "w") as f:
            f.write(make(th))
        print("wrote", path)


if __name__ == "__main__":
    rows = read_rows()
    by = {s: sorted([r for r in rows if r["solver"] == s], key=lambda r: r["nop"])
          for s in ("sem", "ps", "fft")}
    sem = by["sem"]
    nops = [r["nop"] for r in sem]
    sub = "−∇²u = f on [0,2π]², doubly periodic, u = sin 2x cos 3y + sin x cos y"
    uk = [r["unknowns"] for r in sem]
    uticks = [(u, f"{u:,}".replace(",", " ")) for u in uk if u in (1024, 2304, 4096, 9216, 16384)]
    names = [("sem", "SEM (16×16 el.)"), ("ps", "pseudo-spectral"), ("fft", "FFT")]

    # 1. SEM p-convergence (the curve the README section opened with)
    write("sem_pconvergence_periodic_poisson", lambda th: chart(
        dict(th, s=[th["s"][0], th["s"][1]]),
        "SEM error vs polynomial order",
        "−∇²u = f on [0,2π]², doubly periodic, 16×16 elements, direct solve",
        "L-infinity and relative L2 errors of the direct SEM solve versus polynomial order 2 to 8: "
        "both fall exponentially, from about 3e-3 to about 2e-11 and 5e-12.",
        [("‖e‖∞", [(r["nop"], r["linf"]) for r in sem]),
         ("relative ‖e‖₂", [(r["nop"], r["l2rel"]) for r in sem])],
        (nops[0], nops[-1]), False, [(n, str(n)) for n in nops],
        "polynomial order N  (:nop)", "error vs exact solution", (-12, -2)))

    # 2. error vs unknowns, three solvers
    write("periodic_poisson_error_vs_unknowns", lambda th: chart(
        th, "Error vs unknowns: SEM, pseudo-spectral, FFT", sub,
        "L-infinity error versus number of unknowns (1024 to 16384). SEM falls exponentially from "
        "3e-3 to 2e-11; pseudo-spectral and FFT are at round-off (1e-14 to 7e-13, and 2e-15 to 3e-15) "
        "at every size because the exact solution is a trigonometric polynomial.",
        [(lab, [(r["unknowns"], r["linf"]) for r in by[s]]) for s, lab in names],
        (uk[0], uk[-1]), True, uticks, "unknowns (log scale)", "L∞ error vs exact solution", (-16, -2)))

    # 3. solve time vs unknowns, three solvers
    write("periodic_poisson_time_vs_unknowns", lambda th: chart(
        th, "Solve time vs unknowns: SEM, pseudo-spectral, FFT", sub,
        "Solve time versus number of unknowns (1024 to 16384), BenchmarkTools minimum. SEM (sparse LU "
        "factorisation plus solve) grows from 1.2 ms to 88 ms; pseudo-spectral (dense products, "
        "setup excluded) from 5 microseconds to 0.25 ms; FFT (plan excluded) from 2.4 to 53 microseconds.",
        [(lab, [(r["unknowns"], r["solve_s"]) for r in by[s]]) for s, lab in names],
        (uk[0], uk[-1]), True, uticks, "unknowns (log scale)", "solve time [s]", (-6, 0),
        legend="left"))
