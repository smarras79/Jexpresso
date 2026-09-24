#!/usr/bin/env python3
"""Draw the SEM p-convergence figure for the README (light and dark SVG).

    python3 tools/sem_pconvergence/plot.py

Reads sem_periodic_poisson_pconv.csv (written by sweep.jl) and writes
assets/sem_pconvergence_periodic_poisson{,-dark}.svg. No dependencies.
"""
import csv, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

THEMES = {
    "light": dict(surface="#fcfcfb", text1="#0b0b0b", text2="#52514e", grid="#e4e3dd",
                  axis="#b5b3aa", s1="#2a78d6", s2="#eb6834"),
    "dark":  dict(surface="#1a1a19", text1="#ffffff", text2="#c3c2b7", grid="#34342f",
                  axis="#5d5c55", s1="#3987e5", s2="#d95926"),
}

W, H = 680, 420
L, R, T, B = 78, 128, 70, 62          # plot margins
YMIN, YMAX = -12, -2                   # log10 range
SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")

def read_rows():
    with open(os.path.join(HERE, "sem_periodic_poisson_pconv.csv")) as f:
        return [dict(nop=int(r["nop"]), linf=float(r["linf"]), l2rel=float(r["l2rel"]))
                for r in csv.DictReader(f)]

def svg(rows, th):
    x0, x1 = rows[0]["nop"], rows[-1]["nop"]
    px = lambda n: L + (n - x0) / (x1 - x0) * (W - L - R)
    py = lambda v: T + (YMAX - math.log10(v)) / (YMAX - YMIN) * (H - T - B)
    o = []
    a = o.append
    a(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
      f'role="img" aria-labelledby="t d" font-family="-apple-system, Segoe UI, Helvetica, Arial, sans-serif">')
    a('<title id="t">SEM p-convergence, periodic Poisson problem</title>')
    a('<desc id="d">L-infinity and relative L2 errors of the direct SEM solve of -laplacian u = f on '
      '[0,2pi]^2, doubly periodic, 16 by 16 elements, versus polynomial order 2 to 8. Both fall '
      'exponentially, from about 3e-3 at order 2 to about 2e-11 (L-infinity) and 5e-12 (relative L2) '
      'at order 8.</desc>')
    a(f'<rect width="{W}" height="{H}" fill="{th["surface"]}"/>')
    # title
    a(f'<text x="{L}" y="28" font-size="16" font-weight="600" fill="{th["text1"]}">'
      'SEM error vs polynomial order</text>')
    a(f'<text x="{L}" y="48" font-size="12.5" fill="{th["text2"]}">'
      '−∇²u = f on [0,2π]², doubly periodic, 16×16 elements, direct solve</text>')
    # grid + y labels (every decade line, label every 2)
    for e in range(YMIN, YMAX + 1):
        y = py(10.0 ** e)
        a(f'<line x1="{L}" x2="{W-R}" y1="{y:.1f}" y2="{y:.1f}" stroke="{th["grid"]}" stroke-width="1"/>')
        if e % 2 == 0:
            a(f'<text x="{L-10}" y="{y+4:.1f}" font-size="12" text-anchor="end" fill="{th["text2"]}">'
              f'10{str(e).translate(SUP)}</text>')
    # x axis + ticks
    yb = H - B
    a(f'<line x1="{L}" x2="{W-R}" y1="{yb}" y2="{yb}" stroke="{th["axis"]}" stroke-width="1"/>')
    for r in rows:
        x = px(r["nop"])
        a(f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{yb}" y2="{yb+5}" stroke="{th["axis"]}" stroke-width="1"/>')
        a(f'<text x="{x:.1f}" y="{yb+20}" font-size="12" text-anchor="middle" fill="{th["text2"]}">{r["nop"]}</text>')
    a(f'<text x="{(L + W - R)/2:.1f}" y="{H-14}" font-size="12.5" text-anchor="middle" fill="{th["text2"]}">'
      'polynomial order N  (:nop)</text>')
    a(f'<text x="18" y="{(T + H - B)/2:.1f}" font-size="12.5" text-anchor="middle" fill="{th["text2"]}" '
      f'transform="rotate(-90 18 {(T + H - B)/2:.1f})">error vs exact solution</text>')

    series = [("linf", "‖e‖∞", th["s1"], "circle"), ("l2rel", "relative ‖e‖₂", th["s2"], "square")]
    def marker(kind, x, y, col):
        if kind == "circle":
            return (f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4.5" fill="{col}" '
                    f'stroke="{th["surface"]}" stroke-width="2"/>')
        return (f'<rect x="{x-4.5:.1f}" y="{y-4.5:.1f}" width="9" height="9" rx="1.5" fill="{col}" '
                f'stroke="{th["surface"]}" stroke-width="2"/>')
    for key, label, col, kind in series:
        pts = " ".join(f'{px(r["nop"]):.1f},{py(r[key]):.1f}' for r in rows)
        a(f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="2" '
          f'stroke-linejoin="round" stroke-linecap="round"/>')
        for r in rows:
            a(marker(kind, px(r["nop"]), py(r[key]), col))
        # direct label at the end of the line (text in text ink, mark carries the color)
        yl = py(rows[-1][key])
        a(f'<text x="{px(x1)+12:.1f}" y="{yl+4:.1f}" font-size="12.5" fill="{th["text1"]}">{label}</text>')
    # legend (top-right inside the plot, where the lines are low)
    lx, ly = W - R - 150, T + 14
    for i, (key, label, col, kind) in enumerate(series):
        yy = ly + i * 20
        a(f'<line x1="{lx}" x2="{lx+22}" y1="{yy}" y2="{yy}" stroke="{col}" stroke-width="2"/>')
        a(marker(kind, lx + 11, yy, col))
        a(f'<text x="{lx+30}" y="{yy+4}" font-size="12.5" fill="{th["text1"]}">{label}</text>')
    a('</svg>')
    return "\n".join(o) + "\n"

if __name__ == "__main__":
    rows = read_rows()
    for name, th in THEMES.items():
        suffix = "" if name == "light" else "-dark"
        path = os.path.join(ROOT, "assets", f"sem_pconvergence_periodic_poisson{suffix}.svg")
        with open(path, "w") as f:
            f.write(svg(rows, th))
        print("wrote", path)
