#!/usr/bin/env python3
"""Merge the per-solver results of a 2D rerun and compare them with earlier runs.

    python3 tools/periodic_poisson_benchmark/compare2d.py OUTDIR [REFERENCE.csv ...]

1. Merges OUTDIR/parts/*/results.csv (one per solver, written by pipeline.jl)
   into OUTDIR/results.csv, sorted by level, order and solver.
2. For every REFERENCE.csv (a results.csv of an earlier run, e.g. last week's
   ppb_highres/results.csv; default: the committed 16x16 one,
   tools/periodic_poisson_benchmark/results.csv), compares the L-inf and
   relative L2 errors configuration by configuration (solver, level, N) and
   writes OUTDIR/compare.md: per solver the largest |difference|, and every
   configuration whose L-inf error differs by more than 1e-10 (absolute; u
   peaks at 1). Exits with status 1 if any does.

Expected: pseudo-spectral and FFT bit for bit; the four SEM solvers to
round-off (~1e-12). The SEM rows of runs made before the deck switched from
the GMSH mesh to the built-in Cartesian grid (commit b9c76d1) differ by that
much too: the two grids agree to round-off, not bit for bit.
"""
import csv, glob, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ORDER = ["sem", "sem_amg", "sc_direct", "sc_amg", "ps", "fft"]
TOL = 1e-10


def read(path):
    with open(path) as fh:
        rd = csv.DictReader(fh)
        rows = list(rd)
        return rd.fieldnames, rows


def key(r):
    return (r["solver"], int(r.get("level") or 0), int(r["nop"]))


def merge(outdir):
    files = sorted(glob.glob(os.path.join(outdir, "parts", "*", "results.csv")))
    if not files:
        sys.exit(f"compare2d.py: no {outdir}/parts/*/results.csv")
    hdr, byk = None, {}
    for f in files:
        h, rows = read(f)
        hdr = hdr or h
        for r in rows:
            byk[key(r)] = r
    rows = sorted(byk.values(), key=lambda r: (int(r.get("level") or 0), int(r["nop"]),
                                               ORDER.index(r["solver"]) if r["solver"] in ORDER else 99))
    out = os.path.join(outdir, "results.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=hdr)
        w.writeheader()
        w.writerows(rows)
    print(f"merged {len(rows)} configurations from {len(files)} file(s) -> {out}")
    return rows


def compare(new, refpath):
    _, ref = read(refpath)
    refk = {key(r): r for r in ref}
    lines, bad, per = [], [], {}
    for r in new:
        k = key(r)
        if k not in refk:
            continue
        d = abs(float(r["linf"]) - float(refk[k]["linf"]))
        d2 = abs(float(r["l2rel"]) - float(refk[k]["l2rel"]))
        s = per.setdefault(r["solver"], [0, 0.0, 0.0, 0])
        s[0] += 1; s[1] = max(s[1], d); s[2] = max(s[2], d2); s[3] += (d == 0.0)
        if d > TOL:
            bad.append((k, float(refk[k]["linf"]), float(r["linf"]), d))
    missing = sorted(set(refk) - {key(r) for r in new})
    lines.append(f"## against `{refpath}`\n")
    lines.append(f"{sum(s[0] for s in per.values())} configurations in both runs; "
                 f"{len(missing)} in the reference only (not rerun).\n")
    lines.append("| solver | compared | bit-for-bit | max abs diff L∞ | max abs diff rel. L2 |")
    lines.append("|---|---:|---:|---:|---:|")
    for s in ORDER:
        if s in per:
            n, d, d2, same = per[s]
            lines.append(f"| {s} | {n} | {same} | {d:.1e} | {d2:.1e} |")
    if bad:
        lines.append(f"\n**{len(bad)} configuration(s) differ by more than {TOL:g} in L∞:**\n")
        lines.append("| solver | level | N | reference L∞ | new L∞ | abs diff |")
        lines.append("|---|---:|---:|---:|---:|---:|")
        for (s, l, n), a, b, d in bad:
            lines.append(f"| {s} | {l} | {n} | {a:.6e} | {b:.6e} | {d:.1e} |")
    else:
        lines.append(f"\nEvery compared configuration agrees to {TOL:g} in L∞.")
    return "\n".join(lines) + "\n", not bad


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    outdir = sys.argv[1]
    refs = sys.argv[2:] or [os.path.join(HERE, "results.csv")]
    new = merge(outdir)
    report, ok = ["# 2D rerun compared with earlier runs\n"], True
    for ref in refs:
        if not os.path.isfile(ref):
            report.append(f"## `{ref}`: not found, skipped\n"); continue
        txt, good = compare(new, ref)
        report.append(txt); ok &= good
    path = os.path.join(outdir, "compare.md")
    with open(path, "w") as fh:
        fh.write("\n".join(report))
    print("\n".join(report))
    print("wrote", path)
    sys.exit(0 if ok else 1)
