#!/usr/bin/env python3
"""Horizontally averaged near-wall profiles from a Jexpresso .vtu set.

Usage: prof.py <output_dir> <iter_N> [<output_dir> <iter_N> ...]
Reads every rank piece iter_N_*.vtu, bins by z, prints the lowest levels.
"""
import sys, io, re, struct, glob, os
import numpy as np

def read_vtu(path):
    raw = io.open(path, 'rb').read()
    hdr_end = raw.find(b'<AppendedData')
    hdr = raw[:hdr_end].decode('latin1')
    j = raw.find(b'_', hdr_end) + 1
    blob = raw[j:]

    arrays = {}
    for m in re.finditer(
            r'<DataArray type="(\w+)" Name="([^"]*)" NumberOfComponents="(\d+)" format="appended" offset="(\d+)"',
            hdr):
        typ, name, ncomp, off = m.group(1), m.group(2), int(m.group(3)), int(m.group(4))
        if typ != 'Float64':
            continue
        n = struct.unpack('<Q', blob[off:off+8])[0]
        a = np.frombuffer(blob[off+8:off+8+n], dtype='<f8')
        arrays[name] = a.reshape(-1, ncomp) if ncomp > 1 else a
    return arrays

def profiles(outdir, it):
    pieces = sorted(glob.glob(os.path.join(outdir, f'iter_{it}', f'iter_{it}_*.vtu')))
    if not pieces:
        raise SystemExit(f'no pieces under {outdir}/iter_{it}')
    Z, U, V, W = [], [], [], []
    for p in pieces:
        a = read_vtu(p)
        pts = a['Points']
        Z.append(pts[:, 2])
        names = {k.lower(): k for k in a}
        def get(*cands):
            for c in cands:
                if c in names:
                    return a[names[c]]
            return None
        u = get('u'); v = get('v'); w = get('w')
        if u is None:
            raise SystemExit(f'no velocity in {p}; have {sorted(a)}')
        U.append(u); V.append(v if v is not None else np.zeros_like(u))
        W.append(w if w is not None else np.zeros_like(u))
    z = np.concatenate(Z); u = np.concatenate(U)
    v = np.concatenate(V); w = np.concatenate(W)
    lev = np.unique(np.round(z, 4))
    rows = []
    for zz in lev:
        m = np.round(z, 4) == zz
        rows.append((zz, u[m].mean(), v[m].mean(),
                     np.sqrt(u[m]**2 + v[m]**2).mean(),
                     np.abs(w[m]).max(),
                     np.sqrt(u[m]**2 + v[m]**2).max()))
    return rows

if __name__ == '__main__':
    args = sys.argv[1:]
    labels, allrows = [], []
    for k in range(0, len(args), 2):
        labels.append(os.path.basename(os.path.dirname(os.path.dirname(args[k]))) + f'/it{args[k+1]}')
        allrows.append(profiles(args[k], args[k+1]))
    nlev = int(os.environ.get('NLEV', '8'))
    for lbl, rows in zip(labels, allrows):
        print(f'=== {lbl}')
        print(f'{"z":>9} {"<u>":>9} {"<v>":>9} {"<|uh|>":>9} {"max|w|":>9} {"max|uh|":>9}')
        for r in rows[:nlev]:
            print(f'{r[0]:9.3f} {r[1]:9.4f} {r[2]:9.4f} {r[3]:9.4f} {r[4]:9.4f} {r[5]:9.4f}')
        print()
