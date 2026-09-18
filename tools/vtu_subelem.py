#!/usr/bin/env python3
"""Sub-element energy fraction per height: how much of the w (and u) variance
lives at scales below one element. Element-locked noise shows as a fraction
that jumps well above the turbulent-layer value.

Usage: vtu_subelem.py <output_dir> <iter>      Env: HELEM (element size in x,
default 80), ZLEVELS (default "20,100,500,900,1100,1300,1500,1800").
"""
import sys, os, glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vtu_merge import read_vtu

outdir, it = sys.argv[1], sys.argv[2]
helem = float(os.environ.get('HELEM', '80'))
zlev  = [float(s) for s in os.environ.get('ZLEVELS', '20,100,500,900,1100,1300,1500,1800').split(',')]
pieces = sorted(glob.glob(os.path.join(outdir, f'iter_{it}', f'iter_{it}_*.vtu')))
X, Y, Z, U, W = [], [], [], [], []
for p in pieces:
    a = read_vtu(p)
    X.append(a['Points'][:, 0]); Y.append(a['Points'][:, 1]); Z.append(a['Points'][:, 2])
    U.append(a['u']); W.append(a['w'])
x = np.concatenate(X); y = np.concatenate(Y); z = np.concatenate(Z)
u = np.concatenate(U); w = np.concatenate(W)
zr = np.round(z, 3); zl = np.unique(zr)

def frac(f, m):
    ff = f[m] - f[m].mean()
    # element id in x and y from the coordinates (LGL points on element faces
    # belong to two elements; floor puts them in one, which is fine here)
    ex = np.floor(x[m] / helem + 1e-6).astype(int); ey = np.floor(y[m] / helem + 1e-6).astype(int)
    key = ex * 100000 + ey
    order = np.argsort(key); ks = key[order]; fs = ff[order]
    bounds = np.flatnonzero(np.diff(ks)) + 1
    sums = np.add.reduceat(fs, np.concatenate([[0], bounds]))
    cnts = np.diff(np.concatenate([[0], bounds, [len(fs)]]))
    emean = np.repeat(sums / cnts, cnts)
    resid = fs - emean
    tot = (ff**2).mean()
    return (resid**2).mean() / tot if tot > 0 else 0.0, np.sqrt(tot)

print(f'{"z":>8} {"rms w":>8} {"sub-elem frac w":>16} {"rms up":>8} {"sub-elem frac u":>16}')
for h in zlev:
    zz = zl[np.argmin(np.abs(zl - h))]; m = zr == zz
    fw, rw = frac(w, m); fu, ru = frac(u, m)
    print(f'{zz:8.1f} {rw:8.4f} {fw:16.3f} {ru:8.4f} {fu:16.3f}')
print('\nRead it RELATIVE: the fraction above the inversion against the one inside the')
print('PBL. Physical gravity waves are larger than an element, so a clean free')
print('atmosphere sits at or below the PBL value; element-locked noise pushes it up.')
print('(4x4x60 laptop reference, t=4400: PBL 0.4-0.7, above 1.1 km 0.22-0.27.)')
