#!/usr/bin/env python3
"""PNG slices of one Jexpresso snapshot, straight from the rank pieces.

Usage: vtu_slices.py <output_dir> <iter> <prefix>

Writes <prefix>_xy_z<h>.png (horizontal slices of w and θ at the heights in
ZLEVELS, default "20,100,500") and <prefix>_xz.png (vertical slice of w, θ
and u at the y in YCUT, default the domain midline, up to ZTOP, default
2000 m). Nearest LGL level is used, no interpolation. numpy + matplotlib.
"""
import sys, os, glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vtu_merge import read_vtu
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

outdir, it, prefix = sys.argv[1], sys.argv[2], sys.argv[3]
zlevels = [float(s) for s in os.environ.get('ZLEVELS', '20,100,500').split(',')]
ztop    = float(os.environ.get('ZTOP', '2000'))

pieces = sorted(glob.glob(os.path.join(outdir, f'iter_{it}', f'iter_{it}_*.vtu')))
if not pieces:
    raise SystemExit(f'no pieces under {outdir}/iter_{it}')
X, Y, Z, U, W, T = [], [], [], [], [], []
for p in pieces:
    a = read_vtu(p)
    X.append(a['Points'][:, 0]); Y.append(a['Points'][:, 1]); Z.append(a['Points'][:, 2])
    U.append(a['u']); W.append(a['w']); T.append(a['θ'] if 'θ' in a else a['theta'])
x = np.concatenate(X); y = np.concatenate(Y); z = np.concatenate(Z)
u = np.concatenate(U); w = np.concatenate(W); t = np.concatenate(T)
zl = np.unique(np.round(z, 3)); yl = np.unique(np.round(y, 3))
ycut = float(os.environ.get('YCUT', yl[len(yl)//2]))
print(f'{len(x)} points; z levels {len(zl)}; y cut at {ycut:.1f}')

def nearest(vals, target):
    return vals[np.argmin(np.abs(vals - target))]

def scatter(ax, xx, yy, ff, title, cmap, sym, equal=True):
    # tricontourf on the scattered LGL points: filled, no gaps between nodes
    if sym:
        lim = np.percentile(np.abs(ff), 99); lev = np.linspace(-lim, lim, 41)
    else:
        lo, hi = np.percentile(ff, [0.5, 99.5]); lev = np.linspace(lo, hi, 41)
    sc = ax.tricontourf(xx, yy, ff, levels=lev, cmap=cmap, extend='both')
    ax.set_title(title); equal and ax.set_aspect('equal'); plt.colorbar(sc, ax=ax, shrink=0.8)

for h in zlevels:
    zz = nearest(zl, h); m = np.round(z, 3) == zz
    fig, ax = plt.subplots(1, 2, figsize=(14, 6))
    scatter(ax[0], x[m], y[m], w[m], f'w [m/s] at z = {zz:.1f} m, t = iter {it}', 'RdBu_r', True)
    scatter(ax[1], x[m], y[m], t[m], f'θ [K] at z = {zz:.1f} m', 'inferno', False)
    for a_ in ax: a_.set_xlabel('x [m]'); a_.set_ylabel('y [m]')
    f = f'{prefix}_xy_z{int(zz)}.png'; fig.savefig(f, dpi=110, bbox_inches='tight'); plt.close(fig)
    print('wrote', f)

yy = nearest(yl, ycut); m = (np.round(y, 3) == yy) & (z <= ztop)
fig, ax = plt.subplots(3, 1, figsize=(14, 12))
scatter(ax[0], x[m], z[m], w[m], f'w [m/s], y = {yy:.0f} m', 'RdBu_r', True, equal=False)
scatter(ax[1], x[m], z[m], t[m], 'θ [K]', 'inferno', False, equal=False)
scatter(ax[2], x[m], z[m], u[m], 'u [m/s]', 'viridis', False, equal=False)
for a_ in ax: a_.set_xlabel('x [m]'); a_.set_ylabel('z [m]')
f = f'{prefix}_xz.png'; fig.savefig(f, dpi=110, bbox_inches='tight'); plt.close(fig)
print('wrote', f)
