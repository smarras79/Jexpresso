#!/usr/bin/env python3
"""Horizontal-mean profiles and resolved second moments from Jexpresso .vtu dumps.

Usage: vtu_stats.py <output_dir> <iter> [<iter> ...]

For every z level: <u> <v> <theta> <u'u'> <v'v'> <w'w'> <u'w'> <w'theta'>,
fluctuations taken about that snapshot's own horizontal mean, then averaged
over the snapshots given. Ends with a log-law fit of |<u_h>| on the lowest
NFIT levels above the wall (z0 from Z0, default 0.1 m).

Env: NLEV (levels to print, default 40), NFIT (default 4), Z0.
"""
import sys, io, re, struct, glob, os
import numpy as np

def read_vtu(path):
    raw = io.open(path, 'rb').read()
    hdr_end = raw.find(b'<AppendedData')
    hdr = raw[:hdr_end].decode('utf-8', errors='replace')
    j = raw.find(b'_', hdr_end) + 1
    blob = raw[j:]
    out = {}
    for m in re.finditer(r'<DataArray type="(\w+)" Name="([^"]*)" NumberOfComponents="(\d+)" format="appended" offset="(\d+)"', hdr):
        typ, name, nc, off = m.group(1), m.group(2), int(m.group(3)), int(m.group(4))
        if typ != 'Float64':
            continue
        n = struct.unpack('<Q', blob[off:off+8])[0]
        a = np.frombuffer(blob[off+8:off+8+n], dtype='<f8')
        out[name] = a.reshape(-1, nc) if nc > 1 else a
    return out

def pick(a, *names):
    for n in names:
        if n in a:
            return a[n]
    raise SystemExit(f'missing {names}; have {sorted(a)}')

def snapshot(outdir, it):
    pieces = sorted(glob.glob(os.path.join(outdir, f'iter_{it}', f'iter_{it}_*.vtu')))
    if not pieces:
        raise SystemExit(f'no pieces under {outdir}/iter_{it}')
    Z, U, V, W, T = [], [], [], [], []
    for p in pieces:
        a = read_vtu(p)
        Z.append(a['Points'][:, 2])
        U.append(pick(a, 'u')); V.append(pick(a, 'v')); W.append(pick(a, 'w'))
        T.append(pick(a, 'θ', 'theta', 't'))
    z = np.concatenate(Z); u = np.concatenate(U); v = np.concatenate(V)
    w = np.concatenate(W); t = np.concatenate(T)
    zr = np.round(z, 4)
    lev = np.unique(zr)
    rows = np.zeros((len(lev), 9))
    for k, zz in enumerate(lev):
        m = zr == zz
        uu, vv, ww, tt = u[m], v[m], w[m], t[m]
        um, vm, wm, tm = uu.mean(), vv.mean(), ww.mean(), tt.mean()
        up, vp, wp, tp = uu-um, vv-vm, ww-wm, tt-tm
        rows[k] = (zz, um, vm, tm, (up*up).mean(), (vp*vp).mean(), (wp*wp).mean(),
                   (up*wp).mean(), (wp*tp).mean())
    return rows

if __name__ == '__main__':
    outdir, iters = sys.argv[1], sys.argv[2:]
    acc = None
    for it in iters:
        r = snapshot(outdir, it)
        acc = r if acc is None else acc + r
    acc /= len(iters)
    acc[:, 0] = snapshot(outdir, iters[0])[:, 0]   # z is not to be averaged
    nlev = int(os.environ.get('NLEV', '40'))
    print(f'averaged over iters {iters}')
    print(f'{"z":>8} {"<u>":>8} {"<v>":>8} {"<th>":>8} {"<uu>":>8} {"<vv>":>8} {"<ww>":>8} {"<uw>":>9} {"<wth>":>9}')
    for r in acc[:nlev]:
        print(f'{r[0]:8.2f} {r[1]:8.3f} {r[2]:8.3f} {r[3]:8.3f} {r[4]:8.4f} {r[5]:8.4f} {r[6]:8.4f} {r[7]:9.4f} {r[8]:9.4f}')
    # log-law fit on levels 1..NFIT above the wall
    nfit = int(os.environ.get('NFIT', '4')); z0 = float(os.environ.get('Z0', '0.1'))
    zz = acc[1:1+nfit, 0]; uh = np.hypot(acc[1:1+nfit, 1], acc[1:1+nfit, 2])
    A = np.vstack([np.log(zz/z0), np.ones_like(zz)]).T
    slope, icpt = np.linalg.lstsq(A, uh, rcond=None)[0]
    ustar = 0.4*slope
    print(f'\nlog-law fit on z = {zz[0]:.1f}..{zz[-1]:.1f} m (z0 = {z0}):  u* = {ustar:.3f} m/s,'
          f'  intercept {icpt:+.2f} m/s (0 if the profile is a pure log law)')
    pred = slope*np.log(zz/z0) + icpt
    for a, b, c in zip(zz, uh, pred):
        print(f'   z={a:6.2f}  |<uh>|={b:6.3f}  fit={c:6.3f}  ratio={b/c:5.3f}')
    # surface-layer consistency: -<uw> at the first levels vs u*^2
    print(f'\n-<u\'w\'>_res at z={acc[1,0]:.1f}: {-acc[1,7]:.4f}   at z={acc[2,0]:.1f}: {-acc[2,7]:.4f}   u*^2 from fit: {ustar**2:.4f}')
    print(f'<w\'th\'>_res at z={acc[1,0]:.1f}: {acc[1,8]:.4f}   at z={acc[2,0]:.1f}: {acc[2,8]:.4f}   (imposed surface flux 0.12 K m/s; resolved part should be a fraction of it near the wall)')
