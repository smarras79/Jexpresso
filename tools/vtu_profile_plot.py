#!/usr/bin/env python3
"""Vertical-profile figure from Jexpresso snapshots: the standard CBL panels
plus a zoom on the inversion, with the initial sounding overlaid on theta.

Usage: vtu_profile_plot.py <output_dir> <out.png> <iter> [<iter> ...]
Env: ZTOP (default 2000), SOUNDING (path; default the teamx u10 flat file),
     TLABEL (text for the title, e.g. "t = 2500-3500 s").
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vtu_stats import snapshot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

outdir, outpng, iters = sys.argv[1], sys.argv[2], sys.argv[3:]
ztop  = float(os.environ.get('ZTOP', '2000'))
snd   = os.environ.get('SOUNDING', 'data_files/input_sounding_teamx_u10_flat_noheader.dat')
tlab  = os.environ.get('TLABEL', f'iters {",".join(iters)}')

acc = None
for it in iters:
    r = snapshot(outdir, it)
    acc = r if acc is None else acc + r
acc /= len(iters); acc[:, 0] = snapshot(outdir, iters[0])[:, 0]
z, u, v, th, uu, vv, ww, uw, wth = acc.T
m = z <= ztop

# z_i from the maximum <theta> gradient (robust on a small domain, where the
# flux minimum is noisy); the flux minimum is reported alongside.
dthz = np.gradient(th[m], z[m])
zi   = z[m][np.argmax(dthz)]
kf   = np.argmin(wth[m])
band = (z >= zi - 150) & (z <= zi + 150)
# inversion thickness: span where d<theta>/dz exceeds a quarter of its peak
core = dthz > 0.25*dthz.max()
print(f'z_i (max dtheta/dz) = {zi:.0f} m;  inversion layer {z[m][core].min():.0f}-{z[m][core].max():.0f} m '
      f'({z[m][core].max()-z[m][core].min():.0f} m thick);  peak gradient {1000*dthz.max():.1f} K/km')
print(f'min <w theta> = {wth[m][kf]:.4f} at {z[m][kf]:.0f} m  (surface flux 0.12 K m/s, ratio {wth[m][kf]/0.12:+.2f})')
print(f'{"z":>8} {"<theta>":>9} {"<w theta>":>10} {"<ww>":>8}')
for zz, t, f, w2 in zip(z[band], th[band], wth[band], ww[band]):
    print(f'{zz:8.1f} {t:9.3f} {f:10.4f} {w2:8.4f}')

sounding = None
if os.path.isfile(snd):
    s = np.loadtxt(snd); sounding = (s[:, 0], s[:, 1], s[:, 3])

fig, ax = plt.subplots(1, 6, figsize=(22, 7), sharey=False)
ax[0].plot(u[m], z[m], label='<u>'); ax[0].plot(v[m], z[m], label='<v>')
if sounding: ax[0].plot(sounding[2], sounding[0], 'k--', lw=0.8, label='t=0')
ax[0].set_xlabel('m/s'); ax[0].legend(); ax[0].set_title('mean wind')

ax[1].plot(th[m], z[m], label=f'<θ> {tlab}')
if sounding: ax[1].plot(sounding[1], sounding[0], 'k--', lw=0.8, label='t=0 sounding')
ax[1].set_xlim(299, 308); ax[1].set_xlabel('K'); ax[1].legend(); ax[1].set_title('potential temperature')

ax[2].plot(th[band], z[band], 'o-', ms=3, label='<θ>')
if sounding:
    sb = (sounding[0] >= zi - 150) & (sounding[0] <= zi + 150)
    ax[2].plot(sounding[1][sb], sounding[0][sb], 'k--', lw=0.8, label='t=0')
ax[2].axhline(zi, color='r', lw=0.6, ls=':', label=f'z_i = {zi:.0f} m')
ax[2].set_xlabel('K'); ax[2].legend(); ax[2].set_title('inversion zoom (±150 m)')

ax[3].plot(-uw[m], z[m], label="-<u'w'>"); ax[3].plot(wth[m], z[m], label="<w'θ'>")
ax[3].axvline(0, color='k', lw=0.5); ax[3].set_xlabel("m²/s²  |  K m/s"); ax[3].legend(); ax[3].set_title('resolved fluxes')

ax[4].plot(uu[m], z[m], label="<u'u'>"); ax[4].plot(vv[m], z[m], label="<v'v'>"); ax[4].plot(ww[m], z[m], label="<w'w'>")
ax[4].set_xlabel('m²/s²'); ax[4].legend(); ax[4].set_title('resolved variances')

# theta gradient: where the inversion actually is and how thick
dth = np.gradient(th[m], z[m]) * 1000.0
ax[5].plot(dth, z[m]); ax[5].axhline(zi, color='r', lw=0.6, ls=':')
ax[5].set_xlabel('K/km'); ax[5].set_title('d<θ>/dz')
for a in ax: a.set_ylim(0, ztop); a.grid(alpha=0.3)
ax[2].set_ylim(zi - 150, zi + 150)
ax[0].set_ylabel('z [m]')
fig.suptitle(f'{outdir.split("/")[-3] if "/" in outdir else outdir}   {tlab}   (horizontal mean, resolved part)')
fig.tight_layout(); fig.savefig(outpng, dpi=110); print('wrote', outpng)
