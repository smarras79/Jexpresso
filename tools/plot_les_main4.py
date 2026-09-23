#!/usr/bin/env python3
"""The four profiles that carry the intercomparison, from the built-in tavg .dat.

Usage: plot_les_main4.py <les_statistics_tavg.dat> <les_stress_tavg.dat>
                         [les_ustar_tavg.dat] [out.png]

  1  <theta>            with an inversion inset
  2  <w'theta'>         resolved + SFS, normalised by the surface flux
  3  <w'w'>/w*^2        against z/z_i
  4  <u>, <v>

Element-scale zigzag in the vertical is removed by averaging each spectral
element's GLL nodes (ELEM, default 160 m; 0 disables) -- the SGS fluxes are
nu_t times a gradient and the gradient is discontinuous across element faces.
Env: ZTOP (2000), WTHETA_S (0.12), ELEM (160).
"""
import sys, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

def load(p):
    names = open(p).readline().lstrip('#').split()[2:]
    d = np.loadtxt(p); d = d[:, 1:] if d.shape[1] == len(names)+1 else d
    return {n: d[:, i] for i, n in enumerate(names)}

stat, strs = load(sys.argv[1]), load(sys.argv[2])
ust = load(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3].endswith('.dat') else None
out = sys.argv[-1] if sys.argv[-1].endswith('.png') else 'les_main4.png'
ztop = float(os.environ.get('ZTOP','2000')); Q0 = float(os.environ.get('WTHETA_S','0.12'))
elem = float(os.environ.get('ELEM','160'))

z = stat['z']
def esmooth(v):
    """Element average placed at the element centre, then interpolated back onto
    the nodes. Replacing the nodes by their element mean instead would turn the
    zigzag into a staircase, which is not less wrong, only differently wrong."""
    if elem <= 0: return v
    idx = np.floor(z/elem + 1e-9).astype(int)
    e   = np.unique(idx)
    zc  = np.array([z[idx == k].mean() for k in e])
    vc  = np.array([v[idx == k].mean() for k in e])
    return np.interp(z, zc, vc)

th  = stat['t_mean']; u = stat['u_mean']; v = stat['v_mean']
wth_r, wth_s = esmooth(strs['wptp_res']), esmooth(strs['wptp_sfs'])
wth = wth_r + wth_s
ww  = strs['wpwp_res'] + strs['wpwp_sfs']

dth = np.gradient(th, z); m = z <= ztop
zi  = z[np.argmax(dth[m])]
wstar = (9.81/th[np.argmin(abs(z-zi/2))] * Q0 * zi)**(1/3)
ustar = float(np.mean(ust['ustar'])) if ust else np.nan
ml = (z > 200) & (z < 0.8*zi)

fig, ax = plt.subplots(1, 4, figsize=(17, 6))

a = ax[0]
a.plot(th[m], z[m], 'k-', lw=1.8)
a.axhline(zi, color='r', ls=':', lw=0.9)
a.set_xlim(300.8, 306); a.set_xlabel(r'$\langle\theta\rangle$  [K]'); a.set_ylabel('z [m]')
a.set_title(f'potential temperature\nmixed layer {th[ml].mean():.2f} K,  $z_i$ = {zi:.0f} m')
ins = a.inset_axes([0.52, 0.08, 0.45, 0.42])
b = (z > zi-180) & (z < zi+180)
ins.plot(th[b], z[b], 'k.-', ms=3, lw=1.2); ins.axhline(zi, color='r', ls=':', lw=0.8)
ins.set_title('inversion', fontsize=8); ins.tick_params(labelsize=7); ins.grid(alpha=0.3)

a = ax[1]
a.plot(wth[m]/Q0, z[m], 'k-', lw=1.8, label='total')
a.plot(wth_r[m]/Q0, z[m], '--', lw=1.2, label='resolved')
a.plot(wth_s[m]/Q0, z[m], ':', lw=1.4, label='SFS')
a.axvline(0, color='k', lw=0.5); a.axvline(1, color='b', lw=0.6, ls='--')
a.axhline(zi, color='r', ls=':', lw=0.9)
# The entrainment value is the local minimum at the inversion. Above it this
# run has a spurious negative flux that keeps growing with height (gravity
# waves in a layer where f_Ri = 0 leaves no SGS viscosity), so a global
# minimum would report that instead. Search only 0.85-1.15 z_i.
ent = m & (z > 0.85*zi) & (z < 1.15*zi)
k   = np.argmin(np.where(ent, wth, np.inf))
a.set_xlabel(r"$\langle w'\theta'\rangle\,/\,Q_0$"); a.legend(fontsize=8)
a.set_title(f"heat flux\nentrainment {wth[k]/Q0:+.2f} $Q_0$ at {z[k]:.0f} m")

a = ax[2]
kk = np.argmax(ww[m])
a.plot(ww[m]/wstar**2, z[m]/zi, 'k-', lw=1.8)
a.axhline(1.0, color='r', ls=':', lw=0.9)
a.set_xlabel(r"$\langle w'w'\rangle\,/\,w_*^2$"); a.set_ylabel(r'$z/z_i$')
a.set_ylim(0, ztop/zi)
a.set_title(f"vertical velocity variance\npeak {ww[m][kk]/wstar**2:.2f} $w_*^2$ at {z[m][kk]/zi:.2f} $z_i$,  $w_*$ = {wstar:.2f} m/s")

a = ax[3]
a.plot(u[m], z[m], '-', lw=1.8, label=r'$\langle u\rangle$')
a.plot(v[m], z[m], '-', lw=1.8, label=r'$\langle v\rangle$')
a.axhline(zi, color='r', ls=':', lw=0.9); a.axvline(10, color='b', lw=0.6, ls='--')
a.set_xlabel('m/s'); a.legend(fontsize=9)
a.set_title(f'mean wind\nmixed layer {u[ml].mean():.2f} / {v[ml].mean():.2f} m/s' +
            (f',  $u_*$ = {ustar:.3f}' if ust else ''))

for a in (ax[0], ax[1], ax[3]): a.set_ylim(0, ztop); a.grid(alpha=0.3)
ax[2].grid(alpha=0.3)
fig.suptitle('LESICP2-64x64x60   flat, u10   t = 9000-10800 s', fontsize=13)
fig.tight_layout(); fig.savefig(out, dpi=120); print('wrote', out)
