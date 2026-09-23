#!/usr/bin/env python3
"""Twelve-panel CBL profile figure from Jexpresso's built-in time averages.

Usage: plot_les_tavg_dat.py <les_statistics_tavg.dat> <les_stress_tavg.dat>
                            [les_ustar_tavg.dat] [out.png]

Reads the .dat files written by les_finalize! (header line names the columns),
plots resolved + subfilter together where both exist, and prints the numbers
that matter for the LESICP intercomparison. Env: ZTOP (default 2000),
WTHETA_S (surface kinematic heat flux, default 0.12).
"""
import sys, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

def load(path):
    names = open(path).readline().lstrip('#').split()[2:]   # drop time_end, n_samples
    d = np.loadtxt(path)
    d = d[:, 1:] if d.shape[1] == len(names) + 1 else d     # drop the time column
    return {n: d[:, i] for i, n in enumerate(names)}

stat = load(sys.argv[1]); strs = load(sys.argv[2])
ust  = load(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3].endswith('.dat') else None
out  = sys.argv[-1] if sys.argv[-1].endswith('.png') else 'les_tavg_profiles.png'
ztop = float(os.environ.get('ZTOP', '2000')); Q0 = float(os.environ.get('WTHETA_S', '0.12'))

z = stat['z']; m = z <= ztop
def g(d, k): return d[k] if k in d else np.zeros_like(z)

uw  = g(strs,'upwp_res') + g(strs,'upwp_sfs')
vw  = g(strs,'vpwp_res') + g(strs,'vpwp_sfs')
wth = g(strs,'wptp_res') + g(strs,'wptp_sfs')
tke = 0.5*(g(strs,'upup_res') + g(strs,'vpvp_res') + g(strs,'wpwp_res'))
th  = stat['t_mean']

# z_i from the steepest theta gradient; w* from the surface flux
dth = np.gradient(th, z); zi = z[np.argmax(dth[m])]
wstar = (9.81/th[np.argmin(abs(z-zi/2))]*Q0*zi)**(1/3)
ustar = float(np.mean(ust['ustar'])) if ust else np.sqrt(max(-uw[1], 1e-9))

fig, ax = plt.subplots(3, 4, figsize=(19, 12), sharey=True)
A = ax.ravel()
def P(i, series, xlabel, title, norm=1.0):
    for lbl, v, st in series:
        A[i].plot(v[m]/norm, z[m], st, lw=1.6, label=lbl)
    A[i].set_xlabel(xlabel); A[i].set_title(title, fontsize=11)
    A[i].axhline(zi, color='r', ls=':', lw=0.8); A[i].grid(alpha=0.3)
    if len(series) > 1: A[i].legend(fontsize=8)

P(0, [('u','-'),('v','-')] and [('<u>', stat['u_mean'],'-'), ('<v>', stat['v_mean'],'-')], 'm/s', 'mean wind')
P(1, [('<θ>', th, '-')], 'K', 'potential temperature')
P(2, [('total', uw,'k-'), ('resolved', g(strs,'upwp_res'),'--'), ('SFS', g(strs,'upwp_sfs'),':')], "m²/s²", "u'w'")
P(3, [('total', wth,'k-'), ('resolved', g(strs,'wptp_res'),'--'), ('SFS', g(strs,'wptp_sfs'),':')], 'K m/s', "w'θ'")
P(4, [("u'u'", g(strs,'upup_res'),'-'), ("v'v'", g(strs,'vpvp_res'),'-'), ("w'w'", g(strs,'wpwp_res'),'-')], 'm²/s²', 'resolved variances')
P(5, [('TKE', tke,'-')], 'm²/s²', 'resolved TKE')
P(6, [("θ'θ'", g(strs,'tptp_res'),'-')], 'K²', "θ variance")
P(7, [('ε', g(strs,'eps'),'-'), ('ε_θ', g(strs,'eps_t'),'-')], 'm²/s³  |  K²/s', 'SGS dissipation')
P(8, [("w'w'w'", g(strs,'wpwpwp'),'-')], 'm³/s³', "third moment w'w'w'  (plume skewness)")
P(9, [("w'p'", g(strs,'wppp'),'-'), ("u'p'", g(strs,'uppp'),'-')], 'm³/s³', 'pressure transport')
P(10,[('total', wth/Q0,'k-')], "w'θ' / Q₀", 'normalised heat flux', 1.0)
A[10].axvline(0, color='k', lw=0.5); A[10].axvline(1, color='b', lw=0.5, ls='--')
if ust:
    A[11].plot(ust['ustar'], ust['x']/1000.0, '-'); A[11].set_xlabel('u* [m/s]'); A[11].set_ylabel('x [km]')
    A[11].set_title(f"u*(x)   mean {ustar:.4f}, spread {np.ptp(ust['ustar']):.1e}", fontsize=11)
    A[11].grid(alpha=0.3); A[11].set_ylim(0, ust['x'].max()/1000)
else:
    A[11].axis('off')
for a in A[:11]: a.set_ylim(0, ztop)
for r in ax[:, 0]: r.set_ylabel('z [m]')
fig.suptitle(f'LESICP2-64x64x60   t-average, z_i = {zi:.0f} m, u* = {ustar:.3f} m/s, w* = {wstar:.2f} m/s', fontsize=13)
fig.tight_layout(); fig.savefig(out, dpi=110); print('wrote', out)

# ── numbers for the intercomparison ─────────────────────────────────────────
ml = (z > 200) & (z < 0.8*zi); k = lambda zz: int(np.argmin(abs(z - zz)))
core = dth > 0.25*dth[m].max()
print(f"\nz_i = {zi:.0f} m   inversion {z[core].min():.0f}-{z[core].max():.0f} m "
      f"({z[core].max()-z[core].min():.0f} m), peak {1000*dth[m].max():.1f} K/km")
print(f"mixed layer <theta> {th[ml].mean():.3f} K   <u> {stat['u_mean'][ml].mean():.2f}  <v> {stat['v_mean'][ml].mean():.2f} m/s")
print(f"u* (from u'w' at z1) {np.sqrt(-uw[1]):.4f}" + (f"   u*(x) mean {ustar:.4f}" if ust else ""))
print(f"w'theta' at wall {wth[0]:+.4f}, z1 {wth[1]:+.4f}, 100 m {wth[k(100)]:+.4f}  (Q0 = {Q0})")
print(f"   SFS share at z1 {g(strs,'wptp_sfs')[1]/wth[1]*100:.0f}%, at 100 m {g(strs,'wptp_sfs')[k(100)]/wth[k(100)]*100:.0f}%")
print(f"entrainment min w'theta' {wth[m].min():+.4f} at {z[m][np.argmin(wth[m])]:.0f} m  ({wth[m].min()/Q0:+.2f} Q0)")
print(f"<w'w'> peak {g(strs,'wpwp_res')[m].max():.3f} at {z[m][np.argmax(g(strs,'wpwp_res')[m])]:.0f} m "
      f"= {g(strs,'wpwp_res')[m].max()/wstar**2:.2f} w*^2 at {z[m][np.argmax(g(strs,'wpwp_res')[m])]/zi:.2f} z_i")
print(f"u'w' at wall {uw[0]:+.4f}, z1 {uw[1]:+.4f}   (-u*^2 = {-ustar**2:+.4f})")
