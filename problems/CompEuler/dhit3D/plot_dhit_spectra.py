#!/usr/bin/env python3
"""
Plot DHIT kinetic-energy spectra written by Jexpresso's compute_dhit_spectra,
in the style of Flad & Gassner (JCP 350, 2017).

The Jexpresso run writes, at every :diagnostics_at_times, into :output_dir:

    dhit_spectrum_3D_<iout>.dat   columns: k  E(k)  D(k)=2*nu*k^2*E(k)
    dhit_spectrum_1D_<iout>.dat   columns: k  E1_avg E1x E1y E1z
    dhit_integrals.dat            time series: t Ekin(k<=kcut) eps(k<=kcut) Ekin_total eps_total

Usage:
    python3 plot_dhit_spectra.py [output_dir] [--mode 3d|1d|diss] [--save fig.png]

  --mode 3d    3D shell spectrum E(k)            (default)
  --mode 1d    direction-averaged 1D spectrum
  --mode diss  dissipation-rate spectrum D(k)=2 nu k^2 E(k)

With no arguments it looks in ./output-dhit/. Multiple snapshots are overlaid
(faded -> solid with time). A k^(-5/3) reference slope is drawn for energy
spectra.
"""
import argparse, glob, os, re
import numpy as np
import matplotlib.pyplot as plt


def _index(path):
    m = re.search(r"_(\d+)\.dat$", os.path.basename(path))
    return int(m.group(1)) if m else -1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir", nargs="?", default="./output-dhit/")
    ap.add_argument("--mode", choices=["3d", "1d", "diss"], default="3d")
    ap.add_argument("--save", default=None)
    args = ap.parse_args()

    tag = "1D" if args.mode == "1d" else "3D"
    col = 2 if args.mode == "diss" else 1     # 3D file: 1->E(k), 2->D(k)
    files = sorted(glob.glob(os.path.join(args.outdir, f"dhit_spectrum_{tag}_*.dat")),
                   key=_index)
    if not files:
        raise SystemExit(f"No dhit_spectrum_{tag}_*.dat files found in {args.outdir}")

    fig, ax = plt.subplots(figsize=(7, 5))
    n = len(files)
    for i, f in enumerate(files):
        data = np.loadtxt(f, comments="#")
        k = data[:, 0]
        E = data[:, col]                    # E(k), E1_avg, or D(k)
        good = E > 0
        alpha = 0.35 + 0.65 * (i / max(n - 1, 1))
        ax.loglog(k[good], E[good], "-", lw=1.8, alpha=alpha,
                  label=os.path.basename(f))

    # k^(-5/3) inertial-range reference (energy spectra only)
    if args.mode != "diss":
        kk = np.array([k[good][1], k[good][-1]])
        if kk[0] > 0:
            E0 = E[good][1]
            ax.loglog(kk, E0 * (kk / kk[0]) ** (-5.0 / 3.0), "--", color="gray",
                      lw=1.2, label=r"$k^{-5/3}$")

    ax.set_xlabel(r"$k$")
    ax.set_ylabel(r"$D(k)=2\nu k^2 E(k)$" if args.mode == "diss" else r"$E(k)$")
    ax.set_title(f"DHIT {'dissipation' if args.mode=='diss' else 'kinetic-energy'} spectrum ({tag})")
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(True, which="both", ls=":", alpha=0.4)
    fig.tight_layout()
    if args.save:
        fig.savefig(args.save, dpi=150)
        print("saved", args.save)
    else:
        plt.show()


if __name__ == "__main__":
    main()
