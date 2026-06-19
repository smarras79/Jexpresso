#!/usr/bin/env python3
"""
Plot DHIT kinetic-energy spectra written by Jexpresso's compute_dhit_spectra,
in the style of Flad & Gassner (JCP 350, 2017).

The Jexpresso run writes, at every :diagnostics_at_times, two files per output
index <iout> into :output_dir:

    dhit_spectrum_3D_<iout>.dat   columns: k   E(k)            (3D shell spectrum)
    dhit_spectrum_1D_<iout>.dat   columns: k   E1_avg E1x E1y E1z

Usage:
    python3 plot_dhit_spectra.py [output_dir] [--mode 3d|1d] [--save fig.png]

With no arguments it looks in ./output-dhit/. Multiple snapshots are overlaid
(faded -> solid with time). A k^(-5/3) reference slope is drawn.
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
    ap.add_argument("--mode", choices=["3d", "1d"], default="3d")
    ap.add_argument("--save", default=None)
    args = ap.parse_args()

    tag = "3D" if args.mode == "3d" else "1D"
    files = sorted(glob.glob(os.path.join(args.outdir, f"dhit_spectrum_{tag}_*.dat")),
                   key=_index)
    if not files:
        raise SystemExit(f"No dhit_spectrum_{tag}_*.dat files found in {args.outdir}")

    fig, ax = plt.subplots(figsize=(7, 5))
    n = len(files)
    for i, f in enumerate(files):
        data = np.loadtxt(f, comments="#")
        k = data[:, 0]
        E = data[:, 1]                      # E(k) for 3D, E1_avg for 1D
        good = E > 0
        alpha = 0.35 + 0.65 * (i / max(n - 1, 1))
        ax.loglog(k[good], E[good], "-", lw=1.8, alpha=alpha,
                  label=os.path.basename(f))

    # k^(-5/3) inertial-range reference
    kk = np.array([k[good][1], k[good][-1]])
    if kk[0] > 0:
        E0 = E[good][1]
        ax.loglog(kk, E0 * (kk / kk[0]) ** (-5.0 / 3.0), "--", color="gray",
                  lw=1.2, label=r"$k^{-5/3}$")

    ax.set_xlabel(r"$k$")
    ax.set_ylabel(r"$E(k)$")
    ax.set_title(f"DHIT kinetic-energy spectrum ({tag})")
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
