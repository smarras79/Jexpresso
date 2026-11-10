#!/usr/bin/env python3
"""The paper's figures of the Dai & Woodward (1998) shock-cloud cases, all from plot_mhd_matrix.py with fixed options.

    python3 tools/plot_shock_cloud_figures.py --letters output/MHD/shockCloudDaiWoodward1998/<run> \\
                                              --cloud   output/MHD/shockCloudMach15DaiWoodward1998/<run>
    python3 tools/plot_shock_cloud_figures.py --letters <run> --out-dir figures --format png --dry-run

Either case may be left out. Each figure is one plot_mhd_matrix.py call, printed before it runs (--dry-run prints
only), so any figure can be redone by hand. Written to --out-dir (default figures/shockCloud), as --format (pdf):

  letters (problems/MHD/shockCloudDaiWoodward1998, DW Fig. 20), rows = --letters-times:
    letters_fields      log10 rho (gray, DW's Fig. 20) | log10 rho with the field lines over it | log10 p |
                        Mach_fs (the fast Mach number in the shock's frame: 12.0 ahead of the shock, DW's "Mach 12") | nu_DSGS
    letters_paper       log10 rho (gray) | the field lines alone | Mach_fs (0-12) | nu_DSGS
    letters_fieldlines  log10 rho with the field lines over it | the field lines alone (DW Fig. 19 style)
    letters_shockmach   Mach | Mach_s | Mach_fs: the Mach number in the lab frame and in the shock's frame
                        (ahead of the shock: 14.4 sonic, 12.0 fast = DW's "Mach 12")
  cloud (problems/MHD/shockCloudMach15DaiWoodward1998, DW Figs. 23-25), rows = --cloud-times:
    cloud_fields        log10 rho (gray, light where dense, DW's Fig. 23) | log10 rho with the field lines over it |
                        log10 p | Mach_fs (15.0 ahead of the shock, DW's "Mach 15") | nu_DSGS
    cloud_paper         log10 rho (gray) | the field lines alone | Mach_fs (0-15) | nu_DSGS
    cloud_final         rho | Bx | By | nu_DSGS at the last of --cloud-times (DW Figs. 23-25)
    cloud_fieldlines    the field lines alone at the last time, zoomed on the cloud (DW Fig. 19)
    cloud_shockmach     Mach_s | Mach_fs (ahead of the shock: 15.8 sonic, 15.0 fast = DW's "Mach 15")

The times must be output times of the run (the decks write 0.0025 / 0.02 apart). Needs plot_mhd_matrix.py's
dependencies (numpy, matplotlib).
"""
import argparse
import os
import shlex
import subprocess
import sys

MATRIX = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plot_mhd_matrix.py")


def figures(case, run, times):
    """(name, plot_mhd_matrix.py arguments) of every figure of one case."""
    T = ["--times", *times]
    if case == "letters":
        return [
            ("letters_fields", [run, *T, "--vars", "rho", "log10rho", "p", "Mach_fs", "mu_dsgs", "--log", "rho", "p",
                                "--cmaps", "rho", "gray_r", "--fieldlines-on", "log10rho", "--nlevels", "60",
                                "--cmaps", "Mach_fs", "YellowGreyBlue", "--clim", "Mach_fs", "0", "12"]),
            ("letters_paper", [run, *T, "--vars", "rho", "Blines", "Mach_fs", "mu_dsgs", "--log", "rho",
                               "--cmaps", "rho", "gray_r", "--nlevels", "60",
                               "--cmaps", "Mach_fs", "YellowGreyBlue", "--clim", "Mach_fs", "0", "12"]),
            ("letters_fieldlines", [run, *T, "--vars", "rho", "Blines", "--log", "rho", "--cmaps", "rho", "gray_r",
                                    "--fieldlines-on", "rho", "--nlevels", "80"]),
            ("letters_shockmach", [run, *T, "--vars", "Mach", "Mach_s", "Mach_fs",
                                   "--cmaps", "Mach", "YellowGreyBlue", "--clim", "Mach", "0", "3",
                                   "--clim", "Mach_s", "0", "16", "--clim", "Mach_fs", "0", "16"]),
        ]
    last = ["--times", times[-1]]
    return [
        ("cloud_fields", [run, *T, "--vars", "rho", "log10rho", "p", "Mach_fs", "mu_dsgs", "--log", "rho", "p",
                          "--cmaps", "rho", "gray", "--fieldlines-on", "log10rho", "--nlevels", "60",
                          "--cmaps", "Mach_fs", "YellowGreyBlue", "--clim", "Mach_fs", "0", "15"]),
        ("cloud_paper", [run, *T, "--vars", "rho", "Blines", "Mach_fs", "mu_dsgs", "--log", "rho",
                         "--cmaps", "rho", "gray", "--nlevels", "60",
                         "--cmaps", "Mach_fs", "YellowGreyBlue", "--clim", "Mach_fs", "0", "15"]),
        ("cloud_final", [run, *last, "--vars", "rho", "Bx", "By", "mu_dsgs", "--cmaps", "rho", "gray"]),
        ("cloud_fieldlines", [run, *last, "--vars", "Blines", "--nlevels", "100",
                              "--xlim", "0.1", "0.9", "--ylim", "0.1", "0.9"]),
        ("cloud_shockmach", [run, *T, "--vars", "rho", "Mach_s", "Mach_fs", "--log", "rho", "--cmaps", "rho", "gray",
                             "--clim", "Mach_s", "0", "18", "--clim", "Mach_fs", "0", "18"]),
    ]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--letters", help="output directory of a shockCloudDaiWoodward1998 run")
    ap.add_argument("--cloud", help="output directory of a shockCloudMach15DaiWoodward1998 run")
    ap.add_argument("--letters-times", nargs="+", default=["0.025", "0.05", "0.08", "0.105"],
                    help="rows of the letters figures (default: DW's Fig. 20 times as the deck writes them)")
    ap.add_argument("--cloud-times", nargs="+", default=["0.04", "0.08", "0.12", "0.16"],
                    help="rows of the cloud figures; the last one is also the time of cloud_final and cloud_fieldlines")
    ap.add_argument("--out-dir", default=os.path.join("figures", "shockCloud"))
    ap.add_argument("--format", default="pdf")
    ap.add_argument("--dpi", default="300")
    ap.add_argument("--dry-run", action="store_true", help="print the commands only")
    a = ap.parse_args()
    if not (a.letters or a.cloud):
        ap.error("give --letters and/or --cloud")

    jobs = ([("letters", a.letters, a.letters_times)] if a.letters else []) + \
           ([("cloud", a.cloud, a.cloud_times)] if a.cloud else [])
    if not a.dry_run:
        os.makedirs(a.out_dir, exist_ok=True)
    failed = []
    for case, run, times in jobs:
        for name, args in figures(case, run, times):
            out = os.path.join(a.out_dir, f"{name}.{a.format}")
            cmd = [sys.executable, MATRIX, *args, "--format", a.format, "--dpi", a.dpi, "--out", out]
            print("$ " + " ".join(shlex.quote(c) for c in ["python3", os.path.relpath(MATRIX), *cmd[2:]]), flush=True)
            if not a.dry_run and subprocess.run(cmd).returncode != 0:
                failed.append(name)
    if failed:
        sys.exit("failed: " + ", ".join(failed))


if __name__ == "__main__":
    main()
