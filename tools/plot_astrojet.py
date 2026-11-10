#!/usr/bin/env python3
"""astroJetWuShu2018 figure: rows t = 0.0011, 0.0016, 0.002; columns log10 rho, log10 p, M_f, nu_DSGS and the
numerical schlieren; a half-domain run is mirrored about the symmetry axis x = 0 (plot_mhd_matrix.py does it for any
AJ_half_* mesh; --no-mirror keeps the half). Any plot_mhd_matrix.py option given here overrides the default of the same name.
    python3 tools/plot_astrojet.py output/MHD/astroJetWuShu2018/output-<date>
    python3 tools/plot_astrojet.py <run> --times 0.002 --format pdf --out jet.pdf
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_mhd_matrix

DEFAULTS = ["--times", "0.0011", "0.0016", "0.002",
            "--vars", "rho", "p", "Mach_f", "mu_dsgs", "schlieren", "--log", "rho", "p", "--dpi", "200"]

if __name__ == "__main__":
    sys.argv = [sys.argv[0]] + DEFAULTS + sys.argv[1:]
    plot_mhd_matrix.main()
