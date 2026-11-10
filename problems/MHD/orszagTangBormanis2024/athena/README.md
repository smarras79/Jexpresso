# Athena pressure cuts of the Orszag–Tang vortex at t = 0.5

`orszag_tang_athena_0_3125.csv` and `orszag_tang_athena_0_4277.csv` give the pressure
along y = 0.3125 and y = 0.4277 at t = 0.5. They are the red points of Fig. 15 in

> A. Babbar, H. Ranocha, *Compact Runge-Kutta flux reconstruction methods for
> non-conservative hyperbolic equations*, J. Comput. Phys. 563 (2026) 115060.

Columns: `x p`, one point per line, in the order they were digitized (not sorted in x).

## Where they come from

The files are copied unchanged from the paper's reproducibility repository,
<https://github.com/Arpit-Babbar/paper_crk_nonconservative> (commit `19a6900`),
`run/tang/ref_data_tang/`, which is MIT-licensed (see `LICENSE` here). Next to them that
repository keeps WebPlotDigitizer projects, so the points were digitized from a plot rather
than written by a solver.

According to the paper, they are the Athena solution (Stone et al., ApJS 178 (2008) 137)
on 1024 × 1024 points, taken from Rueda-Ramírez et al., J. Comput. Phys. 444 (2021) 110580.
Digitization adds an error of about the size of a marker.

## Normalization

The points are in Athena's normalization: ρ = 25/(36π), p = 5/(12π), B₀ = 1/√(4π),
v = (−sin 2πy, sin 2πx), γ = 5/3.

- The paper's run behind Figs. 14–15 (`run/tang/run_mhd_tang.jl`) uses this normalization,
  even though eq. 82 of the paper prints ρ = 1, p = 1/γ, B/γ.
- `orszagTangBormanis2024` uses it too (`initialize.jl`: ρ = γ²/4π, p = γ/4π, B₀ = 1/√4π).
  Its pressure can therefore be compared with these points without rescaling.

To compare, use `python3 tools/plot_orszagtang_vtu.py --cuts`.
