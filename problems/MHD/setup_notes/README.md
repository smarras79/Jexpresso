# DynSGS and positivity setup of the `problems/MHD` cases

One note per case and a cross-case summary, describing the DynSGS
(`DSGS_MHD`) stabilization and the positivity mechanisms **as coded on
branch `sm/MHD` (commit `dcc3653`)**: deck values, the defaults
`src/io/mod_inputs.jl` fills in, and the kernel paths they reach
(`src/kernel/physics/SGS.jl`, `src/kernel/operators/rhs.jl`,
`src/kernel/positivity/`), with file:line references.

| file | case |
|---|---|
| `summary.pdf` | all nine cases side by side, common items, findings |
| `brioWu1d.pdf` | 1D Brio–Wu shock tube |
| `smoothVortex.pdf` | smooth MHD vortex (accuracy test) |
| `orszagTangBormanis2024.pdf` | Orszag–Tang vortex (physical-form DynSGS) |
| `kelvinHelmholtzChan2022.pdf` | magnetized KHI (Smagorinsky, no DynSGS) |
| `fluxEmergenceSon2025.pdf` | flux emergence, with the case-local floor |
| `fluxEmergenceSon2025DSGS.pdf` | flux emergence, limiter-free (reference-weighted DynSGS) |
| `rotorDaoNazarov2022.pdf` | MHD rotor (conservative element limiter) |
| `blastBalsaraSpicer1999.pdf` | MHD blast (limiter, fast floors, thermal residual) |
| `astroJetWuShu2018.pdf` | Mach-800 magnetized jet (limiter, residual sensor) |
| `dsgsCutoff.pdf` | two pages, in the paper's notation: the legacy and residual sensors, and DynSGS without and with the smoothness cutoff 1e-3 on the Alfvén wave (`dsgsCutoff_fig.py` draws its figure) |
| `dsgsSetups.pdf` | the DynSGS setup of every test in the paper (smooth vortices, Brio–Wu, Orszag–Tang, rotor, blast, jet, flux emergence) in its notation, as the decks stand at `66ebb9c`, with where the paper's text and the decks differ and the dated default changes |

**Underlines.** Solid: identical in every case to which the item applies
(all nine for general settings, the eight `DSGS_MHD` cases for DynSGS, the
three `"conservative"`-limiter cases for the limiter). Dashed: shared by a
majority of those cases, not all.

Shared sources: `preamble.tex` (macros, legend), `codemap.tex` (the
file:line appendix of every note), `cpos_algorithm.tex` (the conservative
limiter, included by the rotor, blast and jet notes).

Build (TeX Live, `pdflatex`):

```sh
cd problems/MHD/setup_notes
for f in *.tex; do
  case $f in preamble.tex|codemap.tex|cpos_algorithm.tex) continue;; esac
  latexmk -pdf -interaction=nonstopmode "$f"
done
latexmk -c
```
