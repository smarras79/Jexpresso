# DHIT kinetic-energy spectra (Flad & Gassner 2017)

The DHIT test case (`problems/CompEuler/dhit3D`) computes 1D and 3D kinetic-energy
spectra of the DGSEM velocity field, following the methodology of

> D. Flad, G. Gassner, *On the use of kinetic energy preserving DG-schemes for
> large eddy simulation*, J. Comput. Phys. 350 (2017) 782–795.

## What is computed

At every `:diagnostics_at_times`, `compute_dhit_spectra` (in
`src/io/dhit_spectra.jl`):

1. **Super-samples** the piecewise-polynomial DG solution from the LGL
   collocation nodes onto a **global equidistant Cartesian grid** of
   `N_uni = nel_per_dir · (N+1)` points per direction — i.e. the paper's
   "DOF per direction". The interpolation is an element-local tensor-product
   Lagrange interpolation; consecutive elements tile a periodic uniform grid
   with no duplicated interface points.
2. Takes a **3D FFT** of the three velocity components.
3. Forms the **3D KE spectrum** by spherical-shell summation,
   `E(k) = Σ_{k-½ ≤ |κ| < k+½} ½(|û|²+|v̂|²+|ŵ|²)`, for integer shells up to the
   Nyquist wavenumber `N_uni/2`.
4. Forms the **1D spectrum** by summing energy over the two perpendicular
   wavenumber planes, averaged over the three homogeneous directions.

Both spectra satisfy `Σ_k E(k) = ½⟨|u|²⟩` (resolved kinetic energy).

5. The **dissipation-rate spectrum** `D(k) = 2ν k² E(k)` and the Fourier-space
   band integrals of Fig. 9,
   `E_kin = Σ_{k≤k_cut} E(k)` and `ε = 2ν Σ_{k≤k_cut} k² E(k)`,
   evaluated with a fixed molecular kinematic viscosity `ν` (= `:dhit_nu`,
   default `mu_molecular`, with `ρ_ref = 1`) so DNS and LES are directly
   comparable, and a cutoff `k_cut = :dhit_kcut` (default 16).

## Matching the paper's resolutions

`DOF/dir = nel · (N+1)`, `k_Ny = DOF/2`:

| Setup            | mesh (`JEXP_MESH`)       | N (`:nop`) | DOF/dir | k_Ny |
|------------------|--------------------------|------------|---------|------|
| coarse LES (32³-equiv) | `cube_periodic_6.msh`  (6³) | 7 | 48  | 24 |
| well-resolved LES      | 18³ periodic cube           | 7 | 144 | 72 |

The reference DNS in the paper is a 512³ pseudo-spectral run; energy and
dissipation are integrated in Fourier space up to `k = 16`.

## Output

Per output index `<iout>`, written to `:output_dir`:

- `dhit_spectrum_3D_<iout>.dat` — columns: `k  E(k)  D(k)=2νk²E(k)`
- `dhit_spectrum_1D_<iout>.dat` — columns: `k  E1_avg  E1_x  E1_y  E1_z`

and a single appended time series:

- `dhit_integrals.dat` — columns: `t  Ekin(k≤kcut)  eps(k≤kcut)  Ekin_total  eps_total`

Each spectrum header records the time, `N_uni`, Nyquist, `dk = 2π/L`, `ν`, and
total resolved KE.

## Plotting

```bash
# energy spectrum E(k)
python3 problems/CompEuler/dhit3D/plot_dhit_spectra.py ./output-dhit/ --mode 3d   --save Ek.png
# dissipation-rate spectrum D(k)=2 nu k^2 E(k)
python3 problems/CompEuler/dhit3D/plot_dhit_spectra.py ./output-dhit/ --mode diss --save Dk.png
```

The energy/dissipation decay over time (Fig. 9) is the `dhit_integrals.dat`
time series: plot columns 2/3 (band-limited) or 4/5 (full) against column 1.

## Enabling / disabling

Controlled by `:dhit_spectra => true` in `user_inputs.jl`. Set to `false` to
skip the spectra computation. The feature is a no-op for every other problem.
