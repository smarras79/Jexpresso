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

- `dhit_spectrum_3D_<iout>.dat` — columns: `k  E(k)`
- `dhit_spectrum_1D_<iout>.dat` — columns: `k  E1_avg  E1_x  E1_y  E1_z`

Each header line records the time, `N_uni`, Nyquist, `dk = 2π/L`, and total
resolved KE.

## Plotting

```bash
python3 problems/CompEuler/dhit3D/plot_dhit_spectra.py ./output-dhit/ --mode 3d --save Ek.png
```

## Enabling / disabling

Controlled by `:dhit_spectra => true` in `user_inputs.jl`. Set to `false` to
skip the spectra computation. The feature is a no-op for every other problem.
