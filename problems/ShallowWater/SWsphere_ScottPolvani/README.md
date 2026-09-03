# SWsphere_ScottPolvani — forced-dissipative shallow-water turbulence on the sphere (giant planets)

The setup of

> Scott, R. K. & Polvani, L. M. (2007), *Forced-dissipative shallow-water
> turbulence on the sphere and the atmospheric circulation of the giant
> planets*, J. Atmos. Sci. **64**, 3158–3176,

on the **same grid, equations and discretisation as `SWsphere`** (the Galewsky
jet): the cubed sphere shipped next to this file, the conservative Cartesian
shallow-water equations of Marras, Kopera & Giraldo (2015), continuous-Galerkin
spectral elements, SSP-RK3 with the Lagrange projection after every stage, the
modal filter and the `ν∇ₛ²` artificial diffusion. Everything the `SWsphere`
README says about the grid, the metrics, the node numbering, the output files
and running in parallel holds here unchanged.

What is different is the *problem*:

| | `SWsphere` | `SWsphere_ScottPolvani` |
|---|---|---|
| planet | Earth | the three of the paper's Fig. 14 — Jupiter, Saturn, Neptune (Uranus/Neptune) — one symbol in the deck or `SP_PLANET` in the environment |
| initial state | balanced mid-latitude jet + bump | **fluid at rest**, uniform depth `H` set by the deformation radius |
| forcing | none (an initial-value problem) | **random, isotropic, small-scale vorticity forcing** injecting energy at a prescribed rate `ε₀` |
| large-scale dissipation | none | **Rayleigh friction** on the momentum and/or **radiative relaxation** of the height |
| what you look at | relative vorticity at day 6 | **zonal-mean zonal velocity vs latitude** as the jets emerge; Ro; the vorticity |
| duration | 6–10 days | hundreds to tens of thousands of rotations |

## Run it

```bash
julia --project=.
julia> using Jexpresso
julia> Jexpresso.run_case("ShallowWater", "SWsphere_ScottPolvani")
```

Pick the planet at the top of `user_inputs.jl` (`planet = …`), and the length
of the run in rotations (`nrot`). Everything else is derived. Both can also be
set from the environment without editing the file, which is how the three
panels of the paper's Fig. 14 are run in a row:

```bash
for p in jupiter saturn neptune; do
  SP_PLANET=$p SP_NROT=2000 julia --project=. -e 'using Jexpresso; Jexpresso.run_case("ShallowWater", "SWsphere_ScottPolvani")'
  mv output/ShallowWater/SWsphere_ScottPolvani/output output/ShallowWater/SWsphere_ScottPolvani/output_$p
done
```

(`:loverwrite_output => true` in the deck, so move each result aside before
the next planet.) In parallel nothing changes; see the `SWsphere` README.

## The model (paper, section 3)

The paper writes the shallow-water equations in vorticity `ζ`, divergence `δ`
and height `h`, Eq. (9), forces the **vorticity** equation with `F` and damps
each field with `D_ξ`:

```
ζ_t + ∇·(u ζ_a) = F + D_ζ
δ_t − k·∇×(u ζ_a) = −∇²(E + gh) + D_δ
h_t + ∇·(u h) = D_h
```

Jexpresso integrates the same equations in the conservative Cartesian form
`q = [φ, φu, φv, φw]`, `φ = gh` (see `user_flux.jl`, `user_source.jl`), so the two
extra terms are translated (`src/kernel/operators/sphere_forcing.jl`):

* **Forcing.** A vorticity forcing `F` is the curl of a velocity forcing:
  with `∇ₛ²ψ_F = F` and `u_F = n̂ × ∇ₛψ_F`, the momentum forcing `φ u_F` forces
  the vorticity with `F` and the divergence with nothing — Eq. (9a)–(9b). `ψ_F`
  is built directly from the paper's spectral band, Eq. (11): real orthonormal
  spherical harmonics of every `(n, m)` with `|n − n_f| ≤ Δn/2`, each `(n, m)`
  pair with a fixed amplitude and a random phase (`e^{iθ}` in Eq. 12). The
  amplitude is set so the **kinetic energy injection rate is `ε₀`**, following
  Eq. (12): each `n` in the band gets `ε₀/n_band`, spread over its `2n+1` zonal
  modes. The forcing is either white in time (Eq. 12) or Markovian with
  decorrelation time `τ` (Eq. 13, the paper's `c_r`); it is drawn once per time
  step and held fixed over the step. `u_F` is obtained from `ψ_F` with the same
  spectral-element surface gradient the RHS uses, so it is exactly tangential
  and discretely non-divergent, and has no pole singularity.

  **The amplitude is renormalised every step**, as in the paper's Eq. (14).
  The fixed amplitude of Eqs. (12)–(13) assumes the velocity the forcing
  creates stays correlated with it for a time `τ`; on a rotating shell it does
  not (Coriolis turns it round in half a rotation, and at `L_f > L_D`
  geostrophic adjustment sheds most of it to gravity waves), and measured on
  Jupiter at `τ = 10` rotations the fixed amplitude injected only `0.19 ε₀`.
  So, every step, the forcing field is scaled by the gain `g` that makes the
  step's kinetic-energy budget `Δt [gA + ½Δt g²B]` equal `Δt ε₀`, with
  `A = ⟨φu·u_F⟩/⟨φ⟩` the forcing–flow correlation and `B = ⟨φ|u_F|²⟩/⟨φ⟩`. The
  quadratic (Itô) term is what makes this regular where Eq. (14) divides by
  zero at a rest start. `:forcing_normalize => false` gives the fixed
  amplitude back; either way the injection actually produced is reported.
* **Rayleigh friction** `D_ζ = −ν_r ζ`, `D_δ = −ν_r δ` ⇒ `∂(φu)/∂t += −ν_r φu`.
* **Radiative relaxation** `D_h = −ν_h h'` ⇒ `∂φ/∂t += −ν_h (φ − φ_ref)`, with
  `φ_ref = gH` the rest state stored in `q.qe`.

Both dissipations are pointwise; both may be on. The paper's third option,
**hypodiffusion** `(−Δ)⁻¹`, is a global inverse Laplacian and is not available
in this element-local RHS — the paper shows (section 4, Figs. 2a vs 4b) that
Rayleigh friction gives equilibria of the same character, and that is what the
deck uses.

The paper's **small-scale** dissipation, `∇⁸` hyperdiffusion damping the highest
wavenumber on a 0.1-rotation time scale, is not reproduced either. The shell
carries the modal filter and `ν∇ₛ²`; both are on. `ν` is kept *small* — a
Laplacian is not scale-selective, and a `ν` large enough to matter at the grid
scale damps the forced band faster than the eddies there can turn over, at
which point the energy is dissipated where it is injected and nothing cascades.
The deck expresses `ν` in units of `a²Ω` and the comment next to it has the
numbers.

## Parameters (paper, sections 3 and 7)

The paper is nondimensional — lengths in planetary radii `a`, times in rotation
periods `T = 2π/Ω` — and characterises a run by

```
Ro  = U/(2aΩ)            Rossby number of the equilibrated flow
L_D = √(gH)/(2Ω a)       polar deformation radius, in planetary radii
```

Its section-7 planetary runs — Fig. 13 (zonal-mean wind) and Fig. 14
(potential vorticity, panels a–c, and vorticity, panels d–f), values from Cho &
Polvani 1996b — and the dimensional constants the deck attaches to them:

| `planet` | Fig. 14 | `a` [m] | `Ω` [1/s] | `g` [m/s²] | `L_D/a` | `Ro` |
|---|---|---|---|---|---|---|
| `:jupiter` | (a)/(d) | 7.1492e7 | 1.7585e-4 | 24.79 | 0.025 | 0.002 |
| `:saturn` | (b)/(e) | 6.0268e7 | 1.6378e-4 | 10.44 | 0.025 | 0.013 |
| `:neptune` | (c)/(f) | 2.4764e7 | 1.0834e-4 | 11.15 | 0.1 | 0.06 |
| `:uranus` | (c)/(f), alias | 2.5559e7 | 1.0124e-4 | 8.87 | 0.1 | 0.06 |

The paper labels the third panel "Uranus/Neptune" and treats the two as one
case; `:neptune` is the option, `:uranus` an accepted alias carrying Uranus's
own constants (the dynamics see only `Ro` and `L_D`, which are the same).

`H = (2Ω L_D a)²/g` follows (≈ 15.9 km for Jupiter; `g` then drops out of the
dynamics, since the equations carry `φ = gh`, and only reappears when `h` is
written). `Ro` is an **outcome**, not an input: with Rayleigh friction `ν_l` the
paper's closure is `E ≈ ε/(2ν_l)`, `ε ≈ 0.4 ε₀` being the fraction of the input
that cascades upscale (section 4a), so the deck turns a target `Ro` into
`ε₀ = 2ν_l · ½(4πRo)² / 0.4` (for Jupiter that is `1.6e-7 a²/T³`, next to the
paper's representative `0.1e-6`). It is a scaling estimate; read the `Ro` the
diagnostics report.

Deck defaults, in the paper's units: `ν_l = 1e-4` per rotation (Fig. 2a),
`c_r = 10` rotations (section 4a's main choice), `Δn = 4`, `n_f = 24`.

### Resolution

The shipped grid (10 elements per panel edge, `nop = 5`) has 200 LGL nodes
around the equator and resolves total wavenumbers up to `N ≈ 50–60` (the paper
runs T170 with `n_f = N/4 ≈ 42`). `n_f = 24` puts 8 nodes per wavelength at the
equator and leaves the jets — energy centroid `n₀ ≈ 5–17` in the paper's
Tables 1–2, Rhines wavenumber `n_Rh ~ √(1/2Ro) ≈ 16` for Jupiter's `Ro` — room
below it, but not much. `initialize.jl` prints `n_Rh`, `α = 2Ro/L_D²` and the
nodes per wavelength, and warns when `n_f` is not well above `n_Rh`. For a
cleaner inverse cascade regenerate the grid finer,

```bash
julia tools/generate_cubed_sphere.jl 20 6.371e6 problems/ShallowWater/SWsphere_ScottPolvani/cubed_sphere.msh equiangular
```

(the radius in the file is immaterial: the deck's `:sphere_radius` rescales it),
and raise `n_f` in step, keeping `N/n_f ≳ 2–4`.

### Time

The paper integrates to `10⁴–10⁵` rotations. At the CFL step (`Δt ≈ 521 s`,
69 steps per rotation on this grid for Jupiter, set by the gravity-wave speed
`√(gH) ≈ 630 m/s`) a rotation costs 1.3 s of wall time in serial, so the deck's
default `nrot = 500` — the spin-up, over which the paper's Fig. 1 shows the
energy growing linearly — takes 11 minutes, and an equilibrated
Rayleigh-friction run (`t ~ 1/(2ν_l) ≈ 5000` rotations) about two hours. Raise
`nrot`.

### What the default deck does, measured

Jupiter, this grid, 500 rotations (34 304 steps, 643 s of wall time, serial):

| rotations | KE/mass [m²/s²] | `U_rms` [m/s] | Ro |
|---|---|---|---|
| 100 | 22.1 | 6.6 | 2.6e-4 |
| 200 | 30.6 | 7.8 | 3.1e-4 |
| 300 | 35.6 | 8.4 | 3.4e-4 |
| 400 | 37.4 | 8.7 | 3.4e-4 |
| 500 | 41.2 | 9.1 | 3.6e-4 |

Mass is conserved to `4e-10` over the run, the off-shell momentum stays at
`1e-9` of the momentum scale, and the injection is `1.000 ε₀` at every print
with a gain of `1.6–2.6`. The kinetic energy rises fast for ~200 rotations and
then slowly: with the forcing at `L_f = a/24 > L_D = 0.025 a`, the flow at the
forcing scale is mostly potential energy and the kinetic-energy cascade is
weak — the paper's own remark for the small-`L_D` cases (section 5b, after
Polvani et al. 1994) — and `2ν n_f(n_f+1)/a² · KE` at these numbers is already
`≈ ε₀`. At 500 rotations the zonal-mean zonal wind is still an alternating
pattern at roughly the forcing scale (`±3 m/s`, ~6° spacing) with a weakly
retrograde equator; the jets of the paper's Fig. 13 are a `10³–10⁴`-rotation
outcome, so raise `nrot` (and consider a smaller `nu_nd`: at `5e-7` the run
holds 26% more kinetic energy at 100 rotations and is stable there, see the
deck comment) to get there. Uranus/Neptune parameters run equally stably and
reach `Ro ≈ 4e-3` within 10 rotations.

## What to look at

Every `:ndiagnostics_prints` steps the time loop prints, in addition to the
`SWsphere` line (mass, energy, drift, `max|ζ|`), a line

```
forced: KE/mass = …  U_rms = … m/s  Ro = …  ε_inj/ε₀ = … (mean)  … (now)
```

with the kinetic energy per unit mass `⟨½φ|u|²⟩/⟨φ⟩`, the velocity scale
`U = √(2KE)`, the Rossby number `U/(2aΩ)` it implies — the number to put on the
paper's Table 1 and Fig. 13 — the **measured** energy injection rate relative
to the `ε₀` asked for (mean since the start, and for the coming step; 1.000
with the renormalisation on, which is the point of it), and the gain `g` the
renormalisation applied to the nominal amplitude. A gain drifting far from 1
says the flow has decorrelated from the forcing; a gain pinned at its cap
(10× the white-noise gain) says the forcing could not inject `ε₀` that step.

Output, in `problems/ShallowWater/SWsphere_ScottPolvani/output/`:

| file | what it is |
|---|---|
| `sphere_grid_ho.vtu` | the initial (rest) state |
| `sphere_NNNN.vtu` | one per output time: `h`, `u` (zonal), `v` (meridional), `w` (radial, ≈ 0), `vorticity` (Fig. 14 d–f), **`pv`** `= (ζ + f)/h` the potential vorticity (Fig. 14 a–c; written when the deck gives `:sphere_Omega` and `:sphere_gravity`), **`vorticity_forcing`** `= ∇ₛ²ψ_F` — the `F` the paper forces with, so its scale and strength can be seen next to the flow |
| `zonal_mean_NNNN.dat` | **the paper's diagnostic**: `lat, ū, v̄, φ̄, nodes` in 2° bands (`:zonal_mean_nbins`), mass-weighted, written with every VTK file. Stack them to get the paper's Hovmöller plots of `ū(φ, t)` (Figs. 2, 4, 6, 7) or plot the last one against Fig. 13. |

In parallel the `.vtu` become `.pvtu`, as for `SWsphere`; the `.dat` is written
by rank 0 from globally reduced sums.

## The switches

`src/kernel/operators/sphere_forcing.jl` reads, all in SI units:

| key | meaning |
|---|---|
| `:lsphere_forcing` | master switch (`false` ⇒ the `SWsphere` equations, untouched) |
| `:forcing_epsilon` | `ε₀` [m²/s³] |
| `:forcing_nf`, `:forcing_dn` | the band `|n − n_f| ≤ Δn/2` |
| `:forcing_tau` | decorrelation time [s]; `0` = white in time |
| `:forcing_seed` | seed of the (private, platform-independent) random draw — same seed, same run |
| `:forcing_normalize` | renormalise every step to inject exactly `ε₀` (Eq. 14); default `true` |
| `:rayleigh_friction` | `ν_r` [1/s] on the momentum; `0` = off |
| `:radiative_relaxation` | `ν_h` [1/s] on `φ − φ_ref`; `0` = off |
| `:sphere_Omega`, `:sphere_gravity` | `Ω` for the `Ro` diagnostic; both for the `pv` output field |
| `:lzonal_mean`, `:zonal_mean_nbins` | the zonal-mean files |

`user_inputs.jl` fills all of them from the planet and the paper's
nondimensional numbers; edit those, not the SI values.

The planet's `Ω` and `g` reach the pointwise callbacks, which receive no
`inputs`, through the typed module globals `SP_OMEGA` / `SP_GRAVITY` declared in
`user_source.jl` and set by `initialize()`.

## Tests

`julia --project=. test/test_sphere_forcing.jl` checks, on this grid: that the
real spherical-harmonic basis is orthonormal under the SEM quadrature; that
`u_F = n̂×∇ₛY` carries `∫|u_F|² = n(n+1)` and is tangential; that the injection
rate is `ε₀` — exactly for white-in-time forcing, statistically for the
Markovian one, by integrating `∂u/∂t = u_F` and reading the energy slope; the
sign and coefficient of both dissipations; the deck switch; reproducibility of
the draw; and the zonal mean of a solid-body rotation.
