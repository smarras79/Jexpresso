# SWsphere_POD — the Galewsky jet, set up so that its POD means something

```julia
julia --project=.
julia> using Jexpresso
julia> Jexpresso.run_case("ShallowWater", "SWsphere_POD")
```

A few minutes on one core. Everything below was **measured on this deck, through
`run_case`**, so a run that reproduces it has a working POD from the grid to the
figure.

The output lands in `output/ShallowWater/SWsphere_POD/output/` (`run.jl` nests
`:output_dir` under the equations and case names).

## What it is

The same physics, grid and discretisation as
[`SWsphere`](../SWsphere/README.md) — the Galewsky, Scott & Polvani (2004)
barotropically unstable jet on the shipped cubed sphere, CG spectral elements,
SSP-RK3 with the Lagrange projection after every stage — configured for the
decomposition rather than for the longest possible integration.

Three differences from `SWsphere`, and each one is about the POD:

| | `SWsphere` | here | why |
|---|---|---|---|
| `:tend` | 20 days | **6 days** | day 6 is where the instability has rolled up — and where the *fluctuation about the temporal mean* is the instability rather than numerical noise |
| `:lfilter` | `false` | **`true`** | `SWsphere` ships `false`, which contradicts its own comments and its README (that configuration goes non-finite at 2.005 d). With the filter off the vorticity fluctuation is 90 % sub-element and the POD correctly returns modes of *that* |
| `:pod_fields` | `[:vorticity, :h]` | **`[:vorticity, "h", :velocity]`** | one derivative field and two primitive ones, so the report's `sub-element content` lines can be compared directly |

## What it should print

```
 #   vorticity: 12 modes from 25 snapshots (svd), Σλ = 2.127231e+04
 #     orthonormality  max|ΦᵀMΦ - I| = 1.22e-15
 #     snapshots       max|q_k - q_1| = 1.868e-04 ,  max|q| = 1.125e-04 ,  ratio = 1.66e+00
 #     field rms       √(Σλ/area) = 6.457728e-06  (area = 5.100997e+14)
 #     mode        λ            E [%]      ΣE [%]      rms_i = √(λ_i/area)
 #        1    1.29364e+04     60.813     60.813     5.035927e-06
 #        2    3.84668e+03     18.083     78.896     2.746093e-06
 #        3    2.15332e+03     10.123     89.019     2.054598e-06
 #        4    1.02591e+03      4.823     93.842     1.418168e-06
 #     90% of the energy is in the first 4 modes
 #     99% of the energy is in the first 8 modes
 #     sub-element content   fluctuation 0.581 , mode 1 0.594   (1 = all grid-scale)
```

with `Δt = 214.2149 s ; 2420 steps to t = 518400.0 s (6.000 days)`, and
`δmass/mass = 2.7e-11`, `δE/E = -5.5e-04` at the end.

## What it should draw

`output/pod_vorticity_modes.png` (and `.pdf`): **the unstable wave train
confined to the 30–60°N band, its wavenumber rising with the mode index.** That
is the barotropic instability, decomposed. `pod_vorticity_coefficients.png`
shows the `(a₁,a₂)` phase portrait as a growing spiral — the amplitude of the
instability against its phase, which is the growth rate read off two numbers per
snapshot.

## The three lines to read, and why

**`snapshots … ratio`** — how far the field moved over the sampling window
against how big it is. At round-off the snapshots are copies of one state and
the modes are of the differences between copies. Here it is 1.66: the field
changed by more than its own size, which is what an instability does.

**`rms_i = √(λ_i/area)`** — mode *i*'s contribution to the RMS of the **field**,
in the field's units. A mode's own magnitude says nothing (Φ is normalised to
`∫Φ² dΩ = 1`, so every mode is ~`1/√area` whatever the flow does); this is the
number to compare with the `max|ζ|` the diagnostics print.

**`sub-element content`** — how much of the fluctuation, and of the leading
mode, lives below the element scale. The two should match; when they do, the
decomposition is reporting what it was fed.

That last line is the one that distinguishes this case from a run whose modes
look like grid noise. Compare the three fields this deck decomposes:

| field | sub-element content (fluctuation / mode 1) | E₁ … E₄ [%] |
|---|---|---|
| `vorticity` (a derivative) | 0.581 / 0.594 | 60.8 / 18.1 / 10.1 / 4.8 |
| `h` | **0.366 / 0.376** | 66.6 / 20.1 / 6.9 / 2.5 |
| `velocity` (u_λ, u_φ jointly) | **0.426 / 0.404** | 61.7 / 21.0 / 9.2 / 4.0 |

The primitive fields are markedly smoother than the derivative, as they should
be, and in every row the mode matches the fluctuation it came from — which is
the statement that the decomposition is faithful.

and compare *those* with the same case at `:tend => 1*24*3600` and
`:lfilter => false`, where the vorticity gives **0.94 / 0.89** and the modes are
an unreadable checkerboard. Nothing is broken in that run either: `ζ = ∇ₛ×u` is
a derivative, its smooth part (the jet) is nearly steady so the temporal mean
absorbs it, and what is left for the modes is the grid-scale velocity noise that
the curl amplifies by the wavenumber. The decomposition is correct; the question
put to it was not. Decomposing `h` or `velocity` over the same window, or
widening the window until the flow really changes, is the fix — and the report
now says so itself.

## Going further

* `:pod_write_snapshots => true` is on, so `output/pod_<f>_snapshots.jld2` holds
  the raw snapshots and their quadrature weights. The decomposition can be
  redone offline over a shorter window or a different rank without re-running:

  ```julia
  d = JLD2.load("output/pod_vorticity_snapshots.jld2")
  P = pod_from_snapshots(d["snapshots"], d["weights"], d["t"]; nmodes = 6)
  ```

* `output/pod_<f>.jld2` is the basis a reduced-order model reads back
  (`pod_load`), with `pod_project` / `pod_reconstruct` as the map to and from
  the reduced coordinates.

* In parallel (`mpiexec -n N`) everything holds: the spectrum and the
  coefficients are global, the modes are written as one `.pvtu` (open that, not
  the pieces in `pod_vorticity/`) plus one `.jld2` per rank, and the maps are
  rendered by reducing the pixel canvas rather than gathering the mesh.

Full documentation: [`docs/POD.md`](../../../docs/POD.md).
