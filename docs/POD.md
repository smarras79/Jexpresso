# Proper Orthogonal Decomposition, and the reduced-order model it is the basis of

`src/kernel/rom/` extracts the dominant modes of a running simulation by Proper
Orthogonal Decomposition and leaves behind everything a reduced-order model
(ROM) needs: the modes themselves, their energies, their temporal coefficients,
and the projection/reconstruction operators that map between the full grid and
the reduced coordinates.

Two shipped cases have it switched on, and they are the two ends of what POD is
used for:

* `problems/ShallowWater/SWsphere` — the barotropically unstable Galewsky,
  Scott & Polvani (2004) jet. A single instability growing on a steady
  background, so the decomposition has something definite to find and the
  leading pair of modes *is* the unstable wave.
* `problems/ShallowWater/SWsphere_ScottPolvani` — forced-dissipative
  shallow-water turbulence with giant-planet parameters (Scott & Polvani 2007).
  Broadband turbulence out of which jets and vortices organise themselves,
  which is the setting POD was invented for.

Nothing in the implementation is specific to either beyond the list of
extractable fields.

---

## 1. What POD computes

Given `K` snapshots `q(x, t_k)` of a field, POD asks for the single spatial
structure `φ(x)` that captures, on average, the most of the signal:

```
maximise   ⟨ |(q', φ)|² ⟩ / (φ, φ)                                          (1)
```

with `q' = q - q̄` the fluctuation about the temporal mean and `⟨·⟩` the average
over snapshots. The maximiser is the leading eigenfunction of the two-point
correlation operator, and the whole ordered family follows from its spectrum:

```
R φ_i = λ_i φ_i ,    R(x,x') = ⟨ q'(x,t) q'(x',t) ⟩                         (2)
```

The result is an orthonormal basis **ordered by energy**. `λ_i` is the mean
square of the projection onto `φ_i`, and no other basis of any given size `r`
captures more of `⟨‖q'‖²⟩` than the first `r` POD modes. That optimality is the
reason POD is the starting point of a ROM: truncating it is the least damaging
truncation available, and the error it leaves behind is known *in advance*:

```
⟨‖q' - Σ_{i≤r} a_i φ_i‖²⟩ / ⟨‖q'‖²⟩  =  Σ_{i>r} λ_i / Σ_i λ_i               (3)
```

The expansion

```
q(x,t) ≈ q̄(x) + Σ_{i≤r} a_i(t) φ_i(x)
```

turns the PDE into `r` ODEs for the `a_i` once the residual is required to be
orthogonal to the retained modes (Galerkin projection). This module supplies
`{φ_i}`, `a_i(t_k)` and the operators such a model is written in terms of.

*References:* Lumley (1967); Sirovich, *Q. Appl. Math.* **45** (1987) 561–571;
Holmes, Lumley, Berkooz & Rowley, *Turbulence, Coherent Structures, Dynamical
Systems and Symmetry*, 2nd ed., CUP (2012).

---

## 2. The inner product is the mass matrix

Equation (1) is posed in `L²(Ω)`, so the discrete inner product must be the
discrete `L²` one:

```
(u, v) = ∫_Ω u v dΩ = Σ_ip M_ip u_ip v_ip                                   (4)
```

with `M` the diagonal SEM mass matrix — on the shell, the surface Jacobian times
the LGL weights (`metrics.M`).

**This is not a detail.** Feeding the raw snapshot matrix to an off-the-shelf
SVD gives the Euclidean inner product `Σ u_ip v_ip`, which weights every node
equally — and LGL nodes are *not* equally spaced. They cluster towards element
edges like `1/nop²`, so the edge nodes of every element would count for several
times their share of the domain, the "modes" would not be orthogonal in `L²`,
and (3) would no longer be the error of anything. `test/test_pod.jl` asserts
both that the modes are orthonormal under (4) and that they are *not* orthonormal
under the Euclidean product, so that this cannot regress silently.

It is also why the decomposition is done **inside the run** rather than
afterwards from the VTK files: the mass matrix belongs to the run that produced
the data and is in no output file. The second reason is sampling — POD wants
snapshots uniformly spaced and dense enough to resolve the dynamics, while the
VTK cadence is chosen to keep a movie small, so `:pod_nsnapshots` carries its own
clock and the snapshots never touch the disk.

---

## 3. Switching it on

Everything is driven from `user_inputs.jl`. Only `:lpod` is required.

| key | default | meaning |
|:--|:--|:--|
| `:lpod`                   | `false`                 | run the decomposition |
| `:pod_fields`             | `[:vorticity, :h]`      | which fields (below) |
| `:pod_nsnapshots`         | `:ndiagnostics_outputs` | sampling **intervals** over the window; one more snapshot than this |
| `:pod_tstart`             | `:tinit`                | start of the POD window — set it past a transient |
| `:pod_tend`               | `:tend`                 | end of the window |
| `:pod_nmodes`             | `0` (all)               | modes kept, saved and written |
| `:pod_nmodes_plot`        | `6`                     | modes drawn in the mode figure |
| `:pod_subtract_mean`      | `true`                  | decompose the fluctuation, not the field |
| `:pod_method`             | `:auto`                 | `:svd` \| `:snapshot` (§6) |
| `:pod_nlon`, `:pod_nlat`  | `720`, `360`            | equirectangular raster size |
| `:pod_write_vtk`          | `true`                  | modes as point data on the sphere |
| `:pod_write_png`          | `true`                  | the three standard figures |
| `:pod_write_data`         | `true`                  | CSV spectrum/coefficients + `.jld2` basis |
| `:pod_time_scale`         | `1/86400`               | multiplies `t` on the coefficient plots |
| `:pod_time_label`         | `"t [days]"`            | its axis label |
| `:pod_cmap`               | `:balance`              | diverging colour map for the modes |
| `:pod_max_memory_gb`      | `4.0`                   | refuse a snapshot set larger than this |

### Fields that can be decomposed

For the shallow water system on the shell (`q = [φ, φu, φv, φw]`, Cartesian
momentum):

| symbol | components | what it is |
|:--|:--|:--|
| `:vorticity` | `ζ`             | relative vorticity `n̂·(∇ₛ×u)` — what the Galewsky test is judged on |
| `:h`         | `h`             | fluid depth `φ/g` |
| `:phi`       | `φ`             | geopotential |
| `:u`, `:v`   | `u_λ`, `u_φ`    | zonal, meridional velocity |
| `:velocity`  | `u_λ, u_φ`      | the horizontal velocity as **one** two-component target |
| `:state`     | `φ, φu, φv, φw` | the conservative state — the basis a Galerkin ROM is projected onto |

A vector-valued target is decomposed **jointly**: its components share one set of
temporal coefficients, so a mode is a velocity field rather than two unrelated
scalars, and the reported energy is the energy of the vector. Decomposing the
components separately answers a different — and for a ROM, wrong — question.

The velocity components are the **zonal and meridional** ones, projected onto
the local tangent basis. The raw Cartesian components of a flow on a sphere are
a property of the frame, not of the flow: their modes would show a rigid zonal
jet as a dipole straddling the prime meridian.

### Cost

One `npoin × ncomp × nsnap` array of `Float64` per field, held for the run. On
the shipped cubed sphere (`npoin = 15 002`), 101 snapshots of `ζ` and `h`
together are 24 MB. The decomposition is `O(N K²)` and runs once, at the end:
sub-second at these sizes. `:pod_max_memory_gb` refuses anything larger, with
the arithmetic in the message.

---

## 4. What it writes

For each field `<f>`, in `:output_dir`:

| file | contents |
|:--|:--|
| `pod_<f>.vtu` (`.pvtu` under MPI) | the mean and every retained mode as point data **on the sphere**, for ParaView |
| `pod_<f>_modes.png`               | the leading modes as equirectangular maps, one panel each |
| `pod_<f>_mode_001.png`, …         | the same modes one per file |
| `pod_<f>_spectrum.png`            | the energy spectrum, and the cumulative energy |
| `pod_<f>_coefficients.png`        | `a_i(t)`, and the `(a₁,a₂)` phase portrait |
| `pod_<f>_mean.png`                | the temporal mean |
| `pod_<f>_spectrum.csv`            | `λ_i`, `E_i`, `ΣE_i`, and the truncation error (3) |
| `pod_<f>_coefficients.csv`        | `a_i(t_k)` |
| `pod_<f>.jld2`                    | the basis itself — what a ROM reads back (§7) |

### Reading the figures

**The spectrum is the plot to read first**, because it says whether a
reduced-order model is possible at all: a spectrum that falls off a cliff after
a handful of modes means a handful of ODEs can carry the flow, and one that
decays slowly means they cannot. The `90 %` and `99 %` lines on the cumulative
panel are annotated with the rank each threshold needs.

**The modes come in pairs when the structure travels.** A travelling wave cannot
be written with one standing pattern, so POD splits it into two of nearly equal
energy, a quarter wavelength apart. `λ₁ ≈ λ₂` in the spectrum and a **circle** in
the `(a₁,a₂)` phase portrait are the signature; a pair of standing structures
traces a line instead. For the Galewsky jet the radius of that circle grows
exponentially and then saturates — that is the barotropic instability, read off
two numbers per snapshot.

**The colour scale of a mode is symmetric about zero** and clipped at the
99.8th percentile of `|φ_i|`. A mode has no preferred sign — (2) fixes it only
up to `±1`, and the code picks "the largest-magnitude entry is positive" purely
so that re-running a case does not invert the colours — so a scale that is not
symmetric would invent structure that is not there.

### The projection

Maps are drawn in the **equirectangular (plate carrée)** projection: longitude
and latitude used directly as the plot axes. It is neither conformal nor
equal-area; it is used because it is the identity map on the coordinates the
data already carries, so nothing in the picture is an artefact of the projection.

The raster is a **rendering of the element tiling**, not an interpolation of
scattered points: every `(ngl-1)²` sub-quad of every element is split into two
triangles and filled by barycentric interpolation
(`src/io/plotting/equirectangular.jl`). Consequences that matter for a POD mode:
there is no smoothing parameter that could round off the small-scale structure
the higher modes consist of; values are convex combinations of nodal values, so
the raster cannot overshoot the data and the colour scale means what it says;
and the panel seams of the cubed sphere are invisible, because the triangles are
drawn from the connectivity and meet exactly there. The dateline is handled by
unwrapping each quad's longitudes and drawing it at every 360° offset that
touches the canvas; the polar caps, where the projection is genuinely
degenerate, are filled from the nearest node **in 3-D**, where there is no
singularity.

`plot_sphere_field(f, mesh, "out.png")` draws any nodal field on the shell the
same way — a vorticity field at one output time, a reconstruction error, the
difference between two runs.

---

## 5. Using it from the REPL

```julia
using Jexpresso

# … after a run, or on snapshots assembled by hand:
#   X :: npoin × ncomp × nsnap,  w :: the diagonal mass matrix
P = pod_from_snapshots(X, w, t; name = "vorticity", nmodes = 20)

P.λ, P.energy, P.cumenergy    # the spectrum
P.Φ[:, 1, 3]                  # mode 3, component 1, as a nodal vector
P.a[:, 1]                     # a₁(t_k)

pod_rank_for_energy(P, 0.99)  # how many modes carry 99 % of the energy
pod_truncation_error(P)       # the a-priori error curve, Eq. (3)

a  = pod_project(P, w, q)     # full state → reduced coordinates
q̂ = pod_reconstruct(P, a)    # reduced coordinates → full state
```

---

## 6. Two algorithms, and MPI

`:pod_method` selects between two routes to the same answer.

* **`:snapshot`** — Sirovich's method of snapshots. Form the `K×K` correlation
  matrix `C_kl = (q'_k, q'_l)/K`, take its eigendecomposition, push the
  eigenvectors back through the snapshots. Cost `O(N K²)`. **The parallel
  path**: `C` is a sum over nodes, so each rank forms its own partial `C` and one
  `Allreduce` of `K²` numbers completes it. The eigenproblem is then solved
  redundantly and identically on every rank, and each rank builds its own slice
  of the modes with no further communication — the field is never gathered.
* **`:svd`** — the SVD of the mass-weighted snapshot matrix `Y = M^{1/2} Q'`.
  Mathematically identical (`σ_i²/K = λ_i`, and the left singular vectors are
  `M^{1/2} φ_i`) but it never forms `C` and so never squares the condition
  number. The difference shows up exactly where it matters for a ROM — in the
  tail of the spectrum, where `:snapshot` loses the modes below `λ₁·ε` while the
  SVD resolves to `λ₁·ε²`. **Serial default** for that reason.

Under MPI two further details keep the result independent of the rank count: a
node shared by several ranks is counted by its **owner** only (`mesh.gip2owner`,
exactly as `sphere_diagnostics` does for the conserved integrals), and the sign
convention is resolved globally. PNG output is skipped under MPI — the
equirectangular raster needs the whole sphere on one rank — and the `.vtu`
carries the same modes.

---

## 7. From the basis to a reduced-order model

`pod_<f>.jld2` is the hand-off. It stores plain arrays rather than a serialised
struct on purpose (the same reason `sem_setup.jl` gives for its metric cache: a
stored struct stops loading the day its definition changes, and a basis is worth
more than the struct that held it).

```julia
P = pod_load("output/pod_state.jld2")   # → St_pod, as the solver had it
a0 = pod_project(P, w, q0)              # the ROM's initial condition
# … march a₁…a_r in time …
q  = pod_reconstruct(P, a)              # back to the grid, to plot or to restart
```

What is deliberately **not** here yet: the Galerkin projection of the shallow
water operator onto `{φ_i}`, which is where a ROM stops being a change of basis
and starts being a model. The pieces it needs — an `L²`-orthonormal basis under
the same mass matrix the solver uses, the mean, the coefficients, and the
projection operator — are all above.

The other classical decompositions build on exactly the same snapshot set, with
the same inner product, and are the natural next additions: **DMD**
(Schmid 2010), which fits `q_{k+1} ≈ A q_k` in the POD subspace and returns
modes with a frequency and a growth rate each — for the Galewsky jet, the growth
rate of the barotropic instability as a number rather than as the slope of a
curve; **balanced POD** (Rowley 2005) when the ROM is meant for control rather
than for description; and **spectral POD** (Towne, Schmidt & Colonius 2018) for
statistically stationary flows.

---

## 8. Tests

| file | what it covers | needs |
|:--|:--|:--|
| `test/test_pod.jl` | the decomposition itself and the rasterizer, against problems with closed-form answers | `Test` only — no package instantiation |
| `test/test_pod_sphere.jl` | the wiring: weights, extractors, recorder, writers, and the raster on the real cubed sphere | `using Jexpresso` |

```bash
julia test/test_pod.jl
julia --project=. test/test_pod_sphere.jl
```
