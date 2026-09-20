# Proper Orthogonal Decomposition, and the reduced-order model it is the basis of

`src/kernel/rom/` extracts the dominant modes of a running simulation by Proper
Orthogonal Decomposition and leaves behind everything a reduced-order model
(ROM) needs: the modes themselves, their energies, their temporal coefficients,
and the projection/reconstruction operators that map between the full grid and
the reduced coordinates.

**It is a property of the framework, not of a case.** Any problem opts in with
one line in its deck,

```julia
:lpod => true,
```

and supplies nothing else: no `user_*.jl` is involved, the fields default to
every variable of the solution vector whatever the case called them, and the
outputs adapt to the geometry the case already declared — 1-D, 2-D, 3-D or a
spherical manifold. Everything else is a `:pod_*` override (§3).

Three shipped cases have it on:

| case | why it is there |
|:--|:--|
| `problems/AdvDiff/PODbenchmark` | **the reference benchmark** (§7): a problem whose POD is known in closed form, so the implementation is checked against arithmetic |
| `problems/ShallowWater/SWsphere` | the barotropically unstable Galewsky jet — one instability growing on a steady background, whose leading pair of modes *is* the unstable wave |
| `problems/ShallowWater/SWsphere_ScottPolvani` | forced-dissipative turbulence with giant-planet parameters — broadband turbulence out of which jets and vortices organise themselves, which is the setting POD was invented for |

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

`:lpod => true` is the whole of the minimum. Everything below overrides a
default that already works.

| key | default | meaning |
|:--|:--|:--|
| `:lpod`                   | `false`                 | run the decomposition |
| `:pod_fields`             | `[:all]`                | which fields — see below |
| `:pod_nsnapshots`         | the output cadence      | sampling **intervals** over the POD window; one more snapshot than this |
| `:pod_tstart`             | `:tinit`                | start of the window — set it past a transient |
| `:pod_tend`               | `:tend`                 | end of the window |
| `:pod_nmodes`             | `0` (all)               | modes kept, saved and written |
| `:pod_nmodes_plot`        | `6`                     | modes drawn in the mode figure |
| `:pod_subtract_mean`      | `true`                  | decompose the fluctuation, not the field |
| `:pod_method`             | `:auto`                 | `:svd` \| `:snapshot` (§6) |
| `:pod_nlon`, `:pod_nlat`  | `720`, `360`            | equirectangular raster (manifold cases) |
| `:pod_nx`, `:pod_ny`      | `400`, `400`            | (x,y) raster (2-D cases) |
| `:pod_write_vtk/png/data` | `true`                  | which outputs to produce |
| `:pod_write_snapshots`    | `false`                 | also dump the raw snapshots, for an offline re-run (§8) |
| `:pod_time_scale`         | `1.0`                   | multiplies `t` on the coefficient plots |
| `:pod_time_label`         | `"t"`                   | its axis label |
| `:pod_cmap`               | `:balance`              | diverging colour map for the modes |
| `:pod_max_memory_gb`      | `4.0`                   | refuse a snapshot set larger than this |

The default snapshot count is `:ndiagnostics_outputs`, or — for the many cases
that ask for `:diagnostics_at_times` instead, which forces that key to zero —
the number of output times they asked for, and failing both, 50.

### Which fields

An entry of `:pod_fields` may be:

| entry | meaning |
|:--|:--|
| `"rho"`, `:theta`, … | a variable **by name**, matched first against the case's solution variables and then against its output variables (the ones `user_uout!` derives). Case-insensitive. |
| `3` | that column of the solution vector, for a case whose variables are unnamed |
| `:all` | **the default** — one scalar target per solution variable |
| `:allout` | one scalar target per output variable instead |
| `:state` | every solution variable, stacked into ONE vector-valued target |

A **vector-valued** target is decomposed jointly: its components share one set of
temporal coefficients, so a mode is a state of the whole system rather than
unrelated scalars, and the reported energy is the energy of the vector. That is
what a Galerkin ROM of a coupled system is projected onto; decomposing the
components separately answers a different — and for that purpose wrong —
question.

On a **spherical shell** four derived fields are available besides, for the
shallow water system `q = [φ, φu, φv, φw]` with Cartesian momentum:

| symbol | components | what it is |
|:--|:--|:--|
| `:vorticity` | `ζ`          | relative vorticity `n̂·(∇ₛ×u)`, supplied by the time loop |
| `:h`         | `h`          | fluid depth `φ/g` |
| `:u`, `:v`   | `u_λ`, `u_φ` | zonal and meridional velocity |
| `:velocity`  | `u_λ, u_φ`   | the horizontal velocity as ONE vector target |

The velocity components are the **tangent-basis** ones. The raw Cartesian
components of a flow on a sphere are a property of the frame, not of the flow:
their modes would show a rigid zonal jet as a dipole straddling the prime
meridian.

### Cost

One `npoin × ncomp × nsnap` array of `Float64` per field, held for the run. On
the shipped cubed sphere (`npoin = 15 002`), 101 snapshots of `ζ` and `h`
together are 24 MB; on the 64×64 one (`npoin ≈ 6.1e5`), 61 snapshots of `ζ` are
300 MB. The decomposition is `O(N K²)` and runs once, at the end.
`:pod_max_memory_gb` refuses anything larger, with the arithmetic in the message.

### What it refuses

`:lamr => true`. Adaptive refinement changes `npoin` between snapshots, so
consecutive snapshots are vectors in different spaces and there is no
correlation matrix to form. (A ROM on an adapting grid needs the snapshots
interpolated onto a common reference mesh first — different machinery.)

## 4. What it writes

For each field `<f>`, in `:output_dir`. The first group is the same in every
dimension — it is a property of the decomposition, not of the grid:

| file | contents |
|:--|:--|
| `pod_<f>_spectrum.png`     | the energy spectrum, and the cumulative energy |
| `pod_<f>_coefficients.png` | `a_i(t)`, and the `(a₁,a₂)` phase portrait |
| `pod_<f>_spectrum.csv`     | `λ_i`, `E_i`, `ΣE_i`, and the truncation error (3) |
| `pod_<f>_coefficients.csv` | `a_i(t_k)` |
| `pod_<f>.jld2`             | the basis itself — what a ROM reads back (§8). **Under MPI: one file per rank**, `pod_<f>_rank0000.jld2`, … (§6) |

The modes themselves are written in the form their geometry asks for:

| geometry | modes |
|:--|:--|
| **1-D** | `pod_<f>_modes.csv` — a table of `x, mean, mode_001, …` sorted by `x` — and `pod_<f>_modes.png`, the modes as curves |
| **2-D** | `pod_<f>.vtu`, the quads of the mesh with one point-data array per mode, and `pod_<f>_modes.png`, filled contours on an (x,y) raster |
| **3-D** | `pod_<f>.vtu`, the hexahedra of the mesh. **No PNG**: a 3-D mode needs a slice or an isosurface, which is what ParaView is for, and any projection this code picked would be the wrong one |
| **manifold** | `pod_<f>.vtu` on the shell itself, so the modes land on the geometry they were computed on, and `pod_<f>_modes.png`, an equirectangular map |

plus `pod_<f>_mode_001.png`, … one per mode, and `pod_<f>_mean.png`, for
everything but 3-D.

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
two numbers per snapshot. §7 makes the same statement exactly.

**The colour scale of a mode is symmetric about zero**, clipped at the 99.8th
percentile of `|φ_i|`, and shared by all components of a vector mode. A mode has
no preferred sign — (2) fixes it only up to `±1`, and the code picks "the
largest-magnitude entry is positive" purely so that re-running a case does not
invert the colours — so a scale that is not symmetric would invent structure
that is not there.

### The rasters

Maps are **rendered from the element tiling**, not interpolated from scattered
points: every `(ngl-1)²` sub-quad of every element is split into two triangles
and filled by barycentric interpolation (`src/io/plotting/mesh_raster.jl`).
Consequences that matter for a POD mode: there is no smoothing parameter that
could round off the small-scale structure the higher modes consist of; values
are convex combinations of nodal values, so the raster cannot overshoot the data
and the colour scale means what it says; and element and panel seams are
invisible, because the triangles are drawn from the connectivity and meet
exactly there.

Spherical fields are drawn in the **equirectangular (plate carrée)** projection:
longitude and latitude used directly as the plot axes. It is neither conformal
nor equal-area; it is used because it is the identity map on the coordinates the
data already carries, so nothing in the picture is an artefact of the
projection. Its two singular places are handled explicitly — the dateline by
unwrapping each quad's longitudes and drawing it at every 360° offset that
touches the canvas, and the polar caps, where the projection is genuinely
degenerate, by filling from the nearest node **in 3-D**, where there is no
singularity.

`plot_mesh_field(f, mesh, rec, "out.png")` draws any nodal field the same way,
and `plot_sphere_field(f, mesh, "out.png")` does it for a shell without a
recorder to hand — a solution at one output time, a reconstruction error, the
difference between two runs.

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
convention is resolved globally — the sign of a mode is fixed by its
largest-magnitude entry, which sits on one rank, and a rank deciding locally
would hand back a mode negated with respect to its neighbours, i.e. a
discontinuity straight through the partition seam.

### What is global, and what is a piece

| quantity | under MPI |
|:--|:--|
| the spectrum `λ`, `E`, `Σλ`, the coefficients `a_i(t)` | **global** — they come out of a correlation matrix summed across ranks, and every rank holds the same numbers. Rank 0 writes the two CSVs once. |
| the modes `Φ` and the mean `q̄` | **partitioned** — each rank holds the slice living on its own nodes, and nothing is gathered. That is what makes the decomposition scale. |
| `pod_<f>.pvtu` | **complete**: every rank writes its piece and ParaView reassembles them, so the picture of the modes is whole. |
| `pod_<f>_rank0000.jld2`, … | **one per rank**, each with its slice plus `ip2gip` (the global node numbers, for stitching) and `rank`/`nparts`. `pod_load` says which piece it got, so a partition can never be mistaken for the whole. A parallel ROM restarted on the same partition reads its own file and stitches nothing. |
| PNGs | **not written**: the raster needs the whole domain on one rank. The `.pvtu` (or, in 1-D, the CSV) carries the same modes. |

The per-step cost is nothing but the snapshot copy; the one collective sequence
is at the end of the run, and it moves `K²` numbers per field, not `N`.

---

## 7. The reference benchmark

`problems/AdvDiff/PODbenchmark` and `test/test_pod_benchmark.jl` are the same
problem: **linear advection of a multi-harmonic wave**, the standard test case of
the transport-dominated model-reduction literature, chosen because its POD can be
written down in closed form. The implementation is therefore checked against
arithmetic rather than against another run.

```
∂u/∂t + c ∂u/∂x = 0 ,  x ∈ [0,L) periodic ,  u(x,0) = Σ_{j=1}^{J} A_j cos(2πj x/L + ϕ_j)
```

Averaging over one period of the translation gives a convolution kernel, so the
eigenfunctions are the Fourier modes and each wavenumber contributes a
**two-dimensional** eigenspace:

```
λ_{2j−1} = λ_{2j} = A_j² L/4 ,   span{ cos(2πjx/L), sin(2πjx/L) }            (B2)
E_j = A_j²/(2 Σ_i A_i²) ,   Σλ = (L/2) Σ_j A_j² ,   ε(2m)² = Σ_{j>m}A_j²/Σ_j A_j²
```

With `A = (1, ½, ¼)` and `L = 2` the computed spectrum is, to every digit
printed:

| mode | λ computed | λ exact | E [%] | E exact [%] |
|---:|---:|---:|---:|---:|
| 1, 2 | 0.5000000000    | 0.5     | 38.0952 | 38.0952 |
| 3, 4 | 0.1250000000    | 0.125   |  9.5238 |  9.5238 |
| 5, 6 | 0.0312500000    | 0.03125 |  2.3810 |  2.3810 |

`Σλ = 1.3125` exactly, and `max|ΦᵀMΦ − I| = 1.3e-15`.

**What the benchmark checks that a spectrum alone would not:**

1. the spectrum against (B2), to `1e-10` relative — a wrong normalisation, a
   missing `1/K`, a mean that was not removed;
2. the **eigenspaces**, not the modes: `λ_{2j−1} = λ_{2j}` exactly, so the two
   members of a pair are defined only up to a rotation between them, and any
   code claiming a particular pair there claims something the problem does not
   determine. What *is* determined is the plane they span;
3. **the inner product**, by running on a deliberately non-uniform
   (Chebyshev–Lobatto) grid and asserting both that the mass-weighted answer is
   right *and* that the unweighted one is measurably wrong. On a uniform grid
   the two agree, so a uniform benchmark cannot tell a correct implementation
   from one that silently dropped the mass matrix;
4. **the ROM quantities**: the a-priori truncation error curve against the
   measured reconstruction error at every rank, and the constant radius of each
   pair's phase portrait;
5. **the Kolmogorov n-width**. With every `A_j` equal the spectrum goes *flat*:
   `2J` modes each carrying `1/(2J)`, and truncation buys nothing. That is the
   known limitation of every linear reduced basis for transport, and it is the
   most consequential thing a POD implementation can get wrong — a decaying
   spectrum there would promise a reduced-order model that cannot exist.

The deck ends its POD window one sampling interval short of `:tend`, which is
the one subtlety worth carrying over to other cases: the degeneracy `λ₁ = λ₂` is
a statement about averaging over a *whole* period, and sampling both ends of one
period repeats the zero phase and splits every pair by `(K/2+1)/(K/2)` — 5 % at
41 snapshots — for a reason that has nothing to do with the decomposition. The
sampling times are added to the integrator's `tstops`, so the snapshots are
taken *at* them rather than at the first step after.

## 8. From the basis to a reduced-order model

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

With `:pod_write_snapshots => true` the raw snapshot set and its quadrature
weights are written out too, which is what lets a decomposition be **redone
offline** — over a shorter window, about a different mean, with a different rank
— without re-running the simulation, the expensive half of building a ROM:

```julia
d = JLD2.load("output/pod_u_snapshots.jld2")
P = pod_from_snapshots(d["snapshots"], d["weights"], d["t"]; nmodes = 8)
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

## 9. Tests

| file | what it covers | needs |
|:--|:--|:--|
| `test/test_pod_benchmark.jl` | **the reference benchmark** (§7): the closed-form POD of an advected multi-harmonic wave, to `1e-10` | `Test` only — no package instantiation |
| `test/test_pod.jl` | the decomposition itself and the rasterizer, against problems with closed-form answers | `Test` only — no package instantiation |
| `test/test_pod_sphere.jl` | the wiring: weights, extractors, recorder, writers, and the raster on the real cubed sphere | `using Jexpresso` |

```bash
julia test/test_pod_benchmark.jl
julia test/test_pod.jl
julia --project=. test/test_pod_sphere.jl
```

The first two need no package instantiation — both files under test are free of
Jexpresso types by design — so they run in seconds in the CI registry job,
before anything heavy is built.
