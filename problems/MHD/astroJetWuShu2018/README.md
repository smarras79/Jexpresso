# `astroJetWuShu2018` — the classical magnetized astrophysical jet

2D ideal GLM-MHD. A Mach 800 dense beam injected into a static, strongly
magnetized ambient medium; the standard closing test of the
positivity-preserving MHD literature.

```
julia --project=. src/Jexpresso.jl MHD astroJetWuShu2018
```

**Read [§6, "What this case is actually up against"](#6-what-this-case-is-actually-up-against)
before you run it.** Every published solution of this test uses an explicit
positivity-preserving limiter. Jexpresso has none: DynSGS is the whole of the
stabilization here, and that is the experiment.

### Two deliberate departures from the paper — state these in any write-up

Both are at the **inlet boundary datum** and nowhere else. Neither touches the
equations, the constants, the domain or the shock physics, and each is one
environment variable from the paper's exact condition. Both were forced by
measurement, not chosen: with either one off, the run aborts inside the first
handful of steps, and §10–§12 give the coordinates and the call numbers.

| | default | the paper | restore it with |
|---|---|---|---|
| nozzle lip, in `x` | smootherstep over `2s = 2h`, centred on `\|x\| = 0.05` | a top hat | `JEXPRESSO_AJ_SMOOTH=0` |
| beam turn-on, in `t` | smootherstep over `τ = 2h/u_jet` (125 steps) | impulsive | `JEXPRESSO_AJ_TRAMP=0` |

The injected flux is preserved: `∫φ dx = x0` to five digits, and the beam is at
full strength from `t = τ = 3.1 %` of `tend` onward. What they buy is that the
boundary datum **agrees with the initial condition at `t = 0` and with its free
neighbours at every later time**, instead of fighting both.

---

## 1. Provenance

The magnetized jet was introduced by

> K. Wu, C.-W. Shu, *A provably positive discontinuous Galerkin method for
> multidimensional ideal magnetohydrodynamics*, SIAM J. Sci. Comput. **40**(5)
> (2018) B1302–B1329, **Example 5.6 ("Astrophysical jets")**,
> [doi:10.1137/18M1168042](https://doi.org/10.1137/18M1168042)

by adding a uniform magnetic field to the *gas dynamical* Mach 800 jet of

> D. S. Balsara, *Self-adjusting, positivity preserving high order schemes for
> hydrodynamics and magnetohydrodynamics*, JCP **231** (2012) 7504–7517,

itself after X. Zhang & C.-W. Shu, JCP **229** (2010) 8918. The same test with
the same constants closes most later structure-preserving MHD papers, e.g.
Peng, Sun & Wu, *Structure-preserving oscillation-eliminating DG schemes for
ideal MHD* ([arXiv:2404.16794](https://arxiv.org/abs/2404.16794)), §4.2.6.

Wu & Shu's statement of the problem, verbatim:

> Initially, the physical domain $[-0.5,0.5]\times[0,1.5]$ is filled with a
> uniform static medium with density $0.1\gamma$ and unit pressure, and the
> adiabatic index $\gamma$ is set as $1.4$. Through the inlet part
> ($|x|<0.05$) on the bottom boundary ($y=0$), a dense jet with speed $800$ is
> injected in the $y$-direction with a density of $\gamma$ and a pressure equal
> to the ambient pressure. The fixed inflow beam condition is specified on the
> nozzle $\{y=0,|x|<0.05\}$, and the others are outflow boundary conditions. We
> initialize the magnetic field with magnitude $B_a$ along the $y$-direction.

and their three configurations:

| | | |
|---|---|---|
| (i)   moderately magnetized         | `B_a = √200`   | `β_a = 1e-2` |
| (ii)  strongly magnetized           | `B_a = √2000`  | `β_a = 1e-3` |
| (iii) extremely strongly magnetized | `B_a = √20000` | `β_a = 1e-4` |

(i) is the default here. Wu & Shu note that the third-order locally
divergence-free conservative DG method with a PP limiter of their earlier paper
"is not able to run this test with `B_a ≥ √200`" — i.e. the *default* is already
past what a PP-limited third-order DG scheme could do at the time.

---

## 2. Every constant, with its value

Nothing here is dimensional; the problem is stated in its own units.

### Fixed by the paper

| symbol | value | where |
|---|---|---|
| `γ` | `1.4` | `user_flux.jl`, `γ_mhd` |
| domain | `[-0.5, 0.5] × [0, 1.5]` | `AJ.geo` |
| nozzle half-width | `0.05` | `user_flux.jl`, `AJ_XNOZZLE` |
| ambient `ρ` | `0.1γ = 0.14` | `AJ_RHO_AMB` |
| ambient `p` | `1` | `AJ_P_AMB` |
| ambient `v` | `(0, 0, 0)` | `initialize.jl` |
| beam `ρ` | `γ = 1.4` | `AJ_RHO_JET` |
| beam `p` | `1` | `AJ_P_JET` |
| beam `v` | `(0, 800, 0)` | `aj_ujet[]` |
| `B` (ambient **and** beam) | `(0, B_a, 0)` | `aj_Ba[]` |
| `B_a` | `√200 ≈ 14.1421356` | default; `JEXPRESSO_AJ_BA2=2000` or `20000` for (ii), (iii) |
| final time | `2e-3` | Wu & Shu draw their figures at `t = 2e-3`; the OEDG paper shows `t = 1e-3, 1.5e-3, 2e-3` |

Note how the numbers are *built* out of `γ`, which is why `γ` is not a free
parameter of this deck:

* beam sound speed `c_j = √(γ p/ρ_j) = √(γ·1/γ) = 1` **exactly**, so "speed 800"
  is exactly **Mach 800**;
* ambient sound speed `c_a = √(1.4/0.14) = √10 ≈ 3.16228`;
* the beam is exactly **10× denser** than what it ploughs into.

### Derived, for `B_a = √200` (the default)

| quantity | ambient | beam |
|---|---|---|
| `ρ` | `0.14` | `1.4` |
| `p` | `1` | `1` |
| sound speed `√(γp/ρ)` | `3.16228` | `1` |
| Alfvén speed `\|B\|/√ρ` | `37.7964` | `11.9523` |
| fast speed bound `√(a²+b²)` | `37.9285` | `11.9940` |
| `\|v\| + c_f` | `37.9285` | `811.994` |
| plasma `β = 2p/\|B\|²` | `1e-2` | `1e-2` |
| `ρE` | `102.5` | `448102.5` |
| of which `p/(γ-1)` | `2.5` (`2.44 %`) | `2.5` (**`5.58e-6`**) |

`c_h`, the GLM divergence-cleaning speed, is the largest of those wave speeds
and is held constant for the run: **`c_h = 811.994`**. `initialize.jl` computes
it from the ambient state *and the prescribed inflow state* and prints it. The
inflow has to be included: the domain is at rest at `t = 0`, so the initial
condition alone would give `c_h ≈ 37.9` while the beam entering from the first
step travels at 800.

For the other two magnetizations, `c_h` barely moves — it is set by the beam
speed, not by the field:

| `B_a` | `β_a` | `c_f` ambient | `c_f` beam | `c_h` |
|---|---|---|---|---|
| `√200`   | `1e-2` | `37.93`  | `11.99`  | `811.99` |
| `√2000`  | `1e-3` | `119.56` | `37.81`  | `837.81` |
| `√20000` | `1e-4` | `377.98` | `119.53` | `919.53` |

which is why one `Δt` per mesh covers all three (§4).

### `∇·B`

`Bx = Bz = 0` and `By = B_a` uniformly, in the initial condition **and** in the
injected beam, so `∇·B = ∂_y B_a = 0` identically — exactly, not to
discretization accuracy. The GLM `ψ` field therefore starts at zero and only
ever carries the divergence error the discretization itself generates. The
inflow prescribes `B` for the same reason: leaving it free at the inlet would
make the nozzle a source of divergence error.

---

## 3. Boundary conditions (`user_bc.jl`)

The four physical-curve tags of the mesh are `bottom`, `right`, `top`, `left`.

* **Nozzle** (`bottom`, `|x| ≤ 0.05`): fixed inflow beam. At Mach 800 every
  characteristic enters, so all nine conserved components are prescribed —
  `ρ = γ`, `v = (0, 800, 0)`, `p = 1`, `B = (0, B_a, 0)`, `ψ = 0`.
* **Everything else** (`bottom` outside the nozzle, `left`, `right`, `top`):
  outflow = *impose nothing*. `build_custom_bcs_dirichlet!` pre-fills `qbdy`
  with a sentinel and copies back only the slots the routine overwrites, so not
  touching `qbdy` leaves the interior solution to convect out and leaves `RHS`
  alone — the same open condition as `problems/CompEuler/ffs_step`'s supersonic
  outflow.

The do-nothing outflow is characteristically correct only where the outgoing
flow is supersonic; on the bottom boundary outside the nozzle, where the cocoon
eventually pushes gas back down subsonically, it is weakly reflective. The time
scale is what makes that acceptable: the ram-pressure balance
`ρ_j(v_j − v_h)² = ρ_a v_h²` gives a head speed

```
v_h = v_j / (1 + √(ρ_a/ρ_j)) = 800/(1 + √0.1) ≈ 608
```

so at `t = 2e-3` the jet head is near `y ≈ 1.2` and has **not** reached the top
of the 1.5-tall box. Nothing has had time to reflect and come back. **That is a
check to run on the output**, not just a remark: a head much slower than
`y ≈ 0.6·800·t` means the stabilization is eating the beam's momentum.

### The nozzle lip

`|x| = 0.05` falls exactly on an element boundary in both shipped meshes
(`0.05/h = 2` and `5`), so there is a node sitting on the discontinuity of the
boundary datum. The closed test `|x| ≤ 0.05` is used, so the inflow patch is
exactly the closed segment of width `0.1` the paper's nozzle is. The tolerance
is an absolute `1e-10` — nine orders below the smallest element — because the
lip comes out of the mesh as `0.050000000000000044`, not `0.05`.

`JEXPRESSO_AJ_SMOOTH=w` blends the injected state into the ambient one across
the lip with `φ(x) = ½(1 − tanh((|x| − 0.05)/w))`. The paper's condition is the
sharp top hat, `w = 0`, which is the default.

**Why the lip is the place to watch.** It is a velocity jump of 800 across one
LGL interval, held open by a Dirichlet condition. That is structurally the same
object that broke `rampCaoEtAl2021_M7`: there a 446 m/s shear across
`2.58e-5 m`, held open by the top Dirichlet, produced the first positivity
repair of the whole run `2.58e-5 m` — exactly one LGL interval — below that
boundary, on step 36, when the fastest signal had travelled 0.07 mm and nothing
could have propagated there. The difference is that on the ramp the jump was a
*bug* in the starting field (a free stream wrongly turned 15°) and could be
removed; here it is the problem statement. So it cannot be fixed, only
represented, and if this case loses positivity somewhere, the nozzle lip is the
first place to look — `x = ±0.05`, `y` within an LGL interval of 0. If you do use it: the blend is
applied to the **primitives** (`ρ, v, p, B`) and `ρE` is rebuilt from them.
Blending the conserved variables instead puts a spurious pressure spike on the
lip — at `φ = ½` the kinetic energy of the mean momentum is not the mean of the
kinetic energies, and `p` comes out ≈ 8000 instead of 1.

---

## 4. Discretization

CG-SEM, LGL nodes, `nop = 4`, Carpenter–Kennedy 2N 5-stage 4th-order explicit.

| mesh | `h` | elements | LGL points | `Δx_min` | default `Δt` | CFL | steps to `2e-3` |
|---|---|---|---|---|---|---|---|
| `AJ_40x60.msh` (default) | `0.025` | 2 400 | 160 × 240 | `4.317e-3` | `5e-7` | `0.09–0.11` | 4 000 |
| `AJ_100x150.msh` | `0.01` | 15 000 | 400 × 600 | `1.727e-3` | `2e-7` | `0.09–0.11` | 10 000 |

`Δx_min = (1 − √(3/7))·h/2 = 0.17267·h` is the smallest LGL gap at `nop = 4`;
the CFL range spans the three magnetizations. `AJ_100x150` is the papers' own
resolution: Wu & Shu compute the right half `[0,0.5] × [0,1.5]` on `200 × 600`
cells (`Δ = 2.5e-3`), and `100 × 150` elements at `nop = 4` give the same nodal
spacing over the **full** width.

The **viscous** limit is not the binding one. DynSGS cannot exceed its own cap
`μ_max = C_max·Δ·(|v| + c_f) = 0.5·(h/5)·c_h`, which is `2.03` on `40x60`, and
`0.5·Δx_min²/μ_max = 4.6e-6` — nine times the advective `Δt`. On `100x150` it is
`1.8e-6`, nine times again. Raising `:dsgs_Cmax` or `:μ` changes that ratio, so
move them and `Δt` together.

The physical **surface must be named `"domain"`** — see the note in `AJ.geo`.
`Geom.jl:163` strips exactly that string from the face labeling, and it has to,
because Gridap propagates a surface group to every edge *interior* to the surface
(4700 of the 4900 edges of the 40×60 mesh). Under any other name every interior
edge is flagged as a boundary edge, which overruns `poin_in_bdy_edge` and — worse
— makes `_dsgs_boundary_pairs!` zero the DynSGS residual over the whole mesh,
turning the shock capturing off silently. Measured on both shipped meshes: with
`"domain"` the flagged set is exactly the 200 / 500 geometric boundary edges.

Both `.msh` files ship with the case, so nothing needs gmsh installed. They were
written directly, with the same entity/node/element layout gmsh produces for
`AJ.geo`, by

```
tools/astro_jet_mesh.py            # regenerates both, byte-identical
tools/astro_jet_mesh.py 60 90      # any nx (multiple of 20) and ny
```

`AJ.geo` is the reference definition and the `gmsh -2` path to the same meshes.
**The one constraint on a mesh of your own**: `0.05·nx` must be an integer
(`nx` a multiple of 20), or the beam edge falls inside an element where the
top-hat datum cannot be represented at all — the script refuses it.

### Full domain, not the half domain

Wu & Shu compute `[0, 0.5] × [0, 1.5]` with a reflecting condition at `x = 0`
and mirror the result for their figures. This case uses the full width instead,
for two reasons: it does not impose a symmetry the beam/cocoon instabilities do
not have, and it puts no artificial wall down the middle of the beam — which is
precisely the kind of boundary that makes a shock case hard. The cost is a
factor 2. Nothing else differs, so a symmetric answer here is evidence that the
run is clean.

---

## 5. Stabilization: DynSGS, in its conserved form

Full description of the model in [`DSGS.md`](../../../DSGS.md); references are
Marras, Nazarov & Giraldo, JCP **301** (2015) 77 and Dao & Nazarov,
J. Sci. Comput. **92**:77 (2022), §4.4. One kinematic coefficient per element,

```
ν|_e = max(0, min( C_max·Δ·(|v| + c_f)|_e ,  C_R·Δ²·max_i ‖R_i‖_∞,e / ‖q_i − ⟨q_i⟩‖_∞,Ω ))
```

with `Δ = h/(N+1)`, `C_R = 1`, `C_max = 0.5` (the paper's own values — the model
is parameter-free; both validated `DSGS_MHD` cases, `brioWu1d` and
`orszagTangBormanis2024`, use exactly these). On `C_max`: the Mach-7 CompEuler
decks run `0.1` on the argument that `Δelem` *is* the nodal spacing, so `0.5`
would overshoot the first-order-upwind bound 5.8×. Every DynSGS kernel in
`SGS.jl` — the MHD ones and the total-energy CompEuler one alike — in fact caps
on `Δ = Δelem/(N+1)`, which at `nop = 4` is `0.2h` against a smallest LGL gap of
`0.17267h`: a factor `1.16`, not `5.8`. So `0.5` is kept here, and their `0.1`
should be read as an empirical 5× loosening of the cap rather than a correction
before being ported. At a discontinuity the normalized residual is `O(10²)`;
smooth flow gives `O(10⁻⁶)`. So the coefficient saturates at the first-order cap
exactly on the Mach shock at the jet head and on the beam/cocoon interface, and
sits at ~0 in the quiescent ambient medium. **That is the property this case is
here to exercise.**

Deck settings that are decisions rather than defaults:

* **`:dsgs_conserved => true`** — the dissipation is `∇·(ν∇q)` on all nine
  *conserved* variables (`user_primitives.jl`), not the physical form (`ν` on
  `ρ`, `ρν` on `u`, `κ` on `T`, `η` on `B`) that
  `orszagTangBormanis2024` uses. The beam is a 10:1 density contact whose
  pressure is `5.6e-6` of its total energy; diffusing `ρ` and `T` separately
  does not keep `p = (γ-1)(ρE − ½ρ|v|² − ½|B|²)` admissible across such a jump —
  that is the mechanism that cost `fluxEmergenceSon2025` its positivity at the
  25× contact of the solar transition region. Diffusing the *conserved* state
  makes a diffused node a convex combination of its neighbours' states, which is
  where any positivity argument for an artificial-viscosity regularization comes
  from. It is also in divergence form on every slot, so the **shock speeds** —
  the thing this benchmark is looked at for — are not altered. Same form as the
  MHD shock tube `brioWu1d`.
* **`:dsgs_sensor => "residual"`** — and here for a reason specific to this
  case, not just because it is the default. It is the only sensor path that
  carries the Dirichlet-boundary treatment of `rhs.jl` (`_dsgs_bdy_zero!`): at a
  constrained node the element's own RHS is the **constraint force**, not an
  under-resolution, and feeding it to the sensor drives `ν` to its cap along the
  whole boundary (measured on the rising bubble, where it blew the run up). This
  case has four non-periodic boundaries and a hard inflow patch. It does not
  blind the model at the nozzle: only the boundary *nodes* are zeroed, and a
  bottom-edge element has 5 of its 25 nodes there, so the other 20 still set its
  coefficient. `JEXPRESSO_AJ_SENSOR=legacy` for the `|∂ₜq|` sensor the
  atmospheric decks use.
* **`:dsgs_Cmin => 0.0`** — *not* `brioWu1d`'s `0.06`, because the floor
  `C_min·Δ·(|v| + c_f)` is proportional to the **local** wave speed, which in the
  beam is 812 and not the `O(1)` of a shock tube. `C_min = 0.06` would put
  `ν = 0.06·(0.025/5)·812 = 0.24` inside the beam, and `√(2νt) = 0.031` over the
  run — 60 % of the beam's own half-width. The floor alone would smear the beam
  away. In the quiescent ambient the same `C_min` is harmless
  (`ν = 0.011`, `√(2νt) = 0.0068`), so if a node-to-node checkerboard mode does
  appear — the one thing a residual sensor is blind to — `JEXPRESSO_AJ_CMIN=0.005`
  to `0.01` is the first lever, and the beam profile is what to check afterwards.
* **`:dsgs_hold_steps => 0`** — the kernel default is 2, and turning it off is a
  decision. The hold exists because on *smooth* data the BDF2 seeded from the
  initial condition makes the residual the whole flux divergence, so the sensor
  reads a fully resolved field as unresolved everywhere (measured on the smooth
  vortex: `ν` at its cap on the very first call). That cannot happen here,
  because this initial condition is **uniform**: `∇·F ≡ 0`, so the residual is
  exactly zero except in the elements the nozzle datum reaches, and `ν` on step
  one is at the cap at the nozzle and zero everywhere else — which is what is
  wanted. Meanwhile the cost of holding is real: the beam is started
  *impulsively* against a medium at rest, the first steps are the most violent
  of the run, and holding `ν` at zero through them integrates exactly the steps
  that need dissipation with none. `ffs_step` turns it off for the same reason
  and states the measurement: with the hold on, that case dies at `t = 1.46e-3`
  in its convex corner instead of reaching `8e-3`.
* **`:dsgs_norms => "domain"`** — the method's own norm over `Ω`, and the only
  choice whose answer does not depend on the MPI partition. On a problem where
  most ranks hold nothing but quiescent ambient gas and a few hold the whole jet,
  `"rank"` is exactly the pathology documented in `SGS.jl`: the quiet ranks
  normalize by their own floor and apply a different viscosity to the same
  solution than their neighbours.
* **`:μ = ones(9)`** — full strength on every slot. In the conserved form `:μ[1]`
  is not optional: it is the mass diffusion `ν∇ρ` that keeps the density jump
  from ringing (on `ffs_step`, zeroing it is the single most destabilizing change
  measured).
* **`:lfilter => false`** — the Boyd–Vandeven filter filters the conservative
  variables *independently*, i.e. it perturbs the three large, nearly cancelling
  terms of `p` by different amounts. At `5.6e-6` of cancellation that is not a
  small error.
* **`:dsgs_nazarov_energy => false`** — it would split the energy primitive into
  a non-thermal part at `ν` and a thermal part at `κ = ρν/Pr`, which breaks the
  single-`ν` convex-combination property above. There is no stratification here
  to conduct heat through.

### One known wrinkle, for the record

`_expansion_visc!` always treats slots 2 and 3 as momentum and builds a
deviatoric stress from their gradients, so under `:dsgs_conserved` those two get
`∇·τ(ρv)` rather than `ν∇²(ρu)`, `ν∇²(ρv)` — for a normal shock, `(4/3)ν`
instead of `ν` on the normal momentum. It is still divergence-form and
dissipative, and it is the established behaviour of the conserved-form MHD path
(`fluxEmergenceSon2025DSGS`); the `4/3` could be taken out of `:μ[2]`, `:μ[3]` if
an exact Lax–Friedrichs form were ever wanted.

### The `B`-normalization caveat

The `B` floor in the kernel is `√ρ̄·c̄ = √0.14·√10 = 1.18`, the field strength of
a `β ≈ 1` plasma. This problem starts at `β_a = 1e-2`, so `B_a = 14.1` is 12×
that (120× at `B_a = √20000`). At `t = 0` the field is uniform, its spread is
zero, and the floor *is* its normalization, so the `B` equations would be read
against a scale well below their own magnitude.

What makes it harmless is that the numerator is zero at the same moment: the
injected field is the **same** `(0, B_a, 0)` as the ambient field, so the nozzle
puts no jump in `B` anywhere, and the induction residual only grows once the beam
starts shearing the field — by which time `By` has a spread of its own order,
which wins the `max`. Worth knowing if the first handful of steps look odd, but
do **not** raise `:dsgs_rel` to "fix" it: that key scales *every* floor and would
desensitize `ρ`, `ρv` and `E` at the same time.

---

## 6. What this case is actually up against

In the beam,

```
ρE = p/(γ-1) + ½ρ|v|² + ½|B|²  =  2.5 + 448 000 + 100  =  448 102.5
```

so the gas pressure is **5.58 parts per million** of the total energy, and

```
p = 0.4·(ρE − 448 100)
```

A **relative** error of `5.6e-6` in `ρE` wipes the pressure out entirely; so does
`2.8e-6` in the momentum `ρv`, since `δ(½(ρv)²/ρ)/KE = 2·δ(ρv)/(ρv)`. That is the
entire difficulty of this benchmark, and it is why:

* Wu & Shu state that with the PP limiter turned off "the simulation will break
  down after several time steps due to nonphysical numerical solutions";
* the OEDG paper says that disabling *any* of its locally-divergence-free
  oscillation-eliminating procedure, its PP limiter, or its upwind discrete
  Godunov–Powell source "the code immediately fails in these challenging tests".

**Jexpresso has no positivity machinery for MHD.** There *is* a node-wise
realizability repair, [`src/kernel/positivity/`](../../../src/kernel/positivity/README.md),
which the Mach-7 CompEuler decks (`rampCaoEtAl2021_M7`, `shock_circle_M7`,
`ffs_step`) run with `:lpositivity => true` — but it is scoped to
`neqs == nsd + 2` exactly and `positivity_validate` **errors** on this
nine-field state rather than guessing, because the magnetic energy in
`p = (γ-1)(ρE − KE − ½|B|² − ½ψ²)` is outside what it knows. Extending it to
the GLM-MHD state (the same momentum rescale θ = √((ρE − e_min − ½|B|² − ½ψ²)/ke),
which still conserves total energy exactly) is the obvious next step for this
case, and it would also give the counted, audited engagement report that is the
whole point of that module.

Until then all there is is a `p_floor_mhd = 1e-6` inside the flux evaluation
(`user_flux.jl`), so that a momentarily inadmissible node still produces a
*finite* flux instead of poisoning the whole RHS. It is not a positivity
guarantee, it is not conservative where it fires, and — unlike the repair — it
does not count itself. **If it fires, the run is already in trouble
and the snapshot is a diagnostic, not a result.** How to tell:

* `log10p` in the VTK output saturates at `-300` where `p ≤ 0` — look for it;
* `soundSpeed.jl` prints the CFL diagnostic each output; note that it uses
  `PhysConst.γ` and does **not** subtract the magnetic energy, so the `p` it
  reports inside the beam is ≈ 41 rather than 1. It is a monitor, not a
  measurement, on an MHD case;
* `JEXPRESSO_AJ_PFLOOR=0` turns the floor off to see the unshielded behaviour.

So: this deck is a faithful setup of the benchmark and an honest test of whether
residual-based viscosity alone can carry it. It may not reach `t = 2e-3` at the
default `B_a = √200` on the first try. That is information, not a bug in the
deck — and the ladder below is how to find out *where* it stops being able to.

### A ramp, if the full problem will not run

Climb it and record where it breaks: the jump that breaks is the one that says
what the model is missing.

Each rung keeps the beam and ambient states, the field and the geometry exactly
as the paper has them and changes only the injection speed (rungs 1–3) or the
field strength (rungs 5–6). `TEND` and `DT` are given per rung so that every run
takes the jet head to the same place, `y ≈ 1.2`, at the same Courant number
(`≈ 0.09` on the default mesh) — otherwise a rung that "works" may simply have
integrated a shorter time or with a smaller step.

| rung | `u_jet` | `B_a²` | `p/(γ-1)` ÷ `ρE` | `c_h` | `v_head` | `TEND` | `DT` |
|---|---|---|---|---|---|---|---|
| **B** | 800 | **2** | `5.6e-6` | `802` | `608` | `2e-3` | `5e-7` |
| 1 | 20  | 200   | `6.5e-3` | `37.9` | `15.2` | `8e-2` | `1e-5` |
| 2 | 80  | 200   | `5.5e-4` | `92.0` | `60.8` | `2e-2` | `4e-6` |
| 3 | 250 | 200   | `5.7e-5` | `262`  | `190`  | `6.3e-3` | `1.5e-6` |
| **4** | **800** | **200** | **`5.6e-6`** | **`812`** | **`608`** | **`2e-3`** | **`5e-7`** |
| 5 | 800 | 2000  | `5.6e-6` | `838` | `608` | `2e-3` | `5e-7` |
| 6 | 800 | 20000 | `5.5e-6` | `920` | `608` | `2e-3` | `5e-7` |

On a cluster, `auxiliary/wulver/submit_astrojet.sh` takes the rung as one word —
`sbatch --export=ALL,AJ_RUNG=B ...` — because rungs 1–3 need `TEND` and `DT` moved
*together* with `u_jet`, and setting three variables by hand and getting one wrong
produces a run that looks like a result.

**Rung B is the one to run first if the failures are in the ambient gas**, and §11
explains why: it keeps the Mach 800 beam exactly but takes `β_a` from `1e-2` to
`1`, which lifts the *ambient* thermal margin from 2.4 % to 71 %. If the run
survives rung B but not rung 4, the difficulty is the **field**, not the beam; if
it dies on rung B too, it is the beam and the shock.

Rung 4 is the deck's default (the paper's case (i)); rungs 5 and 6 are the
paper's cases (ii) and (iii). So, for example, rung 2 is

```
JEXPRESSO_AJ_UJET=80 JEXPRESSO_AJ_TEND=2e-2 JEXPRESSO_AJ_DT=4e-6 \
    julia --project=. src/Jexpresso.jl MHD astroJetWuShu2018
```

Note that `c_h` on rung 1 is set by the **ambient** fast speed (`37.9`) and not
by the beam (`20 + 12 = 32`), which is why it does not keep scaling down with
`u_jet`.

And if the default rung stops short, in order of what to try:

1. `JEXPRESSO_AJ_SMOOTH=0.005` — smooth the nozzle lip over one element. The lip
   is a corner singularity in the boundary datum sitting in the most
   positivity-critical fluid in the domain, and it is the least physical part of
   the setup.
2. `JEXPRESSO_AJ_CMIN=0.005` — a background floor, if the failure looks like a
   node-to-node mode rather than a shock (check `mu_dsgs_*` in the output: a
   checkerboard `ν` says mode, a ridge says shock).
3. `JEXPRESSO_AJ_MESH=100x150` — but not for the reason one might hope. The model
   is *scale-invariant* at a shock: `ν_res = C_R Δ²·R` with `R ∼ λ·δq/Δ` gives
   `ν_res ∼ C_R Δ λ δq`, the same `Δ` scaling as the cap `C_max Δ λ`, so
   refining does *not* buy a better dissipation-to-jump ratio. What it does buy
   is a smaller *absolute* overshoot per element at a given physical feature,
   and it is the absolute size of the wiggle in `ρE` against 2.5 that decides
   positivity here.
4. `JEXPRESSO_AJ_CMAX=1.0` — let the cap go above first-order-upwind strength.
   This is the point at which the run stops being the parameter-free method.
5. `JEXPRESSO_AJ_NODAL=1 JEXPRESSO_AJ_SENSOR=legacy` — `ν` per node instead of
   per element, Dao & Nazarov's own form (their eq. 4.10); a continuous
   coefficient has no jump in the diffusive flux at element interfaces. **Both
   variables together**: the nodal kernel reads the *assembled* residual, which
   the lumped-LGL assembly cancels on an under-resolved solution exactly as it
   does on a resolved one (`DSGS.md` §1.2), so nodal + `residual` is a blind
   sensor — measured on `CompEuler/shock_circle_M7` as a 17x regression (231
   steps against 3969). The deck warns if only one of the two is set.

Diagnostics worth turning on while doing this:
`JEXPRESSO_DSGS_DEBUG=1` prints the per-equation maximum of the normalized
residual, `ν_max`, the cap, and whether the argmax node is on an element edge;
`JEXPRESSO_DSGS_RSPLIT=1` dumps, per element, the time term and the space term
of the residual separately, which is what tells a real under-resolution from a
sensor reading something the other half of the residual is not.

---

## 7. What to look at

The papers plot **schlieren images of `log₁₀ρ` and `log₁₀p`** — Wu & Shu at
`t = 2e-3`, the OEDG paper at `t = 1e-3, 1.5e-3, 2e-3` (its Fig. 14). This deck
writes what is needed for that directly:

* `log10rho`, `log10p` as output fields (`user_primitives.jl`);
* `:lschlieren => true` adds `schlieren_grad_rho` (`|∇ρ|`, quantitative) and
  `schlieren` (`exp(-20|∇ρ|/max|∇ρ|)`, the picture — colour it with a
  **reversed** greyscale in ParaView);
* `beta` and `Mach`, plus the usual primitives and `mu_dsgs_*`.

Structures to look for, in the order they say something:

1. the **Mach shock at the jet head** and the bow shock ahead of it;
2. the **beam/cocoon interface**, with its shear-driven ripples — this is the
   feature that separates a scheme that resolves the jet from one that smears it;
3. the **cocoon**: the low-density, shocked ambient gas that has been swept
   sideways and back;
4. the field: at `β_a = 1e-2` the `y`-aligned field is stiff enough to collimate
   the beam, and the difference between the three magnetizations is mostly what
   the cocoon looks like. Wu & Shu: "the flow structures in different magnetized
   cases are very different."

And the two arithmetic checks from §3 and §6: the head should be near
`y ≈ 0.6·u_jet·t`, and `log10p` should not be anywhere near `-300`.

---

## 8. Files

| file | |
|---|---|
| `user_inputs.jl` | the deck; all the environment overrides are listed at its top |
| `user_flux.jl` | GLM-MHD flux, `γ_mhd = 1.4`, `c_h`, the flux pressure floor, **and every constant of the problem** |
| `initialize.jl` | the uniform magnetized ambient medium; computes and prints `c_h`, the Mach number and the `p/ρE` ratio |
| `user_bc.jl` | nozzle inflow / outflow, the lip treatment |
| `user_source.jl` | Dedner GLM `ψ` damping, `c_r = 0.18`. Nothing else: no gravity, no resistivity |
| `user_primitives.jl` | the conserved-form DynSGS primitives and the 14 output fields |
| `AJ.geo` | the mesh definition, and the `0.05·nx ∈ ℤ` constraint |
| `AJ_40x60.msh`, `AJ_100x150.msh` | the two shipped meshes |
| `../../../tools/astro_jet_mesh.py` | regenerates them without gmsh |
| `../../../auxiliary/wulver/submit_astrojet.sh` | SLURM submit script; `AJ_RUNG` selects a rung of §6 |

The equation set itself is documented in
[`../orszagTangBormanis2024/EQUATIONS.md`](../orszagTangBormanis2024/EQUATIONS.md) —
it is the same nine-field ideal GLM-MHD system, at a different `γ`.

---

## 9. What has been verified, and what has not

The case has **never been integrated in time** — no Jexpresso run exists. What
was checked, with Julia 1.11.9:

| checked | how |
|---|---|
| all six `user_*.jl` / `initialize.jl` parse | `Meta.parseall`, walked for `:error`/`:incomplete` |
| all six **execute**, in `run.jl`'s include order | stub stand-ins for the Jexpresso types; every hook called |
| `user_inputs()` returns the intended 49 keys | `Δt`, `tend`, `nop`, the DynSGS block, the mesh path, the 21 output times incl. `1e-3, 1.5e-3, 2e-3` |
| the deck's `:dsgs_gamma` equals `γ_mhd` | asserted against the constant in `user_flux.jl` |
| `p = 1` in **both** states; beam sound speed exactly `1`; Mach exactly `800` | from the conserved 9-tuples through `pressure_mhd` |
| `c_f` = 37.9285 / 11.9940, `c_h` = 811.9940, `β_a` = 1e-2, beam `ρE` = 448102.5, `p/(γ-1)/ρE` = 5.579e-6 | `aj_wave_speed`, `initialize` |
| `initialize` fills all 9 slots + the pressure slot, sets `qe`, sets `c_h` | on a stub mesh |
| a wrong-domain mesh **warns** | `@test_logs (:warn,)` |
| `user_flux!` values | `G = [1120, 0, 895901, 3.58403e8, 0, 0, 0, 0, 11483.3]`, checked term by term |
| the flux stays **finite** on a `p < 0` state | half the beam energy, `p_floor_mhd` active |
| `user_source!` damps only `ψ`, at `c_h/c_r = 4511` | |
| `user_primitives!` is the conserved form, spare slots untouched | |
| `user_uout!` writes all 14 fields; `log10p` saturates at `-300` on a bad node | |
| `user_bc_dirichlet!` prescribes all 9 slots for `\|x\| ≤ 0.05` **including the lip nodes at the mesh's own ±0.05±4e-17**, and nothing at all on `top`/`left`/`right` or on the bottom outside the nozzle | against the `4325789.0` sentinel |
| both `PERT()` paths error rather than silently misbehave | |
| lip smoothing keeps `p = 1` at `φ = 1` **and** `φ = ½` | the reason it blends primitives, not conserved variables |
| both meshes load through **GridapGmsh**, the real reader | quad cells, extents, `\|Ω\| = 1.5`, tags `bottom/right/top/left/domain` |
| the boundary edge set is **exactly** the geometric boundary (200 / 500 edges) | Jexpresso's own `get_boundary_faces` + label loop, replayed verbatim |
| `41` / `101` nodes on `y = 0`, of which `5` / `11` in the nozzle | `= 2·(0.05·nx) + 1`, so the lip is a node |
| `gmsh -2 AJ.geo` reproduces the shipped `AJ_40x60.msh` | same cells, nodes, tags; node coordinates agree to `2.9e-12` |

**Not** checked, and the honest gaps:

1. **No time integration.** The full dependency set (≈100 packages, MUMPS,
   Pardiso, ONNXRunTime, P4est) was not installable here, so nothing exercised
   `rhs!`, the DSGS kernel, the RK stages or the VTK writer on this case.
2. **This is the first 2D MHD case to pair `:dsgs_conserved => true` with
   genuinely conserved primitives.** `brioWu1d` does it in 1D;
   `fluxEmergenceSon2025DSGS` does it in 2D but with `:dsgs_ref_weight` and the
   split energy on top. (`MHD/smoothVortex` sets `:dsgs_conserved => true` while
   its `user_primitives!` returns the *physical* set `ρ, u, v, T, w, B, ψ` — so
   it is not a precedent for this pairing.) The code path was traced
   branch by branch and every flag it depends on is guarded, but it has not
   been run.
3. **Whether DynSGS alone carries the Mach 800 beam is the open question** —
   that is §6, and it is what the run is for.

---

## 10. Run 1: what happened, and what changed because of it

The first run aborted at **t = 5.45e-5 — step 109 of 4000** — with

```
non-finite solution at t = 5.449999986240073e-5 — 6237 of 6237 local entries (100.0%)
```

on **every** rank, in the same step.

### Reading the signature

100 % of every rank at once is not a local blow-up; a local one shows a handful
of nodes on one rank. It is a **global reduction carrying the damage**. With
`:dsgs_norms => "domain"` every RHS call does `MPI.Allreduce!(avg, MPI.SUM)` and
`MPI.Allreduce!(denom, MPI.MAX)` over the whole domain
(`SGS.jl`, passes 1 and 2). One `Inf` anywhere makes `⟨q⟩` `NaN`, which makes
every `denom` `NaN`, which makes the normalized residual `NaN`, which makes `ν`
`NaN` on **every element of every rank** — and the next stage multiplies that
into the whole field. So:

> the failure was local and became global in one step, and the abort message is
> right that the field can no longer locate it.

Where did the `Inf` come from? The flux divided by `ρ` unguarded: `u = ρu/ρ`. A
single node with `ρ → 0` gives `±Inf`, and `Inf − Inf` in the flux gives `NaN`.

### What changed

1. **`src/kernel/positivity/` now has a GLM-MHD branch** (`positivity_limit_mhd!`).
   The abort message says to "read the positivity report's GLOBAL first repair" —
   and there was none, because the repair errored on this nine-field state. There
   is now: the same energy-conserving momentum rescale, with `½|B|² + ½ψ²` charged
   against `ρE` and `B`/`ψ` never rescaled so `∇·B` is untouched. It is **on** in
   this deck, floors `ρ_min = 1.4e-7` and `p_min = 1e-6` (1e-6 of the ambient
   values), reporting every 50 RHS calls = 10 steps.
2. **The output cadence is front-loaded.** The first frame after the initial
   condition was at `t = 1e-4`; the run died at `5.45e-5`. There was no snapshot
   of the failure at all. There are now ~14 frames before that time, the first at
   `t = 2e-6` (4 steps).
3. **The flux guards its divisions by `ρ`** (`ρ_floor_mhd = 1e-14`), so a local
   defect stays local and locatable instead of being laundered into a global
   `NaN` by the next reduction. Unreachable while `:lpositivity` is on; it is the
   belt to that braces.

### What to run next, and what each outcome means

```
julia --project=. src/Jexpresso.jl MHD astroJetWuShu2018
```

The `POSITIVITY REPAIR ENGAGED` line now answers the question the first run could
not. Read it in this order:

| what the report says | what it means | what to do |
|---|---|---|
| never engages, run still dies | the failure is **not** a realizability one — look at `ν` in the frames before it | `JEXPRESSO_DSGS_DEBUG=1` |
| a few node-visits, run continues | the repair is doing its job; look at the beam profile | keep going |
| **GLOBAL first repair at the nozzle lip** (`x = ±0.05`, `y` within an LGL interval of 0) | the boundary datum's corner is the source, as §3 predicted | `JEXPRESSO_AJ_SMOOTH=0.005` |
| first repair at the jet head or the beam/cocoon interface | genuine under-resolution of the shock | more dissipation — see below |
| `energy-RAISED` dominating | branch 2b: the field energy alone exceeds `ρE − e_min`. At `β_a = 1e-2` this can be the field, not a broken momentum — but a growing count means the state is badly broken, not marginally | lower `B_a` (rung 1–3 of §6) to separate the two |
| engagement growing without bound | no limiter would have saved this; the field is globally wrong | drop to a lower rung of §6 |

Run it with `JEXPRESSO_DSGS_DEBUG=1` the first time. It prints, per call, the
per-equation maximum normalized residual, `ν_max`, **the cap**, and whether the
argmax node sits on an element edge. That last pair settles which knob matters,
and it is not the obvious one:

> `ν = min(C_max·Δ·λ, C_R·Δ²·R)`. Whichever term is smaller is the only one
> that matters: if `ν_max` is **below** the cap then `C_R` is binding and raising
> `C_max` does nothing, and if `ν_max` **is** the cap then the reverse holds.
>
> **Measured, run 2: `ν_max = 2.1013 = C_max·Δ·λ` exactly at `λ = 840.5`
> (`C_max = 0.5`, `Δ = h/5 = 5e-3`). DynSGS is pinned at its cap.** So `C_R` is
> *not* the lever — raising it changes nothing — and `C_max` is. (An earlier
> estimate here said the opposite; the measurement overturned it. This is what the
> debug line is for.)

Also worth knowing while debugging: `JEXPRESSO_AJ_NORMS=element` makes the DynSGS
normalization element-local, which removes the two Allreduce and therefore removes
the mechanism that turned one bad node into a global `NaN`. It changes the model
(the residual is then measured against each element's own spread, not the
domain's), so it is a **diagnostic**, not a fix — but it keeps a failure local and
visible, which is what a first bisection needs.

---

## 11. Run 2: the diagnosis, and the thermal margin nobody warns you about

With the GLM-MHD repair and the front-loaded frames, run 2 got **7× further**:
`t = 3.735e-4`, step 747 of 4000 (18.7 % of the target), and fast. It still
aborted the same way — 100 % non-finite on every rank — but this time the report
said where and why.

```
repaired 103193 node-visits in 1350 RHS calls
  [ρ-floor 191, momentum-scaled 30078, energy-RAISED 72924]
  injected: mass 0.173, energy 2.2158e6
  GLOBAL min ρ -0.005275, GLOBAL min p -6.5139e7
GLOBAL first repair at (x, y) = (-0.075, 0.0)  on RHS call 3, rank 24
```

Three findings, all quantitative.

### (a) The nozzle lip. Confirmed, and it is the clamped/free interface

`(x, y) = (-0.075, 0)` is on the **bottom boundary, one element outside the lip**,
and `RHS call 3` is inside the **first time step**. The fastest signal travels
`c_h·Δt = 4e-4` in a step — 1.6 % of an element — so **nothing propagated there.
The boundary condition put it there**, exactly as a 446 m/s shear held open by a
Dirichlet condition put `rampCaoEtAl2021_M7`'s first repair one LGL interval below
its top boundary.

The mechanism is not the lip coordinate, it is that a top-hat datum imposed
*strongly* pins nodes at `ρE = 4.48e5` next to a **free** node the scheme wants to
leave at `102.5` — a 4400× jump inside one spectral element, `∂(ρE)/∂x ≈ 1.8e7`.
The Gibbs response is the size of the jump, so a node whose `ρE` is 102.5 takes an
excursion of `O(1e5)`.

**Fixed** (§3): the blend is now a compactly supported smootherstep centred on the
lip, with the Dirichlet patch widened to `|x| ≤ x0 + s` so that `φ` reaches
**exactly** zero at the outermost clamped node — which therefore holds precisely
the ambient state, the same thing its free neighbour holds. There is no
clamped/free jump left. The previous `tanh` version could not do this: it only
decays, and 0.02 % of the beam's `ρE` is 100, which is the *entire* ambient `ρE`.
`s` defaults to one element, from `mesh.Δelem_s` — the globally MPI-reduced
smallest element, **not** `mesh.nelem`, which is rank-local and on this 64-rank run
would have given every rank a different width, all ~8× too wide.

Note the patch edge now sits at `|x| = 0.075` — the very node that failed.

### (b) DynSGS is saturated. `C_R` is not the lever; `C_max` is

`max ν = 2.1013`, and `C_max·Δ·λ = 0.5 × 5e-3 × 840.5 = 2.1013`. **Exactly the
cap.** So the residual sensor is asking for more dissipation than its own
first-order-upwind bound allows, and raising `C_R` changes nothing at all. §10's
earlier estimate said the opposite and was wrong.

The viscous CFL was `0.045–0.059`, so there is real headroom: it scales linearly
with `C_max`, so **`C_max` can go to ≈ 3.5 before the viscous limit binds at the
current `Δt`**. `JEXPRESSO_AJ_CMAX=1.0` (twice first-order-upwind) is the honest
first step, and it is the point at which this stops being the parameter-free
method — say so in any write-up.

### (c) The ambient medium has a 2.4 % thermal margin — this is the real difficulty

This is the one that is not in §6, and it changes how to think about the case.
In the **undisturbed ambient gas**:

```
ρE = p/(γ-1) + ½|B|² = 2.5 + 100 = 102.5     at β_a = 1e-2
                       ↑        ↑
                    2.44%    97.56%
```

So an oscillation of **2.4 % in `ρE`, in gas that is doing nothing**, gives `p < 0`.
Not 5.6 ppm as in the beam — but the beam is 0.1 wide and the ambient is the whole
domain. And it gets worse with the paper's other two configurations:

| | `β_a` | ambient `ρE` | thermal margin |
|---|---|---|---|
| (i)   | `1e-2` | 102.5   | **2.44 %** |
| (ii)  | `1e-3` | 1002.5  | 0.249 % |
| (iii) | `1e-4` | 10002.5 | 0.025 % |

That is why **branch 2b dominated 72924 to 30078**: 2b fires when
`ρE − ½|B|² − ½ψ² ≤ e_min`, i.e. when the total energy dips below the magnetic
energy — and in this ambient that is only 2.4 % down. The injected energy
averages `2.2158e6 / 72924 = 30.4` per engagement, so `ρE` had collapsed to about
`70` from `102.5`: a −32 % excursion, not a ripple.

**And 2.2158e6 is 14412× the whole domain's initial total energy** (`102.5 × 1.5 =
153.75`). The repair was not repairing, it was writing the solution. Per its own
README that means the answer is meaningless — which is exactly what the audit
trail is for, and why it is worth having even on a run that fails.

`min p = -6.5e7` says the same thing from the other side: `p = 0.4(ρE − ke − 100)`,
so `ke` exceeded `ρE` by `1.6e8`, i.e. `|v| ≈ 4.8e4` — sixty times the beam speed.
There is a genuine instability, not only ambient ripple.

### What to run next, in order

1. **Just rerun.** The lip fix addresses the *first* cause, at step 1, and nothing
   downstream can be judged until it is gone.
   ```
   JEXPRESSO_DSGS_DEBUG=1 julia --project=. src/Jexpresso.jl MHD astroJetWuShu2018
   ```
   Then read the first-repair coordinate again. If it has moved off the boundary
   and into the jet head or the beam/cocoon interface, (a) is fixed and the
   remaining problem is shock resolution.
2. **`JEXPRESSO_AJ_BA2=2`** (rung B of §6). Mach 800 beam, `β_a = 1`, ambient
   margin 71 %. This separates "the field is the problem" from "the beam is the
   problem" in one run, and finding (c) says it is the single most informative
   knob on this case.
3. **`JEXPRESSO_AJ_CMAX=1.0`**, now that the cap is known to be binding.
4. **`JEXPRESSO_AJ_MESH=100x150`**, which reduces the absolute overshoot at a
   given feature even though it does not change the dissipation-to-jump ratio.

What **not** to do: raise `C_R` (the cap, not the sensor, is binding), or lower
`Δt` alone (CFL 0.09 and viscous 0.06 were both comfortable — `Δt` was not the
constraint).

---

## 12. Run 3: the lip fix worked, and the next discontinuity is in *time*

```
repaired 103384 node-visits in 1800 RHS calls
  [ρ-floor 0, momentum-scaled 21718, energy-RAISED 81666]
  injected: mass 0.0, energy 903446.6
  GLOBAL min ρ 0.11003742612074623, GLOBAL min p -4428.894
GLOBAL first repair at (x, y) = (-0.04999999999999999, 0.025)  on RHS call 3, rank 24
```

### The lip fix worked

| | run 2 | run 3 |
|---|---|---|
| `ρ`-floor engagements | 191 | **0** |
| `GLOBAL min ρ` | **−0.005275** | **+0.110037** (ambient 0.14, so −21 %: healthy) |
| `GLOBAL min p` | −6.514e7 | **−4428.9** (four orders better) |
| injected energy | 2.2158e6 | 9.0345e5 |
| first repair | `(−0.075, 0.000)` — **on** the boundary | `(−0.050, 0.025)` — **off** it |

Density never left the realizable set at all, and the pressure excursion fell by
four orders of magnitude. The clamped/free jump in `x` was real and it is gone.

### The new first repair is the impulsive start

`(x, y) = (−0.05, 0.025)`: the lip abscissa, **one element above the boundary**,
and again on **RHS call 3** — the first time step. A signal travels `c_h·Δt = 4e-4`
in a step, 0.09 of an LGL gap, so once again nothing propagated there.

What is there on step 1 is this: **at `t = 0` the whole domain including `y = 0` is
the ambient medium, and at the first RHS evaluation the boundary condition clamps
the nozzle to the beam.** That is the same 4400× jump in `ρE`, now across the first
LGL gap in `y` (`4.3e-3`), and the lip blend cannot touch it because it only shapes
the datum in `x`. Fixing a discontinuity in space left the one in time.

**Fixed:** the beam is ramped on with the same C² smootherstep,

```
φ = φ_x(x) · φ_t(t),    φ_t(t) = smootherstep(t/τ),  1 for t ≥ τ
```

one multiplicative factor through the same primitive blend. `τ = 2h/u_jet`, which
is **125 time steps on either shipped mesh** because `Δt` scales with `h`, and
which spreads the beam front over `u_jet·τ = 2h` — **exactly the two elements the
lip profile spans in `x`**. The two are matched on purpose: the datum is no steeper
in `y` than in `x`. Measured: the worst change in `ρE` per step at the nozzle
centre drops from 448000 to 9830, 46× smaller, and at `t = 0` the imposed state is
the ambient state *bit for bit* at every clamped node.

Cost: the beam is at full strength from `t = τ = 3.1 %` of `tend`, so the jet head
is delayed by about `τ/2`, 1.6 % of the domain height at the final time.

### What has not changed, and is now the main suspect

`(c)` from §11 is untouched, and it is doing most of the damage:

* **2b still dominates 81666 to 21718** and still injects `9.03e5` — **5876× the
  domain's initial total energy** (`153.75`). The answer is still meaningless.
* `min ρ` is healthy and `min p` is `−4429`, so this is no longer a density
  collapse. It is the energy budget: `ρE` dipping below `½|B|² = 100` in a medium
  whose `ρE` is only `102.5`.
* 287 repairs per step, 0.75 % of the nodes — spread through the domain, not
  concentrated in the beam.

Every one of those is the signature of the **2.4 % ambient thermal margin**, not of
the inlet. So:

### What to run next

1. **Rerun.** The time ramp removes the last step-1 discontinuity. Read the first
   repair coordinate again: if it moves off `RHS call 3` to a later call, both
   inlet problems are gone and whatever remains is physics, not a boundary datum.
   ```
   JEXPRESSO_DSGS_DEBUG=1 julia --project=. src/Jexpresso.jl MHD astroJetWuShu2018
   ```
2. **`JEXPRESSO_AJ_BA2=2`** — rung B of §6, and now the top suspect rather than a
   curiosity. Mach 800 beam unchanged, `β_a = 1`, ambient thermal margin 71 %
   instead of 2.4 %. If this survives to `t = 2e-3`, the difficulty is the field's
   energy budget and the paper's `β_a = 1e-2` needs an invariant-domain-preserving
   scheme, not a repair. If it dies too, the beam and the shock are the problem and
   the field is a bystander.
3. **`JEXPRESSO_AJ_CMAX=1.0`** — still valid from §11(b): `ν` is pinned at the cap,
   and the viscous CFL of 0.055 leaves room to ≈3.5.
4. **`JEXPRESSO_AJ_TRAMP`** = `4h/u_jet` (twice the default) if the first repair is
   still on an early call but has moved off `y = h`.
