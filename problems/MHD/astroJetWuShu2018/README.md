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
| 1 | 20  | 200   | `6.5e-3` | `37.9` | `15.2` | `8e-2` | `1e-5` |
| 2 | 80  | 200   | `5.5e-4` | `92.0` | `60.8` | `2e-2` | `4e-6` |
| 3 | 250 | 200   | `5.7e-5` | `262`  | `190`  | `6.3e-3` | `1.5e-6` |
| **4** | **800** | **200** | **`5.6e-6`** | **`812`** | **`608`** | **`2e-3`** | **`5e-7`** |
| 5 | 800 | 2000  | `5.6e-6` | `838` | `608` | `2e-3` | `5e-7` |
| 6 | 800 | 20000 | `5.5e-6` | `920` | `608` | `2e-3` | `5e-7` |

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

The equation set itself is documented in
[`../orszagTangBormanis2024/EQUATIONS.md`](../orszagTangBormanis2024/EQUATIONS.md) —
it is the same nine-field ideal GLM-MHD system, at a different `γ`.
