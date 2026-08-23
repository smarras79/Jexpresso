# HEVI in Jexpresso

**Horizontally Explicit / Vertically Implicit time integration for the
compressible Euler equations.**

Source: `src/kernel/solvers/hevi/` — `columns.jl`, `operator.jl`,
`factorize.jl`, `ark.jl`, `hevi.jl`, `cfl_diagnostics.jl` (included in that
order).

Worked example throughout: `problems/CompEuler/rtb_hevi`, on
`hevi_10x1x50.msh`.

---

## 1. The problem HEVI solves

An explicit scheme integrating the compressible Euler equations is limited by
the fastest wave on the grid: sound, at `c ≈ 347 m/s`. The stability
constraint is set by the *sum* of the per-direction rates, because the
eigenvalues of a tensor-product operator add along its diagonal:

```
rate_explicit  =  c/h_ξ + c/h_η + c/h_ζ  +  ν/h_ξ² + ν/h_η² + ν/h_ζ²
Δt_max         =  R / rate_explicit
```

where `R` is the scheme's imaginary stability radius (3.95 for
`CarpenterKennedy2N54`).

On an atmospheric mesh the vertical spacing is far smaller than the
horizontal, so `c/h_ζ` dominates that sum and the whole 3D problem is stepped
at a rate set by one direction. On `hevi_10x1x50.msh`:

| direction | LGL spacing | acoustic limit |
|-----------|-------------|----------------|
| ξ (x)     | 172.7 m     | 0.4976 s       |
| η (y)     | 863.4 m     | 2.488 s        |
| **ζ (z)** | **34.53 m** | **0.09953 s**  |

The vertical is 5× finer and costs 5× the step. HEVI removes **that one term**
from the explicit budget by treating it implicitly, and nothing else.

> The gain is therefore bounded by the grid's acoustic anisotropy. On an
> isotropic mesh HEVI buys nothing and is strictly worse, because the implicit
> tableau has a smaller explicit stability radius than a well-tuned explicit
> RK. Run with `:lcfl_report => true` first: if the report says SGS diffusion,
> not vertical acoustics, is the binding term, HEVI as built here will not
> help.

---

## 2. The splitting

### 2.1 Governing equations

In perturbation-capable conservative variables `u = (ρ, ρu, ρv, ρw, ρθ)`:

```
∂ρ/∂t    + ∇·(ρ𝐮)          = 0
∂(ρ𝐮)/∂t + ∇·(ρ𝐮⊗𝐮) + ∇p   = -ρg ẑ
∂(ρθ)/∂t + ∇·(ρθ𝐮)         = 0
p = C₀ (ρθ)^γ
```

### 2.2 What is implicit

Only the terms that carry the **vertical** sound wave, linearised about the
reference state `qe`. Writing `δ· = · - qe`:

```
∂(δρ)/∂t   = -∂/∂z [ δ(ρw) ]
∂(δρw)/∂t  = -∂/∂z [ β δ(ρθ) ]  -  g δρ
∂(δρθ)/∂t  = -∂/∂z [ θ̄ δ(ρw) ]
```

with `β = ∂p/∂(ρθ) = γp/(ρθ) = c²/θ` evaluated on the reference state.
Call this linear operator `A`, so

```
f_imp(u) = A (u - qe)
```

Three of the five equations are involved — the implicit variable set is
`[1, 4, 5]` = `(ρ, ρw, ρθ)`, `nimp = 3`. Horizontal momentum stays fully
explicit.

Vertical **advection** deliberately stays explicit: on these grids `|w|/h_ζ`
is ~1/800 of `c/h_ζ`, so making it implicit would buy nothing and would force
a refactorisation every step (the matrix would depend on the solution).

### 2.3 Two properties that come free

**Affine-free.** `A` acts on the deviation `u - qe`, so `f_imp(qe) = 0`
*exactly*. The reference state is an exact steady state of the implicit part.
Whatever discrete hydrostatic imbalance the code already has stays entirely
inside the explicit part, untouched — the split can neither create nor destroy
balance.

**No lost or double-counted physics.** The explicit part is not written as a
second RHS. It is computed by subtraction:

```
f_exp(u) = rhs!(u) - f_imp(u)
```

That costs one extra (cheap) application of `A` per stage and buys a
guarantee: every source term, boundary condition, SGS closure, wall model and
sponge remains exactly what the explicit code computes. Nothing has to be kept
in sync by hand.

---

## 3. The IMEX Additive Runge–Kutta scheme

### 3.1 Stage equations

An ARK method carries **two** Butcher tableaux — `(aE, bE, cE)` applied to
`f_exp` and `(aI, bI, cI)` applied to `f_imp` — sharing stage values:

```
U_i = uⁿ + Δt Σ_{j<i} aE_ij f_exp(U_j)  +  Δt Σ_{j≤i} aI_ij f_imp(U_j)

uⁿ⁺¹ = uⁿ + Δt Σ_i bE_i f_exp(U_i)  +  Δt Σ_i bI_i f_imp(U_i)
```

`aI` is lower triangular *including* the diagonal (`aI_ii = γ`), so each stage
requires one linear solve:

```
(I - γΔt A)(U_i - qe) = tmp_i - qe
```

with `tmp_i` the accumulated explicit/implicit history. `aE` is strictly lower
triangular, so `f_exp` is never solved for.

### 3.2 The tableaux (ARS family)

Default is `:ARS232`. With `γ = (2-√2)/2`, `δ = -2√2/3`:

```
        0   0   0                    0     0    0
aE =    γ   0   0            aI =     0     γ    0
        δ  1-δ  0                     0   1-γ    γ

bE = bI = [0, 1-γ, γ]
```

Three stages, second order, **stiffly accurate** (`bI = aI[s,:]`, `bE =
aE[s,:]`), so the last stage *is* the answer and its RHS evaluation is reused
as the next step's FSAL. Net cost: **3 `rhs!` evaluations and 2 column solves
per step**, against `CarpenterKennedy2N54`'s 5 `rhs!`.

Available: `:ARS111 :ARS121 :ARS122 :ARS222 :ARS232 :ARS343 :ARS443`.

---

## 4. Joint stability — the part that matters

### 4.1 Why the usual figure of merit is wrong

A tableau is normally chosen by its **explicit imaginary radius**: how far up
the imaginary axis `|R_E(z)| ≤ 1` holds. For HEVI that is the wrong number,
and choosing on it is an active trap.

A HEVI acoustic mode with horizontal wavenumber `k_x` and vertical `k_z`
loads **both halves at once**:

```
z_E = i Δt c k_x      (explicit)
z_I = i Δt c k_z      (implicit)
```

Stability requires the **joint** amplification factor to stay bounded over the
whole rectangle, not just along the `z_I = 0` edge:

```
R(z_E, z_I) = 1 + (z_E bEᵀ + z_I bIᵀ) (I - z_E aE - z_I aI)⁻¹ 𝟙
```

`ark_joint_amplification(tab, zEmax, zImax)` samples `max|R|` over that
rectangle; `ark_joint_dt_max(tab, rate_exp, rate_imp)` inverts it for Δt.

### 4.2 The ARS343 trap, measured

On the anisotropic `rtb` mesh:

| tableau | RHS/step | explicit imag. radius | **joint Δt_max** |
|---------|----------|-----------------------|------------------|
| ARS232  | 3        | 1.732                 | **0.052 s**      |
| ARS443  | 4        | 1.570                 | 0.059 s          |
| ARS343  | 4        | **2.828** (best)      | **0.0004 s**     |

ARS343 has the *largest* explicit imaginary radius in the family and a joint
region ~130× smaller. This case shipped on ARS343 at Δt = 0.03 and amplified
0.79 %/s — slow enough to look like anything but a tableau problem: it
presented as a blow-up at a fixed *model* time, nearly independent of Δt, that
moved with the rank count.

> **A note on ARS343 and the literature.** Giraldo et al. (2023) report
> ARS(3,4,3) among their *fastest* integrators, which looks like a flat
> contradiction of the table above. Both are right. Mapping `|R|` pointwise
> shows ARS343 has a narrow unstable strip near the apex and is heavily damped
> beyond it:
>
> ```
> zE\zI      0.0     0.5     2.0     5.0    20.0
>  0.50   0.9999  1.0040  0.8918  0.5072  0.1066
>  1.70   0.8865  0.9912  1.2649  0.8605  0.2733
> ```
>
> That is the "wedge" shape of their §6, where the explicit region grows with
> the implicit step. They run at a stiffness ratio around 10 **with 4th-order
> hyperdiffusion applied at every stage**, which damps that strip away. This
> code runs at ratio 2–5 with a constant μ (h² damping, not h⁴), so the strip
> survives — and it did, at 0.79 %/s. `ark_joint_amplification` is a
> purely imaginary-axis analysis and is therefore **conservative by
> construction**: it is the undamped answer.

**Judge a tableau by `ark_joint_amplification`, never by the explicit radius.**
A check that samples only `z_I = 0` passes ARS343. `test/hevi/test_joint_stability.jl`
pins exactly that trap.

### 4.3 The startup guard

`build_hevi` computes `max|R|` at the deck's actual Δt on the actual mesh and
**refuses to start** if it exceeds 1, naming the neutral limit, a recommended
step (70 % of it) and the tableaux whose joint region does cover the mesh.
Override with `:hevi_verify => false`.

Measured on `hevi_10x1x50.msh` (unrefined):

```
joint IMEX stability: explicit half 3.7 1/s, implicit half 14.4 1/s
    Δt_max here = 0.3567 s     (recommended: 0.2497 s, i.e. 70%)
```

against an explicit acoustic limit of 0.0995 s — consistent with the observed
behaviour that explicit diverges at Δt = 0.25 while HEVI runs.

---

## 5. The column solve

### 5.1 Why it is cheap

`A` couples a node only to nodes **within the same vertical column**, and only
within `±(ngl-1)` levels of it: inside an element the ζ collocation derivative
reaches all `ngl` nodes of the line, and DSS links an element only to the ones
directly above and below.

So `(I - γΔt A)` is block-banded per column, and a direct banded LU costs
nothing next to one explicit RHS.

For `rtb_hevi` (`nop = 4` → `ngl = 5`, 50 vertical elements):

```
columns          : 205 global   ( = (10·4+1) × (1·4+1) = 41 × 5 )
levels/column    : 201          ( = 50·4 + 1 )
matrix order n   : 603          ( = 201 levels × nimp 3 )
bandwidth kl = ku: 14           ( = nimp·(ngl-1) + 2 = 3·4 + 2 )
storage          : 41.5 MB on the busiest rank
```

### 5.2 The matrix is extracted, not derived

The band is **not** derived on paper and transcribed. It is recovered by
*probing*: colour the levels modulo `S = 2(ngl-1)+1 = 9`. A probe that is 1 on
one colour and 0 elsewhere hits, for each row, exactly one column of the band
— the window a row can see is `S` wide, so it contains exactly one level of
any given colour. `S × nimp = 27` applications of `A` recover every entry
exactly.

Deriving them by hand would mean keeping the ζ collocation weights, metric
terms, DSS sum, inverse mass matrix and wall-flux condition in sync by hand.
A discrepancy would *not* break a convergence test — the split stays
consistent whatever the implicit matrix is, because `f_exp = rhs! - f_imp` —
it would just silently make the solve approximate and bring the vertical
acoustic CFL back with nothing in the output to say why. Probing removes that
failure mode: **the matrix is the operator, by construction.**

### 5.3 One factorisation for the whole run

`γ` is fixed by the tableau and Δt is fixed, so `γΔt` is constant and the LU
is computed **once** at setup.

### 5.4 Under MPI

A column-local solve assumes each column lives on one rank. p4est partitions
by Morton order, so that holds only at tree-column granularity. `columns.jl`
makes it work under any partition:

* columns that are rank-local are solved **in place**, no communication;
* columns split across ranks are **gathered** onto an owner, solved, and
  **scattered** back.

The setup report prints how many columns are split. With `:lxy_partition =>
true` (columnar) it is typically 0 %; under p4est it can be 100 %, which is
correct but communication-heavy. HEVI therefore defaults to the columnar
partition, and `:linitial_refine => true` forces p4est and gives that up.

---

## 6. Worked example: `rtb_hevi`

Mesh `hevi_10x1x50.msh` — 10 × 1 × 50 elements, `nop = 4`, over
`x ∈ [-5000, 5000]`, `y ∈ [0, 5000]`, `z ∈ [0, 10000]` m.

```
 z = 10000 ┌───────┬───────┬───────┬ ... ┬───────┐   ▲
           │       │       │       │     │       │   │  50 elements
           ├───────┼───────┼───────┼ ... ┼───────┤   │  Δz = 200 m
           │       │       │       │     │       │   │  h_ζ(LGL) = 34.53 m
           ├───────┼───────┼───────┼ ... ┼───────┤   │
           │  ...  │  ...  │  ...  │     │  ...  │   │   ← one COLUMN =
           ├───────┼───────┼───────┼ ... ┼───────┤   │     one banded LU
           │       │       │       │     │       │   │     (n = 603)
 z = 0     └───────┴───────┴───────┴ ... ┴───────┘   ▼
          x=-5000                                x=+5000
           └────────── 10 elements, Δx = 1000 m ─────────┘
                       h_ξ(LGL) = 172.7 m

           anisotropy  min(h_ξ,h_η)/h_ζ = 5.0
```

Each of the 205 vertical lines of nodes is one independent tridiagonal-ish
banded system. Horizontally, nothing is coupled implicitly.

### Measured, 1 rank, this mesh

| quantity | explicit (CK2N54) | HEVI (ARS232) |
|---|---|---|
| stability limit on Δt | 0.0995 s (acoustic ζ) | **0.3567 s** (joint) |
| observed | diverges at Δt = 0.25 | runs at Δt = 0.25 |
| cost per step | 0.212 s | **0.156 s** |
| RHS evaluations/step | 5 | 3 (+2 column solves) |

HEVI is cheaper *per step* here and admits a larger step. Note that comparing
the two at the **same** Δt measures only per-step cost and says nothing about
the method — the gain lives in the step size.

### Setup self-check

At startup HEVI verifies the assembled band against the operator it is
supposed to represent, and the solve against the band:

```
self-check: band matches operator to 4.07e-16; solve residual 2.78e-15
```

---

## 7. Reading the CFL report

```
 │  Δt gain available, per scheme (rate-summed, before per-stage cost):
 │    HEVI  (vertical acoustics implicit)              4.82x
 │    HEVI + implicit vertical diffusion               5.24x
 │    HEVI + acoustic substepping (outer step)        58.06x
```

These are **Δt-ratio upper bounds**, computed by summing per-direction
stability rates and removing the terms each scheme would make implicit:

```
gain = rate_full / rate_scheme
```

They are *not* wall-clock speedups (the header says "before per-stage cost"),
and they are idealisations: the realised limit comes from the joint IMEX
region, which also pays for ARS232's smaller explicit radius. On this mesh the
first line reads 4.82× while the joint calculation gives 0.3567/0.0995 = 3.6×.

**Only the first line corresponds to an implemented scheme.** Vertical
diffusion is not in the implicit operator, and acoustic substepping is not
implemented.

---

## 7b. Linearisation: LHEVI-RS vs LHEVI-PS

The implicit operator's coefficients — `beta = dp/d(rho theta)` and `thetabar`
— have to be evaluated somewhere. Two choices, one deck switch:

```julia
:hevi_linearization => :RS,   # default: frozen at the reference state qe
:hevi_linearization => :PS,   # refreshed from the solution, then refactorised
:hevi_update_freq   => 5,     # steps between refreshes (:PS only)
```

**What `:PS` changes, and what it deliberately does not.** The operator's
stiffness-removing power lives entirely in its coefficients: `f_imp(u) =
A(u - qe)` has Jacobian `A` whatever is subtracted, so how well the split
removes the vertical acoustic term depends on `A` alone — the offset cancels in
`f_exp = rhs! - f_imp`. So `:PS` refreshes the coefficients and keeps measuring
the deviation from `qe`.

That is a deliberate departure from Giraldo et al., who linearise about the
previous solution outright. Keeping `qe` as the origin costs nothing in
stiffness and preserves a property their form gives up: `A` acts on `(u - qe)`,
so `f_imp(qe) = 0` **exactly**, `qe` stays an exact steady state of the implicit
half, and whatever discrete hydrostatic imbalance the code has stays entirely
inside the explicit part. Their motivation for the previous-solution form — that
a balanced reference state may not exist for real data, or for whole-atmosphere
runs with large day/night swings — does not apply to a deck that always
defines `qe`.

**Measured on `rtb_hevi`** (1 rank, tend = 100 s):

| mode | s/step | max\|u\| at t=100 |
|---|---|---|
| `:RS` | **0.207** | 2.6760e+00 |
| `:PS` | 0.233 (+12.6%) | 2.6760e+00 |

Identical answer, 12.6% more expensive, and with the stability guard off both
are still stable at Δt = 0.45 — so `:PS` buys no step size either. Expected:
`beta` and `thetabar` barely move for a 2 K bubble on a 300 K background, so
refreshing them costs a refactorisation and recovers a coefficient that was
never stale. **`rtb_hevi` therefore ships `:RS`.**

`:PS` earns its cost when the solution departs far from any fixed reference
state — a 100-day baroclinic instability, real data, or a whole-atmosphere run.

---

## 8. Deck options

```julia
:ode_solver     => HEVI_ARK(:ARS232),   # instead of CarpenterKennedy2N54()
:Δt             => 0.19,

:hevi_verify        => true,     # setup self-check + stability guard (default)
:hevi_linearization => :RS,     # :RS (default) or :PS -- see 7b
:hevi_update_freq   => 5,       # steps between refreshes (:PS only)
:hevi_wall_flux => true,        # zero implicit vertical mass flux at floor/lid
:hevi_vars      => [1,4,5],     # override the implicit variable set
:lcfl_report    => true,        # print the stability table at startup
:lxy_partition  => true,        # columnar partition (HEVI default)
```

`:ladapt => true` is refused: adaptation invalidates the column topology, the
factorised column matrices and the gather/scatter plan.

---

## 9. Limitations

* **Bounded by acoustic anisotropy.** Removes the vertical acoustic term and
  nothing else. On an isotropic mesh HEVI is slower than explicit.
* **Horizontal acoustics stay explicit.** Getting past the anisotropy ratio
  needs acoustic substepping or a full 3D implicit solve. *(Not implemented.)*
* **Vertical diffusion stays explicit.** If the CFL report names SGS diffusion
  as the binding term, this split will not help. *(Not implemented.)*
* **Fixed step.** The tableaux carry no embedded error estimator.
* **No AMR.**

---

## 10. Where the methods come from

Compiled from the standard literature; verify the bibliographic details against
the originals before citing them in a paper.

### IMEX additive Runge–Kutta

* **Ascher, U. M., Ruuth, S. J. & Spiteri, R. J. (1997).** *Implicit–explicit
  Runge–Kutta methods for time-dependent partial differential equations.*
  Applied Numerical Mathematics **25**(2–3), 151–167.
  — The ARS family. `:ARS111 :ARS121 :ARS122 :ARS222 :ARS232 :ARS343 :ARS443`
  in `ark.jl` are these tableaux; the naming ARS(s\_exp, s\_imp, order) is theirs.
* **Kennedy, C. A. & Carpenter, M. H. (2003).** *Additive Runge–Kutta schemes
  for convection–diffusion–reaction equations.* Applied Numerical Mathematics
  **44**(1–2), 139–181.
  — The additive-RK framework: two tableaux sharing stage values, order
  conditions, and stiff accuracy (which is why ARS232's last stage is the
  answer and its RHS evaluation is reused as the FSAL).

### HEVI for the nonhydrostatic atmosphere

* **Giraldo, F. X., Kelly, J. F. & Constantinescu, E. M. (2013).**
  *Implicit–explicit formulations of a three-dimensional nonhydrostatic unified
  model of the atmosphere (NUMA).* SIAM Journal on Scientific Computing
  **35**(5), B1162–B1194.
  — The closest antecedent of this implementation: IMEX/HEVI splitting of the
  compressible Euler equations in a spectral-element atmospheric model,
  including the choice to linearise the implicit operator about the reference
  state so it is affine-free.
* **Giraldo, F. X., de Bragança Alves, F. A. V., Kelly, J. F., Kang, S. &
  Reinecke, P. A. (2023).** *A Performance Study of Horizontally Explicit
  Vertically Implicit (HEVI) Time-Integrators for Non-Hydrostatic Atmospheric
  Models.* arXiv:2311.11425.
  — Times three HEVI variants in NUMA (a spectral element code): NHEVI-GMRES,
  NHEVI-LU and **LHEVI**. LHEVI wins by a wide margin — 5× NHEVI-LU, 10×
  NHEVI-GMRES — with the 2nd/3rd-order ARK pairs fastest. Their Eq. (47),
  `{S(q) - L(q)}_EX + [L(q)]_IM`, is the split implemented here, and their
  banded bandwidth `kl = ku = nvar(Nζ+1) - 1` is the one this code uses. Source
  of the `:PS` variant in §7b, and of the caution that their §6 stability
  analysis is run WITH 4th-order hyperdiffusion at every stage — which is why
  their ARS(3,4,3) result and the ARS343 verdict in §4.2 are both correct (see
  the note there).
* **Weller, H., Lock, S.-J. & Wood, N. (2013).** *Runge–Kutta IMEX schemes for
  the horizontally explicit/vertically implicit (HEVI) solution of wave
  equations.* Journal of Computational Physics **252**, 365–381.
  — The **joint** two-parameter stability analysis of §4: a HEVI mode loads the
  explicit and implicit halves simultaneously, so `R(z_E, z_I)` over a
  rectangle is the figure of merit and the explicit imaginary radius is not.
  `ark_joint_amplification` implements this.
* **Ullrich, P. A. & Jablonowski, C. (2012).** *Operator-split
  Runge–Kutta–Rosenbrock methods for nonhydrostatic atmospheric models.*
  Monthly Weather Review **140**, 1257–1284.
* **Bao, L., Klöfkorn, R. & Nair, R. D. (2015).** *Horizontally explicit and
  vertically implicit (HEVI) time discretization scheme for a discontinuous
  Galerkin nonhydrostatic model.* Monthly Weather Review **143**, 972–990.
* **Durran, D. R. & Blossey, P. N. (2012).** *Implicit–explicit multistep
  methods for fast-wave–slow-wave problems.* Monthly Weather Review **140**,
  1307–1325.
  — The fast-wave/slow-wave framing, and why the slow part must genuinely be
  slow.

### Split-explicit / acoustic substepping

* **Klemp, J. B. & Wilhelmson, R. B. (1978).** *The simulation of
  three-dimensional convective storm dynamics.* Journal of the Atmospheric
  Sciences **35**, 1070–1096.
  — The original time-split scheme with **forward–backward** acoustic
  substeps. The inner scheme in `substep.jl` is this one.
* **Wicker, L. J. & Skamarock, W. C. (2002).** *Time-splitting methods for
  elastic models using forward time schemes.* Monthly Weather Review **130**,
  2088–2097.
  — The **RK3 outer step** (Δt/3, Δt/2, Δt, each starting from uⁿ) used in
  `substep.jl`.
* **Skamarock, W. C. & Klemp, J. B. (1992).** *The stability of time-split
  numerical methods for the hydrostatic and the nonhydrostatic elastic
  equations.* Monthly Weather Review **120**, 2109–2127.
  — **Required reading for the current state of `substep.jl`.** Time-split
  schemes are weakly unstable as constructed: the frozen slow tendency and the
  substepped acoustics interact, and stability needs explicit **divergence
  damping** and/or off-centring of the acoustic substep. Neither is
  implemented yet, which is consistent with the outer step limit measured
  below.
* **Klemp, J. B., Skamarock, W. C. & Dudhia, J. (2007).** *Conservative
  split-explicit time integration methods for the compressible nonhydrostatic
  equations.* Monthly Weather Review **135**, 2897–2913.
  — The formulation used in WRF, including the off-centring parameter of the
  vertically-implicit acoustic step and the 3D divergence damping.

### Explicit reference scheme and the benchmark

* **Carpenter, M. H. & Kennedy, C. A. (1994).** *Fourth-order 2N-storage
  Runge–Kutta schemes.* NASA Technical Memorandum 109112.
  — `CarpenterKennedy2N54`, the explicit scheme every measurement here is
  compared against.
* **Robert, A. (1993).** *Bubble convection experiments with a semi-implicit
  formulation of the Euler equations.* Journal of the Atmospheric Sciences
  **50**, 1865–1873.
* **Giraldo, F. X. & Restelli, M. (2008).** *A study of spectral element and
  discontinuous Galerkin methods for the Navier–Stokes equations in
  nonhydrostatic mesoscale atmospheric modeling.* Journal of Computational
  Physics **227**, 3849–3877.
  — The rising thermal bubble in the form used by `rtb_hevi` and
  `CompEuler/theta`.

### Numerical linear algebra

The per-column solve is LAPACK's banded LU (`dgbtrf` / `dgbtrs`) through
Julia's `LinearAlgebra`. The probing extraction of a banded matrix by
colouring the rows modulo the bandwidth is the standard sparse-Jacobian
colouring argument, e.g. **Curtis, A. R., Powell, M. J. D. & Reid, J. K.
(1974),** *On the estimation of sparse Jacobian matrices*, IMA Journal of
Applied Mathematics **13**, 117–120.

---

## 11. Status of acoustic substepping (`substep.jl`)

Implemented and **not yet delivering its intended benefit.** Recorded here so
the state is not mistaken for a finished feature.

**What works.** The fast operator (§2, `acoustic.jl`) is verified exactly: on a
horizontally uniform field it reproduces the vertical HEVI operator to `0.0`,
bit-identical at 1 and 4 ranks. In isolation, `dV/dt = A V` under
forward–backward substeps is essentially neutral — growth per substep
1.000059 / 1.000225 / 1.000765 at Δτ = 0.02 / 0.04 / 0.08, scaling like Δτ².
The integrator runs and is stable at Δt ≤ ≈0.2 s.

**What does not.** The outer step should be limited by advection and diffusion
(≈4.8 s on this mesh); measured, it is limited to ≈0.2 s, about twice the
explicit acoustic limit of 0.0995 s. Above that it diverges within five or six
outer steps.

**Eliminated by measurement**, each with a dedicated run:

| hypothesis | test | result |
|---|---|---|
| substep Δτ too large | rate-summed limit, stage counts corrected | no change |
| floor/lid flux mismatch between `f_fast` and `rhs!` | `:hevi_wall_flux => false` | no change |
| solution-dependent artificial viscosity in the frozen slow part | `:lvisc => false` | no change |

**Most likely remaining cause**, in order:

1. **No divergence damping or off-centring.** Skamarock & Klemp (1992) show a
   time-split scheme is weakly unstable without them, and Klemp et al. (2007)
   give the forms used in practice. This is the literature-supported next step.
2. **Boundary projection.** `rhs!` applies free-slip by projecting the state at
   boundary nodes, so `S = rhs!(Pu) − f_fast(u)`; nothing re-applies `P` during
   the substeps, so `u` drifts off the boundary condition within a stage, and
   the drift grows with stage length. Split-explicit codes apply the wall
   condition inside the substeps.
