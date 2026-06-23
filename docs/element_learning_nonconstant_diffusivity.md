# Mesh‑induced reference diffusivity `â` for Element Learning

This note documents the **non‑constant diffusivity** capability added to the
Element Learning (EL) solver: what the mesh‑induced reference diffusivity `â`
is, how it is generated/consumed in the code, and how each piece was verified.

It accompanies the code in

- `src/kernel/EL_nonconstant_diffusivity.jl` — the `â` machinery,
- `src/kernel/elementLearningStructs.jl`    — sampling/inference wiring,
- `src/io/write_output.jl`                  — VTU visualisation,
- `problems/Elliptic/elementLearning_hole_nonconstant/` — the geometry test,
- `problems/Elliptic/case1_nonconstant/`    — the physical‑`a(x,y)` reference test,
- `tools/EL_training/train_rfrc.py`         — the matching RF‑RC trainer.

---

## 1. What `â` is

On an element `K` the local stiffness entry of the variable‑coefficient
operator `-∇·(a ∇u)` is pulled back to the reference element `[-1,1]^d`:

```
(A_{vᵒ,v})_{ij} = ∫_K a(x) ∇x φ_i · ∇x φ_j
                = ∫_{[-1,1]^d} (â ∇ξ φ_i) · (∇ξ φ_j),
```

with the **reference diffusivity**

```
â(ξ) = (|K| / |[-1,1]^d|) · a · J_K⁻¹ J_K⁻ᵀ          (a 2×2 SPD tensor in 2-D),
```

where `J_K = ∂x/∂ξ` is the element Jacobian. `â` fuses two ingredients:

1. the **physical** scalar diffusivity `a` of the PDE, and
2. the **geometry** of the element through `J_K`.

The EL surrogate learns the map `â ↦ T^{ie} = (A_{vᵒ,vᵒ})⁻¹ A_{vᵒ,vᵇ}`, i.e. the
local interior‑elimination operator, replacing an explicit per‑element matrix
inverse.

### Two ways `â` becomes non‑constant

| Source of non‑constancy | Physical `a` | Geometry `J_K`            | Test case |
|-------------------------|:------------:|:-------------------------:|-----------|
| **Physical** (`a(x,y)`) | varies       | uniform (structured mesh) | `case1_nonconstant` |
| **Mesh‑induced**        | `a = 1`      | varies (non‑affine elements: general quads / curved) | `elementLearning_hole_nonconstant` |

This document focuses on the **mesh‑induced** case, which is the one the EL
notes prescribe (Option 1: *“use an unstructured mesh, which supplies a set of
random element shapes”*). With constant physical `a = 1`, the only non‑trivial
part of `â` is the Jacobian — exactly the *“simplest case”* of the notes.

### Constant within an element, or varying? (affine vs. bilinear)

`â` is constant **within** an element only when the element map is **affine**,
which for a quadrilateral means a **parallelogram** (opposite sides parallel
*and* equal) — then `J_K` is constant. "Straight sides" is *not* sufficient: a
general straight‑sided quad (trapezoid, kite, or the sheared/fanned quads that
conform to a curved boundary such as the hole) has a **bilinear** map
`x(ξ,η) = Σ Nᵢ(ξ,η) Pᵢ`, so `∂Nᵢ/∂ξ ∝ (1±η)` and `J_K = ∂x/∂ξ` varies
**linearly inside the element**. Hence `â(ξ)` — and the VTU `ahat_eff`,
`ahat_aniso`, `ahat_11/12/22` fields — **vary within each non‑parallelogram
element**, peaking at the most distorted corners. This within‑element gradient
in the plots is expected and correct; only true parallelograms show a single
value per element. (Curved/high‑order geometry is just a further source of the
same `ξ`‑dependence.)

> **Sampling consistency.** The default sampler (`:EL_sample_shape => :affine`)
> draws random **affine** Jacobians (parallelograms ⇒ `â` constant per element).
> To train a surrogate that matches *general* (bilinear) mesh quads, set
> **`:EL_sample_shape => :quad`** — it synthesises random straight-sided bilinear
> quads whose per-node `â` equals `el_ahat_nodes_from_metrics` for the same shape
> (verified exactly), so the training `â` matches what inference feeds. This is
> correct at any solution order `k` (straight sides ⇒ bilinear geometry from 4
> corners; `â` is still evaluated at all `(k+1)²` nodes). (`:warp` is a legacy
> synthetic-warp approximation; `:lEL_xidependent => true` is its alias.)

### A useful identity (area normalisation)

For `a = 1`,

```
det(â) = (|K|/|ref|)² · det(J_K⁻¹ J_K⁻ᵀ) = det(J_K)² · (1/det(J_K))² = 1.
```

So **`det(â) = 1` identically**: the mesh‑induced `â` is *area‑normalised* and
its variation is **purely anisotropic**. Two scalar invariants then fully
describe it (because `λ_max·λ_min = det = 1 ⇒ λ_min = 1/λ_max`):

```
ahat_eff   = tr(â)/2 = (λ_max + λ_min)/2 ≥ 1     (= 1 ⇔ undistorted element)
ahat_aniso = λ_max / λ_min                        (= 1 ⇔ isotropic / square-like)
```

These are the fields written to the VTU (Section 3).

---

## 2. Implementation

### 2.1 Feature representation

The NN input feature is the per‑node `â` tensor, stored as its 3 unique SPD
entries `(â11, â12, â22)` at each of the `(k+1)²` reference nodes:

```
feature length = 3·(k+1)²            (was (k+1)² for the constant-scalar case)
layout (conn node m): [â11_m, â12_m, â22_m]
```

The NN output is the flattened `T^{ie}` of size `nvo·nvb`
(`nvo = (k-1)²` interior nodes, `nvb = (k+1)² - (k-1)²` boundary nodes), unchanged.

### 2.2 Sampling (training‑data generation) — Option 1

`el_nonconstant_sampling!` synthesises random element shapes
(`:EL_sample_shape` = `:affine` parallelogram, `:quad` bilinear straight-sided
quad, or `:warp` synthetic warp) and, per sample,

1. forms the per-node `â = a·det(J)·(JᵀJ)⁻¹` — constant for `:affine`, or varying
   within the element for `:quad`/`:warp` (the `:quad` field equals
   `el_ahat_nodes_from_metrics` for the same shape),
2. assembles the local reference stiffness
   (`el_assemble_local_stiffness` via the precomputed components `K11,K12,K22`,
   or `el_assemble_local_stiffness_field` for the ξ‑dependent case),
3. condenses `T^{ie} = (A_{vᵒvᵒ})⁻¹ A_{vᵒvᵇ}`,
4. writes `(feature, vec(T^{ie}))` in **`mesh.conn` node order**.

### 2.3 Inference

`el_avisc_nonconstant!` builds the per‑element feature from the **real mesh
metrics** `â = Je · a · J_K⁻¹ J_K⁻ᵀ` (with `J_K⁻¹ = [dξdx dξdy; dηdx dηdy]`,
`a = el_diffusivity(x,y)`, default 1), in the same `mesh.conn` order. Because
sampling and inference share that single ordering, `elementLearning_infer!` is
reused **unchanged** — only the feature differs.

Gated by `inputs[:lEL_nonconstant]` (default `false`), so the existing
constant‑amplitude pipeline and its trained model are untouched.

> **Amplitude consistency (training ↔ inference).** The exact map `â ↦ T^{ie}`
> is invariant under a *uniform* rescaling of `â`, so for `a` constant within an
> element `T^{ie}` does not depend on `a`. The **learned** surrogate is *not*
> analytically scale‑invariant, however (it is a standardised NN), so an `â`
> whose magnitude falls outside the training range is extrapolated poorly. The
> sampler therefore draws a physical amplitude `a ~ U(:EL_amin, :EL_amax)` per
> element (default `(1,1)` ⇒ legacy geometry‑only `â`), and `:EL_avar > 0` adds a
> smooth **within‑element** variation of `a` (so `T^{ie}` genuinely varies). To
> infer with a non‑trivial `el_diffusivity` (e.g. `a = 1 + 1(x>5)`), set
> `:EL_amin/:EL_amax` to bracket its range — and `:EL_avar > 0` if `a` varies
> inside an element.

### 2.4 Operator for the physical reference test

For `case1_nonconstant`, `build_laplace_matrix(NSD_2D; afun)` multiplies the
integrand by `a(x,y)` at each node so the discrete operator is genuinely
`-∇·(a∇u)`; `afun === nothing` reproduces the plain Laplacian exactly.

### 2.5 VTU visualisation

`write_vtk(NSD_2D)` writes the geometry‑induced `â` to the VTU when
`inputs[:lEL_nonconstant]` is set and `metrics` are supplied. **Three output
formats** are available via `inputs[:ahat_output]` (default `:cell`):

| `inputs[:ahat_output]` | VTU fields | Location | Description |
|------------------------|------------|----------|-------------|
| `:cell`  (default) | `ahat_eff`, `ahat_aniso` | cell data  | invariants, element‑averaged (piecewise‑constant per element) |
| `:nodal`           | `ahat_eff`, `ahat_aniso` | point data | same invariants, **smoothed nodal** field (`â` averaged over the elements sharing each node, then invariants taken) |
| `:tensor`          | `ahat_11`, `ahat_12`, `ahat_22` | point data | the full smoothed nodal **tensor components** |

with `ahat_eff = tr(â)/2 = (λmax+λmin)/2 ≥ 1` (`1` ⇔ undistorted) and
`ahat_aniso = λmax/λmin` (`1` ⇔ isotropic). The smoothing averages the three
`â` components over the elements incident to each node (count‑weighted).

For the *physical* case, `case1_nonconstant` instead writes the nodal scalar
field `diffusivity = a(x,y)` (gated on `inputs[:diffusivity]`).

---

## 3. Verification

The full Jexpresso solve (GMSH mesh + trained ONNX model) was **not** run in
this environment, so verification targets the mathematics and the data flow
with standalone Julia checks. All tolerances below are machine precision
(`~1e-16`–`1e-10`) unless noted; all checks **passed**.

### 3.1 Manufactured solution (physical reference, `case1_nonconstant`)

For `a(x,y) = 1 + 0.2x`, `u_ex = sin(x)cos(y)`, the analytic source

```
f = -A1·kx·cos(kx x)cos(ky y) + (A0 + A1 x)(kx² + ky²) sin(kx x)cos(ky y)
```

was checked against a central finite‑difference evaluation of `-∇·(a∇u_ex)`:

```
max | -∇·(a∇u_ex)_FD − f | ≈ 2.6e-6     (FD truncation; PASS)
a ∈ [1.2, 2.2] > 0 on the domain
```

### 3.2 Diffusivity‑weighted operator (`build_laplace_matrix`)

Replicating the assembly on a reference element:

- `afun === nothing` reproduces the pure reference Laplacian (matches the EL
  module’s `K11+K22`) → **existing problems unaffected**;
- constant `a = c` scales the operator by exactly `c`;
- non‑constant `a(x,y)` keeps the element matrix **symmetric, PSD, with the
  constant nullspace**.

### 3.3 Local `â` machinery (`EL_nonconstant_diffusivity.jl`)

- the reference Laplacian (`â = I`) is symmetric, PSD, constant nullspace;
- the precomputed `K`‑component assembly equals the full `O(N⁶)` quadrature
  assembly for constant `â`;
- synthesised `â` is **SPD on every draw** (2000 samples);
- `T^{ie}` reproduces the local Dirichlet interior recovery
  `u_vᵒ = -T^{ie} u_vᵇ` vs the direct block solve;
- the metrics‑based `â` (inference) equals the synthesised `â` (sampling) for
  the **same Jacobian** — i.e. **sampling ↔ inference consistency**;
- feature/output dimensions are `3·(k+1)²` and `nvo·nvb` (e.g. `nop=6` →
  `147` and `25×24`).

### 3.4 Inference‑wiring node ordering

- the `mesh.conn ↔ tensor (i,j)` map round‑trips on the element layout;
- the conn‑ordered `T^{ie}` reproduces the direct interior recovery, confirming
  `elementLearning_infer!` can be reused for the non‑constant feature;
- ξ‑dependent (curved‑element) `â` is SPD and genuinely node‑varying, and its
  assembled stiffness is symmetric/PSD with the constant nullspace.

### 3.5 VTU `â` invariants (`write_vtk`)

- **undistorted** element `J = s·I`: `ahat_eff = 1`, `ahat_aniso = 1`,
  `det(â) = 1` at **all scales** `s` (scale‑invariant, as expected);
- **arbitrary distorted** `J` (rotation/stretch/shear sweep): `det(â) = 1`
  always, `â` SPD, `ahat_eff ≥ 1`, `ahat_aniso ≥ 1`;
- the VTU invariants are computed from the **same** `â` entries as
  `el_ahat_nodes_from_metrics`, so *what you see equals what the surrogate
  consumes*.

### 3.6 Trainer

`tools/EL_training/train_rfrc.py` auto‑sizes `N_in` from the CSV (so `3·(k+1)²`
is picked up automatically) and standardises the heterogeneous `â` feature,
baking the standardisation into the exported ONNX so the Julia side keeps
feeding the **raw** `â` feature. The standardise→bake‑in round‑trip
(`export(raw) == train(standardised)`) was verified.

---

## 4. How to run

### Geometry‑induced `â` test (the subject of this note)

```julia
julia --project=.
julia> push!(empty!(ARGS), "Elliptic", "elementLearning_hole_nonconstant");
julia> include("./src/Jexpresso.jl")
```

- **Verify the physical problem and see `â` now** (no model needed): comment
  `:lelementLearning` in `user_inputs.jl`. The direct solver solves `-Δu = f`
  on the hole, prints the MMS `‖e‖_L2`, and writes `ahat_eff` / `ahat_aniso`
  to `output-nonconstant/iter_1.pvtu`. The anisotropy concentrates around the
  circular hole, where elements are most distorted.
- **EL surrogate inference**: keep `:lelementLearning => true` with a model
  trained on the `3·(k+1)²` `â` feature (sample with `:lEL_Sample`, train with
  `tools/EL_training/train_rfrc.py`).

### Physical `a(x,y)` reference test

```julia
julia> push!(empty!(ARGS), "Elliptic", "case1_nonconstant");
julia> include("./src/Jexpresso.jl")
```

Direct solve of `-∇·(a∇u)=f` with `a = 1 + 0.2x`; prints the MMS `‖e‖_L2` and
writes the nodal `diffusivity` field to the VTU.

---

## 5. Known limitations / follow‑ups

- The EL **inference** path needs a model retrained on the `3·(k+1)²` `â`
  feature; running with the old `(k+1)²` model + `:lEL_nonconstant` fails on a
  dimension mismatch (expected).
- For a *physical* non‑constant `a(x,y)` solved **through** the EL surrogate,
  the global skeleton matrix `L` should also be assembled with that `a` (use
  `inputs[:diffusivity]`); the mesh‑induced (`a = 1`) case is already
  self‑consistent because `L` is the standard Laplacian.
- VTU `â` fields and the diffusivity operator are **2‑D** (`NSD_2D`), matching
  the current EL scope.
- Three `â` VTU formats are provided (`inputs[:ahat_output]` = `:cell` /
  `:nodal` / `:tensor`, Section 2.5). The nodal smoothing is a count‑weighted
  average over incident elements; across MPI partitions it is a per‑partition
  average (a visualisation field, not used in the solve).
