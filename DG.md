# Discontinuous Galerkin (`:AD => DiscGal()`)

Jexpresso's default spatial discretization is the continuous spectral element
method (`:AD => ContGal()`). Setting `:AD => DiscGal()` in a case deck selects
the **discontinuous** Galerkin path instead: every element owns its own copy of
the nodes on its faces, and elements are coupled only through a numerical flux
evaluated at the interfaces.

The DG path is **additive**. Every `DiscGal` method sits beside its `ContGal`
twin and the CG/FD paths are untouched, so a case that does not set `:AD` is
unaffected.

---

## 1. What a DG deck sets

```julia
:AD                  => DiscGal(),
:numerical_flux      => rusanov_flux(),   # or upwind_flux()
:lexact_integration  => false,            # collocated LGL
:lvisc               => false,            # REQUIRED -- see §5
```

and the equation set supplies one extra hook next to `user_flux!`:

```julia
# 1D: no normal argument
user_max_wave_speed(q, qe, SD::NSD_1D, ::TOTAL; neqs=1)

# 2D: the face unit normal is supplied by the surface term
user_max_wave_speed(q, qe, SD::NSD_2D, ::TOTAL; nx=1.0, ny=0.0, neqs=1)
```

This is the `λ` of the Rusanov/local-Lax-Friedrichs jump term. For a system
with a wave speed, it is `|u·n| + c` — e.g. `abs(u*nx + v*ny) + sqrt(g*H)` for
the shallow water equations.

## 2. Where the code lives

| file | what it holds |
|---|---|
| `src/kernel/operators/dg_fluxes.jl` | `numerical_flux!` (Rusanov, upwind), `surface_rhs_el!` (1D and 2D), `dg_boundary_ghost!` |
| `src/kernel/mesh/mesh.jl` | `add_high_order_nodes_1D_native_mesh_dg!`, `add_high_order_nodes_2D_gmsh_dg!` (duplicated-DOF numbering), `build_dg_faces_2D!` (face lists) |
| `src/kernel/mesh/meshStructs.jl` | `dg_face_*` (interior + periodic pairs), `dg_bfac_*` (physical boundary faces) |
| `src/kernel/operators/rhs.jl` | `::DiscGal` volume kernels and the `surface_rhs_el!` call site |
| `src/kernel/infrastructure/element_matrices.jl` | `DSS_rhs!`, `divide_by_mass_matrix!`, `matrix_wrapper` under `::DiscGal` |
| `src/kernel/boundaryconditions/BCs.jl` | `::DiscGal` no-ops: DG imposes boundary conditions weakly |

The scheme is the strong form. The volume kernel is the CG one (the weak form
is discretization-agnostic); `surface_rhs_el!` adds the interface correction
`±(F_int − F*)`, and the `1/M` lift comes from the existing
`divide_by_mass_matrix!`.

## 3. Boundary conditions are fluxes, not node values

`build_dg_faces_2D!` builds two lists and prints both:

```
 # build_dg_faces_2D!: 1445 interior + 0 periodic = 1445 faces
 # build_dg_faces_2D!: 110 physical boundary faces (0 periodic facets already paired above, 0 untagged => free/transmissive)
```

* **interior and periodic** faces have two traces and are indistinguishable at
  run time — periodicity under DG is a flux face, never a node merge.
* **physical boundary** faces have one. The second trace is built by
  `dg_boundary_ghost!` from the case's own `user_bc_dirichlet!`:

  ```
  q⁺ = 2·q_bc − q⁻
  ```

  so the mean of the two traces is exactly what the case prescribed. For a
  free-slip wall, where the case zeroes the normal momentum, this is the
  textbook mirror state — no case writes the mirror twice. Components the case
  leaves at the sentinel are copied from the interior trace, i.e. left free to
  leave the domain. An **untagged** boundary facet gets no row, which is
  `F* = F_int`: free/transmissive. The count is printed so a domain meant to
  be closed cannot quietly be open.

The CG strong path (`apply_boundary_conditions_dirichlet!`, which overwrites
`uaux` and zeroes `RHS` at `poin_in_bdy_edge`) is a `DiscGal` no-op. Under DG
that array carries CG point ids while `mesh.x`/`connijk` have been renumbered,
so its writes would land on unrelated nodes.

## 4. Cases

| case | what it is | run with |
|---|---|---|
| `AdvDiff/advection1d_dg` | 1D scalar advection, periodic, `U = 2` on `[-1,1]`; one full traversal returns the IC. **In CI.** | `julia --project=. test/runtests.jl AdvDiff/advection1d_dg` |
| `CompEuler/wave1d_dg` | 1D acoustic wave system, periodic | `Jexpresso.run_case("CompEuler", "wave1d_dg")` |
| `AdvDiff/advection2d_dg` | 2D advection on a doubly periodic 10×20 mesh. **In CI.** | `julia --project=. test/runtests.jl AdvDiff/advection2d_dg` |
| `ShallowWater/SoliWaveIsland_dg` | 2D non-linear shallow water, solitary wave in a closed basin with four free-slip walls — the first DG case with **physical** boundaries. | `julia --project=. -e 'using Jexpresso; Jexpresso.run_case("ShallowWater","SoliWaveIsland_dg")'` |

## 5. What is not implemented

* **Viscous terms.** `_expansion_visc!` has `ContGal` methods only, so
  `:lvisc => true` stops the run with a `MethodError` on the first RHS call —
  loudly, not silently. This holds for `AV()` and for `DSGS_SW()`/`DSGS()`
  alike: DSGS computes its coefficient fine under `DiscGal` and then hands it
  to the CG operator. A DG viscous term needs the element volume part *and* an
  interface term (BR1/BR2 or interior penalty); an element-local Laplacian
  alone has no interface coupling and is not the operator it claims to be.
  This is the single blocker on the full `SoliWaveIsland_dg` run-up — see that
  deck's header.
* **Shock capturing / positivity limiting.** Nothing beyond the interface
  Rusanov jump. Cases with a wet/dry front or a steepening front need more.
  Worth noting for whoever picks this up: `src/kernel/positivity/README.md`
  rules Zhang-Shu out for CG because it needs elements that can be modified
  independently and a cell average kept positive by a positivity-preserving
  interface flux, and CG-SEM has neither. `DiscGal` supplies both. That
  repair module itself is compressible-Euler-specific (it repairs rho and p)
  and off unless `:lpositivity => true`, so it does not apply to the shallow
  water case as it stands.
* **Mesh adaptivity.** Refused at input parsing (`:lamr`, `:ladapt`,
  `:lpreadapt` must be false): the CG mortar projections have no meaning on
  duplicated DOFs.
* **1D physical boundaries.** The 1D face list is interior (+ periodic wrap)
  only; a non-periodic 1D DG run gets no term at the two domain ends. 2D
  carries them.
* **3D.** No `conformity4ncf_q!(::NSD_3D, ::DiscGal)` and no 3D face list.
* **MPI.** `ip2gip` is the identity under DG — serial semantics only.
* **KEP** (`:lkep`) and the **GPU** kernels dispatch on `AD` with no `DiscGal`
  method.
