# Boundary Conditions in JEXPRESSO

**Code structure and user-side workflow**

---

## 1. Overview

In JEXPRESSO the application of boundary conditions (BCs) is split between two layers:

1. A **kernel layer** located in `src/kernel/boundaryconditions/`, which contains the generic machinery that loops over boundary edges (2D) or faces (3D), reads the `tag` associated with every boundary entity in the Gmsh mesh, fetches normals, and pushes the boundary-node values either into `uaux` (for time-dependent problems with explicit time integration) or directly into the global linear system *L x = b* (for elliptic problems).

2. A **user layer** located in `problems/<EQUATION>/<CASE>/user_bc.jl`, where the user prescribes the actual boundary values or boundary fluxes by writing one (or two) Julia functions: `user_bc_dirichlet!` (essential conditions) and `user_bc_neumann` / `user_bc_neumann!` (natural conditions / surface fluxes).

The bridge between the two layers is provided by `src/kernel/boundaryconditions/custom_bcs.jl`, which contains the thin wrappers `dirichlet!` and `neumann` that simply forward the call to the user-defined `user_bc_*` functions.

The high-level driver in `problems/drivers.jl` chooses *which* BC entry point is invoked according to the type of problem (hyperbolic with explicit time integration, elliptic *A x = b*, or hyperbolic with implicit/IMEX time integration that places terms on the left-hand side).

---

## 2. Kernel-side structure

### Files in `src/kernel/boundaryconditions/`

**`BCs.jl`** — Defines the public entry points called by the driver and the time integrators:

- `apply_boundary_conditions_dirichlet!`
- `apply_boundary_conditions_neumann!`
- `apply_boundary_conditions_lin_solve!`
- `apply_periodicity!`

Each public entry point dispatches to a dimension-specific `build_custom_bcs_*!` method using the spatial-dimension type (`NSD_1D`, `NSD_2D`, `NSD_3D`). For the linear-solve case there is in addition a sparse variant (`build_custom_bcs_lin_solve_sparse!`) that mutates a `SparseMatrixCSC`.

**`custom_bcs.jl`** — Wrappers `dirichlet!` and `neumann` that simply forward to the user-defined `user_bc_dirichlet!` / `user_bc_neumann`.

**`surface_integral.jl`** — Helpers that build edge/face mass matrices and assemble surface integrals (`compute_surface_integral!`, `compute_segment_integral!`, `DSS_*_integral!`, `bulk_surface_flux!`). They are used by the Neumann path to add the boundary surface flux contribution to the RHS:

$$\text{RHS}_i \mathrel{+}= \int_{\partial \Omega} F_{\text{surf}}\,\psi_i\,dS$$

### 2.1 Dirichlet machinery

For 2D, `build_custom_bcs_dirichlet!(::NSD_2D, ...)` loops over every boundary edge `iedge`; for each edge it skips `periodicx`, `periodicz`, `periodic1`, `periodic2` and `Laguerre` tags, and for the remaining *ngl* Lagrange-Gauss-Lobatto (LGL) nodes it:

1. Fills a sentinel value `qbdy[ieq] = 4325789.0`;
2. Calls `user_bc_dirichlet!(@view(uaux[ip,:]), ..., bdy_edge_type[iedge], qbdy, nx, ny, @view(qe[ip,:]), inputs[:SOL_VARS_TYPE])`;
3. If the user changed the sentinel and the value differs from the current solution, overwrites `uaux[ip,ieq]` and zeroes `RHS[ip,ieq]` (so the time integrator does not modify the constrained DOF).

At the end the array `u` is repacked from `uaux` via `uaux2u!(u, uaux, neqs, npoin)`.

The 3D version is analogous, looping on `nfaces_bdy` and using `poin_in_bdy_face` together with the 3-component normal *(n_x, n_y, n_z)*.

### 2.2 Neumann machinery

`build_custom_bcs_neumann!` is invoked only when `inputs[:bdy_fluxes] = true`. For each boundary edge/face it:

1. Zeroes `F_surf`;
2. If `inputs[:bulk_fluxes]` is true, calls the analytic `bulk_surface_flux!`; for the special tag `"MOST"` in 3D it calls the Monin-Obukhov routine `CM_MOST!`;
3. Otherwise calls the user routine `user_bc_neumann!(@view(F_surf[...]), uaux[ip,:], uaux[ip1,:], qe[ip,:], qe[ip1,:], bdy_edge_type[iedge], coords, ..., inputs[:SOL_VARS_TYPE])`;
4. Integrates over the boundary segment/face with `compute_segment_integral!` or `compute_surface_integral!`;
5. Performs DSS and adds the result to the global RHS: `RHS[:,ieq] .+= S_flux[:,ieq]`.

### 2.3 Linear-solve machinery

The variant used by elliptic problems is `apply_boundary_conditions_lin_solve!`, which dispatches to either the dense or the sparse kernel:

- **`build_custom_bcs_lin_solve!(::NSD_2D, ...)`** (dense `Matrix{Float64}`). For every Dirichlet boundary node `ip` it:
  1. Calls `user_bc_dirichlet!` *but writes the prescribed value into `RHS[ip,:]`* (and into the helper `qbdy`);
  2. Zeroes the row `L[ip, :] = 0`;
  3. Sets `L[ip, ip] = 1`;
  4. Assigns `RHS[ip, ieq] = qbdy[ieq]`.

  This is the standard "identity row" enforcement of essential BCs in *A x = b*.

- **`build_custom_bcs_lin_solve_sparse!`** does the same on a `SparseMatrixCSC`, calling `apply_dirichlet_bc_inplace!` which collects the unique boundary nodes from `poin_in_bdy_edge` and rewrites `L[ip,:] = 0; L[ip,ip] = 1`.

---

## 3. Driver-side dispatch

The function `driver(...)` in `problems/drivers.jl` performs the high-level dispatch:

```julia
if !inputs[:llinsolve]
    # Hyperbolic / parabolic:  M dq/dt = RHS
    solution = time_loop!(inputs, params, u, partitioned_model)
else
    # Elliptic:  L*q = M*RHS
    ...
    apply_boundary_conditions_lin_solve!(sem.matrix.L, ..., RHS, ...)
    sol = sem.matrix.L \ RHS
    write_output(...)
end
```

Therefore the same `user_bc_dirichlet!` that a user writes in `user_bc.jl` is reused in three different contexts:

1. **Hyperbolic with explicit time integration** — called from inside `time_loop!` via `apply_boundary_conditions_dirichlet!` (and `..._neumann!` when boundary fluxes are activated). The constraint is enforced on `uaux` and the corresponding row of the RHS is zeroed so that the explicit RK step does not drift the boundary values.

2. **Elliptic (*A x = b*)** — called from `apply_boundary_conditions_lin_solve!`. The same user routine is used, but the kernel writes the prescribed value into `RHS` and modifies the system matrix `L` (zeroing the row and putting a 1 on the diagonal).

3. **Hyperbolic with implicit / IMEX time integration** — the implicit operator is assembled in the same Laplace-style matrix `sem.matrix.L` that the elliptic path uses; the same `apply_boundary_conditions_lin_solve!` entry point is invoked at every implicit stage to constrain the rows that correspond to Dirichlet nodes. Neumann contributions are still assembled by `apply_boundary_conditions_neumann!` and added to the right-hand side.

---

## 4. User-side workflow

A user creates a problem under `problems/<EQUATION>/<CASE>/` and is required to provide `user_bc.jl`. The Gmsh `.msh` file declares physical groups whose names become the `tag` string passed to the BC function. Inside `user_bc.jl` the user uses an `if/elseif` ladder on `tag` to set the prescribed values:

```julia
if (tag == "inflow")
    qbdy[1] = 3.0
elseif (tag == "fix_temperature")
    qbdy[2] = 300.0
end
```

Any component that the user does not assign keeps the sentinel `4325789.0`, which the kernel detects and *leaves unchanged* — this is interpreted as a "do nothing" (natural) condition for that component.

### 4.1 Hyperbolic problems with explicit time integration

**Example:** `problems/CompEuler/theta/` (rising thermal bubble).

**`user_inputs.jl`** — The user selects an explicit ODE solver and leaves `:llinsolve` unset (i.e. false, the default), which selects the `time_loop!` branch in the driver:

```julia
:ode_solver    => CarpenterKennedy2N54(),
:Δt            => 0.5,
:tinit         => 0.0,
:tend          => 1000.0,
:SOL_VARS_TYPE => PERT(),
```

**`user_bc.jl`** — A free-slip wall is implemented by projecting out the normal component of the momentum (ρu, ρv). Both the `TOTAL` and the `PERT` variable formulations are provided:

```julia
function user_bc_dirichlet!(q, coords, t, tag, qbdy, nx, ny, qe, ::PERT)
    qnl = nx*(q[2]+qe[2]) + ny*(q[3]+qe[3])
    qbdy[2] = (q[2]+qe[2] - qnl*nx) - qe[2]
    qbdy[3] = (q[3]+qe[3] - qnl*ny) - qe[3]
end
```

The Neumann routine returns a zero flux because no surface flux is needed for this case:

```julia
function user_bc_neumann(q, gradq, coords, t, tag, inputs::Dict)
    return zeros(size(q,2),1)
end
```

**What the kernel does at every RK stage.** The integrator calls `apply_boundary_conditions_dirichlet!`, which for every boundary node `ip` invokes the user routine with the local normal *(n_x, n_y)* taken from `metrics.nx`, `metrics.ny`, copies the returned `qbdy` into `uaux[ip,:]` and zeroes `RHS[ip,:]`. The result is then mapped back into `u` via `uaux2u!`.

### 4.2 Elliptic problems leading to *A x = b*

**Example:** `problems/Elliptic/case1/` (Laplace with non-homogeneous Dirichlet data on a 5 × π rectangle).

**`user_inputs.jl`** — The user activates the linear-solve path:

```julia
:llinsolve  => true,
:ldss_laplace => true,
:rconst     => [0.0],
```

which forces the driver into the `else` branch where the operator `sem.matrix.L` (the Laplacian) is assembled and the system `L * q = M * RHS` is solved by the chosen Krylov method (`:ode_solver => "BICGSTABLE"`) or by the direct backslash.

**`user_bc.jl`** — Only `user_bc_dirichlet!` is needed:

```julia
function user_bc_dirichlet!(q, coords, t, tag, qbdy, nx, ny, qe, ::TOTAL)
    L = 5.0
    if     (tag == "bottom") qbdy[1] = 100.0
    elseif (tag == "right")  qbdy[1] = 100.0
    elseif (tag == "top")    qbdy[1] = 100 + 20.0*sin(pi*coords[1]/L)
    elseif (tag == "left")   qbdy[1] = 100.0
    end
end
```

**What the kernel does once.** Before the linear solve, the driver calls `apply_boundary_conditions_lin_solve!(sem.matrix.L, ..., RHS, ...)`. For every node on a non-periodic boundary, the user routine fills `qbdy`; the kernel then sets

$$L[\text{ip},:] = 0, \qquad L[\text{ip},\text{ip}] = 1, \qquad \text{RHS}[\text{ip}] = \texttt{qbdy[1]}$$

so that the linear system enforces *q_ip = qbdy[1]* exactly. Internally this is the loop in `build_custom_bcs_lin_solve!` (or its sparse counterpart when `:lsparse => true`).

### 4.3 Hyperbolic problems with implicit / IMEX time integration

When part of the operator (typically a stiff diffusion or acoustic operator) is treated implicitly, that operator is assembled into the same matrix `sem.matrix.L` that the elliptic path uses, and the implicit (or IMEX) stage solves

$$(M + \Delta t\,\theta\,L)\,q^{n+1} = M\,q^{n} + \Delta t\,R(q^{n}, q^{n+1})$$

At every implicit stage the BC kernel is therefore called twice:

1. `apply_boundary_conditions_lin_solve!` on the implicit operator and on the right-hand side. This is the same routine described in Section 4.2: rows of the implicit operator corresponding to Dirichlet nodes are replaced by an identity row and the prescribed value is written into the RHS, so that the Krylov solver returns the exact boundary value at each implicit stage.

2. `apply_boundary_conditions_neumann!`, which, when `inputs[:bdy_fluxes] = true`, evaluates the user-defined surface flux through `user_bc_neumann!`, performs the surface quadrature and the DSS, and adds the result to the explicit part of the RHS:

$$\text{RHS}[:,ieq] \mathrel{+}= S_{\text{flux}}[:,ieq]$$

**User responsibility.** The user does *not* have to rewrite the BC code: the same `user_bc_dirichlet!` and `user_bc_neumann!` written for the explicit case are reused verbatim. The user only has to:

- Select an implicit or IMEX integrator in `user_inputs.jl`;
- Set `:llinsolve => true` when the chosen integrator requires the assembly of the implicit operator into `sem.matrix.L`;
- Set `:bdy_fluxes => true` (and optionally `:bulk_fluxes => true` for the built-in bulk surface flux, or `tag = "MOST"` for Monin-Obukhov) if surface fluxes have to be added to the RHS at every stage.

---

## 5. Summary

| Problem class | Driver branch | Kernel entry point | Action |
|---|---|---|---|
| Hyperbolic, explicit | `time_loop!` | `apply_boundary_conditions_dirichlet!` | Overwrite `uaux`, zero RHS row |
| | | `apply_boundary_conditions_neumann!` | Add surface flux to RHS |
| Elliptic (*Ax = b*) | `llinsolve` branch | `apply_boundary_conditions_lin_solve!` | Zero *L* row, *L_ii = 1*, write RHS |
| Hyperbolic, IMEX/impl. | `time_loop!` + implicit stage | `apply_boundary_conditions_lin_solve!` | Same as elliptic on stage matrix |
| | | `apply_boundary_conditions_neumann!` | Add surface flux to RHS |

---

## 6. Files referenced

- `src/kernel/boundaryconditions/BCs.jl`
- `src/kernel/boundaryconditions/custom_bcs.jl`
- `src/kernel/boundaryconditions/surface_integral.jl`
- `problems/drivers.jl`
- `problems/CompEuler/theta/user_bc.jl`, `user_inputs.jl`
- `problems/Elliptic/case1/user_bc.jl`, `user_inputs.jl`
