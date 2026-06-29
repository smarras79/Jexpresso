# SymbolicFD.jl

A small, **stand-alone** "write-the-equation-and-solve-it" engine for Jexpresso.

You write a PDE using Julia's unicode characters (a few common LaTeX spellings
are accepted too), describe the grid with a `user_inputs()`-style `Dict` exactly
as in Jexpresso, and the equation is **discretized and solved automatically** —
time-marched if it has a `∂q/∂t` term, or solved as a steady (elliptic) problem
if it does not.

> **No AI, no PDF parsing, no string-to-case lookup.** This tool is *not*
> `tools/EquationGenerator.jl` (which uses Claude to generate a problem
> directory from a PDF). Here, **each differential operator carries its own
> numerical discretization** and equations are solved by composing those
> operators — see *How it works*.

## The idea (operator algebra, à la Gridap.jl)

`∇`, `∇⋅` and `∇²` are implemented as **standalone, composable operators that
act on discrete `Field`s**, not as patterns that map a whole string to a fixed
equation. This is the same spirit as Gridap.jl's differential operators, but in
*strong form by finite differences* rather than weak form:

| operator | meaning                | rank        | discretization (1D)              |
|----------|------------------------|-------------|----------------------------------|
| `∇f`     | gradient               | r → r+1     | central `∂/∂xᵢ` per direction    |
| `∇⋅F`    | divergence             | r → r−1     | `Σᵢ ∂Fᵢ/∂xᵢ` (central)           |
| `∇²f`,`Δf` | Laplacian            | r → r       | compact `Σᵢ ∂²/∂xᵢ²`             |

The equation is parsed into an **expression tree** whose nodes are these
operators and the field algebra (`+ − * ⋅`); evaluating the tree on the current
state field produces `∂q/∂t`. Because the operators are generic, *any*
composition works — `∇⋅(u q)`, `∇⋅(D∇q)`, `u⋅∇q`, `∇²q`, … — not just
advection-diffusion. Dimensional mistakes (e.g. adding a scalar to a vector) are
caught automatically by the rank bookkeeping.

## Discretization method (`:method`, user-selectable)

The numerical method is a pluggable backend. Everything above the
directional-derivative primitive — the operators, the parser, the RHS / steady
assembly, the time loop and the matrix-free Krylov solve — is method-agnostic;
only `deriv1`/`deriv2` change with the backend, dispatched on the mesh:

| `:method`  | backend                                            | status         |
|------------|----------------------------------------------------|----------------|
| `:fd`      | 2nd-order central finite differences (this file)   | **working** (default) |
| `:sem`     | spectral element, reusing Jexpresso's basis        | **working** (needs Jexpresso) |

The SEM backend does **not** reimplement the spectral basis: it calls Jexpresso's
`basis_structs_ξ_ω!` / `build_Interpolation_basis!`
(`src/kernel/bases/basis_structs.jl`) for the LGL nodes/weights `ξ`,`ω` and the
`dψ` differentiation matrix, builds an element LGL grid, and forms the derivative
with the **same weak-form structure as `rhs.jl`** (`_expansion_inviscid!` +
`DSS_rhs!` + `Minv`):

```
(∂f/∂x)_i = Minv_i · DSS( ω_i · Σ_k dψ[k,i] f_k )
```

i.e. the `dF/dξ|_i = Σ_k dψ[k,i] f_k` kernel, `ω`-weighted, direct-stiffness
summed across elements, and divided by the assembled lumped mass (the affine
metric cancels). The Jexpresso package is loaded lazily, only when
`:method => :sem` is selected, so the FD path keeps its light footprint. Choose
the order with `:nop` and the number of elements with `:nelx`.

## Two dimensions (`:nsd => 2`)

The operator layer is written generically over `nsd`, so 2D needs only a 2D mesh
and the directional-derivative primitive. Three mesh paths are available:

| `:method` | mesh                                    | `deriv1` |
|-----------|-----------------------------------------|----------|
| `:fd`     | structured Cartesian grid               | 2nd-order central in x and y |
| `:sem`    | structured Cartesian grid (no gmsh)     | tensor-product LGL, weak form |
| `:sem` + `:gmsh_filename` | **the existing Jexpresso gmsh grid** | weak form on the read mesh + metric terms |

The gmsh path reads the **same grid the Jexpresso problem uses** — for
`problems/AdvDiff/kopriva` that is `meshes/gmsh_grids/kopriva_periodic.msh` — by
calling Jexpresso's own `sem_setup` (`src/kernel/infrastructure/sem_setup.jl`)
inside a `with_mpi` block. We take `connijk`, the node coordinates, the metric
terms `dξdx,…,Je` and the `dψ`/`ω` basis straight from the returned `sem` bundle
and apply the weak directional derivative exactly as `rhs.jl`'s
`_expansion_inviscid!` (CL, Inexact, ContGal, 2D):

```
∂f/∂ξ|ij = Σ_k dψ[k,i] f[k,j]        ∂f/∂η|ij = Σ_k dψ[k,j] f[i,k]
∂f/∂x    = ∂f/∂ξ·dξdx + ∂f/∂η·dηdx
(∂f/∂x)_ip = Minv_ip · DSS( ω_i ω_j Je_ij · ∂f/∂x|ij )
```

Because `DSS` and `Minv` are linear, `∇⋅(F,G) = ∂F/∂x + ∂G/∂y` reproduces the
fused divergence kernel exactly. `tools/SymbolicFD.jl/run_advdiff_2d.jl` is the
2D analogue of `run_gaussian_1d.jl`, set up like `AdvDiff/kopriva` (gaussian blob
at `(0,3)`, `u=(0.5,1.0)`, `μ=0.1`, periodic, read from the kopriva grid).

Because `DSS` and `Minv` are linear, `∇⋅(F,G) = ∂F/∂x + ∂G/∂y` reproduces the
fused divergence kernel exactly.

## Systems of equations — the Euler θ case (`problems/CompEuler/theta`)

The operators are per-scalar-field; a **system** is just several residual
equations sharing a mesh, each carrying its own `∂t(·)` unknown, with fluxes that
reference the *other* unknowns. The canonical example is the 2D compressible
Euler equations in **potential-temperature (θ) form** — the rising thermal bubble
of `problems/CompEuler/theta`:

```
∂ρ /∂t + ∂x(ρu)        + ∂y(ρv)        = 0
∂ρu/∂t + ∂x(ρu·u + p′) + ∂y(ρu·v)      = 0
∂ρv/∂t + ∂x(ρv·u)      + ∂y(ρv·v + p′) = -ρ′ g
∂ρθ/∂t + ∂x(ρθ·u)      + ∂y(ρθ·v)      = 0
        u = ρu/ρ,  v = ρv/ρ,  p = C₀ (ρθ)^γ   (perfect-gas law for θ)
```

Two small additions to the DSL make this writable verbatim, and they are the
*natural* flux-form primitives (no equation templates, still pure operator
algebra):

| symbol      | node    | meaning                                  |
|-------------|---------|------------------------------------------|
| `∂x`, `∂y`  | `:dx/:dy` | directional first derivative `∂/∂xᵢ` (the same primitive `∇⋅` uses, exposed per component) |
| `a/b`       | `:fdiv` | elementwise field division (`u = ρu/ρ`)  |
| `a^γ`       | `:fpow` | elementwise power by a constant (`(ρθ)^γ`) |

Pass a **vector of residual equations** to `solve`; each must have exactly one
`∂t(·)` term naming its unknown. As in Jexpresso's `theta`
(`:SOL_VARS_TYPE => PERT()`), we evolve the **perturbation** about a hydrostatic
background ρ̄(y), (ρθ)‾(y), p̄(y) (with `∂p̄/∂y = −ρ̄g`); carrying the pressure
perturbation `p′ = p − p̄` and the buoyancy source `−ρ′g` makes the scheme
**discretely well balanced** — the resting background is preserved to machine
precision regardless of the stencil (verified in `runtests.jl`).

```julia
@vars ρ ρu ρv ρθ ρb ρθb pb ν              # unknowns are the perturbations
u  = ρu/(ρb + ρ);  v = ρv/(ρb + ρ)
p′ = C0*(ρθb + ρθ)^γ - pb                  # γ, C0, g are ordinary numbers
eqs = [ ∂t(ρ)  + ∂x(ρu)           + ∂y(ρv),
        ∂t(ρu) + ∂x(ρu*u + p′)    + ∂y(ρu*v)         - ν*Δ(ρu),
        ∂t(ρv) + ∂x(ρv*u)         + ∂y(ρv*v + p′)    - ν*Δ(ρv)  + ρ*g,
        ∂t(ρθ) + ∂x((ρθb+ρθ)*u)   + ∂y((ρθb+ρθ)*v)   - ν*Δ(ρθ) ]
mesh, Q0, Q = SymbolicFD.solve(eqs, inputs)   # ⇒ Vector of component fields
```

`tools/SymbolicFD.jl/run_euler_theta_2d.jl` is the system analogue of
`run_advdiff_2d.jl`: it sets up the rising thermal bubble exactly as
`problems/CompEuler/theta/initialize.jl` (warm anomaly `Δθ = θc(1−r/r0)`,
`θref = 300 K`, `θc = 2 K`, at rest on the hydrostatic background) and time-marches
it with the same RK4 driver. The small constant viscosity `ν` plays the role of
Jexpresso's `:lvisc => true, :μ` artificial viscosity. Extra system inputs:

| key          | meaning                                                       |
|--------------|---------------------------------------------------------------|
| `:q0`        | initial state `(x,y) -> (u₁,…,u_neqs)` (one tuple per node)    |
| `:Δt`        | fixed step; else CFL·Δx / `:wave_speed` (default 350 m/s)      |
| `:diag`      | `(Q, mesh) -> Vector` plotted/printed (default `Q[end]`)       |
| `:diag_name` | label for that diagnostic (default the last unknown's name)   |

```bash
julia --project=. tools/SymbolicFD.jl/run_euler_theta_2d.jl
```

defaults to the self-contained structured-FD backend; flip `:method => :sem`
with a `:gmsh_filename` to run the very same equations on a spectral-element gmsh
grid through Jexpresso's `sem_setup` (the `∂x/∂y` then become the weak-form
derivative on the read mesh). `problems/CompEuler/sod1d`'s **DSGS** residual-based
shock capturing (`src/kernel/physics/SGS.jl`) remains the next follow-up.

## Example — write the equation as live symbols (no string)

```julia
include("src/SymbolicFD.jl"); using .SymbolicFD

@vars q u μ                                  # declare the symbols
equation = ∂t(q) + ∇⋅(u*q) - μ*∇⋅∇(q)        # residual form, " = 0" implied

inputs = Dict(
    :nsd  => 1, :xmin => -1.0, :xmax => 1.0, :npoin => 200, :periodic => true,
    :u    => [1.0],      # advecting velocity (a vector; ∇⋅ contracts it)
    :μ    => 1.0e-3,     # small diffusion coefficient
    :q0   => x -> exp(-(x^2)/(2*0.1^2)),   # gaussian wave
    :tend => 2.0, :CFL => 0.4,
)

mesh, q0, q = SymbolicFD.solve(equation, inputs)
```

`∇`, `Δ`, `∂t` are real objects and `+ − * ⋅` are overloaded on the expression
tree, so `equation` *is* the AST — no parsing of a string. Two Julia-lexing
notes: write **`∂t(q)`** (the literal `∂q/∂t` lexes as two identifiers `∂q`,`∂t`)
and put a **`*`** after a coefficient (`μ*∇⋅∇(q)`, since `μ∇` would lex as one
identifier). The unknown is the symbol carrying `∂t(·)` (or, for steady problems,
the one symbol absent from `inputs`); the rest are resolved from `inputs` as
scalars, vectors, or functions `x->…`.

A plain **string** equation still works too (handy for `=`-form and LaTeX):

```julia
SymbolicFD.solve("∂q/∂t + ∇⋅(\\mathbf{u}q) = \\mu∇⋅∇(q)", inputs)
```

Either way the run banner prints the **discretized RHS it built**, e.g.

```
   discretized RHS : ∂q/∂t = -∇⋅([1]·q) + 0.001·∇²(q)
```

## Running it

`Plots` is the only dependency and it is already a Jexpresso dependency, so run
inside the project environment. From the repository root:

```bash
julia --project=. tools/SymbolicFD.jl/run_gaussian_1d.jl
```

This advects a gaussian wave once around the periodic domain `[-1, 1]` while a
small `μ` diffuses it slightly. The solution is **saved to PNG and plotted on the
fly in exactly the same way as `problems/CompEuler/sod1d`** — using `Plots.jl`
(GR backend), an `INIT-q.png` at `t = 0` and a `fields-it<iout>.png` at every
diagnostic output, with the GR window updating live during the run. Files land
in `tools/SymbolicFD.jl/output/`. A CSV is always written, and a terminal ASCII
plot is available as a display-free fallback (`:outformat => "ascii"`). On a
headless machine set `:plot_live => false` (PNGs are still written).

You can also drive it from the Julia REPL:

```julia
include("tools/SymbolicFD.jl/src/SymbolicFD.jl"); using .SymbolicFD
SymbolicFD.solve("∂q/∂t + ∇⋅(u q) = μ∇²q",
                 Dict(:u => [1.0], :μ => 1e-3, :tend => 2.0))
```

## Time-independent (steady) problems

If the equation has **no `∂q/∂t` term**, it is solved as a steady elliptic
problem rather than time-marched, and the solve is **matrix-free**: no matrix is
ever stored (not even sparse). For a linear PDE the residual is affine,
`residual(q) = A q + c`, so a matrix-vector product is just
`A·x = residual(x) − residual(0)` — i.e. one application of the same finite-
difference operators. That product (with identity rows at the two Dirichlet
boundary nodes) is wrapped as a `LinearOperator` and handed to Krylov.jl's GMRES.
Tune with `:ksp_rtol`, `:ksp_atol`, `:ksp_memory` (default full GMRES).

```bash
julia --project=. tools/SymbolicFD.jl/run_heat_steady_1d.jl
```

solves the steady heat (Poisson) equation `∇⋅∇(q) = f` on `[-1, 1]` with a
localized source and cold ends `q(±1) = 0`. The unknown is inferred as the only
symbol that is not a supplied parameter (override with `:unknown => "q"`). Steady
problems use a non-periodic grid and the Dirichlet keys `:bc_left` / `:bc_right`
(default `0.0`). The banner reports the steady residual it assembled, e.g.

```
   steady residual : ∇²(q) - f(x) = 0   (solve for q)
```

Examples: `∇²q = 0` (Laplace ⇒ linear profile between the end values),
`μ∇²q = f` (Poisson). The direct solve assumes a **linear** operator; nonlinear
steady problems are out of scope for now.

## Tests

Analytic-solution and operator checks live in `runtests.jl`:

```bash
julia --project=. tools/SymbolicFD.jl/runtests.jl
```

They verify the operators against exact derivatives (with 2nd-order
convergence), parsing and rank bookkeeping (dimensional / unknown-symbol
errors), the heat equation against its analytic decay
`q = e^{-μπ²t} sin(πx)`, periodic advection returning to itself after one
period with machine-precision mass conservation, and an end-to-end `solve`
smoke test.

## Changing the equation

Edit **two** things in `run_gaussian_1d.jl` — the `equation` string and the
matching parameters inside `user_inputs()`. Every symbol you use (`μ`, `u`, …)
must have a matching `:symbol => value` entry, otherwise the solver stops with a
clear error naming the missing `:key`. Use a **vector** for a vector quantity
(`:u => [1.0]`) and a scalar for a scalar; in 1D a scalar velocity is also
accepted as a one-component flux.

Worked examples (drop the string into `equation`, adjust inputs):

| What you want                       | `equation =` …                              | extra inputs            |
|-------------------------------------|---------------------------------------------|-------------------------|
| Pure advection (no diffusion)       | `"∂q/∂t + ∇⋅(u q) = 0"`                      | `:u => [1.0]`           |
| Conservative advection–diffusion    | `"∂q/∂t + ∇⋅(\\mathbf{u}q) = \\mu∇⋅∇(q)"`    | `:u => [1.0], :μ`       |
| Non-conservative form               | `"∂q/∂t + u⋅∇q = μ∇²q"`                      | `:u => [1.0], :μ`       |
| Transient diffusion (heat eqn)      | `"∂q/∂t = μ∇²q"`                             | `:μ => 0.01`            |
| Advection–diffusion–decay           | `"∂q/∂t + ∇⋅(u q) = μ∇²q - k q"`             | `:u,:μ,:k`              |
| **Steady** heat / Poisson (no ∂/∂t) | `"∇²q = f"` / `"μ∇²q = f"`                   | `:f, :bc_left,:bc_right`|
| **Steady** Laplace                  | `"∇²q = 0"`                                  | `:bc_left,:bc_right`    |

Change the initial shape with `:q0`, e.g. a narrower gaussian
`x -> exp(-(x^2)/(2*0.05^2))` or a square pulse `x -> abs(x) < 0.3 ? 1.0 : 0.0`.

### Notation accepted

- time derivative `∂q/∂t`;
- gradient `∇q`, divergence `∇⋅(…)`, Laplacian `∇²q` / `Δq` / `∇⋅∇(q)`;
- directional first derivatives `∂x(…)`, `∂y(…)` (symbolic DSL only — the
  flux-form primitive for systems, see the Euler θ case above);
- elementwise field division `a/b` and power `a^γ` (symbolic DSL only — for
  closures like `u = ρu/ρ` and the gas law `(ρθ)^γ`);
- products by juxtaposition (`uq` = `u*q`) or `*`; inner product `⋅`
  (`u⋅∇q`); sums/differences `+ −`;
- single-character variable/parameter names (so `uq` reads unambiguously);
- LaTeX spellings `\nabla`, `\cdot`, `\mu`, `\Delta`, `\mathbf{...}`, etc. are
  translated automatically.

Note: the *contiguous* `∇⋅∇` (and `∇²`, `Δ`) become the **compact** Laplacian
stencil. An explicit `∇⋅(∇q)` is kept as the genuine composition of two
first-difference operators (a wider stencil) — operators are never silently
re-interpreted.

## Inputs reference

| key            | meaning                                         | default          |
|----------------|-------------------------------------------------|------------------|
| `:method`      | `:fd` finite differences / `:sem` spectral elem (next) | `:fd`     |
| `:nop`         | polynomial order per element (`:method => :sem`)| `4`              |
| `:nsd`         | spatial dimension (only `1` implemented so far) | `1`              |
| `:xmin`,`:xmax`| domain                                          | `-1`, `1`        |
| `:npoin`/`:nelx`| number of grid points                          | `100`            |
| `:periodic`    | periodic boundary conditions                    | `true`           |
| `:u`,`:μ`,…    | parameters (vector for a vector, `x->…` for a field) | —           |
| `:q0`          | initial condition `x -> value` (transient)      | gaussian         |
| `:tend`        | final time (transient)                          | `1.0`            |
| `:Δt`          | fixed time step (otherwise CFL-derived)         | auto             |
| `:CFL`         | CFL number for the automatic `Δt`               | `0.5`            |
| `:bc_left`,`:bc_right` | Dirichlet boundary values (steady)      | `0.0`            |
| `:unknown`     | unknown name if it cannot be inferred (steady)  | inferred         |
| `:output_dir`  | directory for figures / `solution.csv`          | `"."`            |
| `:outformat`   | `"png"` (Plots.jl, like sod1d) or `"ascii"`     | `"png"`          |
| `:ndiagnostics_outputs` | number of on-the-fly plot snapshots    | `10`             |
| `:plot_live`   | update an on-screen GR window during the run    | `true`           |

## How it works (pipeline)

1. **Build the tree** — either *symbolic input* (`∇`/`Δ`/`∂t` + overloaded
   `+ − * ⋅` on `Node`s) constructs the expression tree directly, or a *string*
   is normalized (strip `\mathbf{}`, map LaTeX→unicode, collapse `∇⋅∇`/`Δ`→`∇²`)
   and parsed by a recursive-descent parser into the same tree. No equation
   templates either way.
2. **Resolve symbols** — the unknown (the `∂t(·)` argument, or the lone non-input
   symbol) becomes the field; other symbols are resolved from `inputs` as
   scalars / vectors / functions of `x`.
3. **Build the RHS / residual** — if a `∂q/∂t` term is present it is removed and
   the rest moved across `=` to form the tree for `∂q/∂t = …` (transient);
   otherwise the tree is the residual `lhs - rhs` to drive to zero (steady).
4. **Discretize + evaluate** — each `∇/∇⋅/∇²` node calls its finite-difference
   stencil on the field it receives; field algebra composes them.
5. **Solve** — both modes are matrix-free (only the operator is ever applied,
   never stored). Transient: explicit RK4 with an automatic CFL-based `Δt`.
   Steady: wrap `A·x = residual(x) − residual(0)` (plus Dirichlet boundary rows)
   as a `LinearOperator` and solve with Krylov.jl's GMRES.
6. **Output** — PNG + on-the-fly plots like `sod1d`, plus CSV.

## Path to Jexpresso integration

The design intentionally mirrors Jexpresso so it can be folded in later:

- the `Dict`-based `user_inputs()` and the `:outformat`/`:ndiagnostics_outputs`/
  `:output_dir` keys match Jexpresso's problem inputs;
- the PNG / on-the-fly plotting reproduces `src/io/plotting/jeplots.jl`
  (`plot_initial`, `plot_results`, `render_plot_matrix`);
- the **operator layer** (`Field` + `∇/∇⋅/∇²` over `nsd` directions through one
  directional-derivative primitive) is where Jexpresso's own kernels can plug
  in; only a multi-dimensional mesh + that primitive are needed for 2D/3D;
- `:nsd`, `:xmin/:xmax`, `:npoin/:nelx`, `:periodic` reuse Jexpresso's 1D grid
  vocabulary.
