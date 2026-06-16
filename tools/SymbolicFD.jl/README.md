# SymbolicFD.jl

A small, **stand-alone** "write-the-equation-and-solve-it" engine for Jexpresso.

You write a PDE using Julia's unicode characters (a few common LaTeX spellings
are accepted too), describe the grid with a `user_inputs()`-style `Dict` exactly
as in Jexpresso, and the equation is **parsed, finite-difference discretized and
integrated in time automatically** — no external packages required (only the
`Printf` standard library).

## Example

```julia
include("src/SymbolicFD.jl"); using .SymbolicFD

equation = "∂q/∂t + ∇⋅(\\mathbf{u}q) = \\mu∇⋅∇(q)"

inputs = Dict(
    :nsd  => 1, :xmin => -1.0, :xmax => 1.0, :npoin => 200, :periodic => true,
    :u    => 1.0,        # advecting velocity
    :μ    => 1.0e-3,     # small diffusion
    :q0   => x -> exp(-(x^2)/(2*0.1^2)),   # gaussian wave
    :tend => 2.0, :CFL => 0.4,
)

mesh, q0, q = SymbolicFD.solve(equation, inputs)
```

Or just run the bundled example:

```bash
julia run_gaussian_1d.jl
```

This advects a gaussian wave once around the periodic domain `[-1, 1]` while a
small `μ` diffuses it slightly. The solver prints the parsed terms, the chosen
`Δt`, mass/peak diagnostics, writes `output/solution.csv` (`x, q_initial,
q_final`) and draws a terminal ASCII plot of the initial and final fields.

## How it works

1. **Normalization** — `\mathbf{u}`/`\vec{u}` decorations are stripped, LaTeX
   greek/operators (`\mu`, `\nabla`, `\cdot`, `\Delta`, …) are mapped to
   unicode, and every spelling of the Laplacian (`∇⋅∇`, `∇²`, `Δ`) is collapsed
   to a single token `∇²`.
2. **Parsing** — the string is split at `=` into signed additive terms. Each
   term is classified (`∂q/∂t` time derivative, `∇⋅(u q)` advection, `μ∇²q`
   diffusion, linear reaction, constant source) and moved to the right-hand
   side of `∂q/∂t = Σ coeffᵢ Opᵢ(q)` as a typed `PDETerm`.
3. **Discretization** — 2nd-order central finite differences on a uniform
   periodic 1D grid (conservative flux form for advection).
4. **Integration** — explicit RK4 with an automatic CFL-based `Δt`
   (overridable with `:Δt`).

## Recognised equation patterns (1D, single scalar)

| Pattern (unicode / LaTeX)              | Meaning                  |
|----------------------------------------|--------------------------|
| `∂q/∂t`                                 | time derivative          |
| `∇⋅(u q)`  / `∇⋅(\mathbf{u}q)`          | advection (flux div.)    |
| `μ∇⋅∇(q)` / `μ∇²q` / `μΔq`              | diffusion                |
| `k q`                                   | linear reaction          |
| `f` (constant)                          | source                   |

Coefficients may be numbers (`0.1∇²q`) or parameter symbols looked up in the
inputs (`μ∇²q` ⇒ `inputs[:μ]`).

## Inputs reference

| key            | meaning                                         | default          |
|----------------|-------------------------------------------------|------------------|
| `:nsd`         | spatial dimension (only `1` supported for now)  | `1`              |
| `:xmin`,`:xmax`| domain                                          | `-1`, `1`        |
| `:npoin`/`:nelx`| number of grid points                          | `100`            |
| `:periodic`    | periodic boundary conditions                    | `true`           |
| `:u`,`:μ`,…    | parameter symbols used in the equation          | —                |
| `:q0`          | initial condition `x -> value`                  | gaussian         |
| `:tend`        | final time                                      | `1.0`            |
| `:Δt`          | fixed time step (otherwise CFL-derived)         | auto             |
| `:CFL`         | CFL number for the automatic `Δt`               | `0.5`            |
| `:output_dir`  | directory for `solution.csv`                    | `"."`            |
| `:plot`        | terminal ASCII plot                             | `true`           |

## Path to Jexpresso integration

The design intentionally mirrors Jexpresso so it can be folded in later:

- the `Dict`-based `user_inputs()` matches Jexpresso's problem inputs;
- `:nsd`, `:xmin/:xmax`, `:npoin/:nelx`, `:periodic` reuse Jexpresso's 1D grid
  vocabulary;
- the parser emits typed `PDETerm`s, which is the natural place to plug into
  Jexpresso's `user_flux!` / `user_source!` / `CL()`-`PERT()` machinery instead
  of the built-in finite-difference operators;
- extending `FDMesh1D`/operators to 2D–3D and adding upwind/higher-order stencils
  is isolated in sections 5–7 of `src/SymbolicFD.jl`.
