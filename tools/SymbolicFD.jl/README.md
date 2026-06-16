# SymbolicFD.jl

A small, **stand-alone** "write-the-equation-and-solve-it" engine for Jexpresso.

You write a PDE using Julia's unicode characters (a few common LaTeX spellings
are accepted too), describe the grid with a `user_inputs()`-style `Dict` exactly
as in Jexpresso, and the equation is **parsed, finite-difference discretized and
integrated in time automatically**.

The solution is **saved to PNG and plotted on the fly in exactly the same way as
`problems/CompEuler/sod1d`** — using `Plots.jl` (GR backend), an `INIT-<var>.png`
at `t = 0` and a `fields-it<iout>.png` at every diagnostic output, with the GR
window updating live during the run. `Plots` is the only dependency (already part
of Jexpresso's `Project.toml`); a CSV dump is always written, and a terminal
ASCII plot is available as a display-free fallback (`:outformat => "ascii"`).

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
`Δt` and mass/peak diagnostics, writes `output/INIT-q.png`, a
`output/fields-it<iout>.png` per diagnostic output (live-updated on screen) and
`output/solution.csv`.

## Running it

Run it inside the Jexpresso project environment so `Plots` resolves (it is
already a Jexpresso dependency). From the repository root:

```bash
julia --project=. tools/SymbolicFD.jl/run_gaussian_1d.jl
```

PNGs land in `tools/SymbolicFD.jl/output/` (`INIT-q.png`, `fields-it0.png` …
`fields-it20.png`) and, with a display attached, a GR window updates on the fly.
On a headless machine set `:plot_live => false` (PNGs are still written) or
`:outformat => "ascii"` for a terminal plot only.

You can also drive it from the Julia REPL:

```julia
include("tools/SymbolicFD.jl/src/SymbolicFD.jl"); using .SymbolicFD
SymbolicFD.solve("∂q/∂t + ∇⋅(u q) = μ∇²q",
                 Dict(:u => 1.0, :μ => 1e-3, :tend => 2.0))   # other keys take defaults
```

## Changing the equation

Open `run_gaussian_1d.jl` and edit **two** things — the `equation` string and the
matching parameters/grid inside `user_inputs()`:

1. **The equation** (line `equation = "…"`). Write it in unicode or LaTeX. Every
   parameter symbol you use (e.g. `μ`, `u`, `k`) must have a matching entry in
   the inputs Dict, otherwise the solver stops with a clear error telling you
   which `:symbol` to add.

2. **The inputs** — grid (`:npoin`, `:xmin`, `:xmax`), parameter values, the
   initial condition `:q0`, and the run length `:tend`.

Worked examples (drop the string into `equation` and adjust inputs):

| What you want                       | `equation =` …                              | extra inputs            |
|-------------------------------------|---------------------------------------------|-------------------------|
| Pure advection (no diffusion)       | `"∂q/∂t + ∇⋅(u q) = 0"`                      | `:u => 1.0`             |
| Stronger diffusion                  | `"∂q/∂t + ∇⋅(\\mathbf{u}q) = \\mu∇⋅∇(q)"`    | `:μ => 1e-2`            |
| Pure diffusion (heat equation)      | `"∂q/∂t = μ∇²q"`                             | `:μ => 0.01`            |
| Advection–diffusion–decay           | `"∂q/∂t + ∇⋅(u q) = μ∇²q - k q"`             | `:u,:μ,:k`              |
| Faster flow on a finer grid         | `"∂q/∂t + ∇⋅(u q) = μ∇²q"`                   | `:u => 2.0, :npoin => 400` |

Change the initial shape by editing `:q0`, e.g. a narrower gaussian
`x -> exp(-(x^2)/(2*0.05^2))` or a square pulse `x -> abs(x) < 0.3 ? 1.0 : 0.0`.

Accepted notation for each operator: `∂q/∂t`; advection `∇⋅(u q)` /
`∇⋅(\mathbf{u}q)`; diffusion `μ∇⋅∇(q)` / `μ∇²q` / `μΔq`; linear reaction `k q`.
LaTeX spellings (`\nabla`, `\cdot`, `\mu`, `\Delta`, `\mathbf{...}`) are translated
automatically. See the pattern table below.

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
| `:output_dir`  | directory for figures / `solution.csv`          | `"."`            |
| `:outformat`   | `"png"` (Plots.jl, like sod1d) or `"ascii"`     | `"png"`          |
| `:ndiagnostics_outputs` | number of on-the-fly plot snapshots    | `10`             |
| `:plot_live`   | update an on-screen GR window during the run    | `true`           |

## Path to Jexpresso integration

The design intentionally mirrors Jexpresso so it can be folded in later:

- the `Dict`-based `user_inputs()` matches Jexpresso's problem inputs, including
  the `:outformat`, `:ndiagnostics_outputs` and `:output_dir` output keys;
- the PNG / on-the-fly plotting reproduces `src/io/plotting/jeplots.jl`
  (`plot_initial`, `plot_results`, `render_plot_matrix`) so it can be swapped for
  the real `write_output`/`jeplots` calls directly;
- `:nsd`, `:xmin/:xmax`, `:npoin/:nelx`, `:periodic` reuse Jexpresso's 1D grid
  vocabulary;
- the parser emits typed `PDETerm`s, which is the natural place to plug into
  Jexpresso's `user_flux!` / `user_source!` / `CL()`-`PERT()` machinery instead
  of the built-in finite-difference operators;
- extending `FDMesh1D`/operators to 2D–3D and adding upwind/higher-order stencils
  is isolated in sections 5–7 of `src/SymbolicFD.jl`.
