# Adding a New Test (Equation Set + Case)

Jexpresso discovers equation sets and cases purely by **directory structure** —
there is no registration file or config to update. To add a new test you only
need to create a folder and populate it with six Julia files.

## Directory layout

```
problems/
└── YourEquationSet/          ← This name is arbitrary
    └── your_case/            ← This name is arbitrary
        ├── user_inputs.jl     ← solver parameters (time step, mesh, output, …)
        ├── initialize.jl      ← initial conditions
        ├── user_flux.jl       ← flux vectors F and G
        ├── user_source.jl     ← source terms
        ├── user_bc.jl         ← boundary conditions
        └── user_primitives.jl ← conversion from conserved to output variables
```

Create the folder:

```bash
mkdir -p problems/YourEquationSet/your_case
```

Then create each of the six files below. The **2D Euler (`CompEuler/theta`)**
case is used as the reference throughout.

---

## File 1 — `user_inputs.jl`

Returns a `Dict` with all solver parameters. Copy this template and adjust the
values.

```julia
function user_inputs()
    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),  # time integrator
        :Δt                   => 0.5,                     # time step (seconds)
        :tinit                => 0.0,
        :tend                 => 1000.0,                  # end time (seconds)
        :diagnostics_at_times => (0:10:1000),
        :nop                  => 4,                       # polynomial order
        :interpolation_nodes  => "lgl",
        :lread_gmsh           => true,
        :gmsh_filename        => "./meshes/gmsh_grids/hexa_TFI_10x10.msh",
        :lsource              => true,
        :lvisc                => false,
        :SOL_VARS_TYPE        => TOTAL(),                 # or PERT()
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :output_dir           => "./output",
    )
    return inputs
end
```

For 1D problems set `:lread_gmsh => false` and add `:nelx`, `:xmin`, `:xmax`
instead of `:gmsh_filename`.

---

## File 2 — `initialize.jl`

Declares the solution variables and sets their initial values at every mesh
point.

The 2D Euler case (`CompEuler/theta`) solves for `[ρ, ρu, ρv, ρθ]` and
initialises a rising warm bubble:

```julia
function initialize(SD::NSD_2D, PT, mesh::St_mesh, inputs::Dict, OUTPUT_DIR::String, TFloat)

    # Declare the conserved variables (length = number of equations)
    qvars = ["ρ", "ρu", "ρv", "ρθ"]
    q = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, qvars, TFloat, inputs[:backend]; neqs=length(qvars))

    PhysConst = PhysicalConst{Float64}()
    θref = 300.0  # K  — reference potential temperature
    θc   =   2.0  # K  — warm bubble amplitude
    r0   = 2000.0 # m  — bubble radius
    xc   = (maximum(mesh.x) + minimum(mesh.x)) / 2
    yc   = 2500.0 # m

    for ip = 1:mesh.npoin
        x, y = mesh.x[ip], mesh.y[ip]
        r  = sqrt((x - xc)^2 + (y - yc)^2)
        Δθ = r < r0 ? θc * (1.0 - r/r0) : 0.0
        θ  = θref + Δθ

        p    = PhysConst.pref * (1.0 - PhysConst.g*y / (PhysConst.cp*θ))^PhysConst.cpoverR
        ρ    = perfectGasLaw_θPtoρ(PhysConst; θ=θ, Press=p)

        q.qn[ip, 1] = ρ        # density
        q.qn[ip, 2] = 0.0      # x-momentum (fluid at rest)
        q.qn[ip, 3] = 0.0      # y-momentum
        q.qn[ip, 4] = ρ * θ    # potential-temperature density
    end

    return q
end
```

Change `qvars`, the loop body, and `NSD_2D` → `NSD_1D` / `NSD_3D` to match your
equations.

---

## File 3 — `user_flux.jl`

Computes the hyperbolic flux vectors **F** (x-direction) and **G** (y-direction)
at a single quadrature point. The solver calls this inside its element loop.

For 2D Euler the fluxes are:

```julia
function user_flux!(F, G, SD::NSD_2D, q, qe, mesh::St_mesh, ::CL, ::TOTAL; neqs=4, ip=1)

    PhysConst = PhysicalConst{Float64}()
    ρ  = q[1];  ρu = q[2];  ρv = q[3];  ρθ = q[4]
    u  = ρu/ρ;  v  = ρv/ρ;  θ  = ρθ/ρ
    P  = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)

    # x-direction flux
    F[1] = ρu;        F[2] = ρu*u + P;  F[3] = ρv*u;  F[4] = ρθ*u
    # y-direction flux
    G[1] = ρv;        G[2] = ρu*v;      G[3] = ρv*v + P;  G[4] = ρθ*v
end
```

Rules:

- Match the dispatch tags to your dimensionality (`NSD_1D`, `NSD_2D`, `NSD_3D`)
  and variable type (`::TOTAL` or `::PERT`).
- For 1D problems, `G` is unused — just leave it empty.
- The length of `F` and `G` must equal the number of equations.

---

## File 4 — `user_source.jl`

Sets the right-hand-side source term `S` at a single quadrature point. For 2D
Euler the only source is gravity:

```julia
function user_source!(S, q, qe, npoin::TInt, ::CL, ::TOTAL;
                      neqs=4, x=0.0, y=0.0, ymin=0.0, ymax=0.0, xmin=0.0, xmax=0.0)

    PhysConst = PhysicalConst{Float64}()
    ρ = q[1]

    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ * PhysConst.g   # gravitational acceleration
    S[4] = 0.0
end
```

If your equations have no source terms, set every `S[i] = 0.0` and return.

---

## File 5 — `user_bc.jl`

Applies boundary conditions. For 2D Euler, a free-slip wall zeroes the normal
velocity component:

```julia
function user_bc_dirichlet!(q, coords, t::AbstractFloat, tag::String,
                             qbdy::AbstractArray, nx, ny, qe, ::TOTAL)
    # Project out the normal velocity (free-slip)
    qnl    = nx*q[2] + ny*q[3]
    qbdy[2] = q[2] - qnl*nx
    qbdy[3] = q[3] - qnl*ny
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords,
                          t::AbstractFloat, tag::String, inputs::Dict)
    return zeros(size(q, 2), 1)
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, coords,
                          t::AbstractFloat, inputs::Dict)
    return zeros(size(q, 2), 1)
end
```

Both `user_bc_neumann` overloads must be present even if you only use Dirichlet
conditions (the solver calls both signatures).

---

## File 6 — `user_primitives.jl`

Converts the conserved variables stored in `q` to the physical (primitive)
variables written to the output files.

For 2D Euler (`[ρ, ρu, ρv, ρθ]` → `[ρ, u, v, θ]`):

```julia
function user_primitives!(u, qe, uprimitive, ::TOTAL)
    uprimitive[1] = u[1]          # ρ
    uprimitive[2] = u[2] / u[1]  # u  = ρu / ρ
    uprimitive[3] = u[3] / u[1]  # v  = ρv / ρ
    uprimitive[4] = u[4] / u[1]  # θ  = ρθ / ρ
end

function user_uout!(ip, ::TOTAL, uout, u, qe; kwargs...)
    uout[1] = u[1]
    uout[2] = u[2] / u[1]
    uout[3] = u[3] / u[1]
    uout[4] = u[4] / u[1]
end
```

If your variables are already primitive (e.g., a scalar advection problem), just
copy them through unchanged.

---

## Running your new case

```bash
julia --project=. -e 'push!(empty!(ARGS), "YourEquationSet", "your_case"); include("./src/Jexpresso.jl")'
```

Or interactively:

```julia
julia> push!(empty!(ARGS), "YourEquationSet", "your_case")
julia> include("./src/Jexpresso.jl")
```

No other changes to the codebase are needed.

This is **not** required for the solver to find and run your case.

---

## Summary

| File | What it defines |
|---|---|
| `user_inputs.jl` | Solver parameters: time step, mesh, output format, … |
| `initialize.jl` | Variable names and initial values at every mesh point |
| `user_flux.jl` | Flux vectors F(q) and G(q) |
| `user_source.jl` | Source term S(q) |
| `user_bc.jl` | Boundary conditions |
| `user_primitives.jl` | Conversion from conserved to output variables |