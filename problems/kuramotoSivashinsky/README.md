# kuramotoSivashinsky

The Kuramoto-Sivashinsky equation for a front φ(x,t),

```
∂ₜφ + ∂ₓ²φ + ∂ₓ⁴φ + ½(∂ₓφ)² = 0
```

is solved in conservation form, for the slope u = ∂ₓφ:

```
∂ₜu + ∂ₓ(u²/2) + ∂ₓ²u + ∂ₓ⁴u = 0
```

∂ₓ²u drives the long-wave instability and ∂ₓ⁴u damps the short waves. A
mode e^{ikx} grows at σ = k² - k⁴, so every 0 < k < 1 is unstable.

## Discretisation and time integration

The case runs on Jexpresso's continuous spectral elements (LGL, lumped
mass). It is split as du/dt = L u + N(u):

```
L u  = -∂ₓ²u - ∂ₓ⁴u      the stiff linear part
N(u) = -∂ₓ(u²/2)         the flux of user_flux.jl, through the usual rhs!
```

L is declared in the deck (`:lsplit_linear => true`, `:split_d2 => [-1.0]`,
`:split_d4 => [-1.0]`). `kernel/solvers/split_linear.jl` builds it on the
mesh:
- **Second derivative:** ∂ₓ² is the SEM weak Laplacian D₂ = -M⁻¹K. This is
  the same operator the viscous term applies.
- **Fourth derivative:** ∂ₓ⁴ = D₂², the C⁰ mixed (Ciarlet-Raviart) form.
- **Products:** L u uses the sparse matrix, 14 nonzeros per row.

The deck's `ETDRK4()` is OrdinaryDiffEq's exponential Runge-Kutta, the scheme
Kassam & Trefethen use. It integrates L exactly, through its φ-functions.
The time step is therefore set by the nonlinear dynamics, not by the ∂ₓ⁴
eigenvalues (about 5500 on this mesh). The step is h = 1/4, as in Kassam &
Trefethen; an explicit scheme would need Δt ≈ 2e-4.

L is similar to a symmetric matrix: L = S⁻¹(B - B²)S with S = M^{1/2} and
B = S⁻¹KS⁻¹. So one 256×256 symmetric eigendecomposition gives every φₖ(hL)
exactly, in about 10 ms. ExponentialUtilities' general dense algorithm takes
about 10 s for the same matrices.

```bash
julia --project=. src/Jexpresso.jl kuramotoSivashinsky kassamTrefethen2005
```

| case | dim | eqs | grid | boundaries |
|---|---|---|---|---|
| `kassamTrefethen2005` | 1D | 1 | built-in, 64 elements, p = 4, on [0, 32π] | periodic |

## kassamTrefethen2005

The run of Kassam & Trefethen, *SIAM J. Sci. Comput.* 26:1214 (2005),
Fig. 4.1: u(x,0) = cos(x/16)(1 + sin(x/16)) on a periodic domain of length
32π, to t = 150. It takes 600 ETDRK4 steps of Δt = 1/4. Every output time
must be a whole number of steps from :tinit, because the φ operators are
built for that one Δt. The time loop checks this and stops with an error if
it does not hold.

Output (`user_plot.jl`). The solution is sampled every Δt = 0.5:

- `ks-it<n>.png`: u, written every tenth output.
- `ks-spacetime.png`: u(x,t) over the whole run, the picture of their
  Fig. 4.1, drawn once at the end. It is drawn on the native grid, with no
  resampling onto a uniform x. Each LGL node is one column of cells, bounded
  by the midpoints to its neighbours, so the columns keep the non-uniform LGL
  widths. Each output time is one row. The colours are inferno, saturated:
  `:ks_clims` (default (-1.9, 1.0)) clips the range well inside
  max|u| ≈ 2.6, so crests turn yellow, troughs turn black, and u = 0 is
  orange.

`:ks_plot => false` switches the case's plotting hook off. `:outformat =>
"none"` switches off all output.

A few steady cells form early, the first near x ≈ 25 at t ≈ 10. From
t ≈ 45 the cells split and merge, and the run is chaotic.

### Cost

The time loop to t = 150, measured on one core with compilation excluded:

| | time | allocations |
|---|---|---|
| `:outformat => "none"` | 0.20 s | 13 k (3.3 MiB): setup and per-output bookkeeping, none per step |
| default output (301 times, `ks_plot`) | 1.7 s | the rest is the plotting |

The earlier version of this case relaxed a second field onto ∂ₓ²u and
stepped everything explicitly. It needed Δt = 2e-4, which is 750 000 steps
and about 300 s.

### Accuracy

There is no closed-form solution. For checking only, the run was compared
with an independent Fourier-spectral ETDRK4 solution, Kassam & Trefethen's
own method, at N = 512 and h = 1/64. That solution is not part of the case.

| t | 10 | 20 | 40 | 60 | 80 | 100 |
|---|---|---|---|---|---|---|
| max\|u - u_ref\|, Δt = 1/4 | 3.1e-4 | 1.9e-4 | 3.0e-4 | 9.3e-3 | 4.5e-2 | 0.10 |
| max\|u - u_ref\|, Δt = 1/8 | 3.1e-4 | 1.8e-4 | 2.6e-4 | 1.2e-3 | 7.5e-3 | 4.9e-2 |

max|u_ref| is between 1.9 and 2.6 throughout. The Δt = 1/8 row was measured
with the same SEM operators and the same OrdinaryDiffEq ETDRK4, driven
outside Jexpresso. Halving Δt changes nothing to t = 40, so up to there the
error is spatial. After that, chaos amplifies any difference: at t = 150 two
runs of the Fourier solver at different resolutions already differ by 0.12,
so past t ≈ 100 the comparison is statistical, not pointwise.
