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

The case runs on Jexpresso's continuous spectral elements like every other 1D
case. u²/2 is the inviscid flux. The derivatives go through the viscous
operator, which adds `∂ₓ(μ[ieq] ∂ₓ uprimitive[ieq])` to equation `ieq`. To
reach the fourth derivative, a second field w is relaxed onto ∂ₓ²u with a
short time τ, as `../CahnHilliard` does for its chemical potential:

```
∂ₜu + ∂ₓ(u²/2) = -w - ∂ₓ²w
τ ∂ₜw          = ∂ₓ²u - w
```

| equation | `uprimitive` | `:μ` coefficient | source |
|---|---|---|---|
| 1 (u) | w | -1 | -w |
| 2 (w) | u | 1/τ | -w/τ |

As τ → 0, w → ∂ₓ²u and the u row is Kuramoto-Sivashinsky. The linearised
system obeys τs² + s + k⁴ - k² = 0. Its slow root is s ≈ σ(1 - τσ), and its
fast roots have Re s = -1/(2τ) at every k, so the system is stable.

Writing the u row as -∂ₓ²(u + w) instead looks equivalent but is not: once
τk² > 1, w can no longer follow ∂ₓ²u, and the bare anti-diffusion -∂ₓ²u
grows at Re s → k²/2 on the mesh scale.

```bash
julia --project=. src/Jexpresso.jl kuramotoSivashinsky kassamTrefethen2005
```

| case | dim | eqs | grid | boundaries |
|---|---|---|---|---|
| `kassamTrefethen2005` | 1D | 2 | built-in, 64 elements, p = 4, on [0, 32π] | periodic |

## kassamTrefethen2005

The run of Kassam & Trefethen, *SIAM J. Sci. Comput.* 26:1214 (2005),
Fig. 4.1: u(x,0) = cos(x/16)(1 + sin(x/16)) on a periodic domain of length
32π, to t = 150. w starts from the exact second derivative of u
(`ks_initial` in `user_flux.jl`). τ = 1e-4 and Δt = 2e-4. The time-step limit
is set by the fast oscillation, |Im s| ≈ λ/√τ, where λ ≈ 74 is the largest
eigenvalue of -∂ₓ² on this mesh. The run takes about 10 minutes on one core.

Output (`user_plot.jl`). The solution is sampled every Δt = 0.5:

- `ks-it<n>.png`: u and w ≈ ∂ₓ²u, written every tenth output.
- `ks-spacetime.png`: u(x,t) from 0 to the current time, the picture of
  their Fig. 4.1, also redrawn every tenth output and at the end. It is
  drawn on the native grid, with no resampling onto a uniform x. Each LGL
  node is one column of cells, bounded by the midpoints to its neighbours,
  so the columns keep the non-uniform LGL widths. Each output time is one
  row. The colours are inferno, saturated: `:ks_clims` (default
  (-1.9, 1.0)) clips the range well inside max|u| ≈ 2.6, so crests turn
  yellow, troughs turn black, and u = 0 is orange.

A few steady cells form early, the first near x ≈ 25 at t ≈ 10. From
t ≈ 45 the cells split and merge, and the run is chaotic.

There is no closed-form solution. For checking only, this run was compared
with an independent Fourier-spectral ETDRK4 solution, Kassam & Trefethen's
own method, at N = 512 and h = 1/64. That solution is not part of the case.

| t | max\|u - u_ref\| | max\|u_ref\| |
|---|---|---|
| 1 | 3.8e-7 | 1.31 |
| 10 | 3.1e-4 | 2.38 |
| 20 | 1.7e-4 | 1.88 |
| 40 | 6.6e-4 | 2.00 |
| 60 | 5.3e-3 | 2.47 |
| 80 | 9.2e-3 | 2.57 |
| 100 | 4.7e-2 | 2.45 |
| 120 | 0.48 | 2.04 |
| 150 | 2.4 | 2.32 |

The two agree to within 1e-3 through t = 40. After that, the chaos amplifies
any difference: at t = 150, two runs of the Fourier solver itself at
different resolutions already differ by 0.12. Past t ≈ 110 the comparison is
therefore statistical, not pointwise.
