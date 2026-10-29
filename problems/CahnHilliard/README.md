# CahnHilliard

The Cahn-Hilliard equation

```
∂ₜc = ∂ₓ²μ,     μ = c³ - c - ε² ∂ₓ²c
```

is fourth order. These cases write it as two second-order equations by
relaxing the chemical potential μ onto its definition with a short time τ:

```
∂ₜc   = ∂ₓ²μ
τ ∂ₜμ = c³ - c - ε² ∂ₓ²c - μ
```

As τ → 0 this gives μ = c³ - c - ε²∂ₓ²c, which is Cahn-Hilliard.

Both second derivatives go through Jexpresso's viscous operator. That
operator adds `∂ₓ(μ[ieq] ∂ₓ uprimitive[ieq])` to equation `ieq` for whatever
variable the case puts in `uprimitive[ieq]`. Here the variables are crossed
(`user_primitives.jl`):

| equation | `uprimitive` | `:μ` coefficient | term |
|---|---|---|---|
| 1 (c) | μ | 1 | ∂ₓ²μ |
| 2 (μ) | c | -ε²/τ | -(ε²/τ) ∂ₓ²c |

The source carries `(c³ - c - μ)/τ` and the inviscid flux is zero.

**Why not a hyperbolic system like `../KdV`.** The obvious first-order
relaxation of the fourth-order operator, with auxiliary fields for ∂ₓc and
-∂ₓμ, has a flux Jacobian whose eigenvalues satisfy λ⁴ = -ε²/τ³. Those are
complex, so that system is ill-posed. The system above is not. Linearised
about c = 0, a mode e^{ikx+st} obeys

```
τ s² + s + k²(ε²k² - 1) = 0
```

The slow root is s ≈ σ(1 - τσ), with σ = k²(1 - ε²k²) the Cahn-Hilliard
rate. At large k the roots have Re s = -1/(2τ) and |Im s| ≈ εk²/√τ, so every
mode decays except those Cahn-Hilliard itself makes unstable. The explicit
time step is set by that oscillation and shrinks like √τ, not like the
1/(ε²k⁴) of explicit Cahn-Hilliard.

```bash
julia --project=. src/Jexpresso.jl CahnHilliard spinodal1d
```

| case | dim | eqs | grid | boundaries |
|---|---|---|---|---|
| `spinodal1d` | 1D | 2 | built-in, 64 elements, p = 4, on [0, 2π] | periodic |

## spinodal1d

Spinodal decomposition of the unstable mixture c = 0, started from

```
c = δ cos(kx),     μ = c³ - c - ε² ∂ₓ²c
```

The parameters are δ = 1e-3, k = 4, ε = 1/(4√2) and τ = 1e-4, all set in
`ch_parameters()` in `user_flux.jl`. μ is built from the exact second
derivative of c (`ch_initial`), so the run starts on the relaxed manifold.
With this ε, k = 4 is the fastest-growing wavenumber, σ = 1/(4ε²) = 8.
Roundoff in the other modes therefore grows more slowly than the seeded one.

Until the amplitude reaches 0.05 (t ≈ 0.49), every frame overlays the exact
linear solution c = δe^{σt}cos(kx), μ = (ε²k² - 1)c (`user_analytic.jl`).
After that the cubic term takes over, and c saturates into the steady
periodic profile with the same wavelength 2π/k and μ = 0. That profile is a
periodic orbit of ε²c'' = c³ - c. Its amplitude A is fixed by requiring the
orbit's period,

```
P(A) = 4ε ∫₀^{π/2} dθ / √(1 - A²(1 + sin²θ)/2),
```

to equal 2π/k (`ch_equilibrium_amplitude`). For the shipped parameters that
gives A = 0.800781. At this ε the wavelength π/2 is only about six interface
widths, too short for separate c = ±1 plateaus to form.

`ch-it<n>.png` shows c and μ (`user_plot.jl`). Each output time also prints
max|c| and the growth rate it implies, ln(max|c|/δ)/t, next to σ. During the
linear phase it adds max|error|/max|exact| of c and μ, and afterwards the
distance of max|c| from A.

Output of the shipped run, which takes about 10 minutes on one core
(800,000 steps of Δt = 2.5e-6):

| t | max\|c\| | ln(max\|c\|/δ)/t | error against the exact solution |
|---|---|---|---|
| 0.25 | 7.377e-3 | 7.9934 | c: 1.6e-3, μ: 2.5e-3 (linear solution) |
| 0.5  | 5.429e-2 | 7.9888 | — (cubic term starts to act) |
| 1.0  | 0.7707   | 6.647  | max\|c\|/A − 1 = −3.8e-2 |
| 1.5  | 0.80077  | 4.457  | max\|c\|/A − 1 = −1.6e-5 |
| 2.0  | 0.80078  | 3.343  | max\|c\|/A − 1 = −6.7e-9 |

During the linear phase the growth rate is 7.9934 against σ = 8. That gap,
8e-4 relative, is the τσ the relaxation predicts. At t = 2 the solution has
settled on the steady amplitude to 7e-9, with |μ| < 6e-9.
