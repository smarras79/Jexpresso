# PODbenchmark — the POD reference benchmark

A problem whose Proper Orthogonal Decomposition can be **written down in closed
form**, so that what Jexpresso computes is checked against arithmetic and not
against another run. It is the standard test problem of the transport-dominated
model-reduction literature, and it is here for two reasons: to verify the
implementation, and because what it says about POD is worth knowing before
building a reduced-order model on one.

```bash
julia --project=.
julia> using Jexpresso
julia> Jexpresso.run_case("AdvDiff", "PODbenchmark")
```

The same closed form is asserted to `1e-10` by `test/test_pod_benchmark.jl`,
which needs no package instantiation:

```bash
julia test/test_pod_benchmark.jl
```

---

## The problem

Linear advection of a multi-harmonic wave on a periodic line, run for exactly
one revolution:

```
∂u/∂t + c ∂u/∂x = 0 ,   x ∈ [0,L) periodic ,   c = 1 , L = 2

u(x,0) = Σ_{j=1}^{3} A_j cos(2πj x/L + ϕ_j) ,   A = (1, ½, ¼) ,  ϕ = (0, 0.7, −1.3)

⟹  u(x,t) = u(x − ct, 0)
```

## Its POD, exactly

Averaging over one period of the translation gives a convolution kernel

```
R(x,x') = ⟨u'(x,t) u'(x',t)⟩ = Σ_j (A_j²/2) cos(2πj(x−x')/L)                  (B1)
```

whose eigenfunctions are therefore the Fourier modes. Because
`cos(2πj(x−x')/L) = cos·cos + sin·sin`, each wavenumber contributes a
**two-dimensional** eigenspace:

```
λ_{2j−1} = λ_{2j} = A_j² L / 4 ,     span{ cos(2πjx/L), sin(2πjx/L) }         (B2)
E_j      = λ_j / Σλ = A_j² / (2 Σ_i A_i²)                                     (B3)
Σλ       = ⟨‖u'‖²⟩ = (L/2) Σ_j A_j²                                           (B4)
ε(r = 2m)² = Σ_{j>m} A_j² / Σ_j A_j²          (the truncation error)          (B5)
```

With `A = (1, ½, ¼)` and `L = 2` that is:

| pair | λ (exact) | E [%] (exact) | cumulative [%] |
|---:|---:|---:|---:|
| 1, 2 | 0.5    | 38.0952 each | 76.19 |
| 3, 4 | 0.125  |  9.5238 each | 95.24 |
| 5, 6 | 0.03125|  2.3810 each | 100.00 |

`Σλ = 1.3125`, and the temporal mean is zero — a travelling wave averages to
nothing over a period, and the decomposition must find that rather than spend a
mode on it.

**Three things to check in the output**, in `output/`:

* `pod_u_spectrum.csv` — the six eigenvalues against the table above, and the
  truncation-error column against (B5): `1, 0.4880, 0.2182, 0` at
  `r = 0, 2, 4, 6`.
* `pod_u_coefficients.png` — the `(a₁,a₂)` phase portrait is a **circle**.
  That is what a travelling structure looks like in the plane of its own pair;
  a standing one traces a line. The circle closing exactly is the statement that
  the run returned to its initial condition.
* `pod_u_modes.csv`, `pod_u_modes.png` — the modes are `cos(2πjx/L)` and
  `sin(2πjx/L)`. Note that **which** rotation of the pair comes out is not
  determined by the problem: λ₁ = λ₂, so any orthonormal basis of the plane is a
  valid answer, and the phases ϕ_j decide it. What is determined — and what the
  test checks — is the plane.

## Two things the deck does on purpose

**`:pod_tend` is one sampling interval short of `:tend`.** The degeneracy
λ₁ = λ₂ is a statement about averaging over a *whole* period. Sampling `[0,T]`
at both ends repeats the zero phase, and that single repeat splits every pair by
`(K/2+1)/(K/2)` — 5 % at 41 snapshots — for a reason that has nothing to do with
the decomposition. Ending the window one interval early makes the 41 snapshots
tile exactly one period with no repeat. The run still goes to `T`; only the
sampling window stops short. (The sampling times are added to the integrator's
`tstops`, so the snapshots are taken *at* them rather than at the first step
after.)

**The phases ϕ_j are not zero.** They do not appear in (B2)–(B5) — shifting a
harmonic rotates the two modes of its pair within their own plane and changes
nothing observable. They are non-zero so that a run cannot agree with the
reference for the wrong reason: a code that quietly assumed a pure cosine would
still get the spectrum right, and the modes and phase portraits wrong.

## What this benchmark says about POD

Set `A = (1, 1, 1, …)` in `initialize.jl` and the spectrum goes **flat**: by
(B3) every one of the `2J` modes carries `1/(2J)` of the energy, and truncation
buys nothing — keeping half the modes leaves exactly half the energy behind.

That is not a defect of the implementation; it is the Kolmogorov *n*-width of a
transported profile, and it is the known limitation of every *linear* reduced
basis for transport-dominated problems. It is the most consequential thing a POD
implementation can get wrong: reporting a comfortable decaying spectrum here
would promise a reduced-order model that cannot exist. The decaying spectrum in
the table above comes from the amplitudes of the harmonics, not from the
transport.

(Sirovich 1987; Holmes, Lumley, Berkooz & Rowley 2012, §3.3; Ohlberger & Rave
2016; Greif & Urban 2019.)

## The rest of the case

Six files, like every Jexpresso case, and nothing in them knows about POD:

| file | what it carries |
|---|---|
| `initialize.jl` | the multi-harmonic initial condition — **the amplitudes are the answer** |
| `user_flux.jl` | `F = c u`, `c = 1` (must match the deck's `c`) |
| `user_source.jl`, `user_bc.jl`, `user_primitives.jl` | nothing case-specific; periodic, no source |
| `user_inputs.jl` | the grid, the time step, and the three-line POD block |

The POD block reduces, in a case that is not a benchmark, to the single line
`:lpod => true`. See [`docs/POD.md`](../../../docs/POD.md).
