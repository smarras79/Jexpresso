# Elasticity — Timoshenko beams as a first-order hyperbolic system

Two 1D cases solving the **Timoshenko (shear-deformable) beam** equations in
conservation form, which is what makes them a fit for a SEM/DG solver at all.

```bash
julia --project=. src/Jexpresso.jl Elasticity simply_supported
julia --project=. src/Jexpresso.jl Elasticity cantilever
```

or, interactively,

```julia
julia> using Jexpresso
julia> Jexpresso.run_case("Elasticity", "simply_supported")
julia> Jexpresso.run_case("Elasticity", "cantilever")
```

---

## Why Timoshenko and not Euler–Bernoulli

Euler–Bernoulli has no genuine conservation form. Its fourth-order operator
makes it *dispersive*, not hyperbolic: the best available two-field form is

```
∂ₜ(ρA v) + ∂ₓ(∂ₓ(EI κ)) = q,        ∂ₜκ - ∂ₓ(∂ₓ v) = 0
```

— conservative in structure, but with second-order fluxes, so it wants an
LDG-type or mixed treatment rather than a Riemann solver. It is recovered from
the system below as the stiff relaxation limit κₛGA → ∞ (γ → 0).

Timoshenko, by contrast, is a clean **first-order symmetric hyperbolic system**.

## The system

State: the two momentum densities and the two strain measures,

```
U = (ρA v, ρI ω, γ, χ)ᵀ

v = ẇ        transverse velocity            w = transverse displacement
ω = φ̇        cross-section angular rate      φ = cross-section rotation
γ = w_x - φ  shear strain
χ = φ_x      curvature
```

```
∂ₜU + ∂ₓF(U) = S

F(U) = -(κₛGA γ, EI χ, v, ω)ᵀ          S = (q, κₛGA γ, -ω, 0)ᵀ
```

The first two rows are balance of linear and angular momentum, with fluxes
minus the shear force `Q = κₛGA γ` and the bending moment `M = EI χ`. The last
two are the compatibility conditions `γ̇ = v_x - ω` and `χ̇ = ω_x`.

`κₛGA γ` and `-ω` sit in the source because they are **undifferentiated**, not
because they are nonlinear: the whole system is linear with constant Jacobian
`A = ∂F/∂U`.

### Characteristics

The spectrum splits into two decoupled pairs,

```
λ = ±√(κₛG/ρ)     shear wave
λ = ±√(E/ρ)       extensional / rotary wave
```

so the system is strictly hyperbolic away from `κₛG = E`, and symmetrizable by
the energy norm

```
½(ρA v² + ρI ω² + κₛGA γ² + EI χ²)
```

which is the convex entropy; its flux is `-(Qv + Mω)`, the mechanical power.
Standard upwind or Rusanov fluxes apply directly — the cases here run on
Jexpresso's continuous SEM, so no Riemann solver is used, but the structure is
there if a DG variant is wanted.

The two-pair split is also what makes the boundary conditions easy to count:
each pair contributes exactly one incoming characteristic at each end, so each
end takes exactly two conditions.

## Beam and non-dimensionalisation

Both cases use the same beam, defined in one place per case by
`timoshenko_properties()` in `user_flux.jl`:

| | |
|---|---|
| L = 1, ρ = 1, E = 1, ν = 0.3, κₛ = 5/6 | rectangular section |
| h = 0.1, width 1 | so A = 0.1, I = 8.333e-5, **L/h = 10** |
| √(κₛG/ρ) = 0.566 | shear wave |
| √(E/ρ) = 1 | extensional/rotary wave — this one sets the time step |

At L/h = 10 shear deformation is a percent-level correction to
Euler–Bernoulli: large enough to see, small enough that the answer is still
recognisably a beam.

## The cases

### `simply_supported` — verification

Four equations, exactly the system above. The first flexural standing mode of a
pinned–pinned beam, released from rest at maximum deflection:

```
w = W sin(kx) cos(Ωt)      φ = Φ cos(kx) cos(Ωt)      k = nπ/L
```

with Ω the lower (flexural) root of the Timoshenko frequency equation

```
(ρIρA/κₛGA) Ω⁴ - [ρI k² + EI ρA k²/κₛGA + ρA] Ω² + EI k⁴ = 0
```

and `Φ = Wk(1 - ρAΩ²/(κₛGA k²))`. For n = 1 this gives Ω = 0.28023 against the
Euler–Bernoulli 0.28491 — the 1.6% the shear flexibility costs — and a period
T = 22.4215, which is exactly `:tend`.

`user_analytic.jl` supplies this solution, so the shared 1D plotter overlays it
on **every** output frame: the gap between the two curves is entirely the
discretisation's. After one full period every field must land back on its
initial profile.

Boundary conditions (pin: holds the beam down, lets it rotate):

| end | imposed | free |
|---|---|---|
| both | `v = 0` (w = 0), `χ = 0` (M = 0) | ω, γ |

The initial data excites *only* the flexural branch. The shear branch of the
same k sits at Ω ≈ 19.94, seventy times faster — anything of it that shows up
in the solution is discretisation error, and it is fast enough to be obvious.

### `cantilever` — deformation

Six equations: the four above plus **displacement recovery**,

```
U₅ = w      ∂ₜw + ∂ₓ(0) = v
U₆ = φ      ∂ₜφ + ∂ₓ(0) = ω
```

Zero flux, so they add two zero eigenvalues and change neither the
characteristic structure nor the stable time step. They are reconstructions,
not conservation laws — but they are what turns the output into a picture of
the beam actually bending.

Clamped at x = 0, free at x = L, released from the static shape a uniform load
q₀ holds it in:

| end | imposed | free |
|---|---|---|
| x = 0, clamped | `v = 0`, `ω = 0` (and w, φ pinned) | γ, χ |
| x = L, free | `γ = 0` (Q = 0), `χ = 0` (M = 0) | v, ω, w, φ |

The static shape is smooth **and** satisfies every one of those conditions
exactly, so the initial data is compatible and the solution stays smooth.
Suddenly *applying* the load instead would be much worse behaved: v would grow
uniformly everywhere except at the wall, which is a step in x propagating in at
the shear speed, and a non-dissipative spectral element scheme rings on it.

Static tip deflection, for reference:

```
w(L) = q₀L⁴/(8EI)  +  q₀L²/(2κₛGA)  =  1.5000e-2 + 1.560e-4  =  1.5156e-2
        bending         shear
```

`:tend` is the Euler–Bernoulli estimate of the fundamental period (T ≈ 61.90).
There is no closed form for the Timoshenko cantilever, whose fundamental sits
≈0.3% lower here (a fine-grid reference integration puts the half period at
31.06 against the Euler–Bernoulli 30.95), so the beam is deliberately a little
short of closing its cycle at the end of the run — and that gap *is* the shear
correction. The tip also overshoots the static amplitude slightly, to 1.559e-2,
because the released static shape is not exactly the first mode.

## A note on the source term (a bug this equation set flushed out)

Building these cases turned up a real defect in the 1D SEM kernels, fixed in
the same commit. It is worth recording, because the reason nothing caught it
for so long is instructive.

The 1D volume term weighted the flux by `ω[i]` and the source by `ω[i]`, and
the assembled result is divided by the lumped mass matrix `M[i] = Je[i]·ω[i]`:

```
flux:    ω·dFdξ      / (Je·ω)  =  dFdξ/Je = dF/dx    correct
source:  ω·S         / (Je·ω)  =  S/Je               wrong by 2/Δx
```

The flux is right *by cancellation*: in 1D `dF/dx = dFdξ·dξdx` with
`dξdx = 1/Je`, so weighting by `ω·Je` would immediately undo itself and the
kernel simply never wrote `Je` at all. The source has no `dξdx` to cancel
against, so the same (absent) weighting left it scaled by `1/Je`. The 2D and
3D kernels form `dF/dx` explicitly and then weight by `ω·Jac`, which is why
only 1D was affected — in all four 1D paths at once: `_expansion_inviscid!`
and `_expansion_inviscid_laguerre!` on the CPU, and `_build_rhs_gpu_v0!` on
the GPU (used for both the CG and the Laguerre region).

**Why every 1D case still looked fine.** Of the 1D cases in the tree, only
four have a live source, and all four are *sponges* — absorbing-layer damping
whose coefficient is a free tuning knob, so a silent factor `1/Je` was
absorbed into the tuning rather than showing up as a wrong answer. The one 1D
case with a genuinely physical source, `CompEuler/nozzleanderson` (the
quasi-1D area term `p·dA/dx/γ`), is not in CI and has never been checked
against Anderson's solution — at `Je = 0.015` its source was running 67× too
large. No CI case exercised a 1D source that was not a sponge, so there was
nothing to fail.

`CompEuler/wave1d_lag` *is* in CI, and its sponge had been tuned against the
old scaling. There `Je = 0.05` in the CG region and `Je = :yfac_laguerre =
0.05` in the Laguerre region — one value everywhere the sponge is active — so
the compensating `1/0.05` is folded into that deck's coefficient and its
solution is unchanged to round-off. The three `AdvDiff` wave-train sponges
straddle two regions with different `Je`, so no single deck-level constant
restores them; they were left as written and now act at the strength their
coefficients actually state.

For the cases here the fix is not cosmetic: the Timoshenko coupling terms
`κₛGA γ` and `-ω`, and the whole of the cantilever's w/φ recovery, live in
`S`. With the old weighting the cantilever's recovered displacement would have
grown `2/Δx` too fast and the deformed shape would have been meaningless.

## Large deformation

The geometrically exact (Simo–Reissner) beam keeps exactly the same shape,

```
∂ₜ(ρA 𝐯) - ∂ₓ𝐧 = 𝐧̄
∂ₜ(𝐈ω) - ∂ₓ𝐦 = 𝐱,ₓ × 𝐧 + 𝐦̄
```

but `𝐧, 𝐦` become nonlinear functions of the strain measures through the
rotation field, and the angular momentum source picks up the `𝐱,ₓ × 𝐧` term.
That is the version worth having if the beam is a structure to be coupled to a
flow solver.
