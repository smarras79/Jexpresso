# MHD rotor (Balsara & Spicer 1999; Dao & Nazarov 2022, §5.5)

A dense disc spins in a light, magnetized gas at rest. The rotation twists the field lines and launches torsional Alfvén waves; near the disc the gas pressure drops close to zero, which is what makes the test hard (Tóth 2000 reports that many schemes fail on negative pressure).

## Setup

The "first rotor problem" of Tóth (2000, JCP 161:605, §6.6), as used by Dao & Nazarov (2022, J. Sci. Comput. 92:77, §5.5):

| | |
|---|---|
| domain | (0, 1)², doubly periodic (the waves do not reach the boundary by t = 0.15) |
| ambient | ρ = 1, **u** = 0, p = 1, **B** = (5/√(4π), 0) |
| disc, r < r₀ = 0.1 | ρ = 10, **u** = (u₀/r₀)(0.5 − y, x − 0.5) |
| taper, r₀ ≤ r < r₁ = 0.115 | ρ = 1 + 9f, **u** = (f u₀/r)(0.5 − y, x − 0.5), f = (r₁ − r)/(r₁ − r₀) |
| u₀, γ, t_end | 2, 1.4, 0.15 |

Magnetic pressure is ½|**B**|² (Tóth's units), so p_mag = 25/(8π) ≈ 0.995 initially. Dao & Nazarov do not print u₀; u₀ = 2 is Tóth's first rotor.

The equations are the 9-field GLM-MHD system of `problems/MHD/orszagTangBormanis2024` (its `user_flux.jl`, `user_source.jl` and `user_bc.jl` are included, with γ = 1.4), on the same 32×32 periodic mesh refined twice by p4est: 128×128 elements at `:nop => 3`, 385×385 nodes (Dao & Nazarov: P3 on 300×300 nodes).

## Stabilization

DynSGS alone (no background floor, no positivity repair), with the default residual coefficient `:dsgs_CR => 1`:

- **conserved form** (`:dsgs_conserved => true`): one residual-based ν and a Laplacian on (ρ, ρ**v**, E, **B**, ψ). In the spinning ring ½ρ|**v**|² ≈ 20 against p/(γ−1) = 2.5; the physical form (diffusion of u, v and T with mass diffusion) drove p to −10 by t = 0.05.
- **element norms** (`:dsgs_norms => "element"`): the residual is normalized by each element's own spread and scales (floored at ρc, ρc², … with `:dsgs_local_rel = 1`). Inside the disc p ≪ ½ρ|**v**|², so against the domain-wide scale of E the residual of its grid-scale pressure noise is invisible; with domain norms that noise grew to 4% of p (median node-to-node second difference) and Mach spiked to 7 (15 on a 2× finer grid) where p → 0.

## Run

```bash
julia --project=. src/Jexpresso.jl MHD rotorDaoNazarov2022
mpiexec -n 8 julia --project=. src/Jexpresso.jl MHD rotorDaoNazarov2022
```

Outputs at t = 0, 0.05, 0.10, 0.15 (VTK): ρ, u, v, p, **B**, ψ, T = p/ρ, `pmag` = ½|**B**|², `Mach` = |**v**|/√(γp/ρ). For a figure:

```bash
python3 tools/plot_fluxemergence_vtu.py output/MHD/rotorDaoNazarov2022/output-<date> --steps 4 --var rho \
        --xlabel '$x$' --ylabel '$y$' --time-unit ''      # also --var p, pmag, Mach; black lines = field lines
```

## Results at t = 0.15

128×128 elements, `:nop => 3`, 385×385 nodes, Δt = 2e-4, 4 MPI ranks:

| | ρ | p | ½\|**B**\|² | max Mach |
|---|---|---|---|---|
| this case | 0.567 – 10.38 | 0.0347 – 1.971 | 0.074 – 2.541 | 3.52 |
| Tóth (2000) Fig. 18, flux-CT, 400² | 0.483 – 12.95 | 0.0202 – 2.008 | 0.0177 – 2.642 | 8.18 (spurious, at p undershoots) |
| Dao & Nazarov (2022) Fig. 7, P3, 300² nodes | 0.727 – 8.42 | 0.0386 – 1.93 | 0.0551 – 2.30 | 4.82 |

min p = 0.059 at t = 0.10 and 0.035 at t = 0.15; the grid-scale pressure noise inside the ring is 0.03% (median) / 0.2% (95th percentile) of p;
the solution is invariant under the 180° rotation about (0.5, 0.5) to 2e-10 (ρ, p, **B** even, **v** odd); mass and energy change by less than 1e-4.

What each DSGS option did at this resolution (t = 0.15 unless noted):

| DynSGS setting | min p | max Mach | p noise (median) |
|---|---|---|---|
| physical form (u, v, T), domain norms | −9.9 at t = 0.05, blow-up at t = 0.084 | | |
| conserved, domain norms, C_R = 1 | −0.105 | ∞ (p < 0) | |
| conserved, domain norms, C_R = 1, residual sensor | −0.036 | ∞ | |
| conserved, domain norms, C_R = 1, nodal ν | blow-up at t = 0.124 | | |
| conserved, domain norms, C_R = 2 | 0.0058 | 7.3 | 3.8% |
| conserved, domain norms, C_R = 4 | 0.029 | 3.8 | 0.46% |
| conserved, element norms, C_R = 2 | 0.045 | 3.0 | 0.04% (symmetry error 2e-3) |
| **conserved, element norms, C_R = 1** (this deck) | **0.035** | **3.5** | **0.03%** |

## Finer grids

Scale Δt with the smallest LGL spacing, Δx_min ≈ (1 − ξ₁)/2 · 1/(32·2^lvl) with ξ₁ the first interior LGL node, to keep the acoustic CFL ≈ 0.25 (max wave speed ≈ 2.9): `:nop => 3`, level 3: Δt = 1e-4; `:nop => 4`, level 3: Δt = 6e-5.

## References

- D. S. Balsara, D. S. Spicer, J. Comput. Phys. 149 (1999) 270–292.
- G. Tóth, J. Comput. Phys. 161 (2000) 605–652, §6.6 and Fig. 18.
- T. A. Dao, M. Nazarov, J. Sci. Comput. 92 (2022) 77, §5.5 and Fig. 7.
