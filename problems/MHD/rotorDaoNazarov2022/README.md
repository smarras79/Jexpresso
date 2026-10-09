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

DynSGS in its **conserved form** (`:dsgs_conserved => true`): one residual-based ν and a Laplacian on (ρ, ρ**v**, E, **B**, ψ), with the residual coefficient `:dsgs_CR => 2`. In the spinning ring ½ρ|**v**|² ≈ 20 against p/(γ−1) = 2.5, and the physical form (diffusion of u, v and T with mass diffusion) drove p to −10 by t = 0.05; the conserved form keeps p positive. With `:dsgs_nazarov_energy => true` (κ = ρν/Pr on the thermal energy) p dipped to −0.02 on 64² elements.

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

Extrema at t = 0.15 (128×128 elements, :nop => 3, 385×385 nodes, Δt = 2e-4, 4 MPI ranks):

| | ρ | p | ½\|**B**\|² | Mach |
|---|---|---|---|---|
| this case | 0.539 – 10.4 | 0.0058 – 1.965 | 0.066 – 2.561 | 7.31 |
| Tóth (2000) Fig. 18, flux-CT, 400² | 0.483 – 12.95 | 0.0202 – 2.008 | 0.0177 – 2.642 | 8.18 |
| Dao & Nazarov (2022) Fig. 7, P3, 300² nodes | 0.727 – 8.42 | 0.0386 – 1.93 | 0.0551 – 2.30 | 4.82 |

p stays positive at every output (min 0.338, 0.030, 0.0058 at t = 0.05, 0.10, 0.15); the solution is invariant under the
180° rotation about (0.5, 0.5) to 1e-10 (ρ, p, **B** even, **v** odd); total mass and energy change by less than 1e-4.
The Mach maximum sits where p is smallest, at the inner edge of the ring, like Tóth's.

What it took, measured at this resolution with DSGS alone (no floor, no positivity repair):

| DynSGS setting | min p at t = 0.1 / 0.15 |
|---|---|
| physical form (u, v, T), C_R = 1 | −9.9 at t = 0.05, blow-up at t = 0.084 |
| conserved form, C_R = 1 | −0.06 / −0.105 |
| conserved form, C_R = 1, residual sensor | −0.013 / −0.036 |
| conserved form, C_R = 1, nodal ν | −0.36 / blow-up at t = 0.124 |
| **conserved form, C_R = 2** (this deck) | **0.030 / 0.0058** |

## References

- D. S. Balsara, D. S. Spicer, J. Comput. Phys. 149 (1999) 270–292.
- G. Tóth, J. Comput. Phys. 161 (2000) 605–652, §6.6 and Fig. 18.
- T. A. Dao, M. Nazarov, J. Sci. Comput. 92 (2022) 77, §5.5 and Fig. 7.
