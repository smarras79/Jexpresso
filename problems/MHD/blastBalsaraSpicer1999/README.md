# MHD blast wave (Balsara & Spicer 1999b; Balsara 2004, §7.3)

A circle of over-pressured gas expands into a strongly magnetized medium at rest. The ambient plasma β is 2.5·10⁻⁴, so the gas pressure outside the blast is a tiny remainder of the total energy. Schemes without a positivity fix produce negative pressure here, which is why the problem is a standard positivity test (Balsara & Spicer 1999b; Balsara 2004).

## Setup

The non-relativistic problem of Balsara & Spicer (1999b), used in Balsara (2004) §7.3 ("200 × 200 zones … t = 0.01 … ambient β = 0.000251"):

| | |
|---|---|
| domain | (0, 1)², doubly periodic (the fast front stays inside by t = 0.01) |
| everywhere | ρ = 1, **v** = 0, **B** = (100/√(4π), 0, 0), ψ = 0 |
| r < 0.1 about (0.5, 0.5) | p = 1000 |
| r ≥ 0.1 | p = 0.1 |
| γ, t_end | 1.4, 0.01 |

**B** is in Heaviside–Lorentz units (magnetic pressure ½|**B**|²), so |**B**| = 100/√(4π) is Balsara & Spicer's B = 100 in Gaussian units. The ambient β = 2p/|**B**|² = 2.51·10⁻⁴, as stated by Balsara (2004). Other papers use [-0.5, 0.5]² with the same data, which is the same problem shifted.

The equations are the 9-field GLM-MHD system of `problems/MHD/orszagTangBormanis2024`. This case includes its `user_flux.jl`, `user_source.jl` and `user_bc.jl`, with γ = 1.4. The mesh is the same 32×32 periodic mesh, refined once by p4est: 64×64 elements at `:nop => 3`, i.e. 193×193 nodes, close to the paper's 200² zones. The pressure jump is sampled at the nodes as given, with no taper.

## Stabilization

As `problems/MHD/rotorDaoNazarov2022` (see its [algorithm.pdf](../rotorDaoNazarov2022/algorithm.pdf)), with two differences:

- conserved-form DynSGS: one ν per element, ∇·(ν∇q) on every conserved variable, the legacy sensor with element norms, `C_R = 1`, `C_max = 0.5`, and no background viscosity.
  - **Fast-speed floors** (`:dsgs_fast_floors => true`, DSGS.md §4.5). Each residual R_i is divided by D_i = max(spread_i, S_i). The floors S_i use the fast-speed bound c̄_f = √((γp̄ + |**B̄**|²)/ρ̄) of the element-mean state instead of the sound speed c̄ = √(γp̄/ρ̄): S_ρ = ρ̄, S_ρv = ρ̄ c̄_f, S_E = γp̄ + |**B̄**|², S_B = √(γp̄ + |**B̄**|²). In this ambient c_f/c = 75. With sound-speed floors, a resolved fast wave read 75× (ρv, **B**) to 5,700× (E) too strong, and ν sat at its cap over most of the disturbed region (see the comparison below).
  - **No startup hold** (`:dsgs_hold_steps => 0`). The pressure jump at t = 0 is the most violent moment of the run, and outside it the initial data are uniform, so the first residual is nonzero only at the jump.
- the conservative positivity limiter (`:positivity_method => "conservative"`), with ε_ρ = 10⁻⁶ (10⁻⁶ of the ambient ρ) and ε_p = 10⁻⁷ (10⁻⁶ of the ambient p). Elements with a node below ε are scaled toward their mean. Mass, momentum and energy are conserved, and **B** and ψ are untouched.

## Run

```bash
julia --project=. src/Jexpresso.jl MHD blastBalsaraSpicer1999
mpiexec -n 4 julia --project=. src/Jexpresso.jl MHD blastBalsaraSpicer1999
```

Outputs every 0.001 (VTK): ρ, u, v, p, **B**, ψ, T = p/ρ, `pmag` = ½|**B**|², `Mach`, `β` = 2p/|**B**|², and the DynSGS coefficient `mu_dsgs`. Mass and energy are tracked every step in `conservation.dat`. Plot commands:

```bash
D=output/MHD/blastBalsaraSpicer1999/output-<date>
python3 tools/plot_mhd_matrix.py $D --times 0.003 0.007 0.01 --vars rho p pmag Mach mu_dsgs \
        --log p mu_dsgs                                   # rows = times, columns = variables; --fieldlines adds B lines
python3 tools/plot_mhd_matrix.py $D --times 0.01 --vars rho p speed Bmag --log rho p   # Balsara (2004) Fig. 6
python3 tools/plot_fluxemergence_vtu.py $D --steps 11 --var rho --no-fieldlines \
        --xlabel '$x$' --ylabel '$y$' --time-unit ''      # one field per figure; drop --no-fieldlines for B lines
python3 tools/plot_conservation.py $D --out blast_conservation.pdf
```

## Results at t = 0.01

193² nodes, Δt = 10⁻⁵ (acoustic CFL 0.13 at the start, 0.11 at the end), 4 MPI ranks. The second row is the same deck with sound-speed floors:

| floors | ρ | p | ½\|**B**\|² | \|**v**\| | max Mach (p > 10⁻³) | ν > 0.9 cap | mean ν |
|---|---|---|---|---|---|---|---|
| fast (this deck) | 0.2159 – 3.669 | 10⁻⁷ – 247.7 | 223.4 – 604.7 | ≤ 16.6 | 52.8 | 0.5% | 0.0044 |
| sound | 0.2253 – 3.422 | 0.09995 – 242.4 | 234.0 – 574.1 | ≤ 16.3 | 5.41 | 58.7% | 0.0438 |

The dense shells sit along **B**, at x ≈ 0.2 and 0.8. The fast front in magnetic pressure spans y ≈ 0.12 – 0.88, and **B** is expelled from the hot interior. This is the structure of Balsara (2004) Fig. 6. With fast floors, ν is concentrated on the shells and on the thin fast-front ring, and is near zero elsewhere, so the shells are sharper (ρ_max 3.67 against 3.42). ψ stays below 0.11 (|**B**| ≈ 28).

**The cost: a pressure undershoot at the fast front.** In this ambient the thermal energy is 0.06% of E, so p is a small remainder of E and the magnetic energy. With sound floors, the saturated ν damped the dispersive ripples of the fast front, and no node at t = 0.01 was below the ambient p = 0.1. With fast floors, a thin band at the front, deepest where the shells meet it at x ≈ 0.1 and 0.9, falls below ambient:

| p below | 0.0999 | 0.09 | 0.05 | 10⁻³ | 10⁻⁶ |
|---|---|---|---|---|---|
| nodes (of 37,249) | 5,852 | 786 | 220 | 50 | 30 |

There, Mach exceeds 6 at 2,136 nodes. The ~50 nodes at ε are where the limiter holds p.

**Positivity and conservation.**
- **Limiter activity:** the limiter acts in every stage call, with 303,688 element limitings in total (min θ_p = 0.00088).
- **Startup floors:** in steps 1–7, 196 elements next to the initial jump have an inadmissible mean, i.e. p̄ ≤ ε even before limiting. There the non-conservative node floor raises p at 666 nodes, which adds 1.8·10⁻⁵ of the total energy.
- **After step 7:** energy holds to 6.7·10⁻¹⁶ and mass to 2.9·10⁻¹⁵ for the rest of the run.
- **With sound floors:** the limiter stopped by step 200, and there were 554 startup floors (+1.5·10⁻⁵ of energy).

**Why the sound floors saturate here.** The table counts the share of nodes where ν exceeds a given fraction of its local cap C_max Δ (|**v**| + c_f):

| | ν > 0.9 cap | ν > 0.5 cap | ν > 0.1 cap |
|---|---|---|---|
| blast, sound floors, t = 0.01 | 58.7% | 65.8% | 87.6% |
| blast, fast floors, t = 0.01 | 0.5% | 1.8% | 17.8% |
| rotor 128² (sound floors, β ≈ 1), t = 0.15 | 1.6% | 9.5% | 18.3% |

## References

- D. S. Balsara, D. S. Spicer, J. Comput. Phys. 148 (1999) 133–148 (1999b).
- D. S. Balsara, ApJS 151 (2004) 149–184, §7.3.
