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

Same as `problems/MHD/rotorDaoNazarov2022` (see its [algorithm.pdf](../rotorDaoNazarov2022/algorithm.pdf)):

- conserved-form DynSGS: one ν per element, ∇·(ν∇q) on every conserved variable, the legacy sensor with element norms, `C_R = 1`, `C_max = 0.5`, and no background viscosity. One difference from the rotor: there is no startup hold (`:dsgs_hold_steps => 0`). The pressure jump at t = 0 is the most violent moment of the run, and outside it the initial data are uniform, so the first residual is nonzero only at the jump;
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

193² nodes, Δt = 10⁻⁵ (acoustic CFL 0.13 at the start, 0.11 at the end), 4 MPI ranks, 220 s:

| ρ | p | ½\|**B**\|² | \|**v**\| | Mach |
|---|---|---|---|---|
| 0.2253 – 3.422 | 0.09995 – 242.4 | 234.0 – 574.1 | ≤ 16.32 | ≤ 5.41 |

The dense shells sit along **B**, at x ≈ 0.2 and 0.8. The fast front in magnetic pressure spans y ≈ 0.12 – 0.88, and **B** is expelled from the hot interior (|**B**| = 21.6 – 33.9 against 28.2 outside). This is the structure of Balsara (2004) Fig. 6. ψ stays below 0.06 (|**B**| ≈ 28).

**Positivity.** The limiter acts only during the first 200 steps (t < 0.002). Over that window it makes 2,488 element limitings (min θ_p = 0.0016). In steps 1–4, 152 elements also have an inadmissible mean, i.e. p̄ ≤ ε even before limiting. These are elements next to the initial jump, where the thermal energy is 0.06% of the magnetic energy (β = 2.5·10⁻⁴). There the non-conservative node floor raises p at 554 nodes in total, which adds 1.5·10⁻⁵ of the total energy. After step 4 the energy holds to 6.7·10⁻¹⁶ and the mass to 3.3·10⁻¹⁵ for the whole run. With the rotor's two-step hold (ν = 0 on steps 1–2) there were 654 floors and +1.8·10⁻⁵ of energy, and the solution at t = 0.01 is the same to 3–4 digits.

**DynSGS: saturated at low β.** With the same settings as the rotor, ν sits at its first-order cap C_max Δ (|**v**| + c_f) over most of the region the fast front has crossed, not only at the shocks. The table counts the share of nodes where ν exceeds a given fraction of the local cap:

| | ν > 0.9 cap | ν > 0.5 cap | ν > 0.1 cap |
|---|---|---|---|
| blast, t = 0.003 | 12.6% | 14.6% | 22.9% |
| blast, t = 0.01 | 58.7% | 65.8% | 87.6% |
| rotor 128², t = 0.15 | 1.6% | 9.5% | 18.3% |

Ahead of the fast front ν is small but nonzero, ≲ 1.5·10⁻³ at r = 0.4 – 0.5 at t = 0.003. There the solution carries dispersive precursors of the central CG discretization, with |**v**| and |δ**B**| ≈ 3·10⁻⁵.

The likely cause is the floors of the element normalization: ρ_e, ρ_e c_e, ρ_e c_e² and √ρ_e c_e use the sound speed c_e = √(γp_e/ρ_e). That is the right rate at β ≈ 1, as in the rotor, but in this ambient c_e = 0.37 while the fast speed is 28.2. So a magnetosonic disturbance is measured against a scale 75× (momentum, **B**) to 5,700× (E) too small, and reads as unresolved.

`:dsgs_fast_floors => true` builds the floors on c̄_f = √((γp̄ + |**B̄**|²)/ρ̄) instead (DSGS.md §4.5). It is not on in this deck. In a test run to t = 0.01 it changed:
- **ν:** at t = 0.01, ν > 0.9 cap at 0.5% of the nodes instead of 59%, and the mean ν is 10× smaller;
- **ρ_max:** 3.67 instead of 3.42;
- **positivity:** p now reaches ε at about 50 nodes at every output time. They lie on the fast front, at ambient density, where the compressed field leaves 0.06% of E as thermal energy. The limiter therefore acts in every stage call.
- **conservation:** still exact after step 7, with ΔE/E = 1.8·10⁻⁵ from the startup node floors.

## References

- D. S. Balsara, D. S. Spicer, J. Comput. Phys. 148 (1999) 133–148 (1999b).
- D. S. Balsara, ApJS 151 (2004) 149–184, §7.3.
