# Magnetized astrophysical jet (Wu & Shu 2018, Example 5.6)

A Mach 800 beam, ten times denser than its surroundings, is injected into a static and strongly magnetized medium. This is the closing test of most positivity-preserving MHD papers. In the beam the gas pressure is 5.6 parts per million of the total energy, and in the ambient it is 2.4%, so a scheme without a positivity fix produces negative pressure within a few steps. Wu & Shu say exactly this happens without their limiter.

## Setup

Wu & Shu (SISC 40 (2018) B1302, Example 5.6), after Balsara's (2012) gas-dynamical jet:

| | |
|---|---|
| domain | [-0.5, 0.5] × [0, 1.5], full width (the paper computes the right half with a symmetry axis) |
| ambient | ρ = 0.1γ = 0.14, **v** = 0, p = 1, **B** = (0, B_a, 0) |
| nozzle {y = 0, \|x\| ≤ 0.05} | ρ = γ = 1.4, **v** = (0, 800, 0), p = 1, same **B**, ψ = 0 |
| other boundaries | outflow |
| γ, t_end | 1.4, 0.002 |
| B_a | √200 (β_a = 10⁻²), the paper's case (i). `aj_Ba` in `user_flux.jl`: √2000 and √20000 are cases (ii) and (iii) |

The beam sound speed is exactly 1, so u_jet = 800 is exactly Mach 800. The GLM cleaning speed is c_h = 812, the largest of the ambient and beam wave speeds.

**Two departures at the inlet datum**, both in `user_flux.jl` and `user_bc.jl`:

| | default | the paper | restore it with |
|---|---|---|---|
| nozzle lip, in x | C² smootherstep over 2s = 2 elements, centred on \|x\| = 0.05 | top hat | `aj_smooth = 0` |
| beam turn-on, in t | smootherstep over τ = 2h/u_jet (125 steps) | impulsive | `aj_tramp = 0` |

The blend acts on the primitives, so the injected mass and momentum flux are unchanged to five digits (∫φ dx = 0.05), and the beam is at full strength from t = τ = 3% of t_end. Its purpose is that the clamped nodes agree with their free neighbours: a top hat puts a 4400× jump in ρE between a clamped and a free node inside one element.

**Open boundaries.** Nothing is imposed except ψ = 0 (`user_bc.jl`). With every variable left free, the GLM pair (Bₙ, ψ) has no incoming condition. In the first run with the conservative limiter, ψ grew about 10× every 2·10⁻⁵, starting at the first free bottom-boundary node next to the nozzle patch. It reached 458 (B_a = 14) along y = 0 and at the corners, ate the ambient pressure through ½ψ², and the run blew up at the corner (−0.5, 0) at t = 5.0·10⁻⁴.

## Discretization and stabilization

CG spectral elements, LGL nodes, `:nop => 4`, CarpenterKennedy2N54:

| mesh | h | nodes | Δt | steps to 0.002 |
|---|---|---|---|---|
| `AJ_40x60.msh` (default) | 0.025 | 160 × 240 | 5·10⁻⁷ | 4,000 |
| `AJ_100x150.msh` (the paper's spacing) | 0.01 | 400 × 600 | 2·10⁻⁷ | 10,000 |

Both meshes put the nozzle lip |x| = 0.05 on an element boundary. `tools/astro_jet_mesh.py` regenerates them (`AJ.geo` is the gmsh definition; its surface must be named `"domain"`).

- **DynSGS**, conserved form (∇·(ν∇q) on all nine variables), with the residual sensor and domain norms, `C_R = 1`, `C_max = 0.5`, no background viscosity and no startup hold. The residual sensor is the one that excludes the constrained nozzle nodes.
- **Fast-speed floors** (`:dsgs_fast_floors => true`). Each residual R_i is divided by D_i = max(spread_i, S_i), where the floors S_i come from the mean state. With this option the floors use the fast-speed bound c̄_f = √((γp̄ + |**B̄**|²)/ρ̄) instead of the sound speed c̄ = √(γp̄/ρ̄):

  S_ρ = ρ̄,  S_ρv = ρ̄ c̄_f,  S_E = γp̄ + |**B̄**|²,  S_B = √(γp̄ + |**B̄**|²).

  At β_a = 10⁻², c_f/c = 12, so with sound-speed floors a resolved fast wave reads 12× (ρv, **B**) to 144× (E) too strong. With domain norms the floors bind only while the beam is ramped on (t ≲ 10⁻⁴). After that, the domain-mean state, which mixes the beam's kinetic energy into Ē, already gives c̄ ≈ 200.
- **Conservative positivity limiter** (`:positivity_method => "conservative"`, as in `problems/MHD/rotorDaoNazarov2022`), with ε_ρ = 1.4·10⁻⁷ and ε_p = 10⁻⁶ (10⁻⁶ of the ambient values). An element with a node below ε is scaled toward its lumped-mass mean (θ_ρ, then θ_p from the exact root of p = ε), and the increments are assembled like the RHS. Mass, momentum and energy are conserved, and **B** and ψ are untouched. It runs as the stage limiter of the RK, 5 times per step. If an element mean is itself inadmissible, a non-conservative node floor is the fallback, and it is counted in the report.

## Run

```bash
julia --project=. src/Jexpresso.jl MHD astroJetWuShu2018
mpiexec -n 4 julia --project=. src/Jexpresso.jl MHD astroJetWuShu2018
```

Outputs every 2·10⁻⁵ to t = 10⁻⁴, then every 10⁻⁴ (VTK). The fields are ρ, u, v, w, p, **B**, ψ, T, log10rho, log10p, β, Mach, `mu_dsgs`, and the numerical schlieren. For figures:

```bash
D=output/MHD/astroJetWuShu2018/output-<date>
python3 tools/plot_mhd_matrix.py $D --times 0.001 0.0015 0.002 --vars rho p Mach mu_dsgs --log rho p --width 9
```

## Results (`AJ_40x60`, B_a = √200, 4 MPI ranks, 4,000 steps in 12 min)

The run reaches t = 0.002. Earlier setups aborted much sooner: at t = 3.7·10⁻⁴ with the node-wise realizability repair, and at t = 5.0·10⁻⁴ with the conservative limiter but ψ left free at the open boundaries.

| t | ρ | max p | jet head (ρ > 0.3 on \|x\| < 0.05) | bow shock top, half-width | max \|ψ\| | max ν |
|---|---|---|---|---|---|---|
| 0.0011 | 0.062 – 7.71 | 4.70·10⁴ | y = 0.68 | 0.73, 0.30 | 3.5 | 0.60 |
| 0.0016 | 0.032 – 8.83 | 4.79·10⁴ | y = 0.99 | 1.03, 0.40 | 3.7 | 0.70 |
| 0.0020 | 0.020 – 8.79 | 4.83·10⁴ | y = 1.24 | 1.28, 0.47 | 4.0 | 0.74 |

- The head speed is ≈ 620, against the ram-pressure estimate v_h = u_jet/(1 + √(ρ_a/ρ_j)) = 608.
- ρ is symmetric about x = 0 to 3.5·10⁻⁵ at t = 0.002 (1.5·10⁻⁶ at t = 0.0016). With sound-speed floors it was 2.5·10⁻⁸.
- ν is largest at the jet head, moderate along the bow shock, and near zero in the beam and in the ambient.

**Positivity.** The limiter runs in 99% of the stage calls, and each call limits about 57 of the 2,400 elements (min θ_p = 0.047). No element ever has an inadmissible mean, so the non-conservative node floor never fires, and the limiter conserves mass, momentum and energy throughout. With sound-speed floors, 52 elements had an inadmissible mean while the beam was ramped on, giving 222 node floors and min θ_p = 0.0014; the solution at t = 0.002 differs by about 1% in L2. At t = 0.002, 76 nodes (0.2%) sit at p = ε. All of them are at ambient density (0.12 – 0.17) on the foot of the bow shock: the ambient thermal energy is only 2.4% of E, so the pre-shock undershoot of E takes p to ε, and the limiter holds it there. The Mach number at those nodes is meaningless.

## References

- K. Wu, C.-W. Shu, SIAM J. Sci. Comput. 40 (2018) B1302–B1329, Example 5.6.
- D. S. Balsara, J. Comput. Phys. 231 (2012) 7504–7517.
- X. Zhang, C.-W. Shu, J. Comput. Phys. 229 (2010) 8918–8934.
