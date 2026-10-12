# Magnetized astrophysical jet (Wu & Shu 2018, Example 5.6)

A Mach 800 beam, ten times denser than its surroundings, is injected into a static and strongly magnetized medium. This is the closing test of most positivity-preserving MHD papers. In the beam the gas pressure is 5.6 parts per million of the total energy, and in the ambient it is 2.4%, so a scheme without a positivity fix produces negative pressure within a few steps. Wu & Shu say exactly this happens without their limiter.

## Setup

Wu & Shu (SISC 40 (2018) B1302, Example 5.6), after Balsara's (2012) gas-dynamical jet. As in the paper, the right half is computed, with a symmetry axis at x = 0:

| | |
|---|---|
| domain | [0, 0.5] × [0, 1.5], symmetry axis x = 0 (the full width [-0.5, 0.5] × [0, 1.5] runs on `AJ_40x60.msh`, `AJ_100x150.msh`) |
| ambient | ρ = 0.1γ = 0.14, **v** = 0, p = 1, **B** = (0, B_a, 0) |
| nozzle {y = 0, 0 ≤ x ≤ 0.05} | ρ = γ = 1.4, **v** = (0, 800, 0), p = 1, same **B**, ψ = 0 |
| right, top (and left on the full width) | ambient state (far field) |
| bottom outside the nozzle | outflow |
| γ, t_end | 1.4, 0.002 |
| B_a | √200 (β_a = 10⁻²), the paper's case (i). `aj_Ba` in `user_flux.jl`: √2000 and √20000 are cases (ii) and (iii) |

β_a is the ambient plasma beta, gas over magnetic pressure: β = 2p/|**B**|², since the magnetic pressure is ½|**B**|² in the Heaviside–Lorentz units of the code. With p = 1, β_a = 2/B_a². The beam sound speed is exactly 1, so u_jet = 800 is exactly Mach 800. The GLM cleaning speed is c_h = 812, the largest of the ambient and beam wave speeds.

**Two departures at the inlet datum**, both in `user_flux.jl` and `user_bc.jl`:

| | default | the paper | restore it with |
|---|---|---|---|
| nozzle lip, in x | C² smootherstep over 2s = 2 elements, centred on \|x\| = 0.05 | top hat | `aj_smooth = 0` |
| beam turn-on, in t | smootherstep over τ = 2h/u_jet (125 steps) | impulsive | `aj_tramp = 0` |

The blend acts on the primitives, so the injected mass and momentum flux are unchanged to five digits (∫φ dx = 0.05), and the beam is at full strength from t = τ = 3% of t_end. Its purpose is that the clamped nodes agree with their free neighbours: a top hat puts a 4400× jump in ρE between a clamped and a free node inside one element.

**Outflow** (`"bottom"` outside the nozzle). Nothing is imposed except ψ = 0 (`user_bc.jl`). With every variable left free, the GLM pair (Bₙ, ψ) has no incoming condition. In the first run with the conservative limiter, ψ grew about 10× every 2·10⁻⁵, starting at the first free bottom-boundary node next to the nozzle patch. It reached 458 (B_a = 14) along y = 0 and at the corners, ate the ambient pressure through ½ψ², and the run blew up at the corner (−0.5, 0) at t = 5.0·10⁻⁴.

**Far field** (`"right"`, `"top"`, and `"left"` on the full width). The ambient state is imposed there, where the paper has outflow.
- Do-nothing is unstable on these sides. A strong-form CG boundary node advances with the element's one-sided derivative, so the incoming characteristics are not controlled: the fast wave on x = 0.5, sound and Alfvén on y = 1.5. The linearized nop-4 operator then has a growing mode at the corners, Re λ ≈ 1.6 c_f/h, made of u, B_x and B_y. At 256² (h_x ≈ 2·10⁻³) it e-folds every 3·10⁻⁵ and shows as recirculation and field at the right corners.
- The ambient state is exact there until t_end. The bow shock reaches at most x = 0.45 and y = 1.30 at t = 0.002, and nothing runs ahead of it: its lateral speed, ≈ 200, exceeds the ambient fast speed, 38.
- Spurious max |**v**| at lvl 0, t = 1.9·10⁻³:

  | right and top | top-right corner | x = 0.5, y ≈ 1 | bottom-right corner |
  |---|---|---|---|
  | do-nothing + ψ = 0 | 6.0·10⁻² | 0.29 | 1.23 |
  | linearized characteristic condition | 1.1·10⁻² | 0.11 | 1.20 |
  | ambient state (default) | 6.6·10⁻⁴ | 9.9·10⁻³ | 4.5·10⁻² |

- For longer runs, enlarge the box so that the bow shock stays inside.

**Symmetry axis** (`"axis"`, `user_bc.jl`). Under x → −x the solution keeps ρ, ρv, ρw, E, B_y, B_z, ψ even and ρu, B_x odd. The equations preserve this parity, and **B** = (0, B_a) has it at t = 0.
- At x = 0 the odd fields are set to zero, ρu = 0 and B_x = 0; the even fields are left free.
- Their x-fluxes then vanish on the axis: ρu, ρuv − B_xB_y, u(E + p_tot) − B_x(**v**·**B**), uB_y − vB_x, and the GLM flux of ψ, which is ∝ B_x. The diffusive flux ν∂ₓq is zero by the natural boundary condition of the weak form.
- B_x = 0 is the one condition the GLM pair (B_n, ψ) needs.

At t = 2·10⁻⁵ the half-domain run equals the full-width run on x ≥ 0 to 2·10⁻¹¹. Afterwards the two differ by about 1%, because the residual sensor skips every boundary node, the axis column included.

**Why the half domain.** It halves the cost, as in the paper. It also enforces a symmetry that the full-width discretization does not keep with element norms. This is a weakness of the method, not a requirement for stability.
- In the full-width run, round-off differences between mirror elements at the nozzle lip grow about 10¹⁰-fold within 80 steps (t = 2–6·10⁻⁵), far faster than any physical instability.
- The amplifier is the DynSGS sensor. With ν forced mirror-symmetric the full width stays symmetric to 10⁻¹³ (t = 1.2·10⁻⁴), and halving Δt speeds the growth up (×2.8 → ×5.2 per step), which a flow instability would not do.
- By t = 0.002, ρ differs between left and right by 1.7·10⁻¹ at lvl 0 and 2.3·10⁻¹ at lvl 1.
- The full-width runs are nevertheless stable and positive to t = 0.002: no inadmissible element mean at lvl 0 or lvl 1.
- With domain norms the full width stays symmetric to 3.5·10⁻⁵, but fails from lvl 1 on (next section).

## Discretization and stabilization

CG spectral elements, LGL nodes, `:nop => 4`, CarpenterKennedy2N54:

| mesh | domain | h | nodes | Δt |
|---|---|---|---|---|
| `AJ_half_20x60.msh` (default) | right half | 0.025 | 80 × 240 | 5·10⁻⁷ |
| `AJ_half_50x150.msh` | right half | 0.01 (the paper's spacing) | 200 × 600 | 2·10⁻⁷ |
| `AJ_40x60.msh` | full width | 0.025 | 160 × 240 | 5·10⁻⁷ |
| `AJ_100x150.msh` | full width | 0.01 | 400 × 600 | 2·10⁻⁷ |

- **Refinement.** `lvl` at the top of `user_inputs.jl` refines uniformly with p4est and sets Δt = 5·10⁻⁷/2^lvl (CFL ≈ 0.1 against c_h). lvl = 1 and 2 give h = 0.0125 and 0.00625.
- **The h = 0.01 meshes.** On them Δt is 2·10⁻⁷/2^lvl; set it by hand.
- **Mesh files.** Every mesh puts the nozzle lip x = 0.05 on an element boundary. `tools/astro_jet_mesh.py` (with `--half` for the half meshes) regenerates them. `AJ.geo` is the gmsh definition (`half = 1` for the half domain), and its surface must be named `"domain"`.

- **DynSGS**, conserved form (∇·(ν∇q) on all nine variables), residual sensor, **element norms**, `C_R = 1`, `C_max = 0.5`, no background viscosity and no startup hold. The residual sensor is the one that excludes the constrained boundary nodes.
- **Why element norms.** The residual of each equation is divided by a scale D_i before ν_K = min(C_max Δ w_K, C_R Δ² max_i ‖R_i‖_{∞,K}/D_i).
  - With domain norms, D_E = ‖E − ⟨E⟩‖_{∞,Ω} ≈ 4.5·10⁵, the beam's kinetic energy. But in the beam p/(γ−1) is 5.6·10⁻⁶ of E.
  - A grid-scale noise δE carried at v gives ‖R_E‖ ≈ vδE/Δ, so ν ≈ C_R Δ v δE/D_E. A noise δE = 2,000 moves p by 800, yet reads ν ≈ 0.009, which is 1% of the cap.
  - With element norms, D_E = max(‖E − Ē‖_{∞,K}, γp̄ + |**B̄**|²) ≈ max(δE, 200) in the beam, and the same noise puts ν at the cap.

  Measured on the full width:

  | full width | norms | to t | elements with an inadmissible mean | node floors | min θ_p |
  |---|---|---|---|---|---|
  | lvl 0 | domain | 0.002 | 0 | 0 | 0.047 |
  | lvl 1 | domain | 1.5·10⁻⁴ | 635 | 1,901 | 2.5·10⁻⁵ |
  | lvl 0 | element | 0.002 | 0 | 0 | 0.24 |
  | lvl 1 | element | 0.002 | 0 | 0 | 0.24 |
  | lvl 2 | element | 1.5·10⁻⁴ | 0 | 0 | 0.077 |

  At lvl 2 with domain norms the run blew up. The cost of element norms is that ν sits near its cap in the beam: mean ν in the beam core is 0.73 against 0.13 at lvl 1, t = 10⁻⁴. The head is also smoother: at lvl 0, max ρ = 6.0 against 8.8.
- **Fast-speed floors** (`:dsgs_fast_floors => true`). The floors S_i come from the element-mean state, with the fast-speed bound c̄_f = √((γp̄ + |**B̄**|²)/ρ̄) instead of the sound speed c̄ = √(γp̄/ρ̄):

  S_ρ = ρ̄,  S_ρv = ρ̄ c̄_f,  S_E = γp̄ + |**B̄**|²,  S_B = √(γp̄ + |**B̄**|²).

  At β_a = 10⁻², c_f/c = 12, so with sound-speed floors a resolved fast wave reads 12× (ρv, **B**) to 144× (E) too strong.
- **Conservative positivity limiter** (`:positivity_method => "conservative"`, as in `problems/MHD/rotorDaoNazarov2022`), with ε_ρ = 1.4·10⁻⁷ and ε_p = 10⁻⁶ (10⁻⁶ of the ambient values).
  - An element with a node below ε is scaled toward its lumped-mass mean (θ_ρ, then θ_p from the exact root of p = ε), and the increments are assembled like the RHS. Mass, momentum and energy are conserved; **B** and ψ are untouched.
  - It runs as the stage limiter of the RK, 5 times per step.
  - If an element mean is itself inadmissible, a non-conservative node floor is the fallback, and it is counted in the report.

## Run

```bash
julia --project=. src/Jexpresso.jl MHD astroJetWuShu2018
mpiexec -n 4 julia --project=. src/Jexpresso.jl MHD astroJetWuShu2018
```

Outputs every 2·10⁻⁵ to t = 10⁻⁴, then every 10⁻⁴ (VTK). The fields are ρ, u, v, w, p, **B**, ψ, T, log10rho, log10p, β, Mach, Mach_f, `mu_dsgs`, and the numerical schlieren S = exp(−20|∇ρ|/max|∇ρ|). The figure has rows t = 0.0011, 0.0016, 0.002 and columns log₁₀ρ, log₁₀p, M_f, ν and S, with the half domain mirrored about x = 0 (u and B_x change sign):

```bash
D=output/MHD/astroJetWuShu2018/output-<date>
python3 tools/plot_astrojet.py $D                              # -> $D/matrix_rho_p_Mach_f_mu_dsgs_schlieren.png
python3 tools/plot_astrojet.py $D --format pdf --out jet.pdf   # any tools/plot_mhd_matrix.py option; --no-mirror for the full width
```

## Results (half domain, B_a = √200, 4 MPI ranks)

lvl 0 runs 4,000 steps in 3 min; lvl 1 runs 8,000 steps in 28 min.

| lvl | t | ρ | max p | jet head (ρ > 0.3 on x < 0.05) | bow shock top, half-width | max \|ψ\| | max ν (= cap) | nodes at p = ε |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.0011 | 0.080 – 5.47 | 4.29·10⁴ | y = 0.69 | 0.73, 0.30 | 0.34 | 2.03 | 14 |
| 0 | 0.0016 | 0.057 – 6.39 | 4.35·10⁴ | y = 0.99 | 1.03, 0.40 | 0.36 | 2.04 | 17 |
| 0 | 0.0020 | 0.042 – 6.89 | 4.48·10⁴ | y = 1.23 | 1.27, 0.45 | 0.37 | 2.04 | 25 |
| 1 | 0.0011 | 0.042 – 8.13 | 4.80·10⁴ | y = 0.71 | 0.73, 0.29 | 0.44 | 1.02 | 25 |
| 1 | 0.0016 | 0.054 – 11.4 | 4.92·10⁴ | y = 1.03 | 1.04, 0.39 | 0.44 | 1.02 | 29 |
| 1 | 0.0020 | 0.044 – 11.2 | 4.98·10⁴ | y = 1.28 | 1.30, 0.45 | 0.49 | 1.02 | 43 |

- **Head speed.** Between t = 0.0011 and 0.002 it is 600 (lvl 0) and 633 (lvl 1), against the ram-pressure estimate v_h = u_jet/(1 + √(ρ_a/ρ_j)) = 608.
- **Against the full width.** Head and bow shock match the full-width element-norm runs to 0.01.
- **Positivity.** The limiter runs in 99% of the stage calls: 362,853 element limitings at lvl 0 (min θ_p 0.23) and 1,470,098 at lvl 1 (min θ_p 0.12). No element ever has an inadmissible mean, so the non-conservative node floor never fires, and mass, momentum and energy are conserved throughout.
- **Nodes at the floor.** Of the nodes that sit at p = ε (25 and 43 at t = 0.002), all but one or two are at ambient density along the foot of the bow shock: the ambient thermal energy is 2.4% of E, so the pre-shock undershoot of E takes p to ε. The others are in the beam next to the nozzle (y ≤ 0.05).
- **The beam.** At t = 0.002 it keeps the injected state near the nozzle. Median p is 1.2 and Mach 720 for y < 0.01 at lvl 1, and p 4.4, Mach 380 for y = 0.01–0.1 (lvl 0: p 0.8 and 8.9, Mach 910 and 270). It heats further up: median p ≈ 70 at y = 0.1–0.3 at lvl 1.
- **Mach number.** The sonic Mach number M = |v|/√(γp/ρ) is meaningless at the nodes the limiter holds at p = ε_p: M = |v|√(ρ/(γε_p)), which is 8·10⁵ in the beam.
  - `Mach_f` = |v|/c_f with c_f = √((γp + |**B**|²)/ρ) ≥ |**B**|/√ρ stays bounded: 66.7 in the beam, at most 67.3 anywhere.
  - It is the Mach number of MHD shocks at β_a = 10⁻², so plot it instead of M. At t = 0.002 (lvl 0) the floor nodes at ambient density have M_f ≤ 0.17.

## References

- K. Wu, C.-W. Shu, SIAM J. Sci. Comput. 40 (2018) B1302–B1329, Example 5.6.
- D. S. Balsara, J. Comput. Phys. 231 (2012) 7504–7517.
- X. Zhang, C.-W. Shu, J. Comput. Phys. 229 (2010) 8918–8934.
