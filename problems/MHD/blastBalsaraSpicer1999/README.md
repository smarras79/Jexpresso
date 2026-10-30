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

**β** is the plasma beta, the ratio of gas pressure to magnetic pressure: β = p/p_B. **B** is in Heaviside–Lorentz units, where the magnetic pressure is p_B = ½|**B**|², so β = 2p/|**B**|². (In Gaussian units p_B = |**B**|²/(8π), and in SI p_B = |**B**|²/(2μ₀).) Equivalently, β = (2/γ) c_s²/v_A², with sound speed c_s = √(γp/ρ) and Alfvén speed v_A = |**B**|/√ρ. β ≪ 1 means the field dominates the pressure.

|**B**| = 100/√(4π) here is Balsara & Spicer's B = 100 in Gaussian units. The ambient β = 2·0.1/(100²/4π) = 2.51·10⁻⁴, as stated by Balsara (2004); inside the initial circle, β = 2.51. Other papers use [-0.5, 0.5]² with the same data, which is the same problem shifted.

The equations are the 9-field GLM-MHD system of `problems/MHD/orszagTangBormanis2024`. This case includes its `user_flux.jl`, `user_source.jl` and `user_bc.jl`, with γ = 1.4. The mesh is the same 32×32 periodic mesh, refined once by p4est: 64×64 elements at `:nop => 3`, i.e. 193×193 nodes, close to the paper's 200² zones. The pressure jump is sampled at the nodes as given, with no taper.

## Stabilization

As `problems/MHD/rotorDaoNazarov2022` (see its [algorithm.pdf](../rotorDaoNazarov2022/algorithm.pdf)), with two differences:

- conserved-form DynSGS: one ν per element, ∇·(ν∇q) on every conserved variable, the legacy sensor with element norms, `C_R = 1`, `C_max = 0.5`, and no background viscosity.
  - **Fast-speed floors** (`:dsgs_fast_floors => true`, DSGS.md §4.5). Each residual R_i is divided by D_i = max(spread_i, S_i). The floors S_i use the fast-speed bound c̄_f = √((γp̄ + |**B̄**|²)/ρ̄) of the element-mean state instead of the sound speed c̄ = √(γp̄/ρ̄): S_ρ = ρ̄, S_ρv = ρ̄ c̄_f, S_E = γp̄ + |**B̄**|², S_B = √(γp̄ + |**B̄**|²). In this ambient c_f/c = 75. With sound-speed floors, a resolved fast wave read 75× (ρv, **B**) to 5,700× (E) too strong, and ν sat at its cap over most of the disturbed region (see the comparison below).
  - **Thermal-energy residual** (`:dsgs_thermal_residual => true`, DSGS.md §4.5). The residual of ρe = E − ½|**m**|²/ρ − ½|**B**|² − ½ψ² joins the max (below).
  - **No startup hold** (`:dsgs_hold_steps => 0`). The pressure jump at t = 0 is the most violent moment of the run, and outside it the initial data are uniform, so the first residual is nonzero only at the jump.
- the conservative positivity limiter (`:positivity_method => "conservative"`), with ε_ρ = 10⁻⁶ (10⁻⁶ of the ambient ρ) and ε_p = 10⁻⁷ (10⁻⁶ of the ambient p). Elements with a node below ε are scaled toward their mean. Mass, momentum and energy are conserved, and **B** and ψ are untouched.

## Run

```bash
julia --project=. src/Jexpresso.jl MHD blastBalsaraSpicer1999
mpiexec -n 4 julia --project=. src/Jexpresso.jl MHD blastBalsaraSpicer1999
```

Outputs every 0.001 (VTK): ρ, u, v, p, **B**, ψ, T = p/ρ, `pmag` = ½|**B**|², `Mach`, `β` = 2p/|**B**|², `Mach_f` = |**v**|/c_f, and the DynSGS coefficient `mu_dsgs`. Mass and energy are tracked every step in `conservation.dat`. Plot commands:

```bash
D=output/MHD/blastBalsaraSpicer1999/output-<date>
python3 tools/plot_mhd_matrix.py $D --times 0.003 0.007 0.01 --vars rho p pmag Mach_f mu_dsgs \
        --log p mu_dsgs                                   # rows = times, columns = variables; --fieldlines adds B lines
python3 tools/plot_mhd_matrix.py $D --times 0.01 --vars rho p speed Bmag --log rho p   # Balsara (2004) Fig. 6
python3 tools/plot_fluxemergence_vtu.py $D --steps 11 --var rho --no-fieldlines \
        --xlabel '$x$' --ylabel '$y$' --time-unit ''      # one field per figure; drop --no-fieldlines for B lines
python3 tools/plot_conservation.py $D --out blast_conservation.pdf
```

## Results at t = 0.01

193² nodes, Δt = 10⁻⁵ (acoustic CFL 0.13 at the start, 0.11 at the end). The other rows are the same deck without the thermal residual, and with sound-speed floors instead of both:

| sensor | ρ | p | ½\|**B**\|² | \|**v**\| | max Mach (p > 10⁻³) | ν > 0.9 cap | mean ν |
|---|---|---|---|---|---|---|---|
| fast floors + thermal (this deck) | 0.2229 – 3.582 | 0.0959 – 246.1 | 230.6 – 589.8 | ≤ 16.5 | 9.13 (all nodes) | 8.2% | 0.0116 |
| fast floors | 0.2159 – 3.669 | 10⁻⁷ – 247.7 | 223.4 – 604.7 | ≤ 16.6 | 52.8 | 0.5% | 0.0044 |
| sound floors | 0.2253 – 3.422 | 0.09995 – 242.4 | 234.0 – 574.1 | ≤ 16.3 | 5.41 | 58.7% | 0.0438 |

The dense shells sit along **B**, at x ≈ 0.2 and 0.8. The fast front in magnetic pressure spans y ≈ 0.12 – 0.88, and **B** is expelled from the hot interior. This is the structure of Balsara (2004) Fig. 6. With fast floors, ν is concentrated on the shells and on the thin fast-front ring, and is near zero elsewhere, so the shells are sharper (ρ_max 3.67 against 3.42). ψ stays below 0.11 (|**B**| ≈ 28).

**Without the thermal residual: a pressure undershoot at the fast front.** In this ambient the thermal energy is 0.06% of E, so p is a small remainder of E and the magnetic energy. With sound floors, the saturated ν damped the dispersive ripples of the fast front, and no node at t = 0.01 was below the ambient p = 0.1. With fast floors, a thin band at the front, deepest where the shells meet it at x ≈ 0.1 and 0.9, falls below ambient:

| p below | 0.0999 | 0.09 | 0.05 | 10⁻³ | 10⁻⁶ |
|---|---|---|---|---|---|
| nodes (of 37,249) | 5,852 | 786 | 220 | 50 | 30 |

There, Mach exceeds 6 at 2,136 nodes. The ~50 nodes at ε are where the limiter holds p. There M = |**v**|√(ρ/(γε)) is 3·10³ – 1.2·10⁴ at every output. `Mach_f` = |**v**|/c_f stays below 0.9 to t = 0.002, and is 0.04 – 0.06 at those nodes, so plot it instead.

**Thermal-energy residual** (`:dsgs_thermal_residual => true`, on in this deck; DSGS.md §4.5). It adds to the max the residual of ρe = E − ½|**m**|²/ρ − ½|**B**|² − ½ψ², R_ρe = R_E − **v**·R_**m** + ½|**v**|² R_ρ − **B**·R_**B** − ψR_ψ, over max(‖ρe − ⟨ρe⟩_K‖_∞,K, γp̄_K). At fixed E an error δB costs δp/p = −(γ−1)(2/β) δ|**B**|/|**B**| = −3200 δ|**B**|/|**B**| here, so a 0.03% ripple at the fast front empties p. The E residual, measured against |**B**|², does not see it; R_ρe does. 128², same Δt:

| | t | nodes p < 10⁻⁶ | p < 0.0999 | max M | max M_f | ρ_max | mean ν | limiter active (stage calls) |
|---|---|---|---|---|---|---|---|---|
| 64², nop 3, fast floors | 0.01 | 30 | 15.7% | 4,123 | 1.11 | 3.67 | 4.4·10⁻³ | every stage call |
| 64², nop 3, + thermal | 0.01 | 0 | 1.9% | 9.1 | 1.09 | 3.58 | 1.2·10⁻² | 184 (startup) |
| nop 3, fast floors | 0.01 | 88 | 15.5% | 12,830 | 1.23 | 4.37 | 1.3·10⁻³ | 9,982 of 10,000 |
| nop 3, + thermal | 0.01 | 0 | 1.2% | 11.3 | 1.22 | 4.29 | 3.0·10⁻³ | 217 (all before t = 0.003) |
| nop 4, fast floors | 0.001 | 48 | 2.8% | 6,257 | — | 2.90 | 4.9·10⁻⁴ | 2,000 of 2,000 |
| nop 4, + thermal | 0.001 | 0 | 0.2% | 4.3 | — | 2.77 | 9.6·10⁻⁴ | 385 (startup) |
| 256², nop 4, + thermal | 0.002 | 0 | 0.12% | 8.8 (11.6 at 0.0015) | 1.21 | 4.57 | 3.2·10⁻⁴ | 1,025 of 8,000 (mostly startup) |

No node goes below 10⁻³ at any output. At 128² the limiter no longer acts after the startup; at 256², nop 4 (`:init_refine_lvl => 3`, `:nop => 4`, `:Δt => 1.25e-6`, 2 h on 3 ranks to t = 0.002) it acts in a few stage calls later on, and without the thermal residual the same run reaches Mach ≈ 8,100 at the fast front by t = 0.001. The cost is about twice the mean ν and 2–5% off ρ_max. The startup is unchanged: the initial jump still gives 96 (nop 3) and 628 (nop 4) inadmissible element means.

**Resolution.** The same deck with `:init_refine_lvl => 2` and `:Δt => 5.0e-6` gives 128×128 elements (385² nodes). It takes 15 min on 4 ranks. Shares of nodes at t = 0.01:

| | p < 0.0999 | p < 0.09 | p < 0.05 | p < 10⁻³ | elements limited per stage call | startup energy from node floors | ρ_max |
|---|---|---|---|---|---|---|---|
| 64² | 15.7% | 2.11% | 0.59% | 0.13% | 1.5% | 1.8·10⁻⁵ (steps 1–7) | 3.67 |
| 128² | 15.5% | 1.14% | 0.29% | 0.086% | 0.75% | 1.6·10⁻⁶ (step 1) | 4.37 |
| 64², + thermal | 1.9% | 0 | 0 | 0 | 0 after startup | 1.75·10⁻⁵ | 3.58 |
| 128², + thermal | 1.25% | 0 | 0 | 0 | 0 after startup | 1.7·10⁻⁶ | 4.29 |

Without the thermal residual, the deep undershoot is a grid-scale band that halves with h, and the shallow one, less than 1% below the ambient over 15% of the domain, does not change with resolution. With it, no node is 10% below the ambient at either resolution, and the shallow band shrinks to 1–2% of the nodes.

**Positivity and conservation** (64², this deck).
- **Limiter activity:** the limiter acts only during the startup, in 184 of the first 2,000 stage calls (3,000 element limitings). Without the thermal residual it acts in every stage call.
- **Startup floors:** in steps 1–7, 180 elements next to the initial jump have an inadmissible mean, i.e. p̄ ≤ ε even before limiting. There the non-conservative node floor raises p at 574 nodes, which adds 1.75·10⁻⁵ of the total energy.
- **After step 7:** energy holds to 4.8·10⁻¹⁶ and mass to 3·10⁻¹⁵ for the rest of the run.
- **With sound floors:** the limiter stopped by step 200, and there were 554 startup floors (+1.5·10⁻⁵ of energy).

**Why only the blast's floor adds energy.** The node floor adds energy only where a node's kinetic and thermal energy are used up together. Where the gas moves, it removes momentum first, which leaves E unchanged (`cpos_floor!`, `src/kernel/positivity/Positivity.jl`). At t = 0 the blast's gas is at rest. Its only non-magnetic energy is the thermal 6.3·10⁻⁴ of E, the smallest of any initial state among the MHD cases. The pressure jump is also sampled at the nodes without smoothing. The rotor has neither: β ≈ 1, uniform initial pressure and a tapered disc. It never has an inadmissible element mean, and E stays at round-off. The jet is also strongly magnetized (β = 10⁻²), but its beam is ramped in smoothly in space and time, and the floor never fires either.

The table gives the initial states of the MHD cases, sorted by β. Both value columns are the minimum over the domain, with e_th = p/(γ − 1) and E = e_th + ½ρ|**u**|² + ½|**B**|². Unlike β, e_th/E includes the kinetic energy, so it depends on the frame.

| case | β | e_th/E | notes |
|---|---|---|---|
| **MHD blast wave** | **2.5·10⁻⁴** | **6.3·10⁻⁴** | ambient; inside the initial disc β = 2.5, e_th/E = 0.86 |
| magnetized astrophysical jet | 10⁻² | 2.4·10⁻² | ambient; the beam entering through the nozzle has β = 10⁻², e_th/E = 5.6·10⁻⁶ |
| Brio–Wu shock tube | 0.13 | 0.11 | right state; left state β = 1.28, e_th/E = 0.56 |
| circularly polarized Alfvén wave | 0.20 | 0.23 | uniform |
| MHD rotor | 1.0 | 0.11 | β uniform; e_th/E smallest at the edge of the spinning disc, 0.72 in the ambient |
| flux emergence | 1 | 0.95 | at the centre of the field sheet (γ = 1.05); β → ∞ and e_th/E → 1 away from it |
| Orszag–Tang vortex | 5/3 | 0.43 | domain totals β = 10/3, e_th/E = 0.57 |
| shock–cloud, MHD letters | 2.1 | 0.50 | β ahead of the shock (and in the letters), where e_th/E = 0.76; e_th/E behind it, where β = 50 |
| shock–cloud, round cloud | 4.2 | 1.3·10⁻³ | in the cloud, which moves at u = −15.27; ahead of the shock outside the cloud e_th/E = 0.013; behind it β = 222, e_th/E = 0.997 |
| smooth vortex | 78 | 0.54 | in the vortex; far field β → ∞, e_th/E = 0.60 |
| Kelvin–Helmholtz instability | 128 | 0.85 | β uniform; e_th/E smallest in the dense layer |

Notes on the values:

- **Kinetic-dominated cases.** The round cloud and the jet's beam have small thermal fractions too, but from their kinetic energy.
- **Flux-emergence β.** β = 1 at the sheet centre is the deck's inferred `fe_beta_star`.
- **Jet cases.** The jet row is Wu & Shu's case (i). Cases (ii) and (iii) have β = 10⁻³ and 10⁻⁴.

LaTeX for the paper (use `table*` in two-column REVTeX):

```latex
\begin{table}
\caption{Plasma $\beta = 2p/|\mathbf{B}|^2$ and thermal energy fraction
$e_{\mathrm{th}}/E$ of the initial conditions of the MHD test cases, with
$e_{\mathrm{th}} = p/(\gamma - 1)$ and
$E = e_{\mathrm{th}} + \tfrac12\rho|\mathbf{u}|^2 + \tfrac12|\mathbf{B}|^2$
(Heaviside--Lorentz units). Both columns give the minimum over the domain; the
notes give the values elsewhere. Unlike $\beta$, $e_{\mathrm{th}}/E$ includes
the kinetic energy and therefore depends on the frame.}
\label{tab:beta}
\begin{tabular}{lllp{6.2cm}}
\hline
Case & $\beta$ & $e_{\mathrm{th}}/E$ & Notes \\
\hline
MHD blast wave                     & $2.5\times10^{-4}$ & $6.3\times10^{-4}$ & ambient; inside the initial disc $\beta = 2.5$, $e_{\mathrm{th}}/E = 0.86$ \\
Magnetized astrophysical jet       & $10^{-2}$          & $2.4\times10^{-2}$ & ambient; the beam entering through the nozzle has $\beta = 10^{-2}$, $e_{\mathrm{th}}/E = 5.6\times10^{-6}$ \\
Brio--Wu shock tube                & $0.13$             & $0.11$             & right state; left state $\beta = 1.28$, $e_{\mathrm{th}}/E = 0.56$ \\
Circularly polarized Alfv\'en wave & $0.20$             & $0.23$             & uniform \\
MHD rotor                          & $1.0$              & $0.11$             & $\beta$ uniform; $e_{\mathrm{th}}/E$ smallest at the edge of the spinning disc, $0.72$ in the ambient \\
Flux emergence                     & $1$                & $0.95$             & at the centre of the field sheet ($\gamma = 1.05$); $\beta \to \infty$ and $e_{\mathrm{th}}/E \to 1$ away from it \\
Orszag--Tang vortex                & $5/3$              & $0.43$             & domain totals $\beta = 10/3$, $e_{\mathrm{th}}/E = 0.57$ \\
Shock--cloud, MHD letters          & $2.1$              & $0.50$             & $\beta$ ahead of the shock (and in the letters), where $e_{\mathrm{th}}/E = 0.76$; $e_{\mathrm{th}}/E$ behind it, where $\beta = 50$ \\
Shock--cloud, round cloud          & $4.2$              & $1.3\times10^{-3}$ & in the cloud, which moves at $u = -15.27$; ahead of the shock outside the cloud $e_{\mathrm{th}}/E = 0.013$; behind it $\beta = 222$, $e_{\mathrm{th}}/E = 0.997$ \\
Smooth vortex                      & $78$               & $0.54$             & in the vortex; far field $\beta \to \infty$, $e_{\mathrm{th}}/E = 0.60$ \\
Kelvin--Helmholtz instability      & $128$              & $0.85$             & $\beta$ uniform; $e_{\mathrm{th}}/E$ smallest in the dense layer \\
\hline
\end{tabular}
\end{table}
```

**Sensor: legacy, not residual.** The deck keeps `:dsgs_sensor => "legacy"`, the sensor the rotor setup it builds on was tuned with. The residual sensor (the default, used by the smooth vortex and the Alfvén wave) was tried in its place, everything else in the deck unchanged, to t = 0.01 (21 Oct 2026). The thermal residual works with either: it is built from the same per-node residual as the main sensor.

| | 64², legacy | 64², residual | 128², legacy | 128², residual |
|---|---|---|---|---|
| ρ | 0.2229 – 3.582 | 0.2172 – 3.949 | 0.2023 – 4.288 | 0.1942 – 4.648 |
| p min | 0.0959 | 0.0954 | 0.0947 | 4.3·10⁻⁷ |
| nodes with p < 0.0999 | 1.9% | 1.9% | 1.25% | 1.30% |
| nodes with p < 0.09 / 0.05 / 10⁻³ | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 12 / 7 / 2 |
| max M (p > 10⁻³) | 9.13 | 7.66 | 11.3 | 9.9 |
| mean ν | 1.16·10⁻² | 0.98·10⁻² | 3.0·10⁻³ | 2.6·10⁻³ |
| limiter active (stage calls) | 184, startup | 476 | 217, startup | 1,823 of 10,000, to the end |
| element limitings | 3,000 | 4,094 | 5,336 | 10,229 |
| inadmissible means / node floors | 180 / 574 | 88 / 308 | 96 / 288 | 48 / 192 |
| energy added by the node floors | 1.75·10⁻⁵ | 1.42·10⁻⁵ | 1.7·10⁻⁶ | 2.1·10⁻⁶ |

At 64² the residual sensor looks as good or better: sharper shells, a lower Mach peak, half the startup trouble. At 128² it is not. Its ν, 16% below the legacy one, no longer damps the resolved ripple of the fast front: two nodes fall to p ≈ 4·10⁻⁷ and a dozen 10–50% below the ambient, the undershoot the thermal residual was added to remove, and the limiter acts in 18% of the stage calls to the end instead of only at the initial jump. Mass holds to 8·10⁻¹⁵ either way.

**Why the sound floors saturate here.** The table counts the share of nodes where ν exceeds a given fraction of its local cap C_max Δ (|**v**| + c_f):

| | ν > 0.9 cap | ν > 0.5 cap | ν > 0.1 cap |
|---|---|---|---|
| blast, sound floors, t = 0.01 | 58.7% | 65.8% | 87.6% |
| blast, fast floors, t = 0.01 | 0.5% | 1.8% | 17.8% |
| rotor 128² (sound floors, β ≈ 1), t = 0.15 | 1.6% | 9.5% | 18.3% |

## β sweep

Only the ambient field changes: `:blast_B0` in `user_inputs.jl`, with ambient β = 2p/B0² = 0.2/B0². Steps of √10 in B0 give one decade in β per run:

| k | `:blast_B0` | β |
|---|---|---|
| 0 | `100.0/sqrt(4.0*π)` = 28.21 (the paper) | 2.5·10⁻⁴ |
| 1 | `100.0/sqrt(4.0*π)/sqrt(10.0)` = 8.921 | 2.5·10⁻³ |
| 2 | `100.0/sqrt(4.0*π)/10.0` = 2.821 | 2.5·10⁻² |
| 3 | `100.0/sqrt(4.0*π)/sqrt(1000.0)` = 0.8921 | 0.25 (≈ the moderate blast of Gardiner & Stone) |
| 4 | `100.0/sqrt(4.0*π)/100.0` = 0.2821 | 2.5 |

- **Everything else stays fixed:** Δt, ε, the floors and the mesh. Each run writes its own output directory.
- **Δt:** wave speeds only fall as B0 drops, so Δt = 10⁻⁵ holds for the whole table. For B0 > 28.2, scale Δt by 28.2/B0; the acoustic CFL in the log should stay ≲ 0.15.
- **Floors:** with `:dsgs_fast_floors` they reduce to the sound-speed floors at high β, so the sweep varies the physics, not the normalization.
- **Initial log:** each run prints its ambient |B| and β.

One row per β at t = 0.01:

```bash
python3 tools/plot_mhd_matrix.py <run k=4> <run k=3> <run k=2> <run k=1> <run k=0> --times 0.01 \
        --labels '$\beta = 2.5$' '$\beta = 0.25$' '$\beta = 2.5\times10^{-2}$' '$\beta = 2.5\times10^{-3}$' '$\beta = 2.5\times10^{-4}$' \
        --vars rho p pmag mu_dsgs --log p mu_dsgs
```

For each run, the log's `POSITIVITY LIMITER` lines give the limiter activity, and `conservation.dat` gives the energy that the startup node floors add.

**Results at t = 0.01, 64², with fast and with sound floors.** The disturbed region is the set of nodes the blast has reached (|ρ − 1| > 0.01 or |**v**| > 0.1). Each cell gives fast / sound:

| ambient β | mean ν, disturbed | mean ν, ambient | disturbed nodes at the cap | nodes with p < 0.5 p_amb |
|---|---|---|---|---|
| 2.5 | 0.0205 / 0.0206 | 5.0·10⁻⁴ / 5.3·10⁻⁴ | 29% / 29% | 0.25% / 0.25% |
| 0.25 | 0.0204 / 0.0205 | 5.9·10⁻⁴ / 1.1·10⁻³ | 28% / 29% | 0.27% / 0.27% |
| 0.025 | 0.0194 / 0.0204 | 6.7·10⁻⁴ / 2.6·10⁻³ | 25% / 26% | 0.16% / 0.16% |
| 2.5·10⁻³ | 0.0151 / 0.0228 | 7.6·10⁻⁴ / 7.0·10⁻³ | 16% / 32% | 0.19% / 0.22% |
| 2.5·10⁻⁴ | 0.0069 / 0.0504 | 5.2·10⁻⁴ / 3.0·10⁻² | 0.9% / 71% | 0.59% / 0 |

- **High β:** the two floors coincide. At β = 2.5 ν is at the cap on the shell (63% of its nodes, a Mach ≈ 40 shock), about 17% of the cap in the interior, and ~4·10⁻⁴ in the ambient.
- **Falling β:** fast floors keep the ambient ν flat, while sound floors raise it up to 58×.
- **Pressure undershoot:** for β ≥ 2.5·10⁻³ it is the same with either floor (18 – 48 nodes at ε, at the foot of the outgoing shock). Only the saturated sound-floor run at β = 2.5·10⁻⁴ has none.

## References

- D. S. Balsara, D. S. Spicer, J. Comput. Phys. 148 (1999) 133–148 (1999b).
- D. S. Balsara, ApJS 151 (2004) 149–184, §7.3.
