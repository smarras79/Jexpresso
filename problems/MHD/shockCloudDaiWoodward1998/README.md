# Fast shock hitting the clouds "MHD" (Dai & Woodward 1998, Fig. 20)

A Mach 12 MHD fast shock runs to the right into three clouds 50 times denser than the gas ahead of it. The clouds are shaped as the letters **MHD**. Dai & Woodward (1998, JCP 142:331, p. 360) ran it to test the robustness of their scheme. Tóth (2000, JCP 161:605, §6.5) uses the single round cloud of their Fig. 18 instead of the letters.

## Setup

DW's text gives the states, the domain, the grid and the boundaries:

| | (ρ, p, u, v, w, B_x, B_y, B_z) |
|---|---|
| pre-shock, x > x_s | (1, 1, 0, 0, 0, 2, 2, 2) |
| post-shock, x < x_s | (3.91109, 253.724, 13.8388, −0.0500186, −0.0500186, 2, 7.84323, 7.84323) |
| letters | ρ = 50, with the pre-shock p, **v** and **B** |

| | |
|---|---|
| domain | [0, 2] × [0, 1.5], open on all four sides |
| γ | 5/3 |
| DW's grid | 512 × 384 cells; here 128 × 96 elements at N = 4, i.e. 513 × 385 nodes |
| output times | 0.02507, 0.05004, 0.08003, 0.10411 (the panels of Fig. 20) |

**Units of B.** DW write B in Gaussian units, with magnetic pressure |**B**|²/8π. Jexpresso uses Heaviside–Lorentz units, with ½|**B**|², so **B** enters divided by √(4π). Tóth makes the same conversion for DW's single-cloud test: 7.73718/√(4π) = 2.1826. ρ, p and **v** do not change. With this conversion the two states satisfy the Rankine–Hugoniot conditions to 6 digits, with shock speed s = ρ₂u₂/(ρ₂ − ρ₁) = 18.59. The post-shock flow is super-fast, u = 13.84 > c_f = 10.52.

**What the paper does not print.** The paper does not give the initial shock position or the shape of the letters. Both are read off Fig. 20 (scan at 200 dpi, 240 px per unit):

- **Shock position x_s = 0.2.** On the first panel the incident shock is a faint step in the gray level at x = 0.664. That gives x_s = 0.664 − 18.59 × 0.02507 = 0.198.
- **Letters.** The letters are Times Bold. Liberation Serif Bold has the same metrics as Times. Its outlines of M, H and D, scaled to the cap height of the first panel, lie on the figure's H and D to within a pixel. The letters sit between y = 0.575 and 0.970. H starts at x = 0.814, D at x = 1.284, and M ends at x = 0.785. M is placed by its right edge, because in that panel the shock has already swept its left side. The resulting M spans x = 0.243 to 0.785. `mhd_letters.jl` holds the polygons. `tools/shock_cloud_letters.py` writes them; it needs fontTools.

`:sc_shock_x` sets the shock position. `:sc_smooth` sets the width over which the shock and the letters' edges are spread. The spread is a C² ramp (smootherstep) of x − x_s for the shock and of the signed distance to the outline for the letters, applied to the primitive variables. The deck uses one element, 1/64 on `SC_128x96` (use 1/32 on `SC_64x48`). `:sc_smooth => 0` gives DW's sharp data.

**Why the jumps are spread.** On `SC_64x48` with sharp data the run aborted at t = 1.1·10⁻³, after 35 steps:

- the shock reached the thin left stem of the M (x ≈ 0.24–0.34);
- ν sat at its cap there;
- ρ and p fell to the limiter floors;
- then 13% of the solution went non-finite.

With the one-element ramp, the same mesh ran through t = 0.005 and on to t = 0.01. That is as far as it was taken here. In the stem's wake p still touches the 10⁻⁶ floor now and then, and the conservative limiter repairs it. p peaks near 650, behind the shock reflected off the letters (the post-shock p is 254).

**Other DW tests.** The Mach 15 single-cloud comparison of DW's Figs. 23–25 uses a different frame: the pre-shock gas moves at −15.27 and the post-shock gas is at rest. This case does not set it up.

## Boundaries

| tag | condition |
|---|---|
| `left` | the post-shock state, imposed in full. Every characteristic enters at a super-fast inflow; Tóth fixes his inflow the same way. |
| `right`, `top`, `bottom` | open (do-nothing), with ψ = 0, as on the outflow sides of `problems/MHD/astroJetWuShu2018` |

## Stabilization

This case uses the deck of `problems/MHD/blastBalsaraSpicer1999`, with γ = 5/3:

- conserved-form DynSGS with the legacy sensor and element norms;
- fast-speed floors and the thermal-energy residual;
- no startup hold;
- the conservative positivity limiter, with ε_ρ = ε_p = 10⁻⁶.

The equations are the 9-field GLM-MHD system of `problems/MHD/orszagTangBormanis2024`. This case includes that case's flux and its ψ damping. The GLM cleaning speed is c_h = 24.4, the |**v**| + c_f of the post-shock gas.

## Run

```bash
julia --project=. src/Jexpresso.jl MHD shockCloudDaiWoodward1998
```

For a quick look, set `:gmsh_filename` to `SC_64x48.msh`, `:Δt => 3.2e-5` and `:sc_smooth => 1/32`. That run takes about 0.6 s per step on 4 threads, or about 30 min to t = 0.104.

Output is VTK at DW's four times. The fields are ρ, u, v, w, p, **B**, ψ, T, `pmag`, `Mach`, `β`, `Mach_f`, `log10ρ` and `mu_dsgs`. Mass and energy totals are written every 10 steps to `conservation.dat`; they change through the open sides.

To draw Fig. 20, gray with dense regions dark:

```bash
D=output/MHD/shockCloudDaiWoodward1998/output-<date>
python3 tools/plot_mhd_matrix.py $D --times 0.02507 0.05004 0.08003 0.10411 --vars log10rho mu_dsgs \
        --cmaps log10rho gray_r --log mu_dsgs
```

To regenerate the meshes:

```bash
gmsh -2 SC.geo -setnumber ref 1 -o SC_64x48.msh
gmsh -2 SC.geo -setnumber ref 2 -o SC_128x96.msh
```
