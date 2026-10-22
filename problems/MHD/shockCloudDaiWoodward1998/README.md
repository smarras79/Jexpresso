# Fast shock hitting the clouds "MHD" (Dai & Woodward 1998, Fig. 20)

A Mach 12 MHD fast shock runs to the right into three clouds 50 times denser than the gas ahead of it. The clouds are shaped as the letters **MHD**. Dai & Woodward (1998, JCP 142:331, p. 360) ran it to test the robustness of their scheme. Tóth (2000, JCP 161:605, §6.5) uses the single round cloud of their Fig. 18 instead of the letters.

DW's Mach 15 round-cloud test (Figs. 23–25) is `problems/MHD/shockCloudMach15DaiWoodward1998`. It shares this directory's code through `:sc_problem => "cloud"`.

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
- **Letters.** The letters are Times Bold. Liberation Serif Bold has the same metrics as Times. Its outlines of M, H and D, scaled to the cap height of the first panel, lie on the figure's H and D to within a pixel. The letters sit between y = 0.575 and 0.970. H starts at x = 0.814, D at x = 1.284, and M ends at x = 0.785. M is placed by its right edge, because in that panel the shock has already swept its left side. The resulting M spans x = 0.243 to 0.785.

**The letters are part of the initial density, defined in `initialize.jl`.** `SC_LETTERS` holds each letter's outline as a polygon in the domain's coordinates: 25 vertices for M, 28 for H, and 38 + 32 for D and its counter. The curved parts are flattened into 6 straight segments per quadratic arc. Each node's density comes from its coordinates. A point-in-polygon test (odd number of crossings) decides whether the node is inside a letter, and the distance to the nearest edge gives the edge ramp. So the letters are the same on any mesh or polynomial order; nothing is drawn on the grid. The vertices were taken once from the font's outlines with Python's fontTools and are written into `initialize.jl` as numbers, so no tool is needed to run the case.

`:sc_shock_x` sets the shock position. `:sc_smooth` sets the width over which the shock and the letters' edges are spread. The spread is a C² ramp (smootherstep) of x − x_s for the shock and of the signed distance to the outline for the letters, applied to the primitive variables. The deck uses one element, 1/64 on `SC_128x96` (use 1/32 on `SC_64x48`). `:sc_smooth => 0` gives DW's sharp data.

**Why the jumps are spread.** On `SC_64x48` with sharp data the run aborted at t = 1.1·10⁻³, after 35 steps:

- the shock reached the thin left stem of the M (x ≈ 0.24–0.34);
- ν sat at its cap there;
- ρ and p fell to the limiter floors;
- then 13% of the solution went non-finite.

With the one-element ramp the same mesh ran to t = 0.01, the furthest it was taken here. In the stem's wake p still touches the 10⁻⁶ floor now and then, and the conservative limiter repairs it. p peaks near 990 by t = 0.01, behind the shock reflected off the letters (the post-shock p is 254).

## Boundaries

All four sides use one characteristic far-field condition (`user_bc.jl`). At each boundary node the deviation q − q∞ from a far-field state is split into waves. The split uses the eigenvectors of the Jacobian of the normal GLM-MHD flux at q∞. Waves leaving the domain keep their interior values; waves entering it take the far-field values:

q_b = q∞ + Σ_{λ_k > 0} r_k l_k (q − q∞)

- **The far field q∞** is DW's undisturbed planar shock: the post-shock state for x < x_s + s t, the pre-shock state ahead of it. It is exact until waves from the clouds reach the boundary.
- **Super-fast inflow** (the left side here): no MHD wave leaves, so the whole state is imposed. The one exception is the outgoing GLM wave, whose speed c_h exceeds |v|. Tóth fixes his inflow the same way.
- **Super-fast outflow:** every wave leaves, so nothing is imposed.
- **Strong waves:** the condition is linear about q∞, so a strong wave leaving through it, such as a bow shock, passes only approximately.
- **Across the shock,** where it meets the top and bottom, the far field is not known. Nodes within |x − x_s − s t| < 2h + `:sc_smooth` get only ψ = 0 there.

**Why not the plain open condition.** The do-nothing condition imposes nothing and leaves only ψ = 0. That suits a super-fast outflow. In gas at rest it is unstable. In the Mach 15 case (`SC15_48x32`), the post-shock gas at rest on the left side started to drain, uniformly in y:

- u went from −0.007 at t = 0.01 to −2.9 at t = 0.015;
- by t = 0.03 the density within 0.25 of the side had fallen by 30%, before any wave from the cloud had arrived.

With the characteristic condition and the shock band, that side stays at ρ = 3.944–3.947, |u| < 0.005 at t = 0.03. In this case at t = 0.01 on `SC_64x48`, the left and right sides hold their states exactly. The top and bottom stay within the two states, apart from the shock band.

## Stabilization

This case uses the deck of `problems/MHD/blastBalsaraSpicer1999`, with γ = 5/3:

- conserved-form DynSGS with the legacy sensor and element norms;
- fast-speed floors and the thermal-energy residual;
- no startup hold;
- the conservative positivity limiter, with ε_ρ = ε_p = 10⁻⁶.

The equations are the 9-field GLM-MHD system of `problems/MHD/orszagTangBormanis2024`. This case includes that case's flux and its ψ damping. The GLM cleaning speed is c_h = 24.4, the |**v**| + c_f of the post-shock gas.

The files are shared with the Mach 15 case. `user_flux.jl` holds both tests' states (`SC_STATES`), and `:sc_problem` picks one.

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
