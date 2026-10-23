# Mach 15 fast shock hitting a round cloud (Dai & Woodward 1998, Figs. 23–25)

A Mach 15 MHD fast shock hits a round cloud 10 times denser than the gas ahead of it. Dai & Woodward (1998, JCP 142:331, p. 365) use this test to compare their scheme in two forms: one keeps ∇·**B** = 0 exactly, the other only to second order. Figs. 23–25 show ρ, B_x and B_y from both.

The physics, the boundaries and the stabilization are those of `problems/MHD/shockCloudDaiWoodward1998` (the letters "MHD", Fig. 20). This directory holds only its own deck, mesh and output hook. It includes the rest from that case, and `:sc_problem => "cloud"` selects this test. See that README for the units, the boundary condition (with a LaTeX version for the paper, under "Boundaries") and the stabilization.

## Setup

DW's text gives the states, the domain and the grid:

| | (ρ, p, u, v, w, B_x, B_y, B_z) |
|---|---|
| pre-shock, x > x_s | (1, 1, −15.2677, 0, 0, 2, 1, 1) |
| post-shock, x < x_s | (3.94605, 312.063, 0, −0.0229969, −0.0229969, 2, 3.95492, 3.95492) |
| cloud | ρ = 10, with the pre-shock p, **v** and **B** |

| | |
|---|---|
| domain | [0, 1.5] × [0, 1], open on all four sides |
| γ | 5/3 |
| DW's grid | 384 × 256 cells; here 96 × 64 elements at N = 4, i.e. 385 × 257 nodes |

**The frame.** The shocked gas is at rest. The gas ahead of the shock, cloud included, flows in from the right at u = −15.27. In this frame the shock moves right at s = (ρ₂u₂ − ρ₁u₁)/(ρ₂ − ρ₁) = 5.18. Against the gas ahead it moves at 20.45: 14.0 times the fast speed (1.46) and 15.8 times the sound speed (1.29).

**Units of B.** As in the letters case, **B** is DW's Gaussian field divided by √(4π). With that conversion the states satisfy the Rankine–Hugoniot conditions to 5 digits. The gas ahead enters faster than its fast speed (|u| = 15.27 > c_f = 1.46), so the right side is a super-fast inflow. c_h = 16.7.

**What the paper does not print.** DW give neither the cloud's position and radius nor the time of Figs. 23–25 for this test. The deck uses:

- **Cloud centre (0.8, 0.5), radius 0.15.** DW print these for their Mach 10 cloud (Fig. 18), which has the same density ratio and the same frame, in a unit square.
- **Shock at x_s = 0.6.** Tóth (2000, §6.5) uses this for that test.
- **Final time t = 0.16.** In Fig. 23 the main shock is at x ≈ 1.43 (measured on the scan, 299 px per unit). It starts at 0.6 and moves at 5.18, so it reaches 1.43 at t = (1.43 − 0.6)/5.18 = 0.160. The shocked cloud in the figure lies at x = 0.58–0.75, y = 0.29–0.71, which agrees with a cloud hit at x ≈ 0.65.

**Sharp versus spread jumps.** `:sc_smooth => 1/64` (one element) spreads the shock and the cloud edge with the same C² ramp as the letters case. `0` gives DW's sharp data.

## Tests here (`SC15_48x32`, Δt = 4·10⁻⁵, `:sc_smooth => 1/32`, 4 threads)

Run to t = 0.03 (about 300 s):

- the cloud is compressed to ρ ≈ 38;
- p peaks at ≈ 830 behind the reflected shock;
- p touches the 10⁻⁶ floor inside the cloud, and the conservative limiter repairs it;
- the reflected bow shock (ρ ≈ 5.8, u ≈ −5) is at x ≈ 0.3 at t = 0.03 and still moving left.

The open boundaries hold the undisturbed states: the left side stays at ρ = 3.944–3.947 with |u| < 0.005. With the plain do-nothing condition, that side had emptied by then (see the letters README, "Boundaries").

## Run

```bash
julia --project=. src/Jexpresso.jl MHD shockCloudMach15DaiWoodward1998
```

Output is VTK every 0.02 to t = 0.16. The fields are those of the letters case.

**All the paper's figures of both shock-cloud cases in one go:**

```bash
python3 tools/plot_shock_cloud_figures.py --letters output/MHD/shockCloudDaiWoodward1998/<run> \
                                          --cloud   output/MHD/shockCloudMach15DaiWoodward1998/<run>
```

This writes `letters_{fields,paper,fieldlines,shockmach}` and `cloud_{fields,final,fieldlines,shockmach}` to `figures/shockCloud/` as PDF, printing each `plot_mhd_matrix.py` command as it runs it. Useful options:

- `--letters-times` and `--cloud-times` set the rows; the times must be output times of the run;
- `--out-dir` and `--format` set where the figures go and in what format;
- `--dry-run` prints the commands without running them.

The commands below are the individual figures.

To draw Fig. 23 (ρ, light where dense) and the fields of Figs. 24–25:

```bash
D=output/MHD/shockCloudMach15DaiWoodward1998/output-<date>
python3 tools/plot_mhd_matrix.py $D --times 0.16 --vars rho Bx By mu_dsgs --cmaps rho gray --log mu_dsgs
python3 tools/plot_mhd_matrix.py $D --times 0.04 0.08 0.12 0.16 \
        --vars rho log10rho p Mach_fs mu_dsgs --log rho p --cmaps rho gray \
        --fieldlines-on log10rho --nlevels 60 --cmaps Mach_fs YellowGreyBlue --clim Mach_fs 0 15
        # log ρ | log ρ with the field lines | log p | shock-frame fast Mach (15.0 ahead of the shock = DW's "Mach 15") | ν
```

Mach number in the shock's frame (s = 5.18, read from the run): ahead of the shock `Mach_s` = 15.8 (sonic) and `Mach_fs` = 15.0, DW's "Mach 15" (see the letters README):

```bash
python3 tools/plot_mhd_matrix.py $D --times 0.04 0.08 0.12 0.16 \
        --vars rho Mach_s Mach_fs --log rho --cmaps rho gray --clim Mach_s 0 18 --clim Mach_fs 0 18
```

Magnetic field lines as in DW's Fig. 19 (DW show them for their Mach 10 cloud), alone and zoomed on the cloud, or over ρ:

```bash
python3 tools/plot_mhd_matrix.py $D --times 0.16 --vars Blines --nlevels 100 --xlim 0.1 0.9 --ylim 0.1 0.9
python3 tools/plot_mhd_matrix.py $D --times 0.04 0.08 0.12 0.16 --vars rho Blines --log rho --cmaps rho gray --fieldlines-on rho --nlevels 80
```

See the letters README for the options.

To regenerate the meshes:

```bash
gmsh -2 SC15.geo -setnumber ref 1 -o SC15_48x32.msh
gmsh -2 SC15.geo -setnumber ref 2 -o SC15_96x64.msh
```
