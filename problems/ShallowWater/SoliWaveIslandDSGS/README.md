# SoliWaveIslandDSGS — solitary wave on a conical island, stabilized by DynSGS

The `SoliWaveIsland` case (Marras, Kopera, Constantinescu, Suckale, Giraldo,
*Adv. Water Resour.* 114 (2018) 45–63, §5.5: a Synolakis solitary wave
running up and down a conical island in a closed basin, 2D non-linear
shallow water, CG spectral elements of order 4) with the constant artificial
viscosity of the sibling (`AV()`, μ = 0.05 on the three equations) replaced
by the residual-based **DynSGS** of the Euler and MHD cases, in a kernel for
the shallow-water system (`:visc_model => DSGS_SW()`,
`compute_dsgs_viscosity!(::DSGS_SW, ::NSD_2D)` in `kernel/physics/SGS.jl`).
Everything else — initial condition, well-balanced flux/source split,
wet/dry treatment, walls, mesh — is the sibling's, copied file by file.

```bash
julia --project=. src/Jexpresso.jl ShallowWater SoliWaveIslandDSGS
```

## The coefficient

A self-contained description of the shallow-water kernel, equation by
equation, is in [`../SW_DSGS.md`](../SW_DSGS.md).

With $\mathbf q = (H, Hu, Hv)$ and $\mathbf q_e = (H_e, 0, 0)$ the lake at
rest of `initialize.jl`, per element $K$ ($\Delta = \Delta_K/(k+1)$, the
element form; the nodal form of `:ldsgs_nodal` puts the same $\nu$ at every
node, see `DSGS.md` §4.7):

$$
\nu_K = \max\Big(C_{min}\Delta\lambda_K,\ \max\big(0, \min(C_{max}\Delta\lambda_K,\ C_R\Delta^2\mathcal R_K)\big)\Big),
\qquad
\lambda_K = \max_{i\in K}\big(|\mathbf v_i| + \sqrt{gH_i}\big),
$$

$$
\mathcal R_K = \max_{i\in K}\max_{q\in\{H,Hu,Hv\}}\frac{R^q_i}{n^q},\qquad
R^q_i = \Big|\frac{3q_i - 4q_i^{n-1} + q_i^{n-2}}{2\Delta t} - \frac{\mathrm{RHS}^q_i}{m_i}\Big|,\qquad
n^q = \max\big(\lVert\delta q - \langle\delta q\rangle\rVert_{\infty,\Omega},\ 10^{-3}s^q\big),
$$

with $\delta\mathbf q = \mathbf q - \mathbf q_e$ the **departure from the
lake at rest** (the depth itself is cone-shaped over the island, so its
spread would be the island and not the wave), $s^H = \bar H$,
$s^{Hu} = s^{Hv} = \bar H\sqrt{g\bar H}$ the still-water scales of the mean
depth $\bar H$, and the velocity of the wave speed desingularized as
$Hu/\max(H, h_{min})$ with $h_{min}$ the wet/dry threshold (`:dsgs_swe_hmin`
= `_H_WET_SWE` of `user_flux.jl`); $g$ is `:dsgs_swe_g` (= `_G_SWE`).
$C_R$ = `:dsgs_CR` (1), $C_{max}$ = `:dsgs_Cmax` (0.5), $C_{min}$ =
`:dsgs_Cmin` (0): the Dao–Nazarov names used by every DynSGS kernel of the
code.

The **same kinematic $\nu$ goes to the three slots** (times the deck's
`:μ` multipliers), applied on the primitives the sibling already diffuses
(`user_primitives.jl`): $\nabla\cdot(\nu\nabla(H - H_e))$ on the continuity
equation, which vanishes at the lake at rest so the discrete equilibrium of
the well-balanced split survives, and the stress form
$\nabla\cdot\big(\nu(\nabla(H\mathbf v) + \nabla(H\mathbf v)^\top - \tfrac23\nabla\cdot(H\mathbf v)\,I)\big)$
on the momenta.

Where the solution is resolved the residual part of $\nu$ is small; at the
wave front, at the run-up and on the wet/dry ring, where the residual spikes,
it rises towards the first-order cap $0.5\,\Delta(|\mathbf v| + \sqrt{gH})$
(0.14 measured at the run-up on the island, 0.035 in its lee at $t = 25$ s).

**The floor.** With $C_{min} = 0$ the CG solution of the solitary wave grows
an element-scale transverse mode ($Hv$ ripples along the whole front by
$t = 2$ s): a mode of that size has a residual the sensor does not see until
it is large. The run still completes ($t = 25$ s, the residual catches the
mode once it has grown, $\nu$ up to 0.21 on the front at $t = 2$ s, and the
late fields differ little from the floored run's), but the deck keeps a
background floor `:dsgs_Cmin => 0.05` so that the mode never appears,
i.e. $\nu \ge 0.05\,\Delta(|\mathbf v|+\sqrt{gH}) \approx 0.02$ m²/s in the
still basin, 40 % of the sibling's constant 0.05, and the run to $t = 25$ s is
clean. The difference from the sibling is then where the extra dissipation
goes: the residual part is zero at rest and concentrated on the island and the
front.

## Output

`H-it<n>.png`, `Hu-it<n>.png`, `Hv-it<n>.png` as in the sibling, plus
`μ_dsgs_H-it<n>.png`, the coefficient actually applied at that output time
(`:plot_dsgs_vars => ["H"]`; the three slots carry the same value). With
`:outformat => "vtk"` the field is `mu_dsgs` in the `.pvtu` files.

## Files

| file | role |
|---|---|
| `user_inputs.jl` | the deck: `DSGS_SW()`, its constants, the plotting of the coefficient |
| `initialize.jl`, `user_flux.jl`, `user_source.jl`, `user_bc.jl`, `user_primitives.jl` | identical to `SoliWaveIsland` |
| `SoliWaveIsland.geo`, `SoliWaveIsland.msh` | the sibling's mesh |
