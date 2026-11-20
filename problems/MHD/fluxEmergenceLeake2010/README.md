# 2.5D visco-resistive GLM-MHD: flux emergence into a field-free, quadrupolar or dipolar corona (Leake, Linton & Antiochos 2010)

This case runs the 2.5D flux-emergence tests of

> J. E. Leake, M. G. Linton, S. K. Antiochos,
> *Tests of Dynamical Flux Emergence as a Mechanism for Coronal Mass Ejection
> Initiation*, ApJ **722**:550 (2010). https://doi.org/10.1088/0004-637X/722/1/550
> (arXiv:1007.5484)

with the correction of its erratum,

> J. E. Leake, M. G. Linton, S. K. Antiochos, *Erratum: "Tests of Dynamical Flux
> Emergence as a Mechanism for Coronal Mass Ejection Initiation" (2010, ApJ,
> 722, 550)*, ApJ **741**:125 (2011). https://doi.org/10.1088/0004-637X/741/2/125

A buoyant, twisted flux tube 1.8 Mm below the photosphere rises through a
3 Mm deep convection zone and emerges through an isothermal photosphere and
chromosphere, a transition region and a 0.86 MK corona that is field-free, or
holds a quadrupole (the magnetic breakout configuration) or a dipole. The
erratum found that the original runs had an unintended local resistivity,
$\hat\eta = 0.1\exp(-\hat r^2/0.66^2)$ at the origin; without it the axis of the
low-tension tube stays at the surface and the lower coronal flux rope of the
paper does not form. **The default deck is the erratum's Figure 3** (which
replaces Fig. 10 of the paper): the low-tension tube emerging into the medium
("98 G") quadrupole, $\log\rho$ and field lines at $t$ = 920, 1380, 1610,
1955, 2185 and 2530 s.

The numerics are those of
[`fluxEmergenceSon2025DSGS`](../fluxEmergenceSon2025DSGS/README.md): continuous
spectral elements, GLM divergence cleaning, and the reference-weighted,
conserved-form Marras–Nazarov DynSGS operator as the only shock capturing, no
limiter. On top of it come the paper's explicit viscosity and anomalous
resistivity (the 2D `:mhd_dissipation` hook, as in
[`fluxEmergenceLeake2014`](../fluxEmergenceLeake2014/README.md) and
[`fluxEmergenceLeake2022_2D`](../fluxEmergenceLeake2022_2D/README.md), and a new
`:mhd_anomalous_eta`), its energy-relaxation term, an initialization of the
magnetic field that is discretely divergence-free, and a divergence-preserving
(curl) form of the DynSGS resistivity (`:dsgs_curl_resistivity`, new).

## Run

```bash
mpiexec -n 64 julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2010
```

or serially `julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2010`; with
the MPI bundled in Julia (MPICH_jll, no system `mpiexec`):

```bash
julia --project=. -e 'using MPI; run(`$(MPI.mpiexec()) -n 4 $(Base.julia_cmd()) --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2010`)'
```

Output (VTU) goes to `./output/MHD/fluxEmergenceLeake2010/<run>/` every
$5t_0 = 116$ s and at the times of the erratum's figures.
Between two snapshots (an hour of wall time on the fine meshes) the log is
silent but for a line every ten minutes,
`# progress: t = 43.1 t0 = 1000 s of 2530 s; 0.57 h so far, about 9.8 h to go at this rate`
(`JEXPRESSO_LK10_PROGRESS`), and the reports of the positivity floors, which
carry the time too; `conservation.dat` in the output directory gets a line
(`t` in $t_0$, step, mass, energy) every 100 steps.

All the simulations of the paper are selected from the environment, without
editing the deck:

| variable | values | meaning |
|---|---|---|
| `JEXPRESSO_LK10_TUBE` | `1`, `2` (default), `3`, `0` | tube 1 constant twist, 2 increasing twist ("low tension"), 3 decreasing twist; 0 no tube |
| `JEXPRESSO_LK10_FIELD` | `quadrupole` (default), `dipole`, `none` | coronal field |
| `JEXPRESSO_LK10_STRENGTH` | `strong`, `medium` (default), `weak` | quadrupole $q_d$ = 11.7, 5.85, 3.9 $\times10^7$ G cm ("197, 98, 66 G"); dipole $\lvert B(0,0)\rvert$ = 104, 52, 26 G |
| `JEXPRESSO_LK10_REF` | `1`, `2` (default), `3`, `4` | mesh (below); 4 has the core spacing of 2 over the whole figure region |
| `JEXPRESSO_LK10_TEND` | seconds, default 2530 | final time |
| `JEXPRESSO_LK10_DOUT` | seconds, default 116 | snapshot interval (the figure times are always added) |
| `JEXPRESSO_LK10_DT` | $t_0$ units | overrides the time step |
| `JEXPRESSO_LK10_MCZ` | default 1.3 | polytropic index $m$ of the convection zone, $p\propto T^{m+1}$ |
| `JEXPRESSO_LK10_RELAX_YMIN` | $L_0$ units, default $-\infty$ | apply the energy relaxation only above this height |
| `JEXPRESSO_LK10_TUBE_EQ` | `thermal` (default), `relaxed` | tube in thermal equilibrium with its surroundings (the paper), or on the relaxed relation $\rho T_0(\rho) = p$ |
| `JEXPRESSO_LK10_CURLB` | `1` (default), `0` | DynSGS resistivity of the in-plane field in curl form, or the component Laplacian of the parent case |
| `JEXPRESSO_LK10_NU` | `dynamic` (default), `kinematic` | the paper's dynamic viscosity in the convection zone, held at its photospheric kinematic value above; or that kinematic value everywhere |
| `JEXPRESSO_LK10_PARTITION` | `balanced` (default), `uniform` | the ranks hold equal numbers of elements, or equal areas (Meshes, time step and cost) |
| `JEXPRESSO_LK10_SUBNORMALS` | `flush` (default), `keep` | subnormal numbers flushed to zero, or IEEE gradual underflow (Meshes, time step and cost) |
| `JEXPRESSO_LK10_PROGRESS` | minutes, default 10 | interval of the log's progress line (simulated time reached, hours so far and to go); 0: none |
| `JEXPRESSO_LK10_RESTART_DIR` | an output directory, or `latest` | continue that run from a snapshot (Restart, below) |
| `JEXPRESSO_LK10_RESTART_TIME` | seconds | the snapshot nearest this time (default: the last complete one) |
| `JEXPRESSO_LK10_CMIN` | default 0 | DynSGS background floor $C_{min}$ |
| `JEXPRESSO_LK10_FLOORS` | `1` (default), `0` | positivity floors as the Runge–Kutta stage limiter (Positivity, below) |
| `JEXPRESSO_LK10_TFLOOR` | default 0.1 | temperature floor, in $T_{ph}$ = 5700 K |
| `JEXPRESSO_LK10_ETA` | `induction` (default), `all` | DynSGS resistivity from the induction residuals alone, or from the max over all residuals (one coefficient, the parent case) (DynSGS in the corona, below) |
| `JEXPRESSO_LK10_ENERGY` | `thermal` (default), `total` | energy-consistent DynSGS (Positivity, below), or $E$ diffused as a whole (the parent case) |
| `JEXPRESSO_LK10_STABCAP` | default 0.2; 0: none | the DynSGS coefficient never above this viscous number $\mu\Delta t/\Delta x^2$, which the explicit step still integrates (Meshes, time step and cost) |
| `JEXPRESSO_LK10_ETACAP` | `fast` (default), `flow` | DynSGS resistivity capped at $C_{max}\Delta(\lvert v\rvert + c_f)$ like the viscosity, or at $C_{max}\Delta\lvert v\rvert$ (Verification, the REF 4 run) |
| `JEXPRESSO_LK10_ETAFAC` | default 1 | factor on the DynSGS resistivity (the `:μ` multipliers of the B slots), a sensitivity test |
| `JEXPRESSO_LK10_RHOFAC` | default 1 | factor on the DynSGS mass diffusion (the `:μ` multiplier of the ρ slot), a sensitivity test |
| `JEXPRESSO_LK10_HEAT` | `shock` (default), a fraction | what the heat the DynSGS dissipates does (`:dsgs_heating`): `shock`, only the explicit ν, η and the DynSGS viscosity in compression heat, as in Lare (The heat of the DynSGS, below); or that fraction of it into the plasma (1 conserves energy, 0 loses it all) |
| `JEXPRESSO_LK10_SENSOR` | `legacy` (default), `residual` | DynSGS sensor: $\lvert\partial_tq\rvert$ against the assembled RHS (the parent case's), or the element residual once per step (DSGS.md §1.2) |

For example, the erratum's Fig. 1 (tubes 1 and 2, field-free corona) and its
Figs. 6–7 (tubes into the 52 G dipole):

```bash
JEXPRESSO_LK10_FIELD=none JEXPRESSO_LK10_TUBE=1 mpiexec -n 64 julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2010
JEXPRESSO_LK10_FIELD=dipole JEXPRESSO_LK10_TUBE=2 mpiexec -n 64 julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2010
```

## Units and the reference atmosphere

The paper is written in Gaussian units with dimensional values. The code units
are photospheric, as in the authors' Lare runs (Arber, Haynes & Leake 2007,
"normalised ... by photospheric values"; the erratum gives $L_0$ = 150 km,
$v_0$ = 6.5 km/s):

| | value |
|---|---|
| length $L_0$ | 150 km |
| temperature | $T_{ph}$ = 5700 K ($T_{ph}$ = 1) |
| density, pressure | $\rho_{ph} = 2.021\times10^{-7}$ g cm$^{-3}$, $p_{ph} = 8.45\times10^4$ dyn cm$^{-2}$ ($\rho_{ph} = p_{ph} = 1$) |
| velocity $v_0 = (RT_{ph})^{1/2}$ | 6.467 km/s |
| time $t_0 = L_0/v_0$ | 23.19 s (the "23 s" of the paper's Eq. 6) |
| field $B_0 = (4\pi p_{ph})^{1/2}$ | 1030.6 G; magnetic pressure $\tfrac12\lvert\mathbf B\rvert^2$ |
| gravity | $g = g_\odot L_0/v_0^2 = 0.98273$ ($g_\odot$ = 274 m s$^{-2}$) |
| gas law | $p = \rho T$, $\gamma = 5/3$, mean molecular weight $\mu$ = 1.125 (fitted, below) |

**Two constants are not printed in the paper**: the gas constant (equivalently
$\mu$) and the photospheric pressure. Both, and the convection-zone gradient
below, were recovered from the paper's own Fig. 3(a), whose curves are vector
data in the published PDF: a hydrostatic fit of its gas-pressure and
temperature curves gives $\mu = 1.125$ (an isothermal scale height of
152.6 km, the "0.25 times the scale height at the surface" of 0.0375 Mm),
$p(0) = 10^{4.927}$ dyn cm$^{-2}$, and a convection zone with $p\propto T^{2.3}$,
with an rms misfit of 0.006 dex.

**Atmosphere** (Eq. 9), hydrostatic, $T_{cor} = 8.6\times10^5$ K $= 150.9$,
$y_{tr}$ = 3.75 Mm $= 25$, $w_{tr}$ = 0.75 Mm $= 5$:

$$
T(y) = \begin{cases} 1 - \dfrac{g}{m+1}\,y & y\le 0\\[4pt]
1 + \dfrac{T_{cor}-1}{2}\left[\tanh\dfrac{y-y_{tr}}{w_{tr}} + 1\right] & y > 0\end{cases}
$$

The paper sets $m = 1/(\gamma-1)$ and calls the polytrope "marginally unstable
to convection"; its Fig. 3(a) has $p\propto T^{2.3}$, i.e. $m \approx 1.3$, a
slightly super-adiabatic (unstable) gradient. That is the default
(`JEXPRESSO_LK10_MCZ`; 1.5 is the adiabatic polytrope).
The density falls from 18.8 at $y = -3$ Mm to $5.8\times10^{-8}$ at the
transition region and $6.8\times10^{-10}$ at 90 Mm. Above about 12 Mm the
paper's Fig. 3 shows a constant gas pressure, which is not hydrostatic (with
gravity it would have to fall with the 23 Mm coronal scale height); the
hydrostatic corona is kept here.

**Flux tubes** (Eqs. 10–16), axis along $z$ at $y_t = -1.8$ Mm, $B_0$ = 7800 G,
$a$ = 0.3 Mm, $B_z = B_0e^{-r^2/a^2}$, $B_\theta = \Theta(r)\,rB_z$ with

$$
\Theta_1 = c,\qquad \Theta_2 = c\left(1 - e^{-r^2/c_1^2}\right),\qquad \Theta_3 = c\,e^{-r^2/c_2^2},
\qquad c = 1/a,\ c_1 = 0.15\ \text{Mm},\ c_2 = 1.5\ \text{Mm}.
$$

The in-plane field circulates so that $B_x > 0$ at the top of the tube:
parallel to the inner arcade of the quadrupole and antiparallel to the dipole
(the arrows of the paper's Fig. 5). Radial balance (Eq. 16) gives
$p_1 = -\tfrac12B^2 + \int_r^\infty B_\theta^2/r'\,dr'$ in closed form, and the
tube is in thermal equilibrium: $p = p_{bg} + p_1$, $\rho = p/T_{bg}$. At the
axis of tube 2, $\beta = 1.5$ and the density deficit is 33% (the paper's text
says "approximately 10", its Fig. 3(b) shows 1.5–1.9).

**Quadrupole** (Eq. 17): $A_z$ on the source line $y = -7.5$ Mm is
$q_d\,(x/L_0)^2[1 - (x/x_a)^2]$ for $\lvert x\rvert \le x_a$ = 21 Mm, and $A_z$
is harmonic above it with $A_z\to0$ at infinity (its Poisson integral, by
600-point Gauss–Legendre quadrature over the source). **$x$ enters in units of
$L_0$ = 150 km**: with that reading $\beta(0, y)$ of the three quadrupoles
matches the paper's Fig. 3(b) to within 1% at every height, and the null
sits at 3.39 Mm (Fig. 3(b): 3.375 Mm; Fig. 5: "about 3.2 Mm"). The labels
"197, 98, 66 G" are the **maximum** field on $y = 0$ (196, 97.8, 65.2 G at
$x$ = 16.5 Mm); at $x = x_a/2$ the field is 158, 79, 53 G.

**Dipole** (Eq. 18): $A_z = d\,(y - y_d)/r_1^3$, $y_d = -15$ Mm, with
$\lvert B(0,0)\rvert$ = 104, 52, 26 G, the values of Fig. 3(b) and of
Secs. 3.4–3.5 (the text's $d$ = 18, 9, 4.5 $\times10^{28}$ G cm$^3$ give 107, 53, 27 G).

## Equations and dissipation

Ideal GLM-MHD in conservation form (`user_flux.jl`), with gravity, the
paper's dissipation and its energy relaxation:

- **Viscosity** $\nu = 2.6\times10^4$ g cm$^{-1}$ s$^{-1}$ (Sec. 2.1) is a dynamic
  viscosity, 0.01326 in code units. In the convection zone it is applied as
  such: a kinematic viscosity $0.01326/\rho_e(y)$, with $\rho_e = T^m$ the density
  of the polytrope, 0.00126 at the tube's initial depth and 0.0007 at the
  bottom. Above the photosphere the kinematic value is held at its photospheric
  0.01326: taken as dynamic it would be a kinematic viscosity of $10^6$ in the
  corona (the 3D case, `fluxEmergenceLeake2014`, Deviations 1, uses the same
  authors' viscosity as kinematic everywhere; `JEXPRESSO_LK10_NU=kinematic` does
  that here). It enters through the `:mhd_dissipation` hook as a Laplacian of
  the momentum (the paper's stress tensor without its $\nabla\cdot\mathbf v$
  part), and the kinetic energy it removes stays in $E$ as heat. The choice
  hardly matters for the rise (field-free corona, tube 2, REF 1): the DynSGS
  coefficient at the tube, ~0.02, exceeds both forms, and the axis heights of
  the two runs differ by 0.002 Mm at 186 s.
- **Anomalous resistivity** (Eqs. 7–8), $\eta = \eta_0\max(0, \lvert j\rvert/j_{crit} - 1)$,
  $\eta_0 = 1.6\times10^{-11}$ s, $j_{crit} = 2.07\times10^7$ statA cm$^{-2}$: a
  magnetic diffusivity $\eta_0c^2/4\pi = 1.18\times10^{-4}$ and
  $\lvert\nabla\times\mathbf B\rvert_{crit} = 126$ in code units, evaluated per element
  from the largest nodal current (`:mhd_anomalous_eta`, `src/kernel/physics/SGS.jl`).
  Note that $j_{crit}$ is 30 times the current density of the initial tube,
  and in the corona it would need a field reversal of $10^5$ G over one grid
  cell: with the paper's constants the anomalous resistivity can only switch on
  in kilogauss sheets a few km thick. The run prints when it first does.
- **No other resistivity.** The erratum removed the local
  $\hat\eta = 0.1\exp(-\hat r^2/0.66^2)$, $\hat r < 10$, of the original runs.
- **Energy relaxation** (Eqs. 5–6): $D\epsilon/Dt$ gets
  $-(\epsilon - \epsilon_0(\rho))/\tau$, with $\epsilon_0(\rho)$ the internal energy
  of the initial atmosphere at the height where its density is $\rho$
  (tabulated against $\ln\rho$), and $\tau = 23\,\text{s}\,(\rho/\rho_{ph})^{-1.43}$
  with the local density, as in the authors' Lare runs (Arber, Haynes & Leake
  2007, their Eq. 16). It acts in the convection zone ($\tau\approx1$ s) and the
  photosphere (23 s) and is negligible above ~1 Mm. Applied everywhere it also
  cools the buoyant tube within its first second (the tube is in thermal
  equilibrium at a lower density than its surroundings, so $\epsilon_0(\rho)$ is
  that of a higher, cooler layer): the tube contracts and keeps a density
  deficit of ~21% instead of 33%. The paper does not say where the term acts;
  the rise of the tube decides it. Field-free corona, tube 2, REF 1, axis rise
  speed between 139 and 186 s: 1.55 km/s with the relaxation everywhere,
  2.53 km/s with it only above the photosphere (`JEXPRESSO_LK10_RELAX_YMIN=0`),
  against 1.64 km/s in the erratum's Fig. 2(a). It acts everywhere by default.
  (`JEXPRESSO_LK10_TUBE_EQ=relaxed` starts the tube in the contracted state,
  without the transient.)
- **DynSGS** (`DSGS_MHD`, conserved form on the departure from the reference
  state, weighted by $\rho_e$; legacy sensor; element norms; fast-speed floors
  and thermal-energy residual for the low-β corona; $C_R = 1$,
  $C_{max} = 0.5$); the in-plane field is diffused in curl form (next sections).
  It replaces the numerical dissipation of Lare (shock viscosity and remap) and
  is the only shock capturing, with the positivity floors (next item). The background floor
  $C_{min}$ of the parent case is off: $C_{min}\Delta c_f$ with $c_f\approx270$
  in the corona would be a diffusivity of ~4 $L_0v_0$, two orders of magnitude
  above Lare's own estimate of its numerical resistivity.
- **Positivity.** Where the field is strong and β low the pressure is the small
  remainder $p = (\gamma-1)(E - \tfrac12\rho v^2 - \tfrac12B^2 - \tfrac12\psi^2)$,
  and it goes negative: in a run of the default deck without floors from
  ~200 s in the upper chromosphere under the quadrupole's strongest field
  ($\lvert x\rvert$ = 10–19 Mm, $y \approx$ 2.3–2.6 Mm, $\beta \approx 3\times10^{-3}$;
  $p/p_0$ = −2.2 at 232 s, the CFL number 0.11 throughout), then 50 km/s flows
  in the corona and a non-finite solution at 631 s; in the field-free run at
  the top of the emerging tube ($T$ = −0.02 at 1000 s). Two measures:

  *Energy-consistent DynSGS* (`JEXPRESSO_LK10_ENERGY=thermal`, the default;
  `:dsgs_thermal_energy`, `src/kernel/physics/SGS.jl`). The parent operator
  diffuses $E$ as a whole; its pressure change carries
  $(\gamma-1)\mu(2\nabla\mathbf B_e:\nabla\delta\mathbf B + \delta\mathbf B\cdot\nabla^2\mathbf B_e)$,
  first order in the departure $\delta\mathbf B$ from the strong reference field and
  larger than $p$ where $\beta \lesssim 10^{-2}$ (the split energy flux of the parent
  case keeps it). Here the $E$ slot diffuses the internal-energy departure and
  receives, node by node, the work of the other slots' DynSGS terms,
  $\mathbf u\cdot R_{\rho\mathbf u} - \tfrac12\lvert\mathbf u\rvert^2R_\rho + \mathbf B\cdot R_{\mathbf B} + \psi R_\psi$,
  and the heat they dissipate (viscous, Ohmic, all $\ge 0$): the DynSGS changes $p$
  only by diffusing the internal energy and by heating. Test: a 32 km/s kick at
  (17, 2.4) Mm in the quadrupole's chromosphere ($\beta = 3.6\times10^{-3}$), REF 1,
  no floors: with $E$ diffused as a whole the plasma there cools to
  $T$ = 0.71 $T_{ph}$ in 2.3 s and $p/p_0$ falls to 0.28, energy-consistent $T$
  stays at or above 1.01 $T_{ph}$ and $p/p_0 \ge 0.66$, with the same flows
  (14.0 and 13.5 km/s at 11.6 s) and the same total energy (−2.2 and
  −2.3 × 10⁻¹⁰ over 1200 steps). The background stays at rest to the bit;
  the field-free rise to 186 s (REF 1) has the axis of the parent operator to 0.0001 Mm, and $\psi$ stays at round-off (7 × 10⁻¹⁶).

  *Positivity floors* (`JEXPRESSO_LK10_FLOORS`, on by default, as in the other
  Leake cases), the safety net. The stage limiter keeps
  $T = p/\rho \ge T_{min}$ = 0.1 $T_{ph}$ (570 K, below any temperature of the
  model; `JEXPRESSO_LK10_TFLOOR`) by raising $E$ at fixed velocity and field,
  and $\rho \ge 0.2\,p_e(y)/T_{cor}$: $0.2\rho_e$ in the corona, 150 times below
  the density in the chromosphere. Unlike the other cases' $0.2\rho_e(y)$,
  $0.2p_e(y)$, neither acts where the transition region dips (coronal density
  and pressure at chromospheric heights) or on cool plasma lifted into the
  corona. Not conservative where they fire; the log reports the hits.

- **DynSGS in the corona.** The DynSGS has one coefficient per element,
  from the largest normalized residual of all equations, used as viscosity,
  mass and heat diffusivity and resistivity alike (Dao & Nazarov). In the
  low-β corona the waves that reach it from below steepen, and their density
  and thermal residuals set a coefficient of 0.5–3 $L_0v_0$ (REF 2, cluster run),
  0.2–1 % of the first-order value $C_{max}\Delta c_f$ but, as a resistivity, 5–30
  times the $\hat\eta = 0.1$ whose accidental presence at the origin changed the
  paper's results (erratum), and enough to diffuse a 1 Mm flux rope away in a
  few hundred seconds. Lare has no resistivity there: its shock viscosity acts
  on compressions only, and its remap diffuses with the flow speed, not the fast
  speed. Two measures (the first is the default since the cluster run of
  Verification, the second is the mesh of that run): `JEXPRESSO_LK10_ETA=induction`
  (`:dsgs_induction_sensor`) takes the coefficient of the $\mathbf B$ slots from
  the residuals of the induction equation alone, so that an acoustic front,
  which hardly changes $\mathbf B$ at low β, no longer makes the field resistive
  (an Alfvénic or fast-mode structure has $\delta B/B = \delta v/c_f$ and keeps
  $\eta = \nu$; a grid-scale current sheet still gets the cap); and REF 4
  (Meshes), which halves $\Delta$ where the dome and the rope evolve
  (the residual coefficient is $C_R\Delta^2\times$ the residual). The kick test
  of "Positivity" (32 km/s at (17, 2.4) Mm, $\beta = 3.6\times10^{-3}$, REF 1):
  in the disturbed corona (11 < x < 23 Mm, 1.5 < y < 12 Mm) the resistivity is
  0.57 (median) and 1.5 (largest) $L_0v_0$ at 11.6 s with one coefficient, 0.018
  and 0.18 from the induction residuals alone, with the same flows (14.0 and
  13.9 km/s), currents, $\psi$ (round-off) and total energy (−2.2 × 10⁻¹⁰ over
  1200 steps). With the induction sensor the background stays at rest to the
  bit, and the field-free rise to 186 s (REF 1) has the axis of the
  one-coefficient run to 0.004 Mm, with $\psi$ at round-off, no floor hit, and
  the tube's axial field better kept (6216 against 5985 G at 186 s; around the
  tube the resistivity is half the viscosity). The `atmos` figure shows both
  coefficients against height.

- **The heat of the DynSGS.** With the energy-consistent operator every
  erg the DynSGS removes from the flow and the field can become heat, so that
  the total energy is conserved (`JEXPRESSO_LK10_HEAT=1`). Lare does not do that: its remap does not conserve the
  kinetic energy and the loss is not returned as heat unless asked
  (`rke = .FALSE.` in the distributed `control.f90`: "This missing energy can be
  added back into the simulation as a uniform heating"); what heats there is
  the shock viscosity, which acts in compression, and the physical terms of the
  paper's Eq. (5). The difference matters for the cold, dense plasma that the
  emerging field lifts into the corona, which has little internal energy: in
  the REF 4 run (Verification) the DynSGS viscosity could supply all of it in
  100–600 s, the plasma went from below 2 × 10⁴ K to 10⁵ K, stopped draining,
  and the dome stopped rising. `JEXPRESSO_LK10_HEAT=shock` (`:dsgs_heating =>
  "shock"`, `src/kernel/physics/SGS.jl`), the default, heats as Lare does: the explicit ν and
  η in full, the DynSGS viscosity only where the flow is compressed,
  $\rho\mu(\nabla\cdot\mathbf v)_-^2$ and never more than it dissipates there, the
  DynSGS resistivity and ψ diffusion not at all; the rest of what the DynSGS
  dissipates is lost, as in the remap. `JEXPRESSO_LK10_HEAT=0` loses all of it
  (shocks then do not heat either: a test only). `shock` is the default since
  the cluster run of Verification, in which the lifted plasma keeps draining
  and the flux rope of the erratum forms. Here, the kick test of "Positivity" (REF 1, no floors,
  11.6 s) runs the same with the three: flows of 13.9, 13.8 and 13.8 km/s
  (1, `shock`, 0), the same smallest temperature, the pressure near the kick
  1.0 % (`shock`) and 1.3 % (0) below that with all the heat, the total energy change
  −2.98 and −3.22 against −2.21 × 10⁻¹⁰: of the heat the DynSGS dissipates
  there `shock` keeps a quarter, the compressed half of a kick whose velocity
  gradient is half compression. A restart that changes the mode reports it
  and continues.

## What is not in the paper, and what was done

1. **Boundary conditions.** Not stated. As in the authors' 3D Lare runs (Leake,
   Linton & Török 2013): closed walls, all velocity components zero, $\psi = 0$,
   $\rho$, $E$, $\mathbf B$ free (zero normal gradient). Velocity damping
   layers, rising linearly to $1/t_0$ at the boundary, absorb waves for
   $\lvert x\rvert > 37.5$ Mm and $y > 75$ Mm, outside the region of the figures.
2. **Grid.** "Stretched in both x and y", 0.0375 Mm at the surface to 0.15 Mm at
   the top, map not given. `LK10.geo` grades the element size with a two-level
   tanh map (below).
3. **Divergence of B** (next section): the initial in-plane field is a discrete
   curl, and the DynSGS resistivity acts on it in curl form.
4. **The discrete equilibrium** (`:lwell_balanced`): the element RHS of the
   reference state (atmosphere + potential coronal field) is subtracted at every
   stage, so the background is an exact discrete steady state.

## Divergence of B

Lare's staggered field is the discrete curl of $A_z$ and keeps
$\nabla\cdot\mathbf B = 0$ to round-off. Jexpresso's collocated CG field does not,
and GLM carries the error away at $c_h$ (275 $v_0$ here) — into a corona where
$\beta \sim 0.01$–0.1 and a field error $\delta B$ is a Lorentz force
$B\,\delta B/\rho$ on a density of $10^{-8}$. Two sources were found and removed:

- **The initial field.** The nodal interpolant of the analytic tube field has a
  discrete divergence (the strong-form, lumped-mass DSS divergence of
  `rhs.jl`) of 0.2–0.7 $B_0/L_0$ on 2.5 $L_0$ elements. Here the in-plane field
  is instead the lumped-mass average of the element curls of the nodal $A_z$
  (`lk10_discrete_curl` in `initialize.jl`, summed over the ranks through the
  global node numbering), whose discrete divergence is zero to round-off on any
  quadrilateral mesh. The ideal induction term keeps it so (its tendency is the
  same kind of discrete curl, of $E_z = vB_x - uB_y$).
- **The DynSGS resistivity.** The component Laplacian $\nabla\cdot(\eta\nabla\mathbf B)$
  with an element-wise $\eta$ does not keep the discrete divergence: its
  tendency has a divergence as large as itself (measured on stretched
  elements). With `:dsgs_curl_resistivity` (new, default here) the in-plane
  field is diffused as $\partial_t\mathbf B = -\nabla_h\times(\eta\,j_z\hat z)$ in
  strong form, with the nodal $\eta j_z$ summed over the ranks
  (`_dsgs_curl_resistive_2d!`, `src/kernel/physics/SGS.jl`): the same
  operator for a solenoidal field, divergence-free to round-off.

Measured on REF 1 with the quadrupole and tube 2 (relaxation everywhere, so
that the tube contracts violently in its first second, see the relaxation item
above), at $t$ = 1.5 s:

| | max $\lvert\psi\rvert$ | max $\lvert\mathbf v\rvert$ around the null | max $\lvert\mathbf v\rvert$ in the corona ($y$ > 5.25 Mm) |
|---|---|---|---|
| interpolated field, Laplacian | $1.1\times10^{-3}$ | $4.7\times10^{-2}$ | $8.2\times10^{-3}$ |
| discrete-curl field, Laplacian | $1.6\times10^{-3}$ | $5.9\times10^{-3}$ | $9.3\times10^{-4}$ |
| discrete-curl field, curl form (default) | $9\times10^{-16}$ | $1.4\times10^{-14}$ | $1.4\times10^{-15}$ |
| discrete-curl field, no DynSGS on $\mathbf B$ | $1\times10^{-15}$ | $1.4\times10^{-14}$ | $1.4\times10^{-15}$ |

With the Laplacian the disturbance around the quadrupole null kept growing
(0.04 $v_0$ at 3 s), forced by the ψ that the DynSGS made in the contracting
tube; with the curl form the corona is at rest to round-off while the tube
evolves as before. `JEXPRESSO_LK10_CURLB=0` restores the Laplacian of the
parent case.

**The walls.** The resistive $E$ is set to zero at the wall nodes, a
conducting wall for the resistive field, as the ideal $E = -\mathbf u\times\mathbf B$
already is at the no-slip walls. With the one-sided wall values of $E$ kept,
the operator has growing modes, the fastest at the corners of the domain,
where its Poynting flux through the walls feeds them: on 5×5 to 6×6 elements
with walls (uniform or stretched, $\eta$ uniform or varying up to 100× between
elements) the largest growth rate was up to $+2\,\eta/h^2$ with the wall values
kept, and below $10^{-11}$ (round-off) with $E = 0$ on the walls; the weak
Laplacian has none. A run of the default deck blew up at 506.8 s from the top
right corner (41.75, 86.34) Mm, where the largest DynSGS coefficient had grown
from 7 to 113 in 36 s, every other diagnostic steady. The update is still the
curl of a nodal scalar, so $\nabla\cdot\mathbf B$ stays at round-off (quadrupole
and tube, REF 1: $\max\lvert\psi\rvert = 10^{-15}$ at 5.8 s); the log prints the
number of wall nodes. Restarted from 464 s with $E = 0$ on the walls, the same
run went to 600 s with the largest coefficient at 6.3.

## Meshes, time step and cost

`LK10.geo`: $x\in[-45, 45]$ Mm, $y\in[-3, 90]$ Mm; core element size $h_0$
over $\lvert x\rvert < 5$ Mm and $y < 4.5$ Mm (tube, photosphere, chromosphere,
transition region), $2h_0$ over the region of the figures (to 12.75 Mm and
19.5 Mm), up to ~$9h_0$ at the far boundaries.

| `REF` | elements | nodes (`:nop => 4`) | core nodal spacing | file |
|---|---|---|---|---|
| 1 | 90 × 78 | 113 k | 0.62 $L_0$ (94 km) | `LK10_ref1.msh` (smoke tests) |
| 2 | 180 × 156 | 451 k | 0.31 $L_0$ (47 km) | `LK10_ref2.msh` (default; the paper has 37.5 km) |
| 3 | 362 × 312 | 1.8 M | 0.16 $L_0$ | `gmsh -2 LK10.geo -setnumber ref 3 -format msh41 -o LK10_ref3.msh` |
| 4 | 200 × 172 | 552 k | 0.31 $L_0$ (47 km) | `LK10_ref4.msh`: $h_0$ of REF 2 over the whole region of the figures |

REF 2 has its core spacing only where the tube rises and crosses the surface;
where the dome and the coronal flux rope of the erratum's Fig. 3 evolve
($\lvert x\rvert$ to 12.75 Mm, $y$ to 19.5 Mm) its nodes are 0.63 $L_0$ (94 km) apart,
against the paper's 37.5 km at the surface to 150 km at the top. REF 4 keeps
0.31–0.33 $L_0$ (49 km) over that whole region, then grows in one level to ~10 $h_0$
at the far boundaries: 1.23 times the elements of REF 2 at the same time step.

The time step is set by the fastest wave of the initial state, the fast mode
just above the transition region where the quadrupole is strongest
($c_f$ = 275 $v_0$ = 1780 km/s at $x = \pm17$ Mm, $y = 4.7$ Mm, for the medium
quadrupole; 549 and 183 for the strong and weak ones, 258, 130, 66 for the
dipoles, 16 for the field-free corona), held as the GLM speed $c_h$:
$\Delta t = 0.25\,\Delta h/c_h$ on the smallest LGL gap $\Delta h$, i.e.
$3.9\times10^{-4}t_0$ (REF 1) and $2.0\times10^{-4}t_0$ (REF 2) for the
default, $2.8\times10^5$ and $5.6\times10^5$ steps to 2530 s.

**Cost** (measured here, 4 MPI ranks, CPU): 0.5 s per step on REF 1 and about
3 s per step on REF 2, i.e. ~160 core-hours for the default problem on REF 1
and ~1900 core-hours on REF 2 (a day on 64-128 cores, before parallel losses).
The field-free corona allows a ~15 times larger step (`JEXPRESSO_LK10_DT`:
0.006 $t_0$ on REF 1 ran to 1000 s with an acoustic CFL of 0.2, so 0.003 on
REF 2; 1000 s took 80 minutes on REF 1 on 4 cores). A
snapshot takes 84 MB on REF 2 (21 MB on REF 1): about 2.5 GB for the 30
snapshots of the default output.

**Rank partition.** The framework's default (`:lxy_partition`) bins the
elements into a uniform $n_x\times n_y$ grid of the bounding box, equal areas
for every rank. These meshes are graded, so the ranks over the refined region
hold most of the elements and every step waits for them:

| ranks | REF 2, largest rank / ideal | REF 4, largest rank / ideal |
|---|---|---|
| 4 (2 × 2) | 11340 / 7020 = 1.6 | 14700 / 8600 = 1.7 |
| 32 (4 × 8) | 3710 / 878 = 4.2 | 5394 / 1075 = 5.0 |
| 64 (8 × 8) | 2385 / 439 = 5.4 | 3658 / 538 = 6.8 |
| 128 (8 × 16) | 1395 / 219 = 6.4 | 1829 / 269 = 6.8 |

The deck therefore sets `:lxy_partition_balanced` (new; `mesh.jl`,
`_compute_xy_partition`): the same bins with equal *numbers* of elements
($x$-rank columns, $y$-rank blocks within each), 268–269 elements per rank on
REF 4 with 128 ranks. Measured on REF 4, 4 ranks: 1.53 s per step against
2.72 s. The local meshes hold no ghost elements, and the solution is that of
the uniform bins to round-off (REF 1, 333 steps: $\rho$ to $2\times10^{-16}$,
$\mathbf B$ to $6\times10^{-16}$, the same mass and energy to every digit).
`JEXPRESSO_LK10_PARTITION=uniform` restores the old bins, which a restart of a
run made with them needs. The framework's other balanced partition
(`:lxy_partition => false`, METIS) is not used: its local meshes carry a layer
of the neighbours' elements, which the assembly then weights twice one layer
inside each rank, and with it $\psi$ reached $3\times10^{-5}$ in 300 steps where
it stays at round-off here (the discrete curl of `initialize.jl` now counts
every element once whatever the partition).

**The time step and the DynSGS cap.** $\Delta t$ is fixed, $c_h\Delta t/\Delta x_{min} = 0.25$
with $\Delta x_{min}$ the smallest LGL gap and $c_h$ the largest fast speed of
the initial state (274.7 for the medium quadrupole). The fast speed does not
stay there: the cluster's log has a fast CFL number of 0.25 at the start, 0.34
at 928 s and 0.64 at 1624 s (702 $v_0$, 4500 km/s, in the draining dome, where
the density falls and the field is compressed). For `CarpenterKennedy2N54` on
degree-4 elements the linear limits are, in the numbers the log prints
(scratch computation of the spectra of the 1D operators, 2D by tensor product):
a fast CFL number of 1.64 for the advection alone and a viscous number
$\mu\Delta t/\Delta x_{min}^2$ of 0.43 for the diffusion alone. The DynSGS
coefficient at its first-order cap, $C_{max}\Delta(\lvert v\rvert + c_f)$ with
$C_{max} = 0.5$ and $\Delta = \Delta_{elem}/5$, has a viscous number of
$0.58\times$ the fast CFL number, so that with it the scheme is unstable from a
fast CFL number of **0.62**, not 1.64: a run whose fast speed has grown 2.5-fold
blows up at the first strong front in that region (the `JEXPRESSO_LK10_HEAT=shock`
restart did, at 1713 s, in the top of the dome, with a largest coefficient of
107 = a viscous number of 0.46; Verification). `:dsgs_stab_cap`
(`JEXPRESSO_LK10_STABCAP`, default 0.2; `src/kernel/physics/SGS.jl`) therefore
holds the coefficient at or below $0.2\,\Delta x_{min,e}^2/\Delta t$, with
$\Delta x_{min,e}$ the smallest LGL gap of the element: nothing changes where the
coefficient is below it (everywhere at the design CFL number, where the
first-order cap is 0.145), and the scheme stays linearly stable to a fast CFL
number of 1.4. The log names the rank the first time the cap acts. If the fast
CFL number of the log approaches 1, restart with a smaller `JEXPRESSO_LK10_DT`.
Here, the kick test at five times the time step (fast CFL number 1.13) with the
coefficient driven to its cap ($C_R$ = 1000): non-finite at 0.39 $t_0$ without
the stability cap, with a viscous number of 0.71 at the snapshot before; with
it, 0.26 and a clean run to 2 $t_0$ (flows decaying from 32 to 1.2 km/s,
smallest pressure, density and temperature unchanged, $\psi$ at round-off).

**Subnormal numbers** are flushed to zero on every rank
(`set_zero_subnormals(true)` in `initialize.jl`; `JEXPRESSO_LK10_SUBNORMALS=keep`
for IEEE underflow). $\rho\mathbf u$, $\rho w$, $B_z$ and $\psi$ are zero in the
reference state, and what reaches the quiet parts of the domain first is their
numerical precursor, values of $10^{-300}$ and less, on which the processor
takes a path about a hundred times slower. After 121 steps on REF 4 (4 ranks,
induction sensor) a third of the nodes held $\lvert B_z\rvert$ and $\lvert w\rvert$ below
$10^{-200}$ and 2–4 % of them subnormal values; the step took 3.57 s on average
over steps 5–100 against 2.68 s with the flush, and the two solutions differ
by $10^{-291}$ at most ($\rho$, $p$, $B_x$, $B_y$ are equal to the bit). A
128-rank cluster run on REF 4 with the induction sensor, which leaves the far
$B_z$ undiffused, advanced 8 times more slowly than expected before the flush
(5 $t_0$ in 3.5 hours).

## Restart

For a queue with a time limit (the default problem takes about a day on 128
cores), a run continues from its snapshots:

```bash
# job 1
mpiexec -n 128 julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2010
# jobs 2, 3, ...: the most recently written run of the same configuration, from its last complete snapshot
JEXPRESSO_LK10_RESTART_DIR=latest mpiexec -n 128 julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2010
```

With SLURM, for example, `sbatch job1.sh`, then
`sbatch --dependency=afterany:<jobid> job2.sh` (`afterany`: also when job 1
was stopped by the limit); a restart from a run that has reached the end time
stops at once ("nothing to continue"), so a spare job in the chain is harmless.

- **Same switches, mesh, rank count and partition.** Each rank reads its own
  piece of the snapshot (a run made before the balanced partition existed
  continues with `JEXPRESSO_LK10_PARTITION=uniform`). The numerical switches
  may change, and the log says which did: to continue a run made before a
  default changed with its own numerics, set them (`JEXPRESSO_LK10_ETA=all`,
  `JEXPRESSO_LK10_HEAT=1` for runs older than those defaults; its
  `lk10_case.txt` lists what it used). Every run writes its `JEXPRESSO_LK10_*` configuration, mesh and
  rank count to `lk10_case.txt` in its output directory, and a restart from a
  run with another record stops with the list of differences. The time step
  (`JEXPRESSO_LK10_DT`), the end time, the snapshot interval and the numerical
  switches `JEXPRESSO_LK10_FLOORS`, `_CMIN` and `_CURLB` may change (the log
  says which did): that is how a run continues past a failure (below).
- **Which run.** `JEXPRESSO_LK10_RESTART_DIR` is the output directory
  (`output/MHD/fluxEmergenceLeake2010/output-<date>`), or `latest`: the most
  recently written run with the same record under
  `output/MHD/fluxEmergenceLeake2010/`. Give the directory when two runs of
  the same configuration are under way.
- **Which snapshot.** The last complete one (a job stopped by its time limit
  while writing leaves an incomplete one, which is skipped), or the one nearest
  `JEXPRESSO_LK10_RESTART_TIME` (seconds).
- **Output.** The new run writes to a new output directory, starting with the
  snapshot it continues. The figure script reads a restarted run as one from
  its last directory (Figures).
- **What is restored.** The conserved state, rebuilt to round-off from the nine
  primitives of the snapshot (Float64), the discrete $\nabla\cdot\mathbf B$
  included; the reference state, the relaxation table, $c_h$ and $\Delta h$ are
  those of $t = 0$, as in the run that wrote it. The DynSGS residual history
  starts afresh, so the continuation is close to, not bitwise equal to, an
  uninterrupted run (Verification).
- **What a job stopped by the limit loses:** the time since its last snapshot,
  116 s of simulated time by default, ~25,000 steps on REF 2 (an hour or less
  on 128 cores); `JEXPRESSO_LK10_DOUT` shortens it (the figure times are always
  written).

**When a run fails** (`non-finite solution at t = …`, with the place): the
snapshots before it show what grew, and from where (code units; `--dt` is the
run's $\Delta t$ in $t_0$, printed in its log, and turns the rate into the CFL
number with the true fast speed; `--ymin 2` looks at the corona only):

```bash
python3 tools/plot_leake2010.py diag $RUN --dt <the run's Δt>
```

Then continue from a snapshot before the trouble began, one change at a time
(same rank count), e.g. from 464 s with the time step halved, or another
temperature floor (`JEXPRESSO_LK10_TFLOOR`), or the component Laplacian on $B$
(`JEXPRESSO_LK10_CURLB=0`); a run made without the floors gets them by default:

```bash
JEXPRESSO_LK10_RESTART_DIR=$RUN JEXPRESSO_LK10_RESTART_TIME=464 JEXPRESSO_LK10_DT=<half the run's Δt> mpiexec -n <ranks> julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2010
```

A first sign comes before the failure: `anomalous resistivity … active` means
$\lvert\nabla\times\mathbf B\rvert$ passed $j_{crit}$ = 126, a jump of tens of
kilogauss across a cell, which nothing in this problem makes physically.

## Figures

`tools/plot_leake2010.py` draws the paper's figures from the output (numpy and
matplotlib). `RUN` is `output/MHD/fluxEmergenceLeake2010/output-<date>`.

```bash
python3 tools/plot_leake2010.py density  $RUN            # Erratum Fig. 3: log rho, field lines, B_x arrows, 920-2530 s
python3 tools/plot_leake2010.py shear    $RUN            # Erratum Fig. 1: B_z (G) and field lines, 920-2254 s
python3 tools/plot_leake2010.py dsgs     $RUN            # the DynSGS coefficient of the momentum equations at the times of Fig. 3
                                                         #   (--which induction | energy: of those equations)
python3 tools/plot_leake2010.py movie    $RUN            # one image per snapshot and a GIF of them
                                                         #   (--field density | momentum | induction | energy)
python3 tools/plot_leake2010.py dome     tube2=$RUN      # top of the dense dome on x = 0 vs time, with the erratum's Fig. 3 values
python3 tools/plot_leake2010.py axis     tube2=$RUN      # Erratum Fig. 2(a): axis height vs time, with the erratum's curves
python3 tools/plot_leake2010.py eshear   tube2=$RUN      # Erratum Figs. 2(b,c), 5: E_shear above 1.2 and 5 Mm (Eq. 20)
python3 tools/plot_leake2010.py profiles $RUN            # Fig. 3(a,b) at t = 0, with the paper's values
python3 tools/plot_leake2010.py diag     $RUN --dt <Δt>  # per snapshot: largest speed, CFL, current; smallest p/p0, beta
python3 tools/plot_leake2010.py null     base=$RUN       # on x = 0 vs time: the null above the dome, the flux below it, the mass lifted
python3 tools/plot_leake2010.py mass     $RUN            # the plasma above 3 Mm: mass budget, temperature, pressure, draining, heating time
```

Figures are written to `figs_leake2010/` (`--outdir`) as
`flux_emergence_Leake2010_<variable>.png` (`--prefix`, `--format`):
`…_density`, `…_Bz`, `…_mu_dsgs_momentum`, `…_mu_dsgs_induction`,
`…_mu_dsgs_mass_energy`, `…_dome`, `…_axis`, `…_null`, `…_eshear`, `…_atmos`,
`…_profiles`. Every map has a title and a labelled colour bar naming what its
colours show.

**Which DynSGS coefficient.** The solver has one coefficient $\mu$ per element,
from the largest normalized residual of the eight equations ($\rho$, $\rho u$,
$\rho v$, $E$, $\rho w$, $B_x$, $B_y$, $B_z$; and of the thermal energy; $\psi$
is left out), capped (Meshes, time step and cost). What is *applied* differs by equation, and the output has one field per
group of equations with the same coefficient, named after them:

| output field | equations | what it is | `dsgs --which`, `movie --field` |
|---|---|---|---|
| `mu_dsgs_ρu_ρv_ρw` | momentum | $\mu$ plus the paper's explicit viscosity (0.013 above the photosphere) | `momentum` (default; `mu`) |
| `mu_dsgs_Bx_By_Bz` | induction | its own coefficient $\eta$, from the induction residuals alone (`JEXPRESSO_LK10_ETA=induction`), plus the anomalous resistivity where it is active | `induction` (`eta`) |
| `mu_dsgs_ρ_ρE_ψ` | mass, energy, $\psi$ | $\mu$ itself | `energy` (`rho`) |

(With `JEXPRESSO_LK10_ETA=all` the induction equations share $\mu$ and the last
two are one field, `mu_dsgs_ρ_ρE_Bx_By_Bz_ψ`.) All are kinematic, in code units
$L_0v_0$ = 970 km² s⁻¹, the element's value at its nodes (the largest of the
elements sharing a node). `dsgs` and `movie` draw log₁₀ of the chosen one
(`--linear`: the coefficient itself) on the colour scale of the group's DynSGS
figures, dark navy where there is none, through cyan and cream to dark red.
`movie` writes
`frames_<variable>/flux_emergence_Leake2010_<variable>_NNNN_t<seconds>s.png`
for every snapshot (`--tmin`, `--tmax`, `--every`), all with the same window,
colour limits and field-line levels (those of the run's first snapshot, so
that a line is the same flux surface in every frame), and
`flux_emergence_Leake2010_<variable>.gif` (`--fps`, default 4; a frame lasts in
proportion to the time to the next snapshot, which are not equally spaced;
needs Pillow, which matplotlib installs). `--lw` sets the
width of the field lines (0.9 pt); a limit that begins with a minus sign is
written with `=`, `--ylim=-3,19.5`.

The density panels use the paper's window ($\lvert x\rvert < 12.75$ Mm,
$-3 < y < 19.5$ Mm), colour scale ($\log_{10}\rho$ from −16 to −5.3 g cm$^{-3}$, jet),
field lines as equally spaced contours of $A_z$ (integrated from $B_x$, $B_y$
element by element, exactly for the degree-4 interpolant) and arrows along
$B_x$ at $x = 0$, $y$ = 3, 11, 16 Mm. The axis of the tube is where $B_x$ on
$x = 0$ changes sign from − to + (the paper's definition), the crossing nearest
the peak of $B_z$. `axis`, `dome` and `eshear` write one CSV per run.

`RUN` may also be a list, `DIR1,DIR2`, and `DIR@T` uses a directory only up to
`T` seconds: `BASE@920,NEW` is `BASE` to 920 s and `NEW` from there, for a
`NEW` whose records do not lead back to `BASE`.
`RUN` is the last directory of a run continued by restarts: the script reads
the runs it was restarted from through the `lk10_case.txt` files, each up to
the snapshot its successor continued (`--no-chain`: only the directories
given), and lists each directory with its times. A run restarted from the
middle of a finished one, to try other numerics, is thus not mixed with what
the finished one did afterwards, and the two compare as
`dome old=<finished> new=<restarted>`. A figure time with
no snapshot (the run stopped before it) is left out, with a note (`--nearest`:
the nearest snapshot instead).

## Verification here

Short runs on 4 cores; the full runs are for a cluster (above).

- **Initial state against the paper.** At $x = 0$, $t = 0$
  (`plot_leake2010.py profiles`), $T$ and $p$ lie on the curves of the paper's
  Fig. 3(a) and $\beta$ on its Fig. 3(b) (98 G quadrupole). The quadrupole field
  on $x = 0$ is the paper's to within 1% at every height for all three
  strengths, the null at 3.39 Mm against 3.375 Mm, the surface maxima
  195.6, 97.8, 65.2 G against "197, 98, 66 G"; the 52 G dipole within 0.3%.
  The density of the erratum's Fig. 3 away from the tube (read off its colour
  scale at $x = -10$ Mm, 920 s) is this atmosphere's to 0.1 dex between 0.5 and
  3 Mm, and $\log_{10}\rho$ = −14.38 in its corona against −14.36 here.
- **Background at rest** (`JEXPRESSO_LK10_TUBE=0`, quadrupole, REF 1, 4
  ranks): after 220 steps every field equals the initial one to the last bit;
  mass and energy unchanged to all printed digits.
- **Divergence of B**: the table above; with the defaults $\psi$ stays at
  round-off ($10^{-15}$) through the tube's initial contraction, in every
  configuration run (quadrupole, strong quadrupole with tube 3, 52 G dipole with
  tube 1).
- **Rise through the convection zone** against the erratum's Fig. 2(a) (field-free
  corona, tube 2; the paper: "the initial rise of the flux tube is similar to
  the field-free corona case up to about 1000 s" for the quadrupole too):
  REF 1, 4 ranks, $\Delta t = 0.006\,t_0$, the dynamic viscosity (default). Height
  of the axis in Mm:

  | t (s) | 0 | 93 | 186 | 278 | 371 | 417 | 510 | 603 | 696 | 789 | 881 | 1000 |
  |---|---|---|---|---|---|---|---|---|---|---|---|---|
  | REF 1 | −1.79 | −1.70 | −1.57 | −1.41 | −1.22 | −1.14 | −0.95 | −0.77 | −0.60 | −0.45 | −0.34 | −0.22 |
  | erratum | −1.84 | −1.71 | −1.58 | −1.40 | −1.21 | −1.12 | −0.87 | −0.63 | −0.41 | −0.19 | −0.06 | 0.12 |

  To 420 s the axis follows the erratum's curve to 0.02 Mm. (At $t = 0$ it is at
  −1.785 Mm: the zero of $B_x$ of the discrete curl on REF 1 lies 0.1 $L_0$ above
  the analytic axis, −1.8 Mm; the erratum's curve starts at −1.84 Mm.) Then it
  falls behind: 1.9 km/s between 420 and 650 s against the erratum's 2.6, and it
  slows before the photosphere, which the erratum's axis crosses at 920 s; at
  1000 s it is where the erratum's axis was at ~780 s. The kinematic viscosity
  (`JEXPRESSO_LK10_NU=kinematic`) lowers it by less than 0.01 Mm to 417 s. REF 1
  has three nodes across the tube's radius (94 km against the paper's 37.5 km),
  which is the likely cause, untested here: **REF 1 is for smoke tests; the
  comparison with the paper's figures needs REF 2.** The same run on REF 2 is
  the cheapest check of that (field-free, ~50 core-hours, about 2 hours on 32
  cores), before the day-long Fig. 3 run:

  ```bash
  JEXPRESSO_LK10_FIELD=none JEXPRESSO_LK10_TEND=1000 JEXPRESSO_LK10_DOUT=46.389 JEXPRESSO_LK10_DT=0.003 \
      mpiexec -n 32 julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2010
  python3 tools/plot_leake2010.py axis ref2=output/MHD/fluxEmergenceLeake2010/output-<date> --tmax 1000
  ```
- **Restart** (the run above): continued from its 417.5 s snapshot, the state
  is the snapshot's to round-off (4e-16; only the DynSGS coefficient, a
  diagnostic of the residual history, differs). After 46 and 93 s the axis
  agrees with the uninterrupted run to 0.2 km and 2 mm; $\rho$ and $p$ agree to
  2e-5 and $\mathbf B$ to 0.4 % (RMS), the velocity to 1–3 % (RMS, mostly in the
  transition region and the corona), the footprint of the DynSGS coefficient of
  the first steps, whose residual history starts afresh. `latest` then picked the
  restarted run (not the empty directory of a failed one) and continued it to
  557 s, 0.6 km from the uninterrupted axis after two restarts; from a run at its
  end time it stopped with "nothing to continue"; a snapshot cut short (half of
  one piece) was skipped for the one before. `plot_leake2010.py axis` reads the
  three directories as one run.
- **A failure, and the positivity floors.** A cluster run of the default deck,
  then without floors, stopped at $t$ = 631 s ($27.23\,t_0$) with a non-finite
  density in the corona above the tube; the anomalous resistivity had switched on
  just after 580 s. `diag` on its snapshots: the CFL number stayed at 0.11, the
  pressure went negative from ~200 s in the low-β upper chromosphere
  (Positivity, above). The cause is the DynSGS diffusion of $E$ as a whole
  under a strong reference field; the energy-consistent DynSGS removes it and is
  the default since, with the floors as a safety net. With it, the next run had
  no floor hit for 506 s, then blew up from the top right corner of the domain:
  the curl-form resistivity fed a growing mode through the walls ("Divergence of
  B", the walls). $E = 0$ on the walls removes it: a restart of that run from
  464 s with it ran to 600 s (the end it was given), the largest DynSGS
  coefficient 6.3 at 600 s where it had reached 113 at 504 s.
- **The erratum's Fig. 3 on REF 2** (cluster run of the default deck to 2530 s,
  128 ranks; `density`, `dome`, `atmos`). Top of the dense dome on $x = 0$
  (log $\rho$ > −13), in Mm:

  | t (s) | 920 | 1380 | 1610 | 1955 | 2185 | 2530 |
  |---|---|---|---|---|---|---|
  | REF 2 | 3.45 | 7.6 | 9.6 | 11.4 | 9.6 | 7.6 |
  | erratum | 3.5 | 7.5 | 10.4 | 13.05 | 14.7 | 3.65 |

  The rise, the emergence and the dome agree to 1400 s. Then the dome climbs
  too slowly and stops at 11.4 Mm: the $B_x$ arrow at 11 Mm never turns (it
  does at 1955 s in the erratum), the flux rope at the apex forms small at
  (0, 10.5) Mm at 1955 s and is gone at 2185 s, where the erratum has it at
  (0, 13.5) Mm (the 14.7 Mm of the table is its top) and at (4.5, 13.5) Mm at
  2530 s. Beside the dome (8 ≤ $\lvert x\rvert$ ≤ 12.75 Mm) the chromosphere and the
  transition region do not move or heat; a layer up to 100 times the coronal
  density lies at 5–9 Mm from 1160 s (the erratum has it from ~1600 s), and the
  DynSGS coefficient is 0.5–3 $L_0v_0$ between 5 and 12 Mm from 464 s on, before
  anything has emerged. Open; see "DynSGS in the corona" below for the two
  measures under test (REF 4, `JEXPRESSO_LK10_ETA=induction`).
- **REF 4 with `JEXPRESSO_LK10_ETA=induction`** (cluster run to 2530 s, 128 ranks).
  Dome top, in Mm:

  | t (s) | 920 | 1380 | 1610 | 1955 | 2185 | 2530 |
  |---|---|---|---|---|---|---|
  | REF 4, `ETA=induction` | 3.9 | 8.6 | 9.6 | 9.8 | 9.7 | 9.5 |
  | erratum | 3.5 | 7.5 | 10.4 | 13.05 | 14.7 | 3.65 |

  Closer to the figure in its structure than REF 2, and the same stall: ahead
  of the erratum to 1380 s, flat at 9.5–10 Mm from 1610 s, the $B_x$ arrow at
  11 Mm never turns (the null stays below it), a small rope-like bump at
  (−1, 9) Mm from 1955 s, and the side regions ($\lvert x\rvert$ > 8 Mm) full of cool,
  dense plasma at 5–8 Mm, which the erratum's are not at 1380–1610 s. In the
  paper the dome is the expanding shell of the emerged tube, which stretches
  the quadrupole's null into a crescent-shaped current sheet and pushes it up;
  the rope forms in that sheet at its apex.

  *Not reconnection at the null.* `plot_leake2010.py null` on this run: the
  in-plane flux between the tube axis and the null on $x = 0$, which only
  reconnection at the null can lower, is 10.44 × 10¹⁰ Mx cm⁻¹ at $t = 0$, 10.01
  at 920 s, 9.86 at 1380 s and 9.82 at 1955 s, and the overlying flux 4.23,
  4.20, 4.11, 4.04: the dome keeps its flux while it stalls, and the DynSGS
  resistivity at the null is 0.003–0.3 $L_0v_0$ (mostly 0.02–0.1). The tests
  `ETAFAC=0.1` and `ETACAP=flow` (above) were for this hypothesis.

  *The lifted plasma stops draining.* The mass in $3 < y < 20$ Mm,
  $\lvert x\rvert$ < 13 Mm, against $t = 0$: 34 at 1044 s, 606 at 1160 s, 362 at
  1380 s, 258 at 1610 s, then 245–282 to 2530 s. The dome rises while half of
  what it lifted drains (1160–1610 s) and stops when the draining stops; in the
  erratum the dome's density falls about tenfold between 1610 and 1955 s while
  it rises. The dense plasma also reaches 0.5–1.2 Mm above the null (the dome
  top of the table is above the top of the emerged flux), into the overlying
  field, which only a cross-field (numerical) transport can do. The REF 2 run
  shows why it can stay: beside the dome its dense layer at 7.5 Mm has
  log ρ ≈ −12 at 1.3–2 × 10⁵ K: a gas pressure of 10–15 dyn cm⁻², 30–45 times
  the initial coronal one there, and a pressure scale height of 3.5–5.5 Mm, so
  it sits on its own pressure; plasma lifted from the chromosphere and cooled
  by its expansion, as in Lare (no heating there but its shock viscosity), has
  a scale height of a fraction of a Mm (0.16 Mm at 6000 K) and drains.

  `plot_leake2010.py mass` on this run, for the plasma in that box:

  | t (s) | 1160 | 1380 | 1508 | 1610 | 1740 | 1856 | 1955 |
  |---|---|---|---|---|---|---|---|
  | mass / mass(0) | 606 | 362 | 272 | 258 | 279 | 282 | 260 |
  | fraction below 2 × 10⁴ K | 0.99 | 0.88 | 0.52 | 0.23 | 0.08 | 0.02 | 0.00 |
  | $p$ / initial $p$ there (median) | 12 | 2.0 | 21 | 76 | 108 | 113 | 75 |
  | mean $v_y$ (km/s) | +9.3 | −4.0 | −5.0 | −1.9 | +2.0 | +0.7 | −2.2 |
  | DynSGS viscous heating time (s) | 968 | 94 | 352 | 588 | 1060 | 2830 | 11300 |
  | DynSGS Ohmic heating time (s) | 537 | 555 | 2790 | 2930 | 2180 | 4950 | 31000 |

  The plasma arrives cold and cools further as it expands (to twice the
  coronal pressure at 1380 s), draining at 4–5 km/s. Between 1392 and 1740 s it
  is heated to 2 × 10⁴–2 × 10⁵ K, its pressure rises fifty-fold, the draining
  stops, and so does the dome. The heating times are those in which the DynSGS
  dissipation of the snapshot ($\rho\mu\lvert\nabla\mathbf v\rvert^2$ and
  $\eta\lvert\nabla\times\mathbf B\rvert^2$ summed over the plasma denser than
  10⁻¹³ g cm⁻³) supplies that plasma's internal energy: 100–600 s while it
  heats. The mass change is the flow's: after 1392 s the part not carried
  across the box's edges is at most 17 units in an interval and +17 net in
  1140 s, of 260, so the DynSGS mass diffusion does not feed it. See "The heat of the DynSGS" above; the tests are
  restarts of this run from 920 s with `JEXPRESSO_LK10_HEAT=shock` and
  `JEXPRESSO_LK10_HEAT=0`.

  *The restart with `JEXPRESSO_LK10_HEAT=shock`* (from 920 s; cluster). Against
  the run above, in the same units:

  | t (s) | 1380 | 1508 | 1610 | 1624 |
  |---|---|---|---|---|
  | mass / mass(0): all heat, `shock` | 362, 363 | 272, 240 | 258, 201 | 259, 202 |
  | fraction below 2 × 10⁴ K | 0.88, 0.96 | 0.52, 0.70 | 0.23, 0.42 | 0.21, 0.39 |
  | $p$ / initial $p$ (median) | 2.0, 1.5 | 21, 1.8 | 76, 19 | 78, 21 |
  | mean $v_y$ (km/s) | −4.0, −4.4 | −5.0, −7.7 | −1.9, −5.7 | −1.4, −5.0 |

  (The restart's pressure is against its 920 s state, the initial one above
  ~4 Mm.) The heating comes a hundred seconds later and is four times weaker at
  1624 s, and the plasma is still draining at 5 km/s where it had stopped. The
  run then failed at 1713 s (non-finite, in the top of the dome, 1.5 < x < 4.7,
  6.7 < y < 11 Mm), of the time step and not of the heating: the fast CFL
  number had reached 0.64 at 1624 s, above the 0.62 at which the DynSGS
  coefficient on its first-order cap is unstable (Meshes, time step and cost).
  `:dsgs_stab_cap` is the remedy, on by default since. Continued from its
  1508 s snapshot with it (and $\Delta t = 1.5\times10^{-4}$ for margin) the run
  reaches 2530 s, and its density figure has the erratum's sequence: a flux
  rope, a closed field line around a dense core, at the apex of the dome at
  1955 s, at (0.5, 8.3) Mm; at 2185 s it has moved to the right and down, to
  (6, 6.5) Mm; at 2530 s it is gone and the column above the tube has drained
  to coronal density from ~4.5 Mm up (erratum: the rope at (0, 13.7) Mm at
  2185 s, at (4.5, 13.7) Mm at 2530 s, the column coronal above 3.65 Mm). It
  happens some 5 Mm lower and 200 s earlier than in the erratum, and the $B_x$
  arrow at 11 Mm still does not turn. `JEXPRESSO_LK10_HEAT=shock` and
  `JEXPRESSO_LK10_ETA=induction` are the defaults since this run.

  Here: the switches act, are recorded, and a restart that changes them
  reports it and continues; the first 40 s of the run on REF 1 (the corona is
  still quiet then) are the same with `ETACAP=flow` to round-off.
  `JEXPRESSO_LK10_SENSOR=residual` does not lower the coronal coefficient: in
  the kick test of "Positivity" it is 0.45 against 0.59 (median, 11.6 s) and,
  with `ETA=induction`, the resistivity 0.25 against 0.018 — the element
  residual of the induction equation is larger than the legacy sensor's
  assembled one.

## Files

| file | content |
|---|---|
| `user_inputs.jl` | the deck: environment switches, Δt, DynSGS, dissipation, positivity floors |
| `initialize.jl` | atmosphere, tubes, quadrupole/dipole, the discrete curl, $c_h$, $\Delta h$, the $\epsilon_0(\rho)$ table; the restart and `lk10_case.txt` |
| `user_flux.jl` | units and constants of the paper, ideal GLM-MHD flux with $\gamma = 5/3$ |
| `user_source.jl` | gravity, energy relaxation, GLM damping, damping layers |
| `user_bc.jl` | closed walls: $\rho\mathbf v = 0$, $\psi = 0$ |
| `user_primitives.jl` | the reference-weighted DynSGS primitives; output $T$, $\beta$, $\log_{10}\rho$ [g cm⁻³], $B_z$ [G] |
| `LK10.geo`, `LK10_ref{1,2}.msh` | the graded mesh |

Source changes for this case (all opt-in; decks that set none of the keys
are unchanged):

- The 2D `:mhd_dissipation` hook (per-element `νphys`, `ηphys` in
  `compute_dsgs_viscosity!(::DSGS_MHD, ::NSD_2D)`, as in 3D) is shared with
  `fluxEmergenceLeake2022_2D`, which added it in parallel; this case adds
  `:mhd_anomalous_eta` to it: `_mhd_dissipation_2d` (`rhs.jl`) adds
  `_mhd_anomalous_eta_2d!` (`SGS.jl`) of the present state to the static
  $\eta$ at every call.
- `:dsgs_curl_resistivity`: `_dsgs_curl_resistive_2d!` (`SGS.jl`), called from
  `viscous_rhs_el!(…, ::NSD_2D)`, replaces the Laplacian of $B_x$, $B_y$ by
  the strong-form curl of the nodal $\eta j_z$ (section "Divergence of B").
- The keys are documented in `src/io/mod_inputs.jl`.
