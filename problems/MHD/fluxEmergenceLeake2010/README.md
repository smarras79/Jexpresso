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
editing the deck. The log lists what was read
(`Switches set in the environment (JEXPRESSO_LK10_*): ...`), and a
`JEXPRESSO_LK10_*` variable that the code does not read stops the run with the
list of those it does: a misspelling, or a switch newer than the checkout the
job runs from, would otherwise be ignored without a word and the run be the
one it was meant to change (that happened: The viscosity on the velocity,
below).

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
| `JEXPRESSO_LK10_EINT` | `off` (default), `lowbeta`, `all` | the energy equation: the total energy's own (`off`); or, where the internal energy is less than 10⁻³–10⁻² of the total (`lowbeta`) or everywhere (`all`), the changes of the kinetic and magnetic energies plus the internal energy's own equation, as in Lare, so that the pressure of a drained low-β region is not the remainder of two large energies (Positivity, below: The internal energy where it is a small share). Opt-in while the cluster test runs |
| `JEXPRESSO_LK10_STABCAP` | default 0.2; 0: none | the DynSGS coefficient never above this viscous number $\mu\Delta t/\Delta x^2$, which the explicit step still integrates (Meshes, time step and cost) |
| `JEXPRESSO_LK10_DSGS_STAGES` | `all` (default), `1,5`, `1` | the stages of the Runge–Kutta scheme at which the DynSGS terms are evaluated: all five; the first and the last, held in between; or the first alone (the default of `JEXPRESSO_LK10_STABCAP` is then 0.08). The coefficient is that of the first stage, with half the $C_1$ (Meshes, time step and cost: The DynSGS at two stages of the five) |
| `JEXPRESSO_LK10_C1` | default 1 | the DynSGS constant $C_1$ |
| `JEXPRESSO_LK10_ETACAP` | `fast` (default), `flow` | DynSGS resistivity capped at $C_{max}\Delta(\lvert v\rvert + c_f)$ like the viscosity, or at $C_{max}\Delta\lvert v\rvert$ (Verification, the REF 4 run) |
| `JEXPRESSO_LK10_ETAFAC` | default 1 | factor on the DynSGS resistivity (the `:μ` multipliers of the B slots), a sensitivity test; a factor above 1 never carries the coefficient past the stability cap (`JEXPRESSO_LK10_STABCAP`) |
| `JEXPRESSO_LK10_RHOFAC` | default 1 | factor on the DynSGS mass diffusion (the `:μ` multiplier of the ρ slot), a sensitivity test |
| `JEXPRESSO_LK10_VISC` | `momentum` (default), `velocity`, `nodal` | the viscosity (DynSGS and ν) as the diffusion of the momentum $\rho\mathbf v/\rho_e$, or as a stress on the velocity gradient with the momentum of the mass the DynSGS diffuses; that stress weighs every node of an element with the density of its lightest (`velocity`), or each node with its own, up to what the lightest takes at the stability limit (`nodal`) (Equations and dissipation: The viscosity on the velocity). Opt-in while the cluster test runs |
| `JEXPRESSO_LK10_RHOSENS` | `on` with `JEXPRESSO_LK10_VISC=velocity` or `nodal`, `off` with `momentum`; or `a0,a1,C` | the mass slot of an element whose density has a node-to-node part (the highest Legendre mode of $\log_{10}(\rho/\rho_e)$ along any of its lines, from `a0` to `a1` dex) gets at least `C` $\Delta\lvert(u, v)\rvert$; `on` is `0.03,0.15,1` (`:dsgs_rho_sensor`; Equations and dissipation: The density sensor) |
| `JEXPRESSO_LK10_VISCFAC` | default 1 | factor on the DynSGS viscosity (the `:μ` multipliers of the momentum slots, applied after its caps; the case's ν is not scaled), a test of its drag on the plasma that drains from the dome (`drain`, its `visc` column) |
| `JEXPRESSO_LK10_HEAT` | `shock` (default), `viscous`, a fraction | what the heat the DynSGS dissipates does (`:dsgs_heating`): `shock`, only the explicit ν, η and the DynSGS viscosity in compression heat; `viscous`, all the viscous dissipation (DynSGS and ν) and the explicit η heat, the DynSGS resistivity not, as Lare's default tensor shock viscosity heats (The heat of the DynSGS, below); or that fraction of all of it into the plasma (1 conserves energy, 0 loses it all) |
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
  and thermal-energy residual for the low-β corona; $C_1 = 1$,
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

  *The internal energy where it is a small share* (`JEXPRESSO_LK10_EINT=lowbeta`;
  `:mhd_internal_energy`, `src/kernel/operators/rhs.jl`). The velocity-form run of
  22 Nov (REF 4) drains its dome down to the density floor between 2185 and 2530 s,
  $\rho \approx 10^{-15}$ g cm⁻³ in a field of ~50 G: $\beta \sim 10^{-5}$–$10^{-7}$, the
  internal energy that share of the total (the erratum's drains too: its centre,
  $\lvert x\rvert \le 3$ Mm, is at $10^{-12.6}$ g cm⁻³ at 5–7 Mm at 1955 s and at
  $10^{-13.6}$–$10^{-14.8}$ above 5 Mm at 2530 s, `profile` below). The pressure is then what the divergence of
  $E$'s flux leaves after the kinetic and magnetic energies have taken their
  changes, and node by node those agree only to the aliasing error of the
  collocation derivative of products, which repeats in every element and is far
  above $p$. `noise` on that run: at 1955 and 2185 s 10–20 % of the dome's nodes
  (5.8–7.9 Mm) sit at the temperature floor, in rows of nodes one per element in
  height (0.189 Mm apart, the element 0.1875 Mm), while the density there is still
  smooth (node-to-node 0.01–0.04 dex); at 2530 s the density has followed
  (0.06–0.13 dex in 4.3–9.5 Mm, 17–26 % of the nodes at a floor): the horizontal
  stripes of its density maps. With the option, where the internal energy is less
  than $10^{-3}$ of the total (blended up to $10^{-2}$; so off where $\beta \gtrsim 0.01$,
  in the corona and below it), the energy slot's inviscid rate is the chain rule of
  the others plus the internal energy's own equation,
  $\mathbf u\cdot R_{\rho\mathbf u} - \tfrac12\lvert\mathbf u\rvert^2R_\rho + \mathbf B\cdot R_{\mathbf B} + \psi R_\psi + R_e$,
  $R_e$ the rate of $\partial_te + \nabla\cdot(e\mathbf u) + p\nabla\cdot\mathbf u = S_e$: the pressure
  follows that equation, all of whose terms are proportional to $p$, as Lare's
  specific internal energy does (the energy-consistent DynSGS above is the same
  chain rule for the dissipative terms). Total energy is not conserved where it
  acts, by the error it removes. The flux's pressure clamp (`LK10_PFLOOR`, $10^{-9}$,
  above $\rho T_{min}$ of cold plasma at the density floor) is then lowered to
  $10^{-15}$. Tests (local, `noise`-like statistics of $p$): the circularly
  polarized Alfvén wave of `alfvenWaveBabbarRanocha2026` at $\beta = 2\times10^{-6}$
  ($p_0 = 10^{-6}$, $N = 4$, 16² elements, no DynSGS), whose exact pressure is $p_0$
  at every point: after two periods the total-energy form has $p/p_0$ between
  −2.1 and 3.5 (rms error 1.2, 842 of 4225 nodes negative, node-to-node part 1.0),
  the internal-energy form 1 ± 4×10⁻⁶ (node-to-node 3×10⁻⁶); linearly polarized
  ($\lvert\mathbf B\rvert^2$ varies), −3.0 to 5.4 against 0.91–1.11, the wave's own
  compression, smooth to 2×10⁻⁵. At $\beta = 0.2$ (the test's own $p_0$) `lowbeta` is
  off and bit-identical, and the total-energy form is accurate (5×10⁻⁷; `all`,
  2.8×10⁻⁷). In this deck: the background (REF 1, no tube) stays at rest to the
  bit with `lowbeta` and with `all`; and a drained hole in the corona where the
  quadrupole's field is strongest (REF 1, no tube; at (12, 4) Mm, radius 1.5 Mm,
  0.22 of the density and 0.002 of the temperature of the atmosphere, an
  outflow along x of 65 km/s; $\beta$ down to $10^{-6}$), 20 s. Total energy: the
  temperature on a vertical line through the hole is off its neighbours at every
  fourth node, the element faces (1.6, 0.10 and 0.44 $T_{ph}$ among 7–14 at 8 s,
  spikes at 10 s), the floor fires from 1.9 s (10⁴ hits by 17 s), the
  node-to-node part of $\log T$ along y is 0.02–0.10 dex against 0.003–0.01 along
  x. `lowbeta`: the same line is smooth, no floor hit until 11.4 s, 0.002–0.009
  dex up to 10 s. From 11.4 s a spot of 1–4 nodes at the hole's centre, in the
  temperature gradient of the reference's transition region (48 to 134 $T_{ph}$
  between 3.5 and 4.5 Mm), reaches the floor in both runs; it is not rows, and it
  is not traced here.

  The cluster run (23 Nov: the velocity-form run of 22 Nov restarted at 1955 s,
  whose temperature-floor rows had not yet reached the density, with
  `JEXPRESSO_LK10_EINT=lowbeta` and its other switches): `noise` at 2185 and 2530 s
  shows no rows of nodes at the temperature floor. The density stripes remain: in
  the draining legs of the arcade at 2185 s, with almost no floor contact, and
  across the dome at 2530 s, as rows of nodes at the density floor. The pressure
  noise was real, but it is not what makes the stripes (The density sensor,
  below).

  *What the density floor adds.* The mass the floor adds keeps the node's
  momentum and total energy: the node is braked ($v$ by $\rho/\rho_{min}$) and the
  kinetic energy it loses becomes heat there. In the drained dome the kinetic
  energy is 5–50 times the internal energy, so every hit brakes and heats its
  node, and a row of nodes held at the floor drains differently from the rows
  beside it. Evidence (`elements` on the `EINT=lowbeta` run of 23 Nov at 2530 s):
  where the elements on the lines $x = -2, 0, 2$ Mm hold nodes at the density
  floor (1–8 of their 25), the node-to-node part of the vertical velocity is
  9–29 km/s, where they hold none 0.1–8 km/s, while the velocity form's stress
  varies by 0.6–1.5 from one element to the next there. A floor whose added
  mass moves with the node's velocity and has its temperature, as more of the
  same plasma, was tried on restarts of that run from 2185 s and withdrawn: kept
  everywhere without a bound, the run went non-finite in 4 steps, at $x = 0$
  below 3 Mm; with the refill bounded by twice the node's mass, after 236 steps,
  at $x$ = 14–16 Mm, 1–3 Mm high; only where the floor is a fair share of the
  reference density ($\rho_{min} \ge 0.1\rho_e$, the corona above ~3.75 Mm) and for
  up to 1 % of the node's mass per stage, after 1060 steps (2188.7 s), at
  $x$ = 13.6–14.3 Mm, 2.8–5.1 Mm high. The last two under the quadrupole's
  strongest field, outside the dome; the same restart with the floor above ran
  to 2530 s. The erratum's dome drains about as far as this one: at 2530 s its centre
  ($\lvert x\rvert \le 3$ Mm, its colours turned back into log ρ) is at 0.4–0.85 of
  the coronal density of the atmosphere above 7 Mm (−14.2 to −14.8, against
  −14.3 to −14.4), 0.3–0.6 dex above this floor ($0.2\rho_e$ in the corona).
  The run's dome is not thinner than the erratum's (`profile` on the
  `EINT=lowbeta` run of 23 Nov): below its top it is 0.05–0.9 dex denser at
  1955, 2185 and 2530 s, and thinner only at 9–14 Mm, where the erratum's dome
  reaches and this one does not. (The momentum form's dome, 20 Nov, hardly
  drains: its centre at 4–8 Mm is 1.1–1.5 dex denser than the erratum's at
  1955 s, 1.1–1.8 at 2185 s and 1.3–2.5 at 2530 s, with no node at a floor.) Its rows are a minority of its nodes: at
  2530 s 3–11 % of the centre's nodes at 6–10 Mm are at the density floor, its
  10th percentile at the floor (−15.0) under a median of −13.6 to −14.0; and at
  2185 s, when no node of the dome at 4–8 Mm is at a floor, the rows are
  already there. They are made by the advection of the density at the element
  faces, and the floor is where they end, not where they start (The density
  sensor, below).

  *Restarts* hold no DynSGS steps with the legacy sensor (`:dsgs_hold_steps`). A run
  starts with the coefficient at zero for two steps while the history of the
  sensor fills; at a restart that left a developed flow without DynSGS until its
  third step (the first restart with the velocity-keeping floor, above, was at
  its fourth when it failed, the coefficient on its cap). The legacy sensor
  reads $|\partial_tq|$, and on the history seeded with the restart state that is
  its first reading too. Restarting a REF 1 run at
  920 s, both ways end as far from the uninterrupted run at 974 s (velocity rms
  0.081 and 0.086 km/s, $\log\rho$ rms 0.0052 and 0.0050 dex, the same largest
  coefficient): the difference is the restart's, not the hold's.

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
  (the residual coefficient is $C_1\Delta^2\times$ the residual). The kick test
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
  the total energy is conserved (`JEXPRESSO_LK10_HEAT=1`). Lare does not return
  all of it: its remap does not conserve the kinetic energy and the loss is not
  returned as heat unless asked (`rke = .FALSE.` in the distributed
  `control.f90`: "This missing energy can be added back into the simulation as
  a uniform heating"), and what the remap diffuses of the field is lost too.
  What heats there is the paper's Eq. (5) and Lare's shock viscosity, which in
  the default build of the versions of the erratum's time (LareXd 2.6, March
  2011, and 2.9, August 2011: `lagran.F90`, with `#QMONO = -DQ_MONO` commented
  out in the distributed Makefile) is a tensor,
  $\sigma_{ij} = (\nu_1\rho c_fL + \nu_2\rho L^2\lvert s\rvert)(\varepsilon_{ij} - \tfrac13\delta_{ij}\nabla\cdot\mathbf v)$,
  $\nu_1 = 0.1$, $\nu_2 = 0.5$, $L$ the cell, "applied to all cells, not just
  those which are compressed, as significant shear forces may exist across
  expanding cells. The heating term is however always positive" (the LareXd
  manual, "Shock viscosity"); its heat, $\sigma_{ij}\partial_jv_i$, goes into the
  internal energy at every step. That is a first-order viscosity, $0.05\,c_fL$,
  in every cell, all of whose dissipation heats, in shear as in compression.
  `JEXPRESSO_LK10_HEAT=viscous` (`:dsgs_heating => "viscous"`,
  `src/kernel/physics/SGS.jl`) heats so: everything the DynSGS viscosity and ν
  dissipate, and the explicit η; the DynSGS resistivity and ψ diffusion do not
  heat, as the remap's diffusion of the field does not.
  `JEXPRESSO_LK10_HEAT=shock`, the default, was modelled on the scalar option
  (`Q_MONO`, a pressure added where the flow is compressed): the explicit ν and
  η in full, the DynSGS viscosity only where the flow is compressed,
  $\rho\mu(\nabla\cdot\mathbf v)_-^2$ and never more than it dissipates there, the
  DynSGS resistivity and ψ diffusion not at all; the rest of what the DynSGS
  dissipates is lost. In the dome of the velocity-form run of 22 Nov that is
  13–28 % of what the DynSGS viscosity dissipates (`heat`, 1610–1955 s; its
  last column is what Lare's tensor viscosity would heat with the same flow).
  The difference matters for the cold, dense plasma that the emerging field
  lifts into the corona, which has little internal energy: in the REF 4 run
  (Verification; the viscosity on the momentum) all of the heat (`1`) could
  supply it in 100–600 s, the plasma went from below 2 × 10⁴ K to 10⁵ K,
  stopped draining, and the dome stopped rising. `JEXPRESSO_LK10_HEAT=0`
  loses all of it (shocks then do not heat either: a test only). `shock` is the
  default since the cluster run of Verification, in which the lifted plasma
  keeps draining and the flux rope of the erratum forms; `viscous` has not run
  on the cluster. Here, the kick test of "Positivity" (REF 1, no floors,
  11.6 s) runs the same with the three: flows of 13.9, 13.8 and 13.8 km/s
  (1, `shock`, 0), the same smallest temperature, the pressure near the kick
  1.0 % (`shock`) and 1.3 % (0) below that with all the heat, the total energy change
  −2.98 and −3.22 against −2.21 × 10⁻¹⁰: of the heat the DynSGS dissipates
  there `shock` keeps a quarter, the compressed half of a kick whose velocity
  gradient is half compression. With the viscosity on the velocity and
  `JEXPRESSO_LK10_EINT=lowbeta` (24 Nov 2026): flows of 14.02, 14.09 and
  14.10 km/s (`shock`, `viscous`, 1), the same smallest temperature, the mean
  pressure within 2 Mm of the kick 2.26, 2.64 and 2.72 % above the initial one,
  the total energy change −2.85, −2.37 and −2.11 × 10⁻¹⁰: `viscous` returns
  two thirds of what `shock` loses, the rest being the DynSGS resistivity's.
  A restart that changes the mode reports it and continues.

- **The viscosity on the velocity** (`JEXPRESSO_LK10_VISC=velocity`;
  `:dsgs_velocity_viscosity`, `src/kernel/physics/SGS.jl`,
  `_dsgs_velocity_viscosity_2d!`). The conserved form diffuses the momentum
  like every other slot, $\nabla\cdot(\rho_e\mu\nabla(\rho\mathbf v/\rho_e))$.
  With $s = \rho/\rho_e$ that is
  $\nabla\cdot(\rho\mu\nabla\mathbf v) + \nabla\cdot(\mathbf v\,\rho_e\mu\nabla s)$:
  a viscous stress, and the momentum of a mass flux $\rho_e\mu\nabla s$. It
  is a viscosity only where the mass slot diffuses with the same coefficient
  (Guermond & Popov 2014, *SIAM J. Appl. Math.* 74:284: the mass flux of the
  regularization carries its momentum, and only then is the velocity of a
  uniform flow left alone). Here it does not: the REF 4 runs have
  `JEXPRESSO_LK10_RHOFAC=0.1`, and the paper's ν has no mass diffusion at all.
  The difference is momentum taken out of dense plasma that moves through
  thin plasma with no mass to go with it; the operator relaxes
  $\rho\mathbf v/\rho_e$, not $\mathbf v$. In the REF 4 run on the cluster
  (`drain`) that is the force on the cold plasma the emerging field lifts into
  the corona: on its steep field lines at 1160 s the diffusion terms carry
  1.24 times its weight along the field, while the stress on the velocity
  gradient alone is −1.42 of it, the pressure none and the paper's ν 0.03; it
  slides down at 6.5 km/s where free fall gives 50, is compressed (the `heat`
  table: $-p\nabla\cdot\mathbf v$ is 10⁴ K per 100 s on the flanks, ten times
  the DynSGS heat and a hundred times that of ν) until its own pressure holds
  it from 1610 s on, and the dome it loads stops at 9 Mm. That this force is
  what keeps the plasma there, and the plasma what stops the dome, is the
  reading to be tested by the cluster run with `VISC=velocity`; the test
  below shows the force and its removal on a blob. Whether a smaller DynSGS
  viscosity alone would do (`JEXPRESSO_LK10_VISCFAC`) is not known: the two
  restarts from 1044 s made on 21 Nov 2026 to test it, and then this option,
  ran from a checkout older than both switches, which ignored them (no
  `visc_factor` or `visc_form` in their `lk10_case.txt`), and are one run to
  every digit of the `drain` table. What they do show is what a restart alone
  changes: at 1160 s the plasma slides at 5.9 km/s against the uninterrupted
  run's 6.5, with 154.7 against 154.2 units of mass above 3 Mm.

  With `VISC=velocity` the momentum slots get instead

  $$R_{\rho\mathbf v} = \nabla\cdot\big(\rho_e\mu_\rho\nabla(\rho\mathbf v/\rho_e)\big) + \nabla\cdot\big(\rho_e s_{min}(\mu - \mu_\rho)_+\nabla\mathbf v\big),$$

  the diffusion of the momentum with the *mass* slot's coefficient
  $\mu_\rho$ (the momentum of the mass that slot diffuses, and the stress
  $\rho\mu_\rho\nabla\mathbf v$ that goes with it), and the rest of the
  viscosity, DynSGS and ν, as a stress on the velocity gradient. A velocity
  uniform across any density contrast is left alone to round-off, whatever
  the coefficients. $s_{min}$ is the smallest $\rho/\rho_e$ of the element:
  with the density of each quadrature node in the stress, a light node beside
  a dense one is relaxed at $\mu\rho_{dense}/\rho_{light}$, past the explicit
  limit the coefficient is capped for; with $\rho_e s_{min}$ no node is
  diffused faster than by the operator this replaces, and where the
  atmosphere is the reference's it is $\rho\mu\nabla\mathbf v$ exactly. The
  price is a weaker stress where the density jumps inside an element. In the
  drained dome of the `EINT=lowbeta` run of 23 Nov (`elements` on
  $x = -2, 0, 2$ Mm, 3.5–9.5 Mm) the elements that hold a node below the
  reference density have a stress coefficient $s_{min}(\mu - \mu_\rho)$ 65 times
  smaller than the others at 2185 s (median 0.33 against 22 $L_0v_0$), and 45 %
  of them have a node-to-node part of $\log\rho$ above 0.1 dex, against 6 % of
  the others (52 % and 22 % at 2530 s). Those are the elements of the rows of
  its density maps; the weak stress is where the rows are, not what makes them
  (a light node gives its element both). `JEXPRESSO_LK10_VISC=nodal`
  (`:dsgs_velocity_weight => "node"`) weighs each node with its own density,
  up to $h\,s_{min}$, $h = (\mu_{stab} - \mu_\rho)/(\mu - \mu_\rho)$ the margin to the
  stability limit of `JEXPRESSO_LK10_STABCAP`: the stress is $\rho(\mu - \mu_\rho)\nabla\mathbf v$
  wherever the element's contrast is below $h$ (a factor 20–600 in the dome),
  and no node is diffused faster than that limit. It restores the stress and
  leaves the rows (The density sensor). Both
  terms are diffusions with a non-negative coefficient; the mass flux written
  by itself, $\nabla\cdot(\mathbf v F_\rho) + \nabla\cdot(\rho_e s_{min}\mu\nabla\mathbf v)$,
  is not (its stress coefficient is negative at the dense nodes of an element
  whose contrast exceeds $\mu/\mu_\rho$): that first version went non-finite
  at the edge of the test blob below in 33 s. What the two terms dissipate,
  $(\rho\mu_\rho + \rho_e s_{min}(\mu - \mu_\rho)_+)\lvert\nabla\mathbf v\rvert^2 \ge 0$,
  is the heat of the energy slot (`shock`: its ν part in full, its DynSGS
  part in compression). At the first right-hand side with a mass flux the
  kernel checks itself and stops the run otherwise: the mass slot's flux is
  rebuilt from the state and compared with the assembled one (10⁻¹⁵ here),
  and where $\rho/\rho_e$ is uniform in an element the new momentum terms are
  compared with those they replace (10⁻¹³).

  *Test* (a copy of the case, REF 1, no tube; `JEXPRESSO_LK10_RHOFAC=0.1`,
  two-stage DynSGS, fast CFL 1): a blob at rest in pressure balance in the
  corona at (3.0, 5.55) Mm, 10³ times the coronal density, Gaussian radius
  0.3 Mm, on field inclined 35°; free fall along the field is
  $g\,b_y t$.

  | $t$ (s) | free fall (km/s) | `momentum` | `velocity` |
  |---|---|---|---|
  | 23 | 3.50 | 3.46 | 3.49 |
  | 58 | 8.97 | 7.75 | 8.91 |
  | 93 | 14.7 | 11.2 | 14.2 |
  | 116 | 18.6 | 12.4 | 17.7 |

  With `momentum` the blob levels off: from 100 s the diffusion terms carry
  0.8–1.0 of its weight with the shear part opposite (−0.4 to −0.9), the
  signature of the cluster run, and it is heated twofold (mean 6.7 × 10³ K
  against 3.8 × 10³ K). With `velocity` they carry 0.00 ± 0.02 throughout (what
  is left of the lag at 116 s is the pressure the blob builds ahead of it,
  0.10 of its weight), the DynSGS coefficient in the blob is the same in both
  (0.6–0.7 $L_0v_0$ at 93 s), and the run meets the stability cap once where
  `momentum` meets it three times. The default is unchanged to the bit (this case in its five-stage
  and two-stage settings, `fluxEmergenceSon2025DSGS` and
  `orszagTangBormanis2024`, 40–100 steps with an off-grid output, every file
  and the totals of every variable at every step). On the cluster (REF 4, from
  1044 s) `velocity` halves the plasma left above 3 Mm and keeps it cold, and
  the dome is no higher for it ("The erratum's Fig. 3 line by line" in
  Verification): it stays opt-in. `drain` reads the form from `lk10_case.txt`
  and evaluates that operator.

- **The density sensor** (`JEXPRESSO_LK10_RHOSENS=on`; `:dsgs_rho_sensor`,
  `src/kernel/physics/SGS.jl`, `_dsgs_logrho_mode`). The velocity-form runs on
  the cluster (22–23 Nov) have horizontal stripes in the drained dome: rows of
  nodes one element apart, light and hot at the pressure of the plasma around
  them, that end at the density floor (Positivity: What the density floor
  adds). They are made by the advection of the density at the element faces
  where the density is steep, and they are there because the mass slot has a
  tenth of the DynSGS coefficient (`JEXPRESSO_LK10_RHOFAC=0.1`).

  *A stripe in the making* (here, REF 1: the blob test of the previous item with
  100 times the coronal density and a Gaussian radius of 0.75 Mm at
  (3.0, 5.55) Mm, `RHOFAC=0.1`, `VISC=velocity`, `EINT=lowbeta`, and one row of
  13 nodes through it, at an element face, drained to 0.3 of the coronal
  density at $t = 0$). The drained row fills within 12 s and leaves the blob
  with a node-to-node part of $\log\rho$ of 0.04–0.07 dex (0.001–0.004 without
  the row). In the blob's upper edge, where the density falls from 50 to 3
  times the coronal one over 0.8 Mm, the next element face up then empties: on
  $x$ = 2.70 Mm the node at $y$ = 6.18 Mm has $\rho/\rho_e$ = 30, 19, 11, 7.0, 3.6
  and 0.20 (the floor) at 23, 35, 46, 58, 70 and 81 s, while the nodes 0.10 Mm
  below and 0.12 Mm above it stay at 52 ± 1 and 28 ± 1, at 0.86–0.97 of the
  coronal pressure throughout, so that its temperature goes from
  2.8 × 10⁴ K to 3.7 × 10⁶ K; the face at 6.85 Mm follows (4.2 to 0.20 by 93 s),
  then a node inside the element between them. In the density map these are
  the dark horizontal rows of the cluster's figures, with nodes at the density
  floor from 81 s. `massbudget` at 46 s, rates of $\ln\rho$ per 100 s:

  | $y$ (Mm) | $\rho/\rho_e$ | $-\partial_x(\rho u)$ | $-\partial_y(\rho v)$ | DynSGS | sum | observed |
  |---|---|---|---|---|---|---|
  | 6.08 | 53 | 0.5 | −0.3 | −0.2 | −0.1 | −0.1 |
  | 6.18 (face) | 11 | 0.4 | −8.0 | 2.4 | −5.2 | −4.1 |
  | 6.30 | 27 | 0.5 | −0.1 | −0.2 | 0.2 | 0.3 |
  | 6.52 | 10 | 0.5 | −1.7 | 0.1 | −1.0 | −0.8 |
  | 6.74 | 7.3 | 0.5 | 0.4 | −0.2 | 0.6 | 0.8 |
  | 6.85 (face) | 2.4 | 0.3 | −4.5 | 1.3 | −2.9 | −2.4 |

  The mass of a face node leaves through $-\partial_y(\rho v)$ alone. The plasma
  moves down at 2 km/s; the element below the face gives, for the derivative of
  $\rho v$ at its top node, mostly the difference to its next node (the last row
  of the LGL differentiation matrix of degree 4 is 0.5, −1.4, 2.7, −6.8, 5),
  where $\rho v$ is five times larger because that node is five times denser, a
  divergence; the element above gives a third of it with the other sign, and
  the two are averaged. The lighter the face node against the node downstream
  of it, the larger its loss: a dip of a fifth (37 without the drained row, 30
  with it, at 23 s) is at the floor a minute later. The DynSGS returns 2.4 of
  the 8.0: its residual sensor hardly sees a contact, which moves with the
  flow at a pressure that does not change ($\mu$ = 0.3 $L_0v_0$ there, a tenth of
  its first-order cap), and the density gets a tenth of that. `RHOFAC=0.1` was
  set to keep the edge of the dome one element thin (Verification: Sharpness)
  when the viscosity still diffused the momentum $\rho\mathbf v/\rho_e$, and with
  it the mass flux, with the whole coefficient; the velocity form leaves the
  density to its own slot. Neither the weight of the stress (`nodal`: the same
  rows, same nodes at the floor, velocities three times quieter) nor the floor
  (it is where the rows end) is in this.

  *The sensor.* The residual of a contact is small and the mass slot has a
  tenth of it, but the density shows it: a node that drops out of its
  neighbours is the highest Legendre mode of $\log_{10}(\rho/\rho_e)$ along its
  line of the element. `JEXPRESSO_LK10_RHOSENS=on` (`:dsgs_rho_sensor =>
  (0.03, 0.15, 1.0)`) reads it on every line of every element and, where it is
  above 0.03 dex, raises the mass slot's coefficient to at least
  $C\Delta\lvert(u, v)\rvert$, the first-order diffusivity of a contact carried by
  the flow, in full from 0.15 dex, and never above the stability limit of
  `JEXPRESSO_LK10_STABCAP`. Where the plasma is at rest, or its density smooth,
  nothing changes: a contact of 2 dex across one element (10–90 %) has a
  node-to-node part of 0.024 dex, so the sensor keeps an edge one element
  thick, what `RHOFAC=0.1` is for, and diffuses only what is sharper. In the
  test above it gives the stripe elements 0.11 $L_0v_0$ for the mass (0.004
  before) and the blob, with the same drained row:

  | $t$ (s) | 12 | 23 | 46 | 81 | 93 |
  |---|---|---|---|---|---|
  | node-to-node $\log\rho$ (dex), without / with | 0.066 / 0.019 | 0.038 / 0.013 | 0.067 / 0.012 | 0.124 / 0.008 | 0.105 / 0.010 |
  | smallest $\rho/\rho_e$, without / with | 1.00 / 1.00 | 0.98 / 0.98 | 0.85 / 0.87 | 0.20 / 0.78 | 0.20 / 0.78 |
  | nodes at the density floor, without / with | 0 / 0 | 0 / 0 | 0 / 0 | 7 / 0 | 8 / 0 |
  | $u$, node-to-node along $y$ (km/s), without / with | 0.01 / 0.01 | 0.09 / 0.03 | 0.19 / 0.04 | 0.43 / 0.04 | 0.71 / 0.05 |

  The blob falls as without the row (12.25 km/s at 93 s, 12.32 without it,
  14.50 in free fall), and its mass and temperature are those of the run
  without the row to 3 %. The default is unchanged to the bit (the kick test
  of "Positivity", the totals of mass and energy at steps 0, 100 and 200).

  *On the cluster* (REF 4, 24 Nov 2026): the `EINT=lowbeta` run of 23 Nov
  restarted at 1955 s with the sensor, every other switch the same, to
  2530 s. `noise` over $\lvert x\rvert \le 6$, 3–14 Mm, without / with:

  | | 2185 s | 2530 s |
  |---|---|---|
  | largest node-to-node $\log\rho$ along $y$ (dex) | 1.30 / 0.12 | 1.78 / 0.12 |
  | mean along $y$, along $x$ (dex) | 0.053, 0.016 / 0.023, 0.006 | 0.081, 0.032 / 0.029, 0.010 |
  | nodes at the density floor | 192 / 56 | 1150 / 228 |
  | of the nodes at 6.75–8.6 Mm | 0 / 0 | 4.7–10 % / 0–2 % |
  | nodes at the temperature floor | 165 / 256 | 30 / 15 |

  The density itself is the same: `profile`, the median of the centre
  ($\lvert x\rvert \le 3$ Mm) at 4–8 Mm is within 0.17 dex of the run without the
  sensor at 2185 and 2530 s, the flanks' within 0.07 (0.24 at 10 Mm), and the
  10th percentile of the centre at 6–8 Mm is −14.1 to −14.5 at 2530 s where it
  was at the floor (−15.0). What reaches the floor now is a patch at the top of
  the drained dome (7 % of the centre's nodes at 9 Mm at 2530 s), smooth from
  node to node (0.11 dex at most there), not rows. The temperature floor's nodes
  at 2185 s (1–4 % at 5.8–8.6 Mm, filled with mass before its heat) are gone by
  2530 s. The sensor is on by default with the viscosity on the velocity.

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
coefficient driven to its cap ($C_1$ = 1000): non-finite at 0.39 $t_0$ without
the stability cap, with a viscous number of 0.71 at the snapshot before; with
it, 0.26 and a clean run to 2 $t_0$ (flows decaying from 32 to 1.2 km/s,
smallest pressure, density and temperature unchanged, $\psi$ at round-off).

**Where the time goes.** The cost is the number of steps: the fast speed of
the corona (4500 km/s when the dome drains) over the smallest gap between two
nodes (32 km) is 7 ms, the step 3.5–4.6 ms, and 1610 s from the emergence to
the end are 350,000–460,000 steps of five stages. A stage on REF 4 (4 ranks,
`JEXPRESSO_ALLOC_SUMMARY=1`, which also times the pieces of the DynSGS block
since) takes 344 ms: 69 ms the inviscid right-hand side, 229 ms the DynSGS —
its diffusion operator on the nine equations 100 ms, the coefficient 38 ms,
the curl-form resistivity 35 ms, the energy slot 24 ms, the residual of the
legacy sensor 23 ms — and 46 ms the rest (assembly, exchanges, resets). On the
cluster the run makes 0.11–0.12 s a step on 128 ranks of one node, 270
elements a rank: 2.5 times the cost per element and core measured here on 4
ranks. A larger step is the first lever (the stability cap above keeps the
scheme linearly stable to a fast CFL number of 1.4; the log prints it at every
snapshot); the same 128 ranks on two nodes the second, if the node's memory
bandwidth is what the 2.5 is; and the DynSGS block, two thirds of a stage and
evaluated at each of the five, the third.

**The DynSGS at two stages of the five** (`JEXPRESSO_LK10_DSGS_STAGES`;
`:dsgs_stages`, `src/kernel/operators/rhs.jl`). With
`JEXPRESSO_LK10_DSGS_STAGES=1,5` the sensor and the coefficient are computed at
the first stage of the step only; the operators (diffusion, curl-form
resistivity, energy slot) are applied at the first and at the last stage, with
that coefficient; over stages 2–4 the assembled right-hand side of the first
stage is held. `=1` holds it over the whole step. The DynSGS term is then
first order in time, which a regularization can afford; everything else is
integrated as before. Four things make it the same model:

- *The work term, at every stage.* Of the three parts of the energy slot
  (Positivity) — the diffusion of the internal energy, the heat, and the work
  $W = \mathbf u\cdot R_{\rho\mathbf u} - \tfrac12\lvert\mathbf u\rvert^2R_\rho + \mathbf B\cdot R_{\mathbf B} + \psi R_\psi$
  that the other slots' terms $R$ do on the kinetic, magnetic and $\psi$
  energies — the first two are held and $W$ is taken again at every stage from
  the held $R$ and the stage's own $\mathbf u$, $\mathbf B$, $\psi$. It is a
  product of nodal values, so this is exact and costs one loop over the nodes,
  and the pressure sees nothing of the momentum and field diffusion at any
  stage, as before. Held, it would see $(\mathbf u - \mathbf u^n)\cdot R$,
  which is not small against $p$ at $\beta = 10^{-3}$.
- *The coefficient.* The legacy sensor at stage $i$ reads
  $1.5\,(1 - c_i)\lvert\partial_tq\rvert$: 1.5 at the first stage, 0.75 in the mean
  over the five (weights $b_i$). $C_1$ is halved (`:dsgs_C1` is 0.5
  `JEXPRESSO_LK10_C1`).
- *The stability cap.* The diffusion is no longer integrated by the scheme
  itself, and what the step takes changes (linear analysis on degree-4 elements
  in 2D; viscous number $\mu\Delta t/\Delta x_{min}^2$):

  | DynSGS at stages | fast CFL $\le$ 1.2 | 1.4 | 1.6 | `JEXPRESSO_LK10_STABCAP` |
  |---|---|---|---|---|
  | all five | 0.42 | 0.27 | 0.1 | 0.2 |
  | 1 and 5 | 0.6 | 0.6 | 0.6 | 0.2 |
  | 1 | 0.18 to 0.14 | 0.13 | – | 0.08 |

  With stage 1 alone the diffusion is a forward-Euler step, hence the smaller
  cap. With stages 1 and 5 its amplification factor is $1 + z + b_5c_5z^2$,
  $b_5c_5 = 0.147$, close to the two-stage Chebyshev polynomial
  $1 + z + z^2/8$: stable on the real axis to $\lvert z\rvert = 6.8$, where the
  five-stage scheme reaches 4.7. Of the pairs it is the widest (stages 1 and 3:
  0.37; 1 and 4: 0.23).
- *The step clock.* The first stage of a step cannot be told from the stage
  times, because a step is shortened to land on every output time and the steps
  after it start from there. The time loop stamps each step start (a callback,
  added only with this option). The five-stage run has no such stamp: after a
  shortened step its DynSGS history is taken at another stage of every step,
  and its sensor's mean over the stages is 0.78 to 1.11
  $\lvert\partial_tq\rvert$ instead of 0.75, a coefficient up to 1.5 times larger
  until the next output time (DSGS.md §4.4; 1.14 times in the mean over the
  cluster run of Verification from 1508 s, 1.20 from 920 s).

Measured here:

- *Cost* (REF 4, 4 ranks, 74 steps, `JEXPRESSO_ALLOC_SUMMARY=1`). A stage takes
  351 ms with the DynSGS at every stage (233 ms of it the block), 181 ms on
  average with stages 1 and 5, 164 ms with stage 1 alone: the step costs 0.52
  and 0.47 of the five-stage one, 1.94 and 2.14 times faster. With stages 1
  and 5 the block is computed 150 times instead of 372, 76 of them with the
  sensor (one a step, and one where each solve begins), and the held
  right-hand side with its work term is added 222 times at 3.2 ms; with stage
  1, 76 and 296. The counts are the check of the step clock: the 74 steps
  include one shortened to an output time and the ones after it.
- *The same solution.* The kick test of "Positivity" (REF 1, 2 ranks, floors
  off) with a time step that is exact in binary, $3/8192\,t_0$, and a snapshot
  every 128 steps, so that no step is shortened and the five-stage run takes
  its history at the first stage throughout. At 10.9 s, five stages / stages 1
  and 5: flows of 14.25 / 14.44 km/s (32 at the start), smallest pressure 0.653
  / 0.646 of the initial one, smallest density 0.918 / 0.916; pointwise the
  velocity differs by 0.064 km/s rms over the disturbed region (0.5 km/s at
  most), the pressure by 0.17 % rms (5 % at most), the field by 0.02 % rms.
  Stage 1 alone is stages 1 and 5 to four digits (0.0013 km/s rms): at these
  viscous numbers how the diffusion is integrated does not matter, and what
  separates the stage list from the five stages is the coefficient. It is a
  little smaller: at its later stages the five-stage sensor also reads the
  tendency of its own viscous term, which the $1.5\,(1 - c_i)$ of the
  calibration leaves out (its last-stage coefficient is 0.11 of the
  first-stage one in this test, not 0.04). With `JEXPRESSO_LK10_C1` = 1.15
  and 1.3 the stage list ends with flows of 14.31 and 14.18 km/s, on either
  side of the 14.25 of the five stages, and 0.037 and 0.045 km/s rms from
  them: 1.2 is the constant that reproduces a five-stage run whose history is
  in phase, and 1.2 times the factor of the step clock above, 1.4, one like
  the cluster run. For scale: the same five-stage run
  with the default snapshot times instead, which shorten ten of its steps and
  move its history from stage to stage (a mean sensor of 0.98 instead of
  0.75), ends with flows of 14.07 km/s and differs from the first by 0.061 km/s rms
  (0.43 at most; pressure 0.11 % rms, 1.7 % at most). The stage list changes
  the solution as much as the choice of the output times changes the
  five-stage one, in the other direction.
- *Stability.* The kick test at five times the time step (fast CFL number
  1.13) with the coefficient driven to its caps ($C_1\times1000$), to 2 $t_0$:
  five stages with the cap at 0.2 and at 0.5, stages 1 and 5 at 0.2 and 0.5,
  stage 1 at 0.08 and 0.2 all complete, with flows decaying from 32 km/s to
  1.2–1.3 km/s and a smallest pressure of 0.95 of the initial one (2.2 km/s
  and 0.91 with the cap at 0.08, which is what holds the coefficient in this
  test), $\psi$ at round-off. None fails, those beyond their linear limit
  included (five stages at 0.5, stage 1 at 0.2): the defaults are set from the
  linear analysis, not from these runs.
- *The default is untouched.* Without the switch the kick test is bit for bit
  that of the source before the change (every field of its eleven snapshots,
  the DynSGS coefficients and `conservation.dat` included).
- *Restart.* A five-stage run of the default problem (REF 1) continued from
  its 3 s snapshot with `JEXPRESSO_LK10_DSGS_STAGES=1,5`: the log says
  `numerical switches changed: dsgs_stages: all there, 1,5 here` and the run
  goes on.
- *The default problem* (tube 2, medium quadrupole, floors on; REF 1, 2
  ranks) from $t = 0$ to 10 s, five stages / stages 1 and 5: flows of 1.010 /
  1.022 km/s as the tube starts to rise, smallest pressure 0.880 / 0.878 of the
  initial one, 0.005 km/s rms apart (0.018 at most), pressure 0.05 % rms, no
  floor hit in either.
- *The other switches.* The kick test to 2.4 s with stages 1 and 5 and, in
  turn, `JEXPRESSO_LK10_SENSOR=residual` ($C_1$ not halved: that sensor is a
  residual at the first stage already), `_ENERGY=total` (the energy slot held
  whole), `_ETA=all` and `_CURLB=0`: all run, with flows of 18.9–19.0 km/s at
  the end as with the defaults.

Not run here: the emergence itself on REF 4 with the stage list. The check
that costs least on a cluster is a restart of a completed run over its hardest
stretch, 1610 s to 1955 s (the dome drains, the fast CFL number is at its
largest, the stability cap acts), against the snapshot it already has:

```bash
JEXPRESSO_LK10_RESTART_DIR=$RUN JEXPRESSO_LK10_RESTART_TIME=1610 JEXPRESSO_LK10_TEND=1955 JEXPRESSO_LK10_DSGS_STAGES=1,5 mpiexec -n <ranks> julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2010
```

(with the mesh, the rank count and the `JEXPRESSO_LK10_DT` of that run; add
`JEXPRESSO_LK10_C1=1.4` for the coefficient that run had, see *The same
solution* above: without it the stage list runs with about 30 % less DynSGS,
which is a second change).

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
  switches may change (`JEXPRESSO_LK10_FLOORS`, `_CMIN`, `_CURLB` and those of
  the DynSGS: `_ETA`, `_ETACAP`, `_ETAFAC`, `_VISC`, `_VISCFAC`, `_RHOFAC`, `_HEAT`, `_STABCAP`,
  `_SENSOR`, `_ENERGY`, `_DSGS_STAGES`, `_C1`, and `_TFLOOR`; the log says which did): that is
  how a run continues past a failure (below), or with another setting from the
  snapshot where it starts to matter.
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
python3 tools/plot_leake2010.py forces   $RUN --rope     # Erratum Fig. 4: vertical forces through the flux rope (default: x = 0, 2185 s, 10-17 Mm)
python3 tools/plot_leake2010.py diag     $RUN --dt <Δt>  # per snapshot: largest speed, CFL, current; smallest p/p0, beta
python3 tools/plot_leake2010.py null     base=$RUN       # on x = 0 vs time: the null above the dome, the flux below it, the mass lifted
python3 tools/plot_leake2010.py lines    run=$RUN        # Erratum Fig. 3 line by line: the height on x = 0 of each of its field lines, in the run and in the figure
python3 tools/plot_leake2010.py core     run=$RUN        # the plasma at the tube axis vs time: height, density against the atmosphere's, B_z/rho against the initial tube's
python3 tools/plot_leake2010.py noise    $RUN --times 2185,2530   # stripes and grain: the node-to-node part of log rho per element, along y and x, and the nodes at the floors
python3 tools/plot_leake2010.py elements $RUN --times 2185,2530 --x=-2,0,2   # which quantity repeats element by element on vertical lines: s_min, the DynSGS coefficients and the velocity-form stress, node-to-node parts of log rho, log T, u, v, floor nodes
python3 tools/plot_leake2010.py profile  $RUN --times 1955,2185,2530   # the dome's density by height, centre and flanks, beside the erratum's Fig. 3, with the nodes at the floors
python3 tools/plot_leake2010.py massbudget $RUN --time 2185 --x=-2,2 --ylim 3,5.7   # a stripe in the making: node by node on a vertical line, the advection of mass along x and y, the DynSGS mass diffusion, and the change to the next snapshot
python3 tools/plot_leake2010.py --paper-lines density $RUN   # the density figure with exactly the erratum's field lines (any map takes --paper-lines)
python3 tools/plot_leake2010.py mass     $RUN            # the plasma above 3 Mm: mass budget, temperature, pressure, draining, heating time
python3 tools/plot_leake2010.py drain    $RUN            # why that plasma drains or not: temperature, field it sits on, speed and support along the field
python3 tools/plot_leake2010.py heat     $RUN            # what changes its temperature: each term as the code applies it, against the observed rate
python3 tools/plot_leake2010.py column   $RUN --times 1380,1610 --x 0   # the state on a vertical line, as a table
python3 tools/plot_leake2010.py balance  $RUN --times 1380,1610 --x 0   # the vertical forces on that line, and their integrals over the dome
bash tools/lk10_report.sh $RUN --base $RUN0 --log job.out            # all of these tables and the log's lines in one text file
```

Figures are written to `figs_leake2010/` (`--outdir`) as
`flux_emergence_Leake2010_<variable>.png` (`--prefix`, `--format`):
`…_density`, `…_Bz`, `…_mu_dsgs_momentum`, `…_mu_dsgs_induction`,
`…_mu_dsgs_mass_energy`, `…_dome`, `…_axis`, `…_null`, `…_lines`, `…_core`, `…_noise`, `…_eshear`, `…_atmos`,
`…_profiles`. Every map has a title and a labelled colour bar naming what its
colours show.

`axis` draws the height of the tube's axis (where $B_x$ on $x = 0$ changes sign,
the paper's definition) against time, with the two curves of the erratum's
Fig. 2(a). It is the check of the first phase, the rise through the convection
zone to the surface, which the erratum's axis reaches at about 920 s; it then
climbs on, to 0.41 Mm at 1380 s and 0.48 Mm at 1610–1955 s (the same heights
as the O-point of its Fig. 3, `lines` below), and how far the axis gets above
the surface is a measure of how much of the tube has emerged. The erratum's
curves are for a field-free corona, which makes no difference below the
surface. The curve has one point for every snapshot of the
directories given, so a run made in several jobs is `axis LABEL=DIR1,DIR2,...`
(the command prints the snapshots it read from each); `--tmax` is only the
limit of the time axis.

`dome` follows the highest dense plasma on $x = 0$ ($\log_{10}\rho > -13$), whatever
it is. In the erratum's Fig. 3 that is the top of the emerging dome to 1610 s
(3.5, 7.5, 10.4 Mm), then the flux rope that forms at its apex and rises (13.1
Mm at 1955 s, 14.7 Mm at 2185 s), and at 2530 s the base of a column that has
emptied (3.65 Mm), the rope having slid to $x$ = 4.5 Mm. Between 1610 and
1955 s the erratum's dome loses most of its plasma (its density falls about
thirtyfold): the rise to 12–13 Mm is that of a dome that has drained.

`lines` compares the run with the erratum's Fig. 3 field line by field line.
The figure's lines are contours of $A_z$ at a constant interval, and their
levels are known in this case's gauge ($A_z\to0$ at infinity, 820.49 G Mm at the
quadrupole's null): at the 83 places where the lines of panel (a) cross
$x = \pm8,\dots,\pm12$ Mm, far from the tube at 920 s, the case's analytic
quadrupole takes the values $820.0 - 61.35\,k$ G Mm to 2 G Mm (rms) — which
also says that the quadrupole here is the paper's. Line 0 hugs the separatrix;
lines 1, 2, … lie below it in the emerged flux and above it in the overlying
field, and the initial tube holds 17 of them (1044 G Mm between its axis and
the null). A line is a material line until it reconnects, so its height on
$x = 0$ in the figure is a measurement of the paper's run. The tool holds them
for 1380–2530 s (`ERRATUM_FIG3_LINES`; the lines hidden by the figure's $B_x$
arrows are taken from $x = \pm2.3$ Mm) with the O-point of the tube
(`ERRATUM_FIG3_AXIS`: 0.46, 0.52, 0.46 Mm at 1380, 1610, 1955 s, against 0.41,
0.48, 0.48 Mm in the erratum's Fig. 2(a)). The command prints, at the figure's
times, the height of every line in the run beside the figure's, the axis, the
top of the emerged flux and the flux reconnected there ($820.49 - A_z$ at the
top); it writes the heights for every snapshot (`lines_<label>.csv`) and draws
them against time. The run's $A_z$ is the flux function of its field, shifted
so that the top of the column $x = 0$ (on the closed wall, where it cannot
change) has the quadrupole's value: on the initial state it is the analytic
$A_z$ to 0.03 G Mm, on REF 1 and on REF 4. `--paper-lines`, before or after
the name of any map's command, draws exactly these lines instead of 45
equally spaced ones, so that a panel compares with the erratum's line for line.

`noise` measures what a density map shows as grain or stripes: in each
element, the amplitude of the highest Legendre mode of $\log_{10}\rho$ along
$y$ and along $x$ (in dex: the node-to-node part, which the element cannot
tell from noise), mapped over a window (`--xlim`, `--ylim`, default
$\lvert x\rvert \le 6$, $3 \le y \le 14$ Mm) with the nodes held at a positivity
floor (grey: the temperature floor, black: the density floor), and tabulated
per row of elements with the share at each floor. Horizontal stripes that the
discretization makes show as a large amplitude along $y$ in the rows where
they are; plasma at the density floor as floor nodes; real layers as neither.
On the blob test of "The viscosity on the velocity" (REF 1, 58 s) the mean
amplitude along $y$ is 0.006 dex with the momentum form and 0.013 with the
velocity form, the excess in the corona the blob falls through. On the
velocity-form run of 22 Nov (REF 4) it found rows of nodes at a floor one
element apart in height in the drained dome, over a still smooth density at
1955 and 2185 s and with it at 2530 s (Positivity: The internal energy where it
is a small share).

`core` follows the plasma at the tube axis (the O-point of the in-plane field
on $x = 0$): its height, with the erratum's Fig. 2(a) and the O-point of its
Fig. 3; its density over the atmosphere's at that height, which is the tube's
buoyancy (below 1: lighter); its $B_z/\rho$ over that of the initial tube
(1094 G per photospheric density for tube 2); the in-plane flux around the
axis; and over the whole tube $\int B_z^2\,dA$ and the flux-weighted mean of
$B_z/\rho$. In ideal 2.5D MHD the plasma at the O-point keeps its $B_z/\rho$
(there is no in-plane field there to carry $B_z$ away along it) and the flux
around the axis is kept: what a run loses of them, it loses to its
resistivity and its mass diffusion. A table, `core_<label>.csv` and
`…_core.png`.

`drain` says why the lifted plasma drains or does not. For the plasma denser
than $10^{-13}$ g cm⁻³ above 3 Mm, per snapshot: its mass; its median
temperature and the mass fractions at the temperature floor, below and above
$2\times10^4$ K; the fraction in field too weak to guide it ($\beta > 1$) and
the fraction on field lines within 15° of the horizontal, where gravity has
no component along the field; and, on the steeper field lines, its mean speed
down the field and the share of its weight along the field that its pressure
gradient carries (1: hydrostatic, held by its own pressure; 0: nothing in the
way), the share the diffusion terms' force on the velocity carries
($\rho\,dv/dt = R_{\rho v} - vR_\rho$, the DynSGS and the case's $\nu$, with
the code's operators: a drag on plasma that falls), of it the case's $\nu$
alone and the shear part $\nabla\cdot(\rho\mu\nabla v)$ that a viscosity
acting on the velocity would give, the sum of pressure and viscous shares,
whose shortfall from 1 accelerates it down the field, and the median DynSGS
momentum coefficient applied there. The case diffuses the momentum
$\rho v/\rho_e$, not the velocity: what is not shear is momentum carried out
of dense plasma that moves as a block, which the mass equation does not match
(no mass diffusion goes with $\nu$, and `JEXPRESSO_LK10_RHOFAC` times the
DynSGS one). Along a field line only
pressure, gravity and viscosity act, so a dome that keeps its plasma shows as
one of: warm with a support near 1; flat field; $\beta > 1$; or cold on steep
field with no pressure support and yet slow, held by the viscosity. The
derivatives are those of the elements (the LGL differentiation matrix, and the
CG weak form with the lumped mass for the viscous operator): on the initial
atmosphere of REF 1 the support is 1.00 in the corona and 0.99 in the
chromosphere (1.11 with the finite differences used before 21 Nov 2026).

`heat` says what changes the temperature of that plasma. For the plasma denser
than $10^{-13}$ g cm⁻³ above 3 Mm (`--ymin`, `--ymax`, `--thr`; `--thr -99` for
all of it), split into the center ($|x| < 2.5$ Mm, where the dome's field lies
flat) and the flanks, each term of its temperature equation as the code applies
it with the run's `JEXPRESSO_LK10_HEAT` (`lk10_case.txt`), as a mass-weighted mean $DT/Dt$ in K per
100 s: compression and expansion ($-p\nabla\cdot\mathbf v$, ideal), the case's
viscosity (in full), the DynSGS heat (`DSGS heat`: its viscosity where the flow
is compressed with `shock`, all of it with `viscous`), the
case's anomalous resistivity (in full; zero while
$|\nabla\times\mathbf B| < j_{crit}$, a 40 G reversal thinner than about 50 m),
the DynSGS diffusion of the internal energy's departure from the reference
atmosphere, the DynSGS mass diffusion (which leaves the internal energy where it
is, so that $T = p/\rho$ rises where it takes mass away and falls where it
brings it), the energy relaxation, and the GLM damping of $\psi$, whose energy
$\tfrac12\psi^2$ the pressure gains ($\alpha_p c_h/\Delta h\,\psi^2$; `--ch`,
the run's log). Each is evaluated from the snapshot with the code's operators
(the reference-weighted $\nabla\cdot(\rho_e\mu\nabla((q - q_e)/\rho_e))$, the
reference atmosphere rebuilt from its formula: on the $t = 0$ snapshot every
term is below $10^{-4}$ K per 100 s), with the elements' derivatives and the CG
weak form on the mesh's lattice of rectangles, where the 2D operator is exactly
the sum of the 1D ones along the rows and the columns (60 times closer to the
exact operator than finite differences on a graded test lattice). Then what the
DynSGS dissipates and the run does not heat, and last (`Lare visc`) what Lare's
default tensor shock viscosity (The heat of the DynSGS, above) would heat with
the snapshot's flow, $\sigma_{ij}\partial_jv_i$ with its cell $L$ = 0.048 Mm
(`--lare-dx`; the paper's grid in the figure window, to a factor ~1.3). Where two
snapshots are at most 20 s apart (the off-grid outputs, 12–17 s apart) the
observed $DT/Dt = \partial_tT + \mathbf v\cdot\nabla T$, node by node, is set
against the budget on the nodes that stay in the set, below $2\times10^5$ K,
and change by less than 30 % (a front sweeping over a node is not an advection a
finite difference follows): the difference is what the budget misses, the
temperature floor and its own error. On REF 1 at 920–928 s it is 1 % on the
flanks and 5 % in the center above 3 Mm.

`profile` puts the dome's density beside the erratum's, height by height: for
bands 0.5 Mm thick around $y$ = 2, 3, …, 15 Mm, the run's median log ρ over the
centre ($\lvert x\rvert \le 3$ Mm) and over the flanks ($3 < \lvert x\rvert \le 8$ Mm),
the erratum's at that height (its Fig. 3 turned back into log ρ through its
colour bar, medians over the same ranges of $x$; at 920 s, before the tube
reaches 5 Mm, it gives −14.32 for the corona at 5–15 Mm where the atmosphere
here has −14.18 to −14.48), their difference, the run's 10th percentile, and
the share of its nodes at the density and temperature floors. The erratum's
centre: −12.3 to −12.6 at 4–8 Mm at 1955 s, −12.5 to −13.3 at 2185 s, and at
2530 s −13.1 at 4 Mm, −13.6 at 5 Mm and −14.2 to −14.8 above 7 Mm, at or below
the atmosphere's coronal density; its flanks are denser below 7 Mm (−11.2 to
−12.9 at 2530 s).

`massbudget` is for the stripes (The density sensor): node by node on a
vertical line, between two consecutive snapshots, $\rho/\rho_e$ now and at the
next snapshot, and as rates of $\ln\rho$ per 100 s the advection of mass along
$x$ and along $y$ (the code's element derivatives, averaged at the faces), the
DynSGS mass diffusion with the mass slot's coefficient of the snapshot, their
sum and the change observed, then the node-to-node part of $\log_{10}(\rho/\rho_e)$
of the node's element along $y$ (what the sensor reads). A stripe in the making
is one row of nodes, at an element face (marked `|`), that loses its mass
through $-\partial_y(\rho v)$ alone while the rows 50 km above and below keep
theirs, and to which the diffusion returns a fraction. Use two snapshots that
are close (the off-grid outputs, 12–18 s after 1380, 1610, 1955, 2070 and
2185 s).

`balance` is the erratum's Fig. 4 as a table, for the question of what holds
the dome: on a vertical line, at each time, gravity, the magnetic tension, the
magnetic pressure force and the gas pressure force against height (10⁻⁵ dyn
cm⁻³, the elements' derivatives), their sum, and each integrated over three
layers (dyn cm⁻²): from 3 Mm to the top of the dense plasma, from there to
where $B_x$ changes sign (the top of the emerged field on that line), and from
there to 14 Mm (the field above). A dome at rest has a sum near zero in each;
the terms say what carries the weight of the lifted plasma and what presses on
it from above. In the initial potential field the tension and the magnetic
pressure force cancel to six digits.

`tools/lk10_report.sh RUN [--base RUN0] [--log FILE] [--py PYTHON]` writes
`lk10_report_<run>.txt`: the configuration of each directory of the run, the
`drain`, `heat`, `mass`, `null` and `dome` tables (the last two with the base
run beside it), `balance` at $x$ = 0 and 5 Mm, `column` at $x$ = 0, 5 and 10 Mm at
1380, 1610, 1955 and 2185 s, and from the job's log the set-up, the reports of
the positivity floors and of the stability cap, the CFL numbers at each
snapshot and the end of the run. It reads the snapshots only. `--base` is the
run that started at $t = 0$, whose first snapshot is the unit of mass, so that
the tables of two restarts can be compared.

`forces` draws the erratum's Fig. 4: on a vertical line at one time, gravity
$-\rho g$, the magnetic tension $(\mathbf B\cdot\nabla)B_y/4\pi$, the magnetic
pressure force $-\partial_y(B^2/8\pi)$ and the gas pressure force
$-\partial_yP$ in 10⁻⁵ dyn cm⁻³ against height, and their sum (and a CSV of
them). The defaults are the paper's ($x = 0$, 2185 s, 10–17 Mm); `--rope` puts
the line through the run's own flux rope at `--time` (the deepest O-point of
the field above 3 Mm, from the flux function; it lists what it finds) with
7 Mm of height around it, and `--x`, `--ylim` set them by hand. Checked on the
initial state, where they must cancel: gravity against the gas pressure to
3 % at 2 Mm (differences on the LGL nodes across a 150 km scale height) and
0.5 % above 3 Mm, tension against magnetic pressure to 0.1 %.

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

**Which stage.** A snapshot holds the coefficient of the last evaluation
before it. With the DynSGS at every stage (the default) that is the last stage
of the last step, and the legacy sensor is not the same at every stage: to
leading order it reads $1.5\,(1 - \tau/\Delta t)\,\lvert\partial_tq\rvert$, with
$\tau$ the time since the DynSGS history was last taken, anything up to a step
(DSGS.md §4.4). The mean over the stages, which is what acts, is 0.75 to 1.11
$\lvert\partial_tq\rvert$ depending on where in the step the history is taken, and
that changes at every snapshot, because the step before a snapshot is
shortened to land on it. The factor of the value written is anywhere between
0 and 1.5: for a run with $\Delta t = 1.5\times10^{-4}$ continued from 1508 s
the gate's logic gives 0.69 at 1610 s, 0.76 at 1955 s, 1.39 at 2185 s and 0.02
at 2530 s, against acting means of 0.75, 0.92, 0.75 and 0.92. Maps of the
coefficient at different times of such a run are therefore not on one scale:
the one at 2185 s is 1.8 times too high, the one at 2530 s several times too
low (not forty: where the factor is small the sensor reads what is left, the
tendency of its own viscous term; in the kick test of "Positivity" a factor of
0.06 against a mean of 0.75 left the snapshot's coefficient four to five
times below the acting one, and one of 1.05 against 0.92, 1.14 times above it,
as predicted). Values on a cap, the explicit viscosity and the anomalous
resistivity are not affected. With `JEXPRESSO_LK10_DSGS_STAGES=1,5` or `1`
(Meshes, time step and cost) the snapshot holds the coefficient of the step,
the one that acts, at every time.

*The first snapshot of a restart* is written with the coefficient of the
first evaluation, made before the sensor has a history: on REF 1 it is 0.05 of
what the run it continues wrote at the same instant, with the state equal to
round-off, and at the next snapshot the two runs' coefficients agree to four
digits (a restart at 3.75 s with every switch unchanged, compared at 4.5 s;
the states then differ by 10⁻⁶). In the REF 4 run the dense lifted plasma has
0.017 $L_0v_0$ in the first snapshot of a restart at 1044 s against 0.075 in
the run's own: that is this, not a setting.

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
  by its expansion, unless something heats it (in Lare its shock viscosity,
  which in the default build heats in shear too: "The heat of the DynSGS"), has
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

  *Sharpness: the same nodes as the paper, five times its edge.* In the
  figure window REF 4 has a node every 47–49 km on average (elements of
  187–197 km, degree 4, gaps of 32–64 km; 544 × 480 nodes across the window);
  the paper has 37.5 km at the surface and 150 km at the top of its 90 Mm box,
  some 38–62 km in the window if its stretching is linear (it does not say).
  Outside the window this mesh is coarser than the paper's coarsest (a node
  every ~470 km beyond 30 Mm from the axis and above 60 Mm). Yet the paper's
  dome has an edge, ripples along its top and curls at its flanks, and this
  one a gradient. Measured on the two published images (colour turned back
  into log ρ along vertical lines at $x = \pm2, \pm3$ Mm): the top edge of the
  dome, log ρ from −12 down to −14, is 0.06–0.17 Mm thick in the erratum at
  1380 and 1610 s (one to three pixels of the image: at most three or four of
  its cells) and 0.6–0.8 Mm here, 13–17 nodes. That is the DynSGS mass
  diffusion: the ρ slot carries the element's $\mu$, 0.2–0.5 $L_0v_0$ at the top
  of the dome (`null`), and an edge that has diffused for 300 s with 0.16
  $L_0v_0$ is 0.7 Mm thick (two decades of an error-function profile,
  $3.2\sqrt{\mu t}$). Lare does not diffuse mass at all: its Lagrangian step
  carries the density with the flow and its remap is limited, so a contact
  stays three cells wide. `JEXPRESSO_LK10_RHOFAC=0.1` would bring the edge to
  ~0.2 Mm, one element, which is what this mesh can hold; the kick test of
  "Positivity" runs the same with 1, 0.1 and 0 (smallest density 0.91 of the
  initial one, smallest pressure 0.62, flows of 13.8–14.0 km/s, no floor
  needed), but it has no such contact, and the test is a restart of the run
  from 920 s. Beyond that, half the element size in the window is 2.4 times
  the elements and twice the steps.

  Here: the switches act, are recorded, and a restart that changes them
  reports it and continues; the first 40 s of the run on REF 1 (the corona is
  still quiet then) are the same with `ETACAP=flow` to round-off.
  `JEXPRESSO_LK10_SENSOR=residual` does not lower the coronal coefficient: in
  the kick test of "Positivity" it is 0.45 against 0.59 (median, 11.6 s) and,
  with `ETA=induction`, the resistivity 0.25 against 0.018 — the element
  residual of the induction equation is larger than the legacy sensor's
  assembled one.

- **The erratum's Fig. 3 line by line: less of the tube emerges here**
  (22 Nov 2026; `lines`, `core`). Two cluster runs on REF 4 continue the same
  run from 920 s: with the DynSGS viscosity on the momentum (20 Nov, to 2530 s)
  and, from 1044 s, on the velocity (`JEXPRESSO_LK10_VISC=velocity`, 22 Nov, to
  2530 s).

  *The velocity form drains the dome and does not raise it.* Against the
  momentum form, the plasma above 3 Mm (in units of 7.64 × 10⁴ g cm⁻¹) is 110
  against 166 at 1380 s, 75 against 130 at 1610 s and 66 against 133 at 1955 s,
  its flanks stay at 0.3–3 × 10⁴ K where they had reached 3.5–6.6 × 10⁴ K, and
  the interior of the dome falls to log ρ = −12.3 to −12.8 by 1955 s, the
  erratum's level (−12.6). The top of the dome is at 8.0, 8.4, 8.3, 9.0 Mm at
  1380, 1610, 1955, 2185 s (momentum form: 7.7, 8.9, 9.6, 8.5; erratum: 7.5,
  10.4, 13.05, 14.7): the drained dome does not rise either.

  *Weight is not what holds it.* The erratum's density on the centre line
  (its colours turned back into log ρ, median over $\lvert x\rvert \le 3$ Mm; the
  figure's atmosphere far from the tube is the case's to 0.07 Mm in height) is
  −10.5 inside the dome and −10.1 in a shell under its top at 1380 s, −11.4 and
  −11.0 at 1610 s, −12.6 and −12.0 at 1955 s: the plasma above 3 Mm weighs
  440, 120, 23 and 7 dyn cm⁻² at 1380, 1610, 1955 and 2185 s. In the runs it
  weighs 272, 250 and 9 dyn cm⁻² (velocity form, `balance`) and about 185, 35
  and 87 (momentum form, from the 0.5 Mm table of its report) at the first
  three. The erratum's dome is the heaviest of the
  three at 1380 s and rises 3 Mm in the next 230 s; the momentum-form dome has
  a third of its weight at 1610 s and has stopped.

  *The field lines.* Heights on $x = 0$ of the erratum's lines (`lines` above;
  line 1 is next to the null, 17 next to the axis), read off its figure, and
  the same lines in the two runs (from the $x = 0$ tables of their reports,
  $B_x$ every 0.5 Mm; `lines` on the snapshots is exact):

  | line | 2 | 4 | 6 | 8 | 10 | 12 | 14 | axis |
  |---|---|---|---|---|---|---|---|---|
  | 1380 s, erratum | 5.37 | 4.28 | 3.39 | 2.70 | 2.10 | 1.64 | 1.16 | 0.46 |
  | velocity form | 5.36 | 3.88 | 2.87 | 2.14 | 1.59 | 1.12 | 0.69 | 0.08 |
  | momentum form | 5.16 | 3.74 | 2.79 | 2.09 | 1.56 | 1.12 | 0.70 | 0.00 |
  | 1610 s, erratum | 7.67 | 5.89 | 4.54 | 3.48 | 2.70 | 2.01 | 1.52 | 0.52 |
  | velocity form | 5.69 | 3.96 | 2.85 | 2.08 | 1.53 | 1.10 | 0.68 | 0.07 |
  | momentum form | 5.99 | 4.06 | 2.90 | 2.12 | 1.55 | 1.12 | 0.71 | 0.07 |
  | 1955 s, erratum | 9.40 | 6.75 | 4.91 | 3.65 | 2.73 | 1.98 | 1.44 | 0.46 |
  | velocity form | 6.69 | 3.56 | 2.30 | 1.58 | 1.13 | 0.79 | 0.47 | 0.04 |
  | momentum form | 5.34 | 3.24 | 2.18 | 1.52 | 1.11 | 0.78 | 0.45 | 0.01 |

  At 1380 s, when the tops of the three domes agree, every line from the
  fourth down already lies 0.4–0.6 Mm lower in the runs, and so does the axis
  (0.0–0.1 Mm against 0.46 Mm; the erratum's Fig. 2(a) has 0.41 Mm). From then
  on the erratum's lines all rise, the low ones too (line 13 from 1.38 to
  1.75 Mm by 1610 s), and its axis reaches 0.5 Mm; the runs' lines stay where
  they are to 1610 s and sink afterwards. The flux of the emerged system above
  3 Mm (up to the null) is 318–341 G Mm in the runs at 1380 s and 321–332 at
  1610 s, where the erratum has 433 and 561 less what has reconnected (less
  than one line, 61 G Mm, through 1955 s: its line 1 is still there); above
  5 Mm at 1955 s, 55–61 against 361. The dome of the runs is the top of a tube
  that has come up 0.4 Mm less and then stops emerging; its field at 3–6 Mm is
  0.7 of the erratum's at 1380 s (133 G against ~180 G at 3 Mm, from the line
  spacing), and at 1955 s half of it at 5 Mm and a fifth at 7–8 Mm. Reconnection is not ahead of the
  erratum's until then: the runs have reconnected 7–18 G Mm at 1380 s, 15–33
  at 1610 s and 79–108 at 1955 s. The shear field is not what lifts the
  erratum's dome either: the erratum has 0.3–0.5% of the magnetic energy above
  5 Mm in $B_z$.

  *Why: the tube loses its buoyancy on the way up* (REF 1, field-free corona,
  4 ranks; `core`). At the O-point of the tube there is no in-plane field, so
  in ideal 2.5D MHD the plasma there keeps its $B_z/\rho$, and with it the
  magnetic pressure that makes it lighter than its surroundings. In the REF 1
  run of "Rise through the convection zone" above it has 0.77 of the initial
  value at 46 s, 0.61 at 420 s and 0.51 at 1000 s; the in-plane flux around
  the axis falls to 0.74; and the density at the axis, 0.69 of the atmosphere's
  at $t = 0$, is 0.91 of it at 420 s — where the axis leaves the erratum's
  curve — 0.98 at 742 s and 1.05 at 1000 s: the tube reaches the surface with
  no buoyancy left. The rate is that of the DynSGS resistivity at the tube
  (0.018–0.03 $L_0v_0$, 17–29 km² s⁻¹, with either sensor): the peak of a
  Gaussian of radius $a$ decays at $4\eta/a^2$, which summed over the run gives
  a fall of $\ln(B_z/\rho)$ of 0.80 against the 0.62 measured. With today's
  deck the first 93 s are

  | at 93 s | $B_z/\rho$ at the axis / initial | in-plane flux / initial | $\rho/\rho_e$ at the axis |
  |---|---|---|---|
  | `RHOFAC=1` | 0.76 | 0.981 | 0.831 |
  | `RHOFAC=0.1` | 0.78 | 0.968 | 0.819 |
  | `RHOFAC=0.1`, `ETAFAC=0.1` | 0.87 | 1.018 | 0.793 |

  (0.96 at $t = 0$ is the discrete axis of REF 1). The first seconds cost the
  most: the tube contracts when the relaxation cools it, the DynSGS
  coefficients are 0.15–0.18 while it does, and a REF 1 run of 17 Nov had lost
  7.5% after 6 s; at 46 s the three variants have lost 17, 13 and 9%. After
  that the loss is 1.4 × 10⁻³ s⁻¹ with the full resistivity and nothing
  measurable with a tenth of it (0.87 at 46 and at 93 s; the density deficit
  stays at 20.7%, the value the relaxation leaves, where it shrinks to 18% and
  17% in the other two). Lare has no explicit resistivity there (the
  anomalous one never switches on), and the erratum's figure shows that its
  numerical one hardly touches the tube's core: the flux function of the
  analytic tube is −223.1 G Mm at its axis, its line 17 is the contour
  −223.0 (to the 2.4 G Mm of the calibration), and that line is there, around
  the axis, in every panel to 2530 s. The flux function at the paper's axis
  has therefore changed by less than about 3 G Mm in 2530 s; at ours it has
  risen by 230 G Mm by 920 s on REF 1 (field-free) and by 43 G Mm on REF 4
  (1044 → 1001 G Mm between axis and null; 52 G Mm by 1380 s, where only
  line 17 is gone). The B_z profile, peaked where the flux function is flat
  (tube 2 has no current on its axis), loses more: what the REF 4 axis has
  kept of $B_z/\rho$ is to be read with `core` on its snapshots (one node of
  the 20 Nov run's table, at the axis at 1380 s, has 0.65 of the initial
  value and 1.19 times the atmosphere's density there).

  *A tenth of the resistivity is not the fix.* The same REF 1 rise with
  `ETAFAC=0.1` (and `RHOFAC=0.1`) to 1400 s: the axis keeps 0.79–0.87 of the
  initial $B_z/\rho$ through 920 s (the 13 Nov run: 0.52) and reaches the
  surface 0.15 Mm behind the erratum's curve (0.27 Mm in the 13 Nov run, and
  0.15 Mm on REF 4 with the full resistivity): the lag from ~500 s on, at a
  rise speed of 2.1 against 2.45 km/s, is not the resistivity's. Then the run
  went non-finite at 1062 s, where the tube flattens against the surface
  (first non-finite nodes on the axis at −0.37 Mm and at (−6.4, 0.75) Mm),
  after the first-order cap had run to 10⁵⁸–10¹⁰⁰: the resistivity is needed
  there.

  *What the DynSGS sees in the rising tube.* The coefficient there is
  0.02–0.07 $L_0v_0$, 10–40 times the paper's kinematic viscosity at those
  depths (0.0013–0.008), so that the tube rises at a Reynolds number of ~10
  against several hundred in the paper. The legacy sensor is the time
  derivative, ~$1.5\lvert\partial_t q\rvert$ against its mean (DSGS.md, "Selected
  stages"), so it fires on any structure that moves, resolved or not. A
  count over the tube's elements ($\lvert B_z\rvert > 1$ kG, REF 1, first 140 s,
  an instrumented copy of the kernel): the largest normalized residual is the
  energy's (at 14 s) and then the thermal residual's in every one; the
  induction residuals alone would give 0.01–0.03, the coefficient applied is
  0.03–0.09.

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
- `:dsgs_stages`: the DynSGS block at selected stages of the Runge–Kutta
  scheme (`rhs.jl`, `_visc_lag`; the step-start callback in
  `TimeIntegrators.jl`; section "Meshes, time step and cost").
- The keys are documented in `src/io/mod_inputs.jl`.
