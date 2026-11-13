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

All the simulations of the paper are selected from the environment, without
editing the deck:

| variable | values | meaning |
|---|---|---|
| `JEXPRESSO_LK10_TUBE` | `1`, `2` (default), `3`, `0` | tube 1 constant twist, 2 increasing twist ("low tension"), 3 decreasing twist; 0 no tube |
| `JEXPRESSO_LK10_FIELD` | `quadrupole` (default), `dipole`, `none` | coronal field |
| `JEXPRESSO_LK10_STRENGTH` | `strong`, `medium` (default), `weak` | quadrupole $q_d$ = 11.7, 5.85, 3.9 $\times10^7$ G cm ("197, 98, 66 G"); dipole $\lvert B(0,0)\rvert$ = 104, 52, 26 G |
| `JEXPRESSO_LK10_REF` | `1`, `2` (default), `3` | mesh (below) |
| `JEXPRESSO_LK10_TEND` | seconds, default 2530 | final time |
| `JEXPRESSO_LK10_DOUT` | seconds, default 116 | snapshot interval (the figure times are always added) |
| `JEXPRESSO_LK10_DT` | $t_0$ units | overrides the time step |
| `JEXPRESSO_LK10_MCZ` | default 1.3 | polytropic index $m$ of the convection zone, $p\propto T^{m+1}$ |
| `JEXPRESSO_LK10_RELAX_YMIN` | $L_0$ units, default $-\infty$ | apply the energy relaxation only above this height |
| `JEXPRESSO_LK10_TUBE_EQ` | `thermal` (default), `relaxed` | tube in thermal equilibrium with its surroundings (the paper), or on the relaxed relation $\rho T_0(\rho) = p$ |
| `JEXPRESSO_LK10_CURLB` | `1` (default), `0` | DynSGS resistivity of the in-plane field in curl form, or the component Laplacian of the parent case |
| `JEXPRESSO_LK10_NU` | `dynamic` (default), `kinematic` | the paper's dynamic viscosity in the convection zone, held at its photospheric kinematic value above; or that kinematic value everywhere |
| `JEXPRESSO_LK10_CMIN` | default 0 | DynSGS background floor $C_{min}$ |
| `JEXPRESSO_LK10_FLOORS` | `0` (default), `1` | positivity floors as the Runge–Kutta stage limiter |

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
  is the only shock capturing. The background floor
  $C_{min}$ of the parent case is off: $C_{min}\Delta c_f$ with $c_f\approx270$
  in the corona would be a diffusivity of ~4 $L_0v_0$, two orders of magnitude
  above Lare's own estimate of its numerical resistivity.

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
The field-free corona allows a ~15 times larger step (`JEXPRESSO_LK10_DT`;
0.006 $t_0$ is safe on REF 1 until the tube reaches the photosphere). A
snapshot takes 84 MB on REF 2 (21 MB on REF 1): about 2.5 GB for the 30
snapshots of the default output.

## Figures

`tools/plot_leake2010.py` draws the paper's figures from the output (numpy and
matplotlib). `RUN` is `output/MHD/fluxEmergenceLeake2010/output-<date>`.

```bash
python3 tools/plot_leake2010.py density  $RUN            # Erratum Fig. 3: log rho, field lines, B_x arrows, 920-2530 s
python3 tools/plot_leake2010.py shear    $RUN            # Erratum Fig. 1: B_z (G) and field lines, 920-2254 s
python3 tools/plot_leake2010.py dome     tube2=$RUN      # top of the dense dome on x = 0 vs time, with the erratum's Fig. 3 values
python3 tools/plot_leake2010.py axis     tube2=$RUN      # Erratum Fig. 2(a): axis height vs time, with the erratum's curves
python3 tools/plot_leake2010.py eshear   tube2=$RUN      # Erratum Figs. 2(b,c), 5: E_shear above 1.2 and 5 Mm (Eq. 20)
python3 tools/plot_leake2010.py profiles $RUN            # Fig. 3(a,b) at t = 0, with the paper's values
```

The density panels use the paper's window ($\lvert x\rvert < 12.75$ Mm,
$-3 < y < 19.5$ Mm), colour scale ($\log_{10}\rho$ from −16 to −5.3 g cm$^{-3}$, jet),
field lines as equally spaced contours of $A_z$ (integrated from $B_x$, $B_y$
element by element, exactly for the degree-4 interpolant) and arrows along
$B_x$ at $x = 0$, $y$ = 3, 11, 16 Mm. The axis of the tube is where $B_x$ on
$x = 0$ changes sign from − to + (the paper's definition), the crossing nearest
the peak of $B_z$. `axis`, `dome` and `eshear` write one CSV per run.

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
  RISE_PLACEHOLDER
- **Not done here:** the emergence itself (after ~1000 s) and the
  quadrupole run to 2530 s, i.e. the comparison with the erratum's Fig. 3, which
  needs the cluster (above). `plot_leake2010.py density` and `dome` produce
  the comparison: the erratum's dense dome reaches 3.5, 7.5, 10.4, 13.05 and
  14.7 Mm on $x = 0$ at 920, 1380, 1610, 1955 and 2185 s (log $\rho$ > −13), and
  its flux rope sits at $(x, y) \approx (0, 13.5)$ Mm at 2185 s and $(4.5, 13.5)$ Mm
  at 2530 s.

## Files

| file | content |
|---|---|
| `user_inputs.jl` | the deck: environment switches, Δt, DynSGS, dissipation, optional floors |
| `initialize.jl` | atmosphere, tubes, quadrupole/dipole, the discrete curl, $c_h$, $\Delta h$, the $\epsilon_0(\rho)$ table |
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
