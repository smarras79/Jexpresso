function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        #
        # 2D CompEuler + Navier-Stokes: Mach-7.7 hypersonic flow over a
        # 15-degree compression ramp with laminar separation.
        #
        # THIS DECK IS problems/CompEuler/rampCaoEtAl2021 PLUS THE THREE
        # THINGS THAT CAME OUT OF THE Mach-7 DEBUGGING (ffs_step_M7,
        # ffs_step_M7_round, shock_circle_M7).  The original is untouched
        # and is still the paper-faithful reference; run both if you want
        # the comparison.  What is different here, and nothing else:
        #
        #   (a) :lkep => true with :volume_flux => ranocha().  The inviscid
        #       volume term becomes entropy-conservative flux differencing
        #       instead of the aliased collocation flux.  This is the
        #       stabilisation aimed at the actual M-squared mechanism --
        #       see the block at :lkep below and the header of user_flux.jl.
        #
        #   (b) :lpositivity => true.  A last-resort repair of rho and of
        #       the internal energy, applied to uaux once per RHS call
        #       (src/kernel/positivity/).  It is NOT a positivity-preserving
        #       scheme -- Zhang-Shu cannot be built on a continuous
        #       Galerkin space -- it is a floor that keeps one bad node from
        #       taking the run down, and it AUDITS ITSELF so you find out
        #       whether it ever fired and how much mass and energy it
        #       injected.  Read src/kernel/positivity/README.md before
        #       trusting a run it engaged on.
        #
        #   (c) ONE Mach number.  initialize.jl used to write 7.7 twice --
        #       once for the free stream and once, hardcoded, inside the
        #       recovery temperature of the starting profile.  It is
        #       const RAMP_MACH now, so a stepping-stone run at M = 7.0 is
        #       genuinely one line.
        #
        # And one deprecation fix: initialize.jl reads mesh.coords[dim,ip]
        # instead of the deprecated per-axis node arrays.
        #
        # WHAT IS DELIBERATELY UNCHANGED: every DynSGS setting.  :dsgs_sensor
        # stays "legacy", :dsgs_Cmax stays 0.1, :mu stays [1,1,1,1],
        # :ldsgs_nodal stays off.  Changing the sensor at the same time as
        # adding (a) and (b) would make the result unattributable, and that
        # mistake has already been made once in this debugging campaign.
        # If this deck survives longer than the original, (a) is why.
        #
        #   S. Cao, J. Hao, I. Klioutchnikov, H. Olivier, C.-Y. Wen,
        #   "Unsteady effects in a hypersonic compression ramp flow with
        #    laminar separation", J. Fluid Mech. 912, A3 (2021),
        #   doi:10.1017/jfm.2020.1093
        #
        # Flat plate L = 100 mm with a sharp leading edge, then a 15-degree
        # ramp also 100 mm long; M = 7.7, Re_L = 4.2e5, T_inf = 125 K,
        # p_inf = 760 Pa, isothermal wall at 293 K (Table 1, shock tunnel
        # TH2).  See initialize.jl for the conditions and ramp15.geo for
        # the grid.
        #
        # This is the x-y flow of the paper: the base flow its 3D DNS is
        # started from (Section 2.3) and the one its global stability
        # analysis linearises about (Section 3.2).  The ramp-induced
        # pressure rise separates the laminar boundary layer at x/L = 0.59
        # and it reattaches at x/L = 1.26, with a separation shock, a
        # reattachment shock and their interaction above the bubble.  The
        # streaks and the low-frequency unsteadiness of the paper are the
        # 3D instability of this flow and cannot appear here.
        #
        # THREE DECISIONS, all forced by that description:
        #
        # (1) :energy_equation => "energy".  Slot 4 is rho*E, not rho*theta.
        #     theta is an entropy variable: conserved across a contact, NOT
        #     across a shock, so the Euler-theta system carries the wrong
        #     speed for the separation and reattachment shocks no matter
        #     how it is stabilized.
        #
        # (2) :visc_model => DSGS().  Residual-based artificial viscosity as
        #     shock capturing (Nazarov & Hoffman, IJNMF 71:339-357, 2013,
        #     eq. 3.4-3.7): the coefficient is proportional to the local
        #     residual of the conservation laws, so it appears at the shocks
        #     and stays near zero over the smooth 90% of the field.  A
        #     constant coefficient large enough to hold a Mach-7.7 shock
        #     would drown the boundary layer this case exists to resolve.
        #     This is the same choice, and the same total-energy branch, as
        #     problems/CompEuler/ffs_step.
        #
        # (3) :lsutherland => true.  DynSGS is a SENSOR, so it vanishes
        #     exactly where the physics of this case lives: the laminar
        #     boundary layer is smooth and well resolved.  Without a
        #     molecular viscosity there is no boundary layer, no separation
        #     and no bubble — the case degenerates into inviscid flow over a
        #     ramp.  The flag adds mu(T) of Sutherland's law at every node
        #     on top of DynSGS (rhs.jl, _viscous_rhs_el_2d_dsgs!); the 2D
        #     assembly it feeds is the real Navier-Stokes viscous operator,
        #     deviatoric stress tensor + viscous work + Fourier conduction.
        #     The air constants below reproduce the paper's Reynolds number
        #     from its own p, T and u to within 0.5% — see the note in
        #     initialize.jl.
        #---------------------------------------------------------------------------
        #
        # TIME INTEGRATOR.  Explicit, and it should stay explicit.
        #
        # WHY NOT IMEX / HEVI.  The grid is anisotropic -- hx/hy reaches 23
        # at the wall -- so a directionally-split integrator looks
        # attractive.  The CFL is NOT anisotropic, and that is what
        # matters.  Swept over all 16,140 elements against a modelled
        # compressible laminar profile (Crocco-Busemann temperature,
        # Pohlhausen velocity, calibrated to the paper's delta = 1.38 mm
        # at separation):
        #
        #   streamwise only   (|u|+c)/dx      ->  dt = 1.77e-8 s
        #   wall-normal only  (|v|+c)/dy      ->  dt = 1.79e-8 s      ratio 0.99
        #   fully explicit    both            ->  dt = 1.27e-8 s
        #   y taken implicit  (HEVI-like)     ->  dt = 1.77e-8 s      1.40x
        #
        # The two anisotropies cancel.  Where dy is smallest (7.98e-6 m at
        # the wall) the wave speed is smallest too -- no slip pins u = 0 and
        # the isothermal wall pins T = 293 K, so |u|+c = 343 m/s.  Where the
        # wave speed is largest (|u|+c ~ 2000 m/s at the boundary-layer edge
        # and in the free stream) the wall-normal mesh has already stretched
        # past 2.5e-5 m.  A vertical implicit solve would buy 1.40x and cost
        # more than that per stage.
        #
        # The second reason is independent of the grid: at M = 7.7 the
        # ACOUSTICS ARE THE SLOW WAVES.  c/u = 0.13, so advection is the
        # fast part.  HEVI and IMEX in atmospheric codes take the pressure
        # terms implicitly precisely because M << 1 there makes sound the
        # fast, uninteresting wave; here that split removes the small
        # eigenvalue and leaves the large one explicit.  (src/kernel/
        # operators/imex.jl splits exactly that pressure Jacobian, has no
        # directional form, and is not wired to any live case.)
        #
        # WHY CarpenterKennedy2N54.  It is what the total-energy DynSGS path
        # is calibrated against -- ffs_step's :μ sweep was measured with this
        # stage layout, and with :dsgs_sensor => "legacy" the sensor takes a
        # BDF2 of the stage state at EVERY stage, so the stage layout is part
        # of the effective viscosity.  SSPRK54 trades imaginary-axis
        # stability for the SSP property, which is worth nothing here: DynSGS
        # is the shock capturing, not a TVD limiter.  If you want to
        # experiment, RDPK3SpFSAL49 is the one to try (built for compressible
        # DG, better stability per RHS call) -- but re-check the shocks, not
        # just that it runs.
        #
        # DO NOT set :ode_adaptive_solver => true.  params.Δt is the deck
        # constant (params_setup.jl:299), and DynSGS reads it twice: the
        # BDF2 weights are built on h = params.Δt and the once-per-step
        # history gate fires on time - thist >= 0.999*params.Δt (rhs.jl).
        # An integrator changing dt underneath that desynchronises the
        # sensor from the step it is supposed to measure, silently.
        #
        :ode_solver           => CarpenterKennedy2N54(),
        :tinit                => 0.0,
        # 2.0e-3 s is t*u_inf/L = 34.5, about 17 flow-throughs of the
        # 197 mm domain.  The paper's 3D field is fully developed by
        # t*u_inf/L = 60 (Section 2.4), but that clock starts AFTER the 2D
        # solution has converged; the 2D bubble itself settles in well
        # under half of it.  Watch the separation and reattachment points
        # in the output and extend with :lrestart if they are still moving.
        :tend                 => 2.0e-3,
        :lrestart             => false,
        :restart_time         => 0.0,
        #
        # TIME STEP.  Explicit, and the binding constraint is VISCOUS, not
        # advective -- the same lesson ffs_step's deck records ("DynSGS
        # saturates its own mu_max bound ... mu*dt/(rho*dx^2) is already
        # 0.22"), and getting it backwards here is what made the first
        # version of this deck blow up.
        #
        # THE MECHANISM.  DynSGS bounds its coefficient by the first-order
        # upwind viscosity, nu_cap = Cmax * Delem * (|u| + c).  Delem is the
        # SHORTEST element side (mesh.jl, compute_element_size!), but the
        # diffusion is resolved on the LGL spacing, which at nop = 4 is only
        # 0.17267 of it.  So wherever the sensor saturates, the diffusion
        # number it implies is
        #
        #     nu_cap*dt/dy_LGL^2 = (Cmax/0.17267^2) * (|u|+c)*dt/Delem
        #
        # i.e. 5.8*Cmax times the advective CFL at the same node.  At the
        # default Cmax = 0.5 that is 2.9x, and the viscous limit binds by a
        # factor of a few everywhere the shocks live.  THIS is where the
        # grid anisotropy bites -- not in the advective CFL, which the note
        # below the integrator shows is isotropic to 1%.
        #
        # :dsgs_Cmax => 0.1 (set under the DynSGS block) is the answer to
        # the 5.8: Nazarov's cap is meant to be the first-order upwind
        # viscosity at the NODAL spacing, and using Delem with Cmax = 0.5
        # overestimates it by exactly that factor on an LGL grid.  ffs_step
        # keeps 0.5 because its grid is isotropic and it simply accepts the
        # small dt; this grid cannot afford to.
        #
        # WHAT BINDS NOW.  With Cmax = 0.1 and the upstream strip removed
        # (see the mesh block), at a diffusion number of 0.22:
        #
        #   where                   Delem      |u|+c    nu_cap     dt
        #   leading-edge column     4.62e-5   1950     9.0e-3    1.6e-9
        #   mid-boundary-layer      1.48e-4   2094     3.1e-2    4.6e-9
        #   first wall element      4.62e-5    343     1.6e-3    8.9e-9
        #   top of the domain       4.53e-3   1950     8.8e-1    1.5e-7
        #
        # so ~1.6e-9, set by the leading edge, where the inflow holds the
        # free stream right down onto the wall-clustered mesh.  The
        # advective limit is 2.0e-8 (CFL = 1), so it is NOT what sets dt.
        #
        # 1.0e-9 is that with room to spare.  It is MEASURED, not derived:
        # 5.0e-9 died at step 3 and 1.0e-9 ran 304 steps before dying at the
        # strip, which is the cell this deck no longer has.
        :Δt                   => 1.0e-9,
        :diagnostics_at_times => (0:2.5e-5:2.0e-3),
        # Wall-clock note, not a setting: 2.0e-3 s at 1e-9 is 2,000,000
        # steps on 16,140 elements, about 380 core-hours (see README).
        # A long silence between the CFL/VTK lines is the run working, not a
        # hang; JEXPRESSO_STEP_HEARTBEAT=1 turns on a per-step trace.
        :lsource              => false,
        :SOL_VARS_TYPE        => TOTAL(),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,               # ramp15.geo is sized for this
        #---------------------------------------------------------------------------
        # Physical parameters / constants
        #---------------------------------------------------------------------------
        :energy_equation      => "energy",        # slot 4 is rho*E — see note (1)
        :lvisc                => true,
        #
        # Molecular viscosity — note (3).  Standard air Sutherland
        # constants; mu(125 K) = 8.656e-6 Pa.s, which with the paper's
        # p_inf, T_inf, u_inf and L gives Re_L = 4.22e5 against its 4.2e5.
        #
        :lsutherland          => true,
        :sutherland_muref     => 1.716e-5,        # Pa.s
        :sutherland_Tref      => 273.15,          # K
        :sutherland_S         => 110.4,           # K
        :Pr_lam               => 0.71,            # molecular Pr (Section 2.1)
        #---------------------------------------------------------------------------
        # (a) KINETIC-ENERGY / ENTROPY-CONSERVATIVE FLUX DIFFERENCING.
        #
        # With :lkep the inviscid volume term is not assembled from the
        # pointwise flux in user_flux.jl.  _expansion_inviscid_KEP! (rhs.jl)
        # builds it from SYMMETRIC TWO-POINT fluxes instead, and user_flux.jl
        # supplies the per-node auxiliary state (user_fluxaux!) and the
        # two-point fluxes themselves (flux_turbo) that it needs.  A deck
        # that sets :lkep without those methods dies with
        # `UndefVarError: flux_turbo not defined` inside the first RHS.
        #
        # WHY THIS AND NOT MORE VISCOSITY.  A collocation CG/SEM integrates
        # the nonlinear flux with the same LGL rule it interpolates on, so
        # the flux is ALIASED: energy goes into the grid-scale modes and CG
        # has nothing to take it back out.  At hypersonic Mach that is not a
        # cosmetic error, because
        #
        #     p = (gamma-1)(rhoE - rho|u|^2/2)
        #
        # is a difference of two nearly equal numbers: a relative error in
        # rhoE or rho u comes out of that subtraction amplified by
        #
        #     (gamma-1) rhoE/p = 1 + gamma(gamma-1) M^2/2
        #
        # which is 3.5 at M = 3 and 17.5 at M = 7.7.  The same aliasing
        # error buys five times the pressure error here, the pressure drives
        # the momentum flux, and the loop closes.  Flux differencing with an
        # EC two-point flux removes that transfer BY CONSTRUCTION, which is
        # the term this failure is made of -- artificial viscosity only
        # damps the symptom, and enough of it to hold a Mach-7.7 shock also
        # drowns the boundary layer this case exists to resolve.
        #
        # ranocha() is entropy conservative and is the default.
        # kennedy_gruber() is kinetic-energy preserving only.
        # central_euler() is the plain central flux written in
        # flux-differencing form: it is the CONTROL -- if the case behaves
        # the same under it as under ranocha(), the two-point machinery is
        # not what is helping and the diagnosis is somewhere else.
        #---------------------------------------------------------------------------
        :lkep                 => true,
        :volume_flux          => ranocha(),
        #---------------------------------------------------------------------------
        # (b) POSITIVITY REPAIR.  src/kernel/positivity/, applied to uaux
        # once per RHS call, after the filter and before the fluxes are
        # taken.  What it does, exactly:
        #
        #   rho < rho_min            -> rho := rho_min
        #   rhoE - ke < e_min        -> scale the MOMENTUM by
        #                               theta = sqrt(max((rhoE - e_min)/ke, 0))
        #                               leaving rhoE untouched, so energy is
        #                               conserved and only kinetic energy is
        #                               moved into internal energy.  Only if
        #                               ke is already ~0 does it instead
        #                               raise rhoE, which is the one branch
        #                               that injects energy.
        #
        # It is a REPAIR, not a preserving scheme.  The CG-native answer is
        # invariant-domain preservation with convex limiting
        # (Guermond-Popov-Tomas); this is the floor that keeps a single bad
        # node from ending a 2,000,000-step run, and its value is as much
        # diagnostic as numerical -- it prints how many node-visits it
        # repaired out of how many RHS calls, how much mass and energy that
        # injected, the smallest rho and p it saw, and where and when the
        # FIRST repair happened.  That report is MPI-collective (summed over
        # all ranks) and goes to stdout.
        #
        # Floors are 1e-6 of the free stream: rho_inf = 760/(287*125) =
        # 2.1185e-2 kg/m3 and p_inf = 760 Pa.  Both are far below anything
        # physical in this flow -- the coldest state in the domain is the
        # isothermal wall, rho = p_inf/(R*293) = 9.04e-3, and every shock
        # here RAISES p -- so the floors engage only on a state that has
        # already left the realizable set.
        #
        # IF THE REPORT SHOWS REPAIRS AWAY FROM THE FIRST FEW STARTUP STEPS,
        # THE RUN IS WRONG.  Do not raise the floors to make the message go
        # away; the message is the finding.
        #---------------------------------------------------------------------------
        :lpositivity          => true,
        :positivity_rho_min   => 2.1e-8,          # 1e-6 * rho_inf
        :positivity_p_min     => 7.6e-4,          # 1e-6 * p_inf
        :positivity_report    => true,
        # Every 1000 RHS calls the ranks reduce their counters and rank 0
        # prints if anything fired.  It is a collective, so every rank must
        # reach it: the trigger is the CALL COUNT, which is identical on all
        # ranks, never the repair count, which is not.
        :positivity_report_every => 1000,
        #
        # Shock capturing — note (2).
        #
        # :dsgs_sensor "legacy" is the sensor the total-energy DynSGS path
        # was validated with on ffs_step (the assembled RHS against a fixed
        # BDF2 of the stage state, in effect a |dq/dt| sensor); "residual",
        # the default, is the element-wise strong residual, DSGS.md 1.2.
        :visc_model           => DSGS(),
        :dsgs_sensor          => "legacy",
        # Per-equation multiplier on the DynSGS coefficient.  The method is
        # parameter-free, so 1.0 is the papers' own setting and it is what
        # this case starts from — unlike ffs_step, which needed x4 on the
        # momentum and energy slots, here those slots also carry the
        # physical mu and an over-large artificial part would show up
        # directly in the wall heat flux, the quantity being measured.
        #
        # Slot 1 is NOT zero: the total-energy DynSGS carries the density
        # diffusion beta*grad(rho) of eq. (3.3), and it is what keeps the
        # density jump across the leading-edge and separation shocks from
        # ringing.  Turning it off is the single most destabilising change
        # measured on ffs_step.
        #
        # If the shocks still ring, raise slots 2-4 together (ffs_step's
        # note applies: dissipation only helps when :Δt is cut to match),
        # and check the wall Stanton number afterwards — an artificial
        # conductivity that reaches the wall falsifies figure 3(c).
        :μ                    => [1.0, 1.0, 1.0, 1.0],
        # Cap on the DynSGS coefficient, nu_cap = Cmax*Delem*(|u|+c).  The
        # default 0.5 is Nazarov's, for a grid whose Delem IS the nodal
        # spacing; on an LGL grid at nop = 4 the nodal spacing is 0.17267
        # of Delem, so 0.5 overestimates the intended first-order-upwind
        # bound by 5.8x and the viscous CFL it implies is 2.9x tighter than
        # the advective one. 0.1 restores the intended magnitude and buys
        # 5x in :Δt. It lowers the CAP only -- the residual coefficient
        # itself is untouched, so this is a bound on how much artificial
        # dissipation the sensor may apply, not a change to the sensor.
        :dsgs_Cmax            => 0.1,
        # Made EXPLICIT here, not because the values change -- both are the
        # kernel defaults -- but because they are the two DynSGS knobs that
        # were got wrong during the Mach-7 debugging and a reader of this
        # deck should not have to go look them up.
        #
        # :dsgs_norms => "domain" normalises the residual by norms taken over
        # the WHOLE domain, so the viscosity does not depend on how many MPI
        # ranks the run uses.  "rank" makes the answer partition-dependent.
        #
        # :ldsgs_nodal => false is the elementwise DynSGS kernel.  The
        # nodal kernel reads the ASSEMBLED residual, which vanishes on an
        # under-resolved solution exactly as it does on a resolved one
        # (DSGS.md 1.2) -- so :ldsgs_nodal => true together with
        # :dsgs_sensor => "residual" is a BLIND sensor, and on
        # shock_circle_M7 that pairing was measured as a 17x regression
        # (231 steps against 3969).  Paired with "legacy", which is what
        # this deck runs, the nodal form is legitimate and it is the knob
        # worth trying next if this deck still needs more -- change it
        # ALONE, and only after (a) has been scored against the original.
        :dsgs_norms           => "domain",
        :ldsgs_nodal          => false,
        # Artificial Prandtl number P of eq. (3.7): kappa = P/(gamma-1)*mu.
        # Nazarov & Hoffman use P ~ 0.1.  This is NOT :Pr_lam above; the
        # two coefficients are added on the same slot and are separately
        # meaningful.
        :Pr                   => 0.1,
        #---------------------------------------------------------------------------
        # Mesh
        #
        # ramp15.msh is a two-block transfinite quad mesh of the ramp,
        # 269 x 60 elements, which at :nop => 4 is 1077 x 241 LGL points —
        # case G1 of Section 2.2 (1080 x 240), the coarser of the paper's
        # two grids in x and y.  Wall spacing 7.98e-6 m against its 8e-6 m.
        #
        # NO UPSTREAM STRIP.  Section 2.3 puts 20 points in 1 mm ahead of
        # the leading edge; this deck does not, and the 269 streamwise
        # elements are spent on the body instead (137 plate + 132 ramp).
        # The strip was the worst cell on the grid: conforming blocks force
        # it to carry the wall-clustered dy = 8e-6 m while its lower
        # boundary is a symmetry line, so it has no boundary layer and runs
        # the full 1726 m/s at the finest wall-normal spacing in the domain.
        # Measured: with the strip the case died at step 304 at exactly that
        # node. The cost is that the leading-edge singularity now sits on
        # the inflow plane, which is what Section 2.3 used the strip to
        # avoid; user_bc.jl gives that node to the wall.
        #
        # Regenerate with `python3 generate_mesh.py` (no gmsh needed) or
        # `gmsh -2 ramp15.geo -o ramp15.msh`; both files carry the same
        # numbers.  To go to case G2 (1600 x 320) set NX_PLATE/NX_RAMP/NY
        # in generate_mesh.py — but the mesh refiner below is the cheaper
        # way to add resolution uniformly.
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        # The meshes are NOT duplicated into this directory -- ramp15.msh and
        # ramp15_uniform.msh are 9.6 MB together and identical to the
        # originals, so this deck reads them from the case it was forked
        # from.  Regenerate them there (`python3 generate_mesh.py`).
        :gmsh_filename        => "./problems/CompEuler/rampCaoEtAl2021/ramp15_uniform.msh",
        #
        # DIAGNOSTIC GRID, no vertical stretching.  ramp15_uniform.msh is the
        # same geometry with a UNIFORM wall-normal distribution: 401 elements
        # over H, dy_wall = 2.58e-5 m, y+ = 1.00 (against 0.3 on G1).  Swap
        # the line above for
        #
        #   "./problems/CompEuler/rampCaoEtAl2021/ramp15_uniform.msh"
        #
        # It exists because the stretching is what breaks this case, and that
        # is measured, not suspected: a five-rung ladder from ffs_step showed
        # that ffs_step passes; this geometry and these Mach-7.7 conditions
        # pass on ffs_step's isotropic grid; this case passes on an
        # unstretched grid; and it fails with the free-slip wall, the
        # adiabatic wall and the isothermal wall alike on the stretched one.
        # Conditions, geometry and boundary conditions are therefore all
        # cleared, and dy_wall = 7.98e-6 m is the variable left standing.
        #
        # 269 x 401 = 107,869 elements, 6.7x the production grid, because a
        # uniform grid fine enough at the wall has to carry that spacing all
        # the way to H = 60 mm.  :Δt => 1.0e-9 runs it unchanged (the limits
        # relax by ~5x on this grid); up to ~4e-9 should hold.
        #---------------------------------------------------------------------------
        # Plotting
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        :output_dir           => "./output",
        :loutput_pert         => false,           # plot the total state
        # Numerical schlieren from rho, computed at output times only
        # (kernel/physics/schlieren.jl).  Adds two point-data fields to the
        # VTU on top of :outvars —
        #   schlieren_grad_rho  |grad rho| [kg/m^4], quantitative
        #   schlieren           exp(-k|grad rho|/max|grad rho|), the picture
        # Colour "schlieren" with a REVERSED greyscale in ParaView for the
        # familiar dark-shock look.  This is the field that makes the
        # leading-edge shock, the separation shock, the reattachment shock
        # and the shear layer over the bubble visible at once, and it is
        # the direct counterpart of the experimental schlieren of Roghelia
        # et al. (2017b) that the paper is validated against.
        :lschlieren           => true,
        :schlieren_k          => 20.0,            # contrast; Hadjadj uses 10-100
        #---------------------------------------------------------------------------
        # Refinement.
        #
        # OFF: ramp15.msh already IS the paper's coarse grid (case G1), so
        # a level of initial refinement puts this case between G1 and G2 in
        # x-y at four times the cost.  Set :linitial_refine => true with
        # :init_refine_lvl => 1 to do exactly that — and remember to halve
        # :Δt with it, the advective limit scales with the element size.
        #---------------------------------------------------------------------------
        :linitial_refine      => false,
        :init_refine_lvl      => 1,
        :ladapt               => false,
        :lxy_partition => false,
    ) #Dict

    return inputs

end
