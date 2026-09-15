function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        #
        # 2D CompEuler + Navier-Stokes: Mach-7.7 hypersonic flow over a
        # 15-degree compression ramp with laminar separation.
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
        # TIME STEP.  Explicit, on a hypersonic viscous grid, so it is small
        # and there is no way around it.  At :nop => 4 the tightest LGL
        # interval of an element is 0.17267 h, and the advective limit
        # dx/(|u|+c) over the grid is
        #
        #   where                element h   tightest dx  |u|+c     dt
        #   x, inflow strip      2.0e-4 m    3.45e-5 m    1950 m/s  1.8e-8
        #   x, leading edge      4.88e-4 m   8.43e-5 m    1950 m/s  4.3e-8
        #   y, first element     4.62e-5 m   7.98e-6 m     343 m/s  2.3e-8
        #   y, mid-boundary-layer 1.48e-4 m  2.56e-5 m    2100 m/s  1.2e-8
        #   y, top element       4.53e-3 m   7.82e-4 m    1950 m/s  4.0e-7
        #
        # The binding one is NOT the wall: there u -> 0 and the gas is at
        # the 293 K wall temperature, so the wave speed is only 343 m/s.
        # It is the middle of the boundary layer, where the mesh has grown
        # to 2.6e-5 m but the gas is at 1200 K and 1400 m/s -- 1.2e-8 s.
        #
        # The viscous limit rho dy^2/mu at the wall is ~1e-7 s and does NOT
        # bind, which is the opposite of ffs_step: there the viscosity was
        # the artificial one, saturating at the step corner; here it is the
        # physical one and the wall is where the density is highest.
        #
        # 5.0e-9 s is CFL ~ 0.4 against that 1.2e-8 s, which is what the
        # impulsive start needs.  1.0e-8 s is usually fine once the
        # leading-edge transient has washed out (the first ~2000 steps).
        #
        # WHERE THE REMAINING dt IS, if you want it.  The binding node is
        # not on the ramp at all: it is at x = -1 mm, y = 2.3e-5 m, in the
        # free-stream strip ahead of the leading edge.  That strip carries
        # the wall-clustered y-spacing (the blocks are conforming, so it
        # must) but has NO boundary layer -- its lower boundary is a
        # symmetry line -- so it runs the full 1726 m/s at the finest
        # wall-normal spacing on the grid, resolving nothing.  Cutting
        # NX_UP in generate_mesh.py from 5 to 3 moves the limit from
        # 1.27e-8 to 1.57e-8 s, +24%, and then it saturates: the next
        # bottleneck is the leading-edge region at x = +1 mm.  The cost is
        # 20 -> 12 points ahead of the leading edge, which is a deviation
        # from Section 2.3.  Left at 5 here; it is the only free dt in the
        # setup and it is worth exactly 24%.
        :Δt                   => 5.0e-9,
        :diagnostics_at_times => (0:2.5e-5:2.0e-3),
        # Wall-clock note, not a setting: 2.0e-3 s at 5e-9 is 400,000 steps
        # on 16,140 elements.  Run it on several ranks and expect hours,
        # not minutes; a long silence between the CFL/VTK lines is the run
        # working, not a hang.  JEXPRESSO_STEP_HEARTBEAT=1 turns on a
        # per-step trace without editing this deck.
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
        # Artificial Prandtl number P of eq. (3.7): kappa = P/(gamma-1)*mu.
        # Nazarov & Hoffman use P ~ 0.1.  This is NOT :Pr_lam above; the
        # two coefficients are added on the same slot and are separately
        # meaningful.
        :Pr                   => 0.1,
        #---------------------------------------------------------------------------
        # Mesh
        #
        # ramp15.msh is a three-block transfinite quad mesh of the ramp,
        # 269 x 60 elements, which at :nop => 4 is 1077 x 241 LGL points —
        # case G1 of Section 2.2 (1080 x 240), the coarser of the paper's
        # two grids in x and y.  Wall spacing 7.98e-6 m against its 8e-6 m.
        #
        # Regenerate with `python3 generate_mesh.py` (no gmsh needed) or
        # `gmsh -2 ramp15.geo -o ramp15.msh`; both files carry the same
        # numbers.  To go to case G2 (1600 x 320) set NX_PLATE/NX_RAMP/NY
        # in generate_mesh.py — but the mesh refiner below is the cheaper
        # way to add resolution uniformly.
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/CompEuler/rampCaoEtAl2021/ramp15.msh",
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
    ) #Dict

    return inputs

end
