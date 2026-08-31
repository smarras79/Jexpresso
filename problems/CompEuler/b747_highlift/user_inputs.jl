function user_inputs()

    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),
        :tinit                => 0.0,
        :tend                 => 2.0e-4,          # 500k steps
        :lrestart             => false,
        :restart_time         => 0.0,
        # TIME STEP. Read this before changing :μ, the mesh, or Δt — they are
        # not independent, and getting it wrong does not look like a CFL
        # failure. It looks like "Instability detected" on step 3.
        #
        # The binding limit here is NOT the advective CFL. It is the DIFFUSIVE
        # limit of the DynSGS artificial viscosity, and it is 16x tighter:
        #
        #   advective   Δt ≲ ½·Δx/(|u|+c)              = 7.4e-8 s
        #   diffusive   Δt ≲ ½·Δx²/ν,  ν = μ⃗·C2·Δ·(|u|+c)
        #
        # with (all printed by the run itself, under "ELEMENT SIZES")
        #   Δx = smallest LGL node spacing  (the CFL length scale)
        #   Δ  = smallest element size      (what DynSGS caps on)
        #   |u| + c ≈ 1372 m/s,  C2 = 0.5 (SGS.jl),  μ⃗ = :μ below.
        #
        # ν is the DynSGS CAP μ_max = C2·Δ·(|u|+c), which the residual actually
        # reaches inside the elements that hold the bow shock — so the cap is
        # the operative value there, not a bound nobody attains.
        #
        # Note ν ∝ Δ ∝ Δx, so Δt_diffusive ∝ Δx: refining the nose costs time
        # step LINEARLY, and it is the nose that sets Δx. That is the whole
        # economics of this case.
        #
        # Measured, on an earlier grid (Δ = 1.6e-3, Δx = 2.0e-4) with
        # :μ => [1,4,4,4], for which the formula gives 4.7e-9:
        #     Δt = 3.0e-8  ->  NaN on step 3
        #     Δt = 5.0e-9  ->  ran clean
        # The formula is trustworthy; it is not a rule of thumb.
        #
        # FOUR ELEMENTS COST EIGHT ROUNDED TIPS, and the binding one is not a
        # trailing edge — it is a FLAP LEADING EDGE. The 64A210 nose radius is
        # 0.0056c, so a 0.15 m flap has a 0.84 mm nose, and an element that does
        # not resolve it gets badly distorted when the curving pulls its edge
        # onto the arc. Measured, by the Jacobian ratio the curving reports:
        #
        #     H_TIP = 0.0020  ->  min|J|/max|J| = 0.013, unstable at 6.0e-10
        #                         (i.e. 2.4x under the diffusive limit — the
        #                          time step was NOT the problem)
        #     H_TIP = 0.0012  ->  min|J|/max|J| = 0.211, runs
        #
        # So if this case misbehaves, look at that ratio in the run's output
        # before reaching for :Δt. On this grid the diffusive limit is 9.5e-10
        # and the value below sits 2.4x under it.
        :Δt                   => 4.0e-10,
        :diagnostics_at_times => (0:5.0e-6:2.0e-4),
        :lsource              => false,
        :SOL_VARS_TYPE        => TOTAL(),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,               # polynomial order
        #---------------------------------------------------------------------------
        # Physical parameters / constants
        #---------------------------------------------------------------------------
        :energy_equation      => "energy",        # slot 4 is ρE
        :lvisc                => true,
        :visc_model           => DSGS(),          # residual-based shock capturing
        # DynSGS per-equation multiplier — shock_circle's value, and do not
        # lower it here without checking what comes out.
        #
        # An earlier version of this deck used the unamplified [1, 1, 1, 1] to
        # buy a 4x larger time step. It ran, in the sense that it never tripped
        # the NaN check, but it produced UNPHYSICAL STATES: grid-scale
        # oscillations along the aft surface grew into the wake and drove the
        # pressure negative there — p < 0, hence T < 0 (min -2462 K) and an
        # imaginary sound speed, so Mach came out as 58 or undefined. A run can
        # be "stable" and still be worthless; the NaN check is not a physics
        # check. If you change this, look at min(p) and min(T) in the output,
        # not just at whether the run finished.
        #
        # The cost is real and linear: ν ∝ μ⃗, so this 4 divides the time step
        # by 4 against the diffusive limit below.
        :μ                    => [1.0, 4.0, 4.0, 4.0],
        # Artificial Prandtl number P of Nazarov & Hoffman eq. (3.7):
        # κ = P/(γ-1)·μ, with P ≈ 0.1.
        :Pr                   => 0.1,
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/CompEuler/b747_highlift/b747_highlift.msh",
        # EXACT GEOMETRY — the point of this case.
        #
        # naca64A210.msh is straight-sided: its vertices sit on the section, and
        # every high-order node in between sits on the CHORD. Around the nose,
        # where the radius of curvature is 0.0056c = 3.4e-3 m and a boundary
        # edge is 3.4e-3 m long, that puts the mid-edge nodes ~3.2e-4 m inside
        # the aerofoil — a ring of small forward-facing steps through the
        # stagnation region of a Mach-3 flow, and no value of :nop removes them.
        #
        # This key NAMES the boundary; ./user_exactGeo.jl says what shape it is.
        # It has to, because a NACA 6A-series section has no closed form: it is
        # the interpolating spline through a published table of ordinates, which
        # is exactly the kind of definition that cannot live in the kernel. The
        # kernel does the shape-independent half — put the boundary nodes on the
        # curve, blend the correction into the touching elements with the
        # Gordon-Hall transfinite map, check that the grid still conforms and
        # that nothing folded (src/kernel/mesh/exact_geometry.jl).
        #
        # The tuple is a SPEC: Jexpresso never looks inside it, it hands it back
        # to user_exactGeo. Here it is the placement of the section —
        #
        #     (:naca64A210, x_LE, y_LE, chord, incidence [deg], r_TE/chord)
        #
        # — and make_mesh.jl reads this very line to write the .geo, so the mesh
        # and the exact geometry cannot drift apart. Change it, re-run
        # make_mesh.jl, and both follow.
        #
        # α = 0 on purpose: the tunnel walls below are free-slip, so a free
        # stream at incidence would fight them. The 64A210 is cambered (design
        # c_l = 0.2), so the flow is asymmetric anyway.
        :exact_geometry       => Dict(
            # A four-element high-lift section, in deck order from nose to tail.
            # Each entry is one BODY, and each is its own gmsh Physical Curve —
            # the kernel already curves any number of named boundaries, so a
            # multi-element aerofoil needs no new machinery, only more tags.
            #
            #      (:naca64A210, x_LE,  y_LE,   chord, deflection [deg], r_TE/chord)
            #
            # Deflection is nose-up positive, so a drooped slat is negative and a
            # lowered flap is positive. r_TE is a FRACTION OF CHORD, and it is
            # chosen per element to give roughly the same ABSOLUTE radius (~3 mm,
            # 6 mm on the main): what has to be meshed is the arc in metres, not
            # a fraction of somebody's chord.
            "slat"  => (:naca64A210, 0.860, -0.085, 0.160, -25.0, 0.0188),
            "main"  => (:naca64A210, 1.000,  0.000, 0.600,   0.0, 0.0100),
            "flap1" => (:naca64A210, 1.470, -0.055, 0.240,  30.0, 0.0125),
            "flap2" => (:naca64A210, 1.690, -0.185, 0.150,  55.0, 0.0200)),
        #---------------------------------------------------------------------------
        # Plotting
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        :output_dir           => "./output",
        :loutput_pert         => false,           # plot the total state
        # Numerical schlieren from ρ, computed at output times only
        # (kernel/physics/schlieren.jl). Adds two point-data fields to the VTU
        # on top of :outvars —
        #   schlieren_grad_rho  |∇ρ| [kg/m⁴], quantitative
        #   schlieren           exp(-k|∇ρ|/max|∇ρ|), the picture
        # For the familiar dark-shock look, colour "schlieren" with a REVERSED
        # greyscale in ParaView. On this case it shows the detached bow shock at
        # the nose, the expansion over the upper surface, the recompression
        # shock off the trailing edge, and the reflections of all three off the
        # tunnel walls.
        :lschlieren           => true,
        :schlieren_k          => 20.0,            # contrast; Hadjadj uses 10-100
        #---------------------------------------------------------------------------
        # AMR off: make_mesh.jl already grades the grid onto the aerofoil.
        #
        # :linitial_refine => true DOES work here, unlike on a fitted geometry:
        # refinement puts each new boundary vertex on the midpoint of a chord, a
        # sagitta inside the section (3.2e-4 m at the nose), but the section is
        # STATED in user_exactGeo.jl rather than fitted from the grid, and the
        # snap moves vertices as well as high-order nodes — so the wall comes
        # back onto the aerofoil either way.
        #---------------------------------------------------------------------------
        :linitial_refine      => false,
        :ladapt               => false,
    ) #Dict

    return inputs

end
