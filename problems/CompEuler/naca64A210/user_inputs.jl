function user_inputs()

    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),
        :tinit                => 0.0,
        :tend                 => 6.0e-3,          # ≈ 2 tunnel flow-throughs
        :lrestart             => false,
        :restart_time         => 0.0,
        # CFL. The smallest element is at the nose (h = 2.3e-3 m, measured by
        # ./check_mesh.jl); at :nop => 4 the tightest LGL spacing inside it is
        # 0.173h ≈ 4.1e-4 m, and the fastest wave is |u| + c ≈ 1372 m/s, so the
        # explicit stability limit is ~3.0e-7 s. This runs at a tenth of it,
        # which is what the DynSGS residual — evaluated on the previous step —
        # wants across a shock. 200k steps.
        :Δt                   => 3.0e-8,
        :diagnostics_at_times => (0:2.0e-5:6.0e-3),
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
        :μ                    => [1.0, 4.0, 4.0, 4.0],
        # Artificial Prandtl number P of Nazarov & Hoffman eq. (3.7):
        # κ = P/(γ-1)·μ, with P ≈ 0.1.
        :Pr                   => 0.1,
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/CompEuler/naca64A210/naca64A210.msh",
        # EXACT GEOMETRY — the point of this case.
        #
        # naca64A210.msh is straight-sided: its vertices sit on the section, and
        # every high-order node in between sits on the CHORD. Around the nose,
        # where the radius of curvature is 0.0056c = 3.4e-3 m and a boundary
        # edge is 2.3e-3 m long, that puts the mid-edge nodes ~1.7e-4 m inside
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
        #     (:naca64A210, x_LE, y_LE, chord, incidence in degrees)
        #
        # — and make_mesh.jl reads this very line to write the .geo, so the mesh
        # and the exact geometry cannot drift apart. Change it, re-run
        # make_mesh.jl, and both follow.
        #
        # α = 0 on purpose: the tunnel walls below are free-slip, so a free
        # stream at incidence would fight them. The 64A210 is cambered (design
        # c_l = 0.2), so the flow is asymmetric anyway.
        :exact_geometry       => Dict("airfoil" => (:naca64A210, 1.0, 0.0, 0.6, 0.0)),
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
        # sagitta inside the section (1.7e-4 m at the nose), but the section is
        # STATED in user_exactGeo.jl rather than fitted from the grid, and the
        # snap moves vertices as well as high-order nodes — so the wall comes
        # back onto the aerofoil either way.
        #---------------------------------------------------------------------------
        :linitial_refine      => false,
        :ladapt               => false,
    ) #Dict

    return inputs

end
