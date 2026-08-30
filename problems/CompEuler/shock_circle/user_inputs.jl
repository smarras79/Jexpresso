function user_inputs()

    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),
        :tinit                => 0.0,
        :tend                 => 2.0e-2,          # ≈ 2.7 tunnel flow-throughs
        :lrestart             => false,
        :restart_time         => 0.0,
        :Δt                   => 0.5e-7,
        :diagnostics_at_times => (0:5.0e-5:2.0e-2),
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
        :energy_equation      => "energy",        # slot 4 is ρE — see note (1)
        :lvisc                => true,
        :visc_model           => DSGS(),          # residual-based shock capturing     
        :μ                    => [1.0, 4.0, 4.0, 4.0],
        # Artificial Prandtl number P of eq. (3.7): κ = P/(γ-1)·μ. Nazarov &
        # Hoffman use P ≈ 0.1.
        :Pr                   => 0.1,
        #:ldsgs_global_norms   => true,
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./meshes/gmsh_grids/ffs_step_transfinite.msh",
        :gmsh_filename        => "./problems/CompEuler/shock_circle/plate_hole_circle_unit.msh",
        # EXACT GEOMETRY. plate_hole_circle_unit.geo defines the hole as a gmsh
        # `Circle`, but the .msh is straight-sided: the 16 boundary vertices sit
        # on the circle of radius 0.2 at (1, 0) and every high-order node in
        # between sits on the CHORD, up to 3.8e-3 inside the wall. That polygon
        # is what the solver would otherwise see however large :nop is, and its
        # sixteen corners are sixteen places for a slip wall to shed spurious
        # vorticity into a Mach-3 flow.
        #
        # This snaps the boundary nodes onto the circle and blends the
        # correction into the touching elements (Gordon-Hall, see
        # src/kernel/mesh/exact_geometry.jl). The wall becomes the degree-:nop
        # isoparametric interpolant of the true circle — 4.4e-8 from it at
        # :nop => 4 instead of 3.8e-3, a factor of 10^5.
        #
        # THIS KEY ONLY NAMES THE BOUNDARY. What shape it is — the analytic
        # definition, and the fit — is stated in ./user_exactGeo.jl, which is
        # where you put an ellipse, an aerofoil or a spline instead. The value
        # below is a SPEC: Jexpresso does not look inside it, it hands it back to
        # user_exactGeo, so the numbers stay next to the mesh they describe.
        #
        # THE CIRCLE IS STATED, NOT FITTED, because :linitial_refine is on below.
        # `:circle` on its own asks user_exactGeo_setup to fit the centre and
        # radius from the boundary vertices the mesh file carries, and on the
        # UNREFINED grid that recovers (1.0, 0.0) and 0.2 exactly. Refinement
        # breaks it: each new boundary vertex lands on the midpoint of a chord, a
        # sagitta inside the circle, so the fit sees half its points 3.8e-3 in,
        # returns r = 0.19809, and is refused — correctly, but the net effect is
        # that nothing gets curved at all. The values below are the ones in
        # plate_hole_circle_unit.geo (xc, yc, radius).
        :exact_geometry       => Dict("circle_boundary" => (:circle, 1.0, 0.0, 0.2)),
        #---------------------------------------------------------------------------
        # Plotting
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        #:output_dir           => "/scratch/smarras/smarras/output/shock_circle/",
        :output_dir           => "./output",
        :loutput_pert         => false,           # plot the total state
        # Numerical schlieren from ρ, computed at output times only
        # (kernel/physics/schlieren.jl). Adds two point-data fields to the
        # VTU on top of :outvars —
        #   schlieren_grad_rho  |∇ρ| [kg/m⁴], quantitative
        #   schlieren           exp(-k|∇ρ|/max|∇ρ|), the picture
        # For the familiar dark-shock look, colour "schlieren" with a
        # REVERSED greyscale in ParaView. This is the field to look at on
        # this case: the bow shock, its reflection off the roof, the Mach
        # stem and the slip line downstream of the triple point are all
        # density features, and the exponential map keeps the weak ones
        # visible next to the strong bow shock.
        :lschlieren           => true,
        :schlieren_k          => 20.0,            # contrast; Hadjadj uses 10-100
        #---------------------------------------------------------------------------
        # AMR off: the mesh already resolves the shocks at h/nop = 1/80, and
        # DynSGS is what handles what is left under-resolved.
        #---------------------------------------------------------------------------
        :linitial_refine      => true,
        :init_refine_lvl      => 1,
        :ladapt               => false,
    ) #Dict

    return inputs

end
