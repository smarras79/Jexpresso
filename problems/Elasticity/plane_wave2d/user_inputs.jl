function user_inputs()

    #---------------------------------------------------------------------------
    # 2D linear elastodynamics — oblique plane P-wave on a doubly periodic
    # square. The verification case for the equation set.
    #
    #   ∂ₜ(ρu, ρv, σxx, σyy, σxy)ᵀ + ∂ₓF + ∂_y G = 0
    #
    # The wave travels along n = (1,1)/√2, so every component of the state and
    # both flux directions are active — this is a real 2D test, not a 1D
    # solution embedded in a 2D mesh. Because the domain is periodic in both
    # directions there are no boundaries at all: what is being checked is the
    # flux operator, the metrics, and the periodic node matching, against an
    # exact solution.
    #
    # THE CHECK: the wave is exactly periodic in time with period
    # T = 2π/(cp|k|). :diagnostics_at_times below lands on whole multiples of
    # T, so every output frame must reproduce the t = 0 frame. Ten periods in,
    # whatever has drifted is dispersion and time-integration error and
    # nothing else.
    #
    # Material and wave are defined in user_flux.jl (elastic_properties /
    # elastic_plane_wave) and read from there, so the deck and the initial
    # condition cannot drift apart.
    #---------------------------------------------------------------------------
    W    = elastic_plane_wave()
    Tper = 2π/W.ω                       # ≈ 0.6094494002

    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),
        :tinit                => 0.0,
        :tend                 => 10*Tper,
        :Δt                   => 2.0e-3,
        :diagnostics_at_times => collect(range(0.0, 10*Tper; length=11)),  # one per period
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # No source: every term of 2D elastodynamics is a divergence.
        #---------------------------------------------------------------------------
        :lsource              => false,
        #---------------------------------------------------------------------------
        # DELIBERATELY UNSTABILISED — no viscosity and no filter.
        #
        # This is the verification case, and both would corrupt what it
        # measures. The check is that each output frame reproduces the t = 0
        # frame exactly; a Boyd-Vandeven filter deletes the top modal
        # coefficient every RK stage, so any of the wave living in that mode
        # would be removed and read as amplitude loss that is not the scheme's.
        #
        # It is also the DIAGNOSTIC. The element-scale instability that killed
        # an earlier bounded case in this directory grew at ≈2.6% per step; over
        # this run's ~3050 steps that would take round-off to ≈1e18. So if this
        # case blows up too, the problem is in the volume operator or the time
        # step and not in the boundary treatment — there are no boundaries here.
        # If it stays clean for ten periods, the volume operator is sound and
        # the boundaries are where to look.
        #
        # Once that question is settled, the filter can be switched on here to
        # match beam2d — see the long note in that deck for what the parameters
        # mean:
        #
        #     :lfilter     => true,
        #     :filter_type => "erf",
        #     :mu_x        => 0.05,
        #     :mu_y        => 0.05,
        #---------------------------------------------------------------------------
        :lvisc                => false,
        :lfilter              => false,
        #---------------------------------------------------------------------------
        # Mesh: 16x16 quads on the unit square, both directions periodic.
        # Regenerate with problems/Elasticity/make_meshes.py.
        #
        # Stability is set by the FASTER family, cp = 1.1602, against the
        # smallest LGL spacing (≈ 0.0108): Δt = 2e-3 sits at CFL ≈ 0.21.
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/Elasticity/plane_wave2d/square_periodic_16x16.msh",
        #---------------------------------------------------------------------------
        # Output
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :output_dir           => "./output",
    ) #Dict

    return inputs
end
