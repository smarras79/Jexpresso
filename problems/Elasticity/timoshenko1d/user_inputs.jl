function user_inputs()

    #---------------------------------------------------------------------------
    # 1D Timoshenko beam — first flexural standing mode of a simply supported
    # (pinned-pinned) beam, released from rest at maximum deflection.
    #
    #     ∂ₜ(ρA v, ρI ω, γ, χ)ᵀ + ∂ₓ(-κₛGA γ, -EI χ, -v, -ω)ᵀ = (q, κₛGA γ, -ω, 0)ᵀ
    #
    # The exact solution is in user_analytic.jl and is overlaid on every output
    # frame, so what this deck produces is a convergence/verification picture
    # and not just a movie: after one full period every field must be back on
    # top of its initial profile.
    #
    # Beam properties and the mode itself are in user_flux.jl
    # (timoshenko_properties / timoshenko_mode); the period below is read from
    # there rather than hardcoded so the two can never drift apart.
    #---------------------------------------------------------------------------
    L    = timoshenko_properties().L    # 1.0
    Tper = 2π/timoshenko_mode().Ω       # ≈ 22.4214714536 for the shipped beam

    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),
        :tinit                => 0.0,
        :tend                 => Tper,                       # exactly one period
        :Δt                   => 2.0e-3,
        :diagnostics_at_times => collect(range(0.0, Tper; length=9)),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # Sources: the Timoshenko coupling terms κₛGA γ and -ω are
        # undifferentiated, so they live in S rather than in the flux.
        #---------------------------------------------------------------------------
        :lsource              => true,
        :lperiodic_1d         => false,
        #---------------------------------------------------------------------------
        # Physical parameters / constants
        #
        # No artificial viscosity: the system is linear, the mode is smooth and
        # fully resolved, and the point of the case is that the scheme carries
        # it without help. Any dissipation here would only hide phase error.
        #
        # AND NO FILTER IS AVAILABLE: `filter!` (src/kernel/operators/filter.jl)
        # has methods for NSD_2D and NSD_3D only, so setting :lfilter => true in
        # a 1D deck does not filter — it raises a MethodError on the first RHS
        # evaluation. The Boyd-Vandeven filter the 2D cases use (see
        # problems/Elasticity/beam2d/user_inputs.jl) has no 1D counterpart. If a
        # 1D case ever needs modal stabilisation, that gap has to be closed in
        # filter.jl first.
        #---------------------------------------------------------------------------
        :lvisc                => false,
        #---------------------------------------------------------------------------
        # Mesh.
        #
        # 32 elements at p = 4 over a unit beam. The stability limit is set by
        # the FASTER of the two wave families, √(E/ρ) = 1 here, against the
        # smallest LGL spacing (≈ 0.0054): Δt = 2e-3 sits at CFL ≈ 0.37.
        #---------------------------------------------------------------------------
        :lread_gmsh           => false,
        :xmin                 => 0.0,
        :xmax                 => L,
        :nelx                 => 32,
        #---------------------------------------------------------------------------
        # Output
        #---------------------------------------------------------------------------
        :outformat            => "png",
        :loverwrite_output    => true,
        :output_dir           => "./output",
    ) #Dict

    return inputs
end
