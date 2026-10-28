function user_inputs()

    #---------------------------------------------------------------------------
    # KdV soliton through the hyperbolic approximation KdVH (Biswas, Ketcheson,
    # Ranocha & Schütz, J. Sci. Comput. 103:90, 2025):
    #
    #     ∂ₜ(u, v, w)ᵀ + ∂ₓ(u²/2 + w, -v/τ, u/τ)ᵀ = (0, -w/τ, v/τ)ᵀ
    #
    # with v → ∂ₓu and w → ∂ₓ²u as τ → 0, so u solves ∂ₜu + u∂ₓu + ∂ₓ³u = 0.
    # τ, the soliton speed and its start position are in kdv_parameters()
    # (user_flux.jl). The exact KdV soliton and its first and second derivatives
    # are overlaid on every output frame (user_analytic.jl).
    #---------------------------------------------------------------------------
    p    = kdv_parameters()
    tend = 20.0                 # the soliton travels c·tend = 20 to x = +10

    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),
        :tinit                => 0.0,
        :tend                 => tend,
        #---------------------------------------------------------------------------
        # Time step. Explicit, so bounded by the FASTEST characteristic, the
        # -1/τ of the v row (1000 for τ = 1e-3), not by the soliton speed:
        # with the smallest LGL spacing ≈ 0.086 below, Δt = 2.5e-5 is a CFL of
        # ≈ 0.29. Halve Δt with τ.
        #---------------------------------------------------------------------------
        :Δt                   => 2.5e-5,
        :diagnostics_at_times => collect(range(0.0, tend; length=5)),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # The relaxation terms -w/τ and v/τ are undifferentiated: they are the
        # source.
        #---------------------------------------------------------------------------
        :lsource              => true,
        :lperiodic_1d         => true,
        #---------------------------------------------------------------------------
        # No viscosity and no filter (1D has none, see the note in
        # problems/Elasticity/timoshenko1d/user_inputs.jl): the soliton is
        # smooth and resolved, and dissipation would only hide the dispersion
        # error this case exists to measure.
        #---------------------------------------------------------------------------
        :lvisc                => false,
        #---------------------------------------------------------------------------
        # Mesh: 120 elements of width 0.5 at p = 4 on [-30, 30]; the soliton
        # (c = 1) has half-width ≈ 1.8, so ≈ 30 nodes across it.
        #---------------------------------------------------------------------------
        :lread_gmsh           => false,
        :xmin                 => -30.0,
        :xmax                 =>  30.0,
        :nelx                 => 120,
        #---------------------------------------------------------------------------
        # Output
        #---------------------------------------------------------------------------
        :outformat            => "png",
        :loverwrite_output    => true,
        :output_dir           => "./output",
    ) #Dict

    return inputs
end
