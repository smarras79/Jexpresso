function user_inputs()

    #---------------------------------------------------------------------------
    # 1D Timoshenko beam — clamped-free cantilever, released from the static
    # shape a uniform load q₀ holds it in.
    #
    #   ∂ₜ(ρA v, ρI ω, γ, χ, w, φ)ᵀ
    #     + ∂ₓ(-κₛGA γ, -EI χ, -v, -ω, 0, 0)ᵀ = (q, κₛGA γ, -ω, 0, v, ω)ᵀ
    #
    # Rows 1-4 are the Timoshenko conservation laws; rows 5-6 recover the
    # displacement and rotation so the output is the deformed beam itself
    # (see user_flux.jl).
    #
    # Beam properties, the load, and the static shape are in user_flux.jl
    # (timoshenko_properties / timoshenko_static_cantilever).
    #---------------------------------------------------------------------------
    p = timoshenko_properties()

    #
    # Fundamental period, Euler-Bernoulli estimate: Ω₁ = (β₁L)² √(EI/ρA L⁴)
    # with β₁L = 1.8751… the first root of cos(βL)cosh(βL) + 1 = 0.
    #
    # It is an ESTIMATE on purpose — there is no closed form for the Timoshenko
    # cantilever, whose fundamental sits ≈0.3% lower at L/h = 10 (a reference
    # integration puts the half period at 31.06 against the 30.95 below). So at
    # :tend the beam is deliberately a little short of closing its cycle, and
    # that gap IS the shear-deformation correction this equation set exists to
    # carry. Compare against Elasticity/simply_supported, where the exact
    # Timoshenko frequency IS known and the case does close.
    #
    β1L = 1.8751040687119611
    Ω1  = (β1L/p.L)^2*sqrt(p.EI/p.ρA)
    T1  = 2π/Ω1                          # ≈ 61.90 for the shipped beam

    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),
        :tinit                => 0.0,
        :tend                 => T1,
        :Δt                   => 2.0e-3,
        :diagnostics_at_times => collect(range(0.0, T1; length=9)),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # Sources: the Timoshenko coupling terms κₛGA γ and -ω are
        # undifferentiated, and the whole of the w/φ recovery is source.
        #---------------------------------------------------------------------------
        :lsource              => true,
        :lperiodic_1d         => false,
        #---------------------------------------------------------------------------
        # Physical parameters / constants
        #
        # No artificial viscosity: the system is linear and the released static
        # shape is smooth and fully resolved. Dissipation here would damp the
        # very oscillation the case is about.
        #---------------------------------------------------------------------------
        :lvisc                => false,
        #---------------------------------------------------------------------------
        # Mesh.
        #
        # 32 elements at p = 4 over a unit beam. The stability limit is set by
        # the FASTER of the two wave families, √(E/ρ) = 1 here, against the
        # smallest LGL spacing (≈ 0.0054): Δt = 2e-3 sits at CFL ≈ 0.37. The
        # flexural period is four orders of magnitude longer than that limit,
        # which is what makes this run ~31 000 steps.
        #---------------------------------------------------------------------------
        :lread_gmsh           => false,
        :xmin                 => 0.0,
        :xmax                 => p.L,
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
