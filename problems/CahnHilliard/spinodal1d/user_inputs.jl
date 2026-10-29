function user_inputs()

    #---------------------------------------------------------------------------
    # Cahn-Hilliard spinodal decomposition, with the chemical potential relaxed
    # onto its definition:
    #
    #     ∂ₜc   = ∂ₓ²μ
    #     τ ∂ₜμ = c³ - c - ε² ∂ₓ²c - μ
    #
    # which is ∂ₜc = ∂ₓ²(c³ - c - ε²∂ₓ²c) as τ → 0. ε, τ and the initial mode
    # are in ch_parameters() (user_flux.jl). A perturbation δ cos(4x) of the
    # mixture c = 0 grows at the Cahn-Hilliard rate σ = 8 until t ≈ 0.5, where
    # it is overlaid on the exact linear solution (user_analytic.jl), and then
    # saturates into a periodic array of c = ±1 domains.
    #---------------------------------------------------------------------------
    p    = ch_parameters()
    tend = 2.0

    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),
        :tinit                => 0.0,
        :tend                 => tend,
        #---------------------------------------------------------------------------
        # Time step. Explicit, so bounded by the fastest mode of the relaxed
        # system, which oscillates at εk²/√τ for the largest k the mesh
        # carries; Δt shrinks like √τ. (Plain explicit Cahn-Hilliard would
        # need Δt ~ 1/(ε²k⁴) instead.)
        #---------------------------------------------------------------------------
        :Δt                   => 2.5e-6,
        :diagnostics_at_times => [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0],
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        :lsource              => true,
        :lperiodic_1d         => true,
        #---------------------------------------------------------------------------
        # The second derivatives. The viscous operator adds
        # ∂ₓ(μ[ieq] ∂ₓ uprimitive[ieq]) to equation ieq, and user_primitives.jl
        # crosses the variables: ∂ₓ²μ into the c row, -ε²/τ ∂ₓ²c into the μ row.
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :μ                    => [1.0, -p.ε^2/p.τ],
        #---------------------------------------------------------------------------
        # Mesh: 64 elements at p = 4 on [0, 2π]. The interfaces that form
        # have width ≈ √2 ε = 0.25, against an element width of 0.098.
        #---------------------------------------------------------------------------
        :lread_gmsh           => false,
        :xmin                 => 0.0,
        :xmax                 => 2π,
        :nelx                 => 64,
        #---------------------------------------------------------------------------
        # Output
        #---------------------------------------------------------------------------
        :outformat            => "png",
        :loverwrite_output    => true,
        :output_dir           => "./output",
    ) #Dict

    return inputs
end
