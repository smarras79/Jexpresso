function user_inputs()

    #---------------------------------------------------------------------------
    # Kuramoto-Sivashinsky equation in conservation form,
    #
    #     ∂ₜu + ∂ₓ(u²/2) + ∂ₓ²u + ∂ₓ⁴u = 0,
    #
    # (u = ∂ₓφ of ∂ₜφ + ∂ₓ²φ + ∂ₓ⁴φ + ½(∂ₓφ)² = 0), on a periodic domain of
    # length 32π from u(x,0) = cos(x/16)(1 + sin(x/16)): the run of Kassam &
    # Trefethen (2005), Fig. 4.1. The fourth derivative is reached through a
    # relaxed w ≈ ∂ₓ²u; see user_flux.jl, where L and τ are set.
    #---------------------------------------------------------------------------
    p    = ks_parameters()
    tend = 150.0

    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),
        :tinit                => 0.0,
        :tend                 => tend,
        #---------------------------------------------------------------------------
        # Time step. Explicit, so bounded by the fastest mode of the relaxed
        # system: Re s = -1/(2τ), |Im s| ≈ λ/√τ with λ ≈ 74 the largest
        # eigenvalue of -∂ₓ² on this mesh. That puts the limit near 3e-4 for
        # τ = 1e-4; Δt shrinks like √τ.
        #---------------------------------------------------------------------------
        :Δt                   => 2.0e-4,
        :diagnostics_at_times => collect(0.0:1.0:tend),
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
        # crosses the variables: -∂ₓ²w into the u row, ∂ₓ²u/τ into the w row.
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :μ                    => [-1.0, 1.0/p.τ],
        #---------------------------------------------------------------------------
        # Mesh: 64 elements at p = 4 on [0, 32π], 256 nodes. Kassam &
        # Trefethen used 128 Fourier modes for this run.
        #---------------------------------------------------------------------------
        :lread_gmsh           => false,
        :xmin                 => 0.0,
        :xmax                 => p.L,
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
