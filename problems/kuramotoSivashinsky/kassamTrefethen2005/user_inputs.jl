function user_inputs()

    #---------------------------------------------------------------------------
    # Kuramoto-Sivashinsky equation in conservation form,
    #
    #     ∂ₜu + ∂ₓ(u²/2) + ∂ₓ²u + ∂ₓ⁴u = 0,
    #
    # (u = ∂ₓφ of ∂ₜφ + ∂ₓ²φ + ∂ₓ⁴φ + ½(∂ₓφ)² = 0), on a periodic domain of
    # length 32π from u(x,0) = cos(x/16)(1 + sin(x/16)): the run of Kassam &
    # Trefethen (2005), Fig. 4.1, with their scheme (ETDRK4) and their step.
    #
    # du/dt = L u + N(u): N = -∂ₓ(u²/2) is the flux of user_flux.jl, and
    # L = -∂ₓ² - ∂ₓ⁴ the split linear operator below, built on the spectral
    # element mesh and integrated exactly (kernel/solvers/split_linear.jl).
    #---------------------------------------------------------------------------
    p    = ks_parameters()
    tend = 150.0

    inputs = Dict(
        :ode_solver           => ETDRK4(),
        :tinit                => 0.0,
        :tend                 => tend,
        #---------------------------------------------------------------------------
        # Time step: Kassam & Trefethen's h = 1/4. ETDRK4 treats L exactly, so
        # nothing ties Δt to the mesh's ∂ₓ⁴ eigenvalues (≈ 5500 here: an
        # explicit scheme would need Δt ≈ 2e-4). Every output time has to be
        # a whole number of steps (the φ operators are built for this one Δt).
        #---------------------------------------------------------------------------
        :Δt                   => 0.25,
        :diagnostics_at_times => collect(0.0:0.5:tend),   # rows of ks-spacetime.png
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        :lsource              => false,
        :lperiodic_1d         => true,
        #---------------------------------------------------------------------------
        # The linear part, L = c2 ∂ₓ² + c4 ∂ₓ⁴ = -∂ₓ² - ∂ₓ⁴ (one coefficient per
        # equation). No viscous operator: the derivatives are all in L.
        #---------------------------------------------------------------------------
        :lsplit_linear        => true,
        :split_d2             => [-1.0],
        :split_d4             => [-1.0],
        :lvisc                => false,
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
        # The case's figures (user_plot.jl): ks-it<n>.png and ks-spacetime.png.
        # false switches the plotting hook off altogether; "none" above
        # switches off all output.
        :ks_plot              => true,
        :ks_clims             => (-1.9, 1.0),   # saturated inferno range of ks-spacetime.png
        :loverwrite_output    => true,
        :output_dir           => "./output",
    ) #Dict

    return inputs
end
