function user_inputs()

    #---------------------------------------------------------------------------
    # 2D plane-strain cantilever — a real elastic beam, not a beam model.
    #
    #   ∂ₜ(ρu, ρv, σxx, σyy, σxy, uₓ, u_y)ᵀ + ∂ₓF + ∂_y G = (0,0,0,0,0, u, v)ᵀ
    #
    # Clamped on the x = 0 face, traction-free on the other three, struck at
    # t = 0 with a transverse velocity in the shape of the first bending mode
    # and then left alone. Rows 6-7 recover the displacement so the output can
    # be read as a bending beam (warp by (uₓ, u_y) in ParaView).
    #
    # Same material and same beam as the 1D Timoshenko case — ν = 0.3,
    # L/h = 10 — so the frequency that comes out can be read against it. In
    # plane strain the bending modulus is E/(1-ν²), which is why Ω₁ here is not
    # the Ω₁ of the 1D deck.
    #
    # Material, geometry and the mode shape are in user_flux.jl
    # (elastic_properties / cantilever_mode1 / cantilever_omega1) and read from
    # there, so the deck and the initial condition cannot drift apart.
    #---------------------------------------------------------------------------
    Tper = 2π/cantilever_omega1()      # ≈ 59.0528, Euler-Bernoulli estimate

    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),
        :tinit                => 0.0,
        :tend                 => Tper/2,                     # half a swing
        :Δt                   => 1.0e-3,
        :diagnostics_at_times => collect(range(0.0, Tper/2; length=9)),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # Source: the elastodynamic rows have none — every term of the system is
        # a divergence. S carries only the u̇ = v of the displacement recovery.
        #---------------------------------------------------------------------------
        :lsource              => true,
        #---------------------------------------------------------------------------
        # No artificial viscosity: linear system, smooth initial state, and the
        # motion this case is about is the one dissipation would damp first.
        #---------------------------------------------------------------------------
        :lvisc                => false,
        #---------------------------------------------------------------------------
        # Mesh: 40x4 quads over [0,1] x [-0.05,0.05] — square elements of side
        # 0.025, four of them (so 16 LGL points at p = 4) through the depth,
        # which is what resolves the bending stress profile.
        # Regenerate with problems/Elasticity/make_meshes.py.
        #
        # Stability is set by the FASTER family, cp = 1.1602, against the
        # smallest LGL spacing (≈ 0.0043): Δt = 1e-3 sits at CFL ≈ 0.27.
        #
        # RUNTIME: a slender beam's bending period is four orders of magnitude
        # longer than its elastic transit time, which is what makes even half a
        # swing ~30 000 steps. That ratio is physics, not a deck choice.
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/Elasticity/cantilever2d/beam_40x4.msh",
        #---------------------------------------------------------------------------
        # Output
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :output_dir           => "./output",
    ) #Dict

    return inputs
end
