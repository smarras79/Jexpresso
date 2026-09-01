function user_inputs()

    #---------------------------------------------------------------------------
    # 2D plane-strain beam on two hinged supports, loaded in the middle.
    #
    #   ∂ₜ(ρvx, ρvy, σxx, σyy, σxy, ux, uy)ᵀ + ∂ₓF + ∂_y G = (0,0,0,0,0, vx, vy)ᵀ
    #
    # Starts undeformed, unstressed and at rest. A downward pressure is raised
    # smoothly onto a patch of the top surface at midspan; the beam sags,
    # overshoots, and rings about the static deflection.
    #
    # WHAT TO EXPECT: midspan uy settles around -9.70e-3 (beam theory: -8.53e-3
    # of bending plus -1.17e-3 of shear), having first swung to roughly twice
    # that. σxx is the bending stress — tension on the bottom, compression on
    # the top, largest at midspan, ±O(1e-2). σyy stays small away from the
    # loaded patch. Everything is printed by initialize.jl at startup.
    #
    # Material, geometry and load live in user_flux.jl (elastic_properties /
    # beam_load_pressure / beam_static_midspan) and are read from there, so the
    # deck and the physics cannot drift apart.
    #---------------------------------------------------------------------------
    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),
        :tinit                => 0.0,
        :tend                 => 22.0,                        # ≈ 2 bending periods
        :Δt                   => 1.0e-3,
        :diagnostics_at_times => collect(range(0.0, 22.0; length=12)),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # Source: the elastodynamic rows have none — every term is a divergence.
        # S carries only the u̇ = v of the displacement recovery.
        #---------------------------------------------------------------------------
        :lsource              => true,
        #---------------------------------------------------------------------------
        # STABILISATION. This is not optional and it is not cosmetic.
        #
        # A continuous-Galerkin spectral element discretisation of a first-order
        # hyperbolic system, run undamped for tens of thousands of steps with
        # strongly imposed boundary values, grows an element-scale checkerboard
        # mode. Every other 2D case in this repository carries a viscosity for
        # the same reason (CompEuler/theta, AdvDiff/kopriva, ShallowWater/…).
        #
        # μ = 2e-4 damps the grid mode at ≈ 68 per unit time while costing the
        # physical bending mode ≈ 4% of its amplitude over the whole run. The
        # trade-off, at this resolution:
        #
        #     μ       grid mode      mode 1 loss over t = 22
        #     1e-4    34 /time       2%      may be marginal
        #     2e-4    68 /time       4%      the value used here
        #     1e-3   339 /time      19%      safe, visibly over-damped
        #
        # If the run still goes unstable, raise it; if the ringing damps out too
        # fast to see, lower it. Rows 6-7 (ux, uy) get zero: they are ODEs, and
        # diffusing a recovered displacement would be meaningless.
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :visc_model           => AV(),
        :μ                    => [2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 0.0, 0.0],
        #---------------------------------------------------------------------------
        # Mesh: 32x6 quads over [0,1] x [-0.1,0.1] — near-square elements of
        # 0.03125 x 0.03333, six of them (25 LGL points at p = 4) through the
        # depth, which is what resolves the bending stress profile.
        # Regenerate with problems/Elasticity/make_meshes.py.
        #
        # Stability is set by the FASTER family, cp = 1.1602, against the
        # smallest LGL spacing (≈ 0.0054): Δt = 1e-3 sits at CFL ≈ 0.22.
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/Elasticity/beam2d/beam_32x6.msh",
        #---------------------------------------------------------------------------
        # Output
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :output_dir           => "./output",
    ) #Dict

    return inputs
end
