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
        :tend                 => 15.0,                        # ≈ 2 bending periods
        :Δt                   => 1.0e-3,
        :diagnostics_at_times => collect(range(0.0, 15.0; length=100)),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 5,
        #---------------------------------------------------------------------------
        # Source: the elastodynamic rows have none — every term is a divergence.
        # S carries only the u̇ = v of the displacement recovery.
        #---------------------------------------------------------------------------
        :lsource              => true,
        #---------------------------------------------------------------------------
        # STABILISATION. Not optional, and not cosmetic.
        #
        # A continuous-Galerkin spectral element discretisation of a first-order
        # hyperbolic system, run undamped for tens of thousands of steps with
        # strongly imposed boundary values, grows an element-scale checkerboard
        # mode. An earlier version of this case shipped without stabilisation
        # and reached σyy ≈ 1e+276 by its eighth output frame — about 2.6%
        # growth per step, out of round-off. That number is the target: whatever
        # is switched on here has to remove the top mode faster than that.
        #
        # (1) BOYD-VANDEVEN MODAL FILTER — the primary mechanism.
        #
        # :filter_type => "erf" selects the erf-log (Boyd-Vandeven) transfer
        # function of order 12 in init_filter (src/kernel/operators/filter.jl).
        # It leaves every mode up to 2p/3 EXACTLY alone and rolls off above it,
        # so at p = 4 the transfer weights are
        #
        #     mode   0      1      2      3        4
        #     w      1      1      1      0.9957   0
        #
        # i.e. it is a top-mode-only filter here: it deletes the one mode the
        # checkerboard lives in and does not touch the resolved solution at all.
        # That is exactly the right shape for this failure, and it is why the
        # filter is preferable to a Laplacian, which damps every wavenumber
        # including the beam's own bending mode.
        #
        # :mu_x / :mu_y are BLEND FACTORS in [0,1], one per reference direction:
        # init_filter builds F = μ·F_BV + (1-μ)·I. They are not viscosities and
        # they do not scale with Δx. The filter is applied once per RK STAGE,
        # so with CarpenterKennedy2N54 there are five applications per step:
        #
        #     μ       removal per time step, mode 3 / mode 4
        #     0.01     0.0%  /  4.9%      ~2x margin over the growth above
        #     0.02     0.0%  /  9.6%      ~4x
        #     0.05     0.1%  / 22.6%      ~9x  ← used here
        #     0.10     0.2%  / 41.0%      heavy-handed, no extra benefit
        #
        # (2) A RESIDUAL ARTIFICIAL VISCOSITY — a backstop, deliberately weak.
        #
        # Kept at a TENTH of what this deck used to carry (2e-4 → 2e-5). The
        # reason it is not simply zero: at p = 4 the Vandeven filter really only
        # bites on the single top mode, and the failure that motivated all this
        # has never been re-run, so which mode it actually lived in is not
        # established. At 2e-5 the cost to the physical bending mode is ≈0.4%
        # of its amplitude over the whole run — cheap insurance.
        #
        # ONCE A RUN CONFIRMS THE FILTER HOLDS ON ITS OWN, set :μ to zeros and
        # :lvisc => false. That is the preferred end state: ∇²σ has no physical
        # meaning, and the filter does the job without inventing one.
        #
        # :energy_equation is pinned because the shared 2D AV kernel carries one
        # equation-specific branch: at ieq == 4 it augments the flux with the
        # Navier-Stokes viscous-work term τ_ij·u_j, built from uprimitive[2] and
        # [3] on the assumption that those are velocities. Here ieq == 4 is σyy
        # and those slots hold ρvy and σxx, so the term would be nonsense. It is
        # gated off when :energy_equation == "theta", which is also the default;
        # pinned here so a stray setting elsewhere cannot switch it on.
        #---------------------------------------------------------------------------
        :lfilter              => true,
        :filter_type          => "erf",      # Boyd-Vandeven
        :mu_x                 => 0.05,
        :mu_y                 => 0.05,
        :lvisc                => true,
        :visc_model           => AV(),
        :energy_equation      => "theta",
        :μ                    => [2.0e-5, 2.0e-5, 2.0e-5, 2.0e-5, 2.0e-5, 0.0, 0.0],
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
