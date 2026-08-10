function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # SWsphere — shallow water equations on a spherical shell.
        #
        # STATUS: GRID + INITIAL CONDITION.
        #
        # This deck builds the high-order (LGL) spectral-element grid on the
        # closed spherical shell read from :gmsh_filename, verifies it, builds
        # the Galewsky et al. (2004) barotropically unstable jet on it, writes
        # both to VTK, and STOPS. The equations and the metric terms on the
        # manifold are not written yet.
        #
        #   :lspherical_shell => true   use src/kernel/mesh/sphere_mesh.jl (a
        #                               2D manifold embedded in 3D: x, y AND z)
        #                               instead of the flat 2D gmsh reader in
        #                               mesh.jl, which keeps only (x,y) and
        #                               would collapse the sphere onto its
        #                               equatorial disc.
        #
        #   :linit_only       => true   <<<< stop after the grid AND the initial
        #                               condition are built and written. Set it
        #                               to false once the flux/source kernels
        #                               and the manifold metrics exist.
        #
        #   :lgrid_only       => true   stop one step earlier, right after the
        #                               grid, without touching initialize.jl.
        #                               Takes precedence over :linit_only.
        #---------------------------------------------------------------------------
        :lspherical_shell     => true,
        :lgrid_only           => false,
        :linit_only           => true,
        #---------------------------------------------------------------------------
        # Grid checks (see check_sphere_mesh in src/kernel/mesh/sphere_mesh.jl).
        #
        # A spherical shell is CLOSED: no boundary, every edge shared by exactly
        # two elements, and gmsh files built panel-by-panel routinely carry
        # duplicated nodes along the panel seams. These switches control what is
        # done about that.
        #
        #   :lmerge_coincident_nodes  fuse geometrically coincident vertices
        #                             before the topology is built (this is what
        #                             turns "six panels that look joined" into
        #                             one watertight shell)
        #   :node_merge_tol           merge tolerance RELATIVE to the radius
        #   :lcheck_grid              run the T1..T9 consistency tests
        #   :lstop_on_bad_grid        abort if any of them fails
        #   :lproject_to_sphere       snap every node radially onto the shell
        #---------------------------------------------------------------------------
        :lmerge_coincident_nodes => true,
        :node_merge_tol          => 1.0e-8,
        :lcheck_grid             => true,
        :lstop_on_bad_grid       => true,
        :lproject_to_sphere      => true,
        # Radius of the shell. Setting it explicitly RESCALES the grid onto a
        # sphere of that radius (a unit-sphere .msh is blown up to the Earth);
        # comment it out to take the radius from the gmsh file instead, as the
        # mean |x| over its nodes.
        #
        # Pinned here to the Galewsky Earth radius so the test is the published
        # one whatever the .msh happens to carry.
        :sphere_radius        => 6.37122e6,
        #---------------------------------------------------------------------------
        # Initial condition: Galewsky, Scott & Polvani (2004), Tellus 56A,
        # 429-440 — barotropically unstable mid-latitude jet.
        #
        # All of these are optional: initialize.jl falls back to the published
        # values shown. They are spelled out so a parameter sweep never needs to
        # touch the source.
        #
        #   the jet          u_max, and the band φ ∈ [φ0, φ1]
        #   the balance      Ω, g, and the global mean depth h_mean that fixes
        #                    the constant of integration h₀
        #   the trigger      ĥ cos φ exp[-(λ/α)²] exp[-((φ2-φ)/β)²]
        #   :galewsky_nquad  Simpson intervals for the balance integral. 512
        #                    gives ~5e-8 m on a 10 km depth; the answer is
        #                    unchanged to 10 digits from 128.
        #
        # Set :lgalewsky_perturbation => false for the BALANCED control run: an
        # exact steady solution, which is the natural first test of the scheme.
        #---------------------------------------------------------------------------
        :lgalewsky_perturbation => true,
        :galewsky_umax          => 80.0,
        :galewsky_phi0          => π/7,
        :galewsky_phi1          => π/2 - π/7,
        :galewsky_Omega         => 7.292e-5,
        :galewsky_gravity       => 9.80616,
        :galewsky_hmean         => 10000.0,
        :galewsky_hhat          => 120.0,
        :galewsky_alpha         => 1.0/3.0,
        :galewsky_beta          => 1.0/15.0,
        :galewsky_phi2          => π/4,
        :galewsky_nquad         => 512,
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        #---------------------------------------------------------------------------
        # Mesh parameters and files:
        #
        # A CLOSED all-quadrilateral surface grid of the sphere. Generate one
        # without needing gmsh:
        #
        #   julia tools/generate_cubed_sphere.jl 10 6.371e6 ./meshes/gmsh_grids/cubed_sphere.msh
        #
        # or with gmsh, from the .geo shipped next to this file:
        #
        #   gmsh -2 problems/ShallowWater/SWsphere/cubed_sphere.geo \
        #        -format msh2 -o ./meshes/gmsh_grids/cubed_sphere.msh
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./meshes/gmsh_grids/cubed_sphere.msh",
        #---------------------------------------------------------------------------
        # Time integration — placeholders, unused while :lgrid_only => true.
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK54(),
        :Δt                   => 1.0,
        :tinit                => 0.0,
        :tend                 => 1.0,
        :case                 => "swsphere",
        :SOL_VARS_TYPE        => TOTAL(),
        :lsource              => true,
        :lvisc                => false,
        #---------------------------------------------------------------------------
        # Plotting parameters
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        :output_dir           => "./output",
        :loutput_pert         => false,
        #---------------------------------------------------------------------------
    ) #Dict

    return inputs

end
