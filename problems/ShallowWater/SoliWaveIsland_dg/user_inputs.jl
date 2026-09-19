function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # SoliWaveIsland_dg: 2D non-linear shallow water equations, solved with
        # the DISCONTINUOUS Galerkin discretization (:AD => DiscGal()).
        #
        # Same problem as problems/ShallowWater/SoliWaveIsland (Section 5.5 of
        # Marras, Kopera, Constantinescu, Suckale, Giraldo, "A residual-based
        # shock capturing scheme for the continuous/discontinuous spectral
        # element solution of the 2D shallow water equations", Advances in
        # Water Resources 114 (2018) 45-63): a right-propagating solitary wave
        # runs up and runs down a circular conical island in a closed basin.
        # Run the two side by side to compare the CG and DG solutions of the
        # same equations on the same mesh.
        #
        # What differs from the CG deck, and why:
        #
        #   :AD => DiscGal()        every element owns its own nodes; elements
        #                           are coupled only through the interface flux
        #                           (surface_rhs_el!, src/kernel/operators/).
        #   :numerical_flux         Rusanov (local Lax-Friedrichs). The shallow
        #                           water system is genuinely multi-speed and
        #                           non-linear, so upwind_flux() -- which in
        #                           this code is Rusanov with a single speed --
        #                           has no advantage here and the name would
        #                           misdescribe the scheme.
        #   :lvisc => false         REQUIRED, not a preference: the viscous
        #                           kernel _expansion_visc! has no DiscGal
        #                           method, so :lvisc => true stops the run
        #                           with a MethodError on the first RHS call.
        #                           The dissipation here is the Rusanov jump
        #                           term at the interfaces.
        #   :lexact_integration     collocated LGL. Over-integration is a
        #             => false      separate quadrature on the volume term and
        #                           is not wired through the DG surface term.
        #
        # Boundary conditions. All four walls are free-slip (user_bc.jl, as in
        # the CG case), but DG imposes them WEAKLY: build_dg_faces_2D! lists
        # the physical boundary faces, and surface_rhs_el! pairs the interior
        # trace with the mirror state built from user_bc_dirichlet! and sends
        # the pair through the same numerical flux as an interior face. No
        # nodal value is overwritten.
        #
        # Conservative variables:
        #   q = [H, Hu, Hv]
        # where H is the local water depth, Hb(x,y) is the (time-independent)
        # bathymetry, and the still-water depth far from the island is h0.
        #---------------------------------------------------------------------------
        # tend = 3.0, NOT the CG deck's 25.0. What stops it there is the
        # missing shock capturing, measured on this mesh:
        #
        #   t <= 3.00   max|u| holds at 0.33 m/s -- the solitary wave's own
        #               speed -- while the water away from it stays at rest.
        #   t  = 3.25   the wave front reaches the toe of the cone (x = 8.9)
        #               and max|u| jumps to 1.2, then 2.4 by t = 3.75.
        #   t  = 5.15   non-finite; the run aborts.
        #
        # This is the wet/dry front, not the DG operator. Three measurements
        # separate them: (1) a lake-at-rest run of this same deck holds
        # max|u| at 3e-11 m/s, so the interface flux, the wall flux and the
        # well-balanced source are exact on the equilibrium; (2) the same
        # deck with the cone removed runs 20 s clean, solitary wave and all,
        # including the reflection off the x = 25 wall at t ~ 12.7; (3) the
        # CG deck with :lvisc => false -- the same equations, same mesh, no
        # stabilization -- aborts at t = 4.85, EARLIER than this one.
        #
        # What carries the CG side through the run-up is a viscosity: the
        # constant mu = 0.05 of problems/ShallowWater/SoliWaveIsland, or,
        # properly, the residual-based shock capturing of the source paper
        # in problems/ShallowWater/SoliWaveIslandDSGS (:visc_model =>
        # DSGS_SW(), which runs the island case correctly).
        #
        # Neither reaches this deck, and both fail at the SAME single point.
        # Measured by running this case with the DSGS_SW settings of that
        # sibling deck: the sensor and the coefficient are computed fine
        # under DiscGal, and the run then stops on the first RHS call with
        #
        #   MethodError: no method matching _expansion_visc!(..., ::DSGS_SW,
        #                                    ::NSD_2D, ::DiscGal)
        #
        # DSGS supplies nu; the operator that applies it is the shared CG
        # viscous kernel, which has ContGal methods only. So the missing
        # piece is one thing, not two: a DiscGal viscous term (the element
        # volume part plus an interface term -- BR1/BR2 or interior penalty;
        # a bare element-local Laplacian would have no interface coupling
        # and would not be the operator it claims to be). Once that lands,
        # DSGS_SW drops in unchanged, and :tend can go back to the CG deck's
        # 25.0 so this case covers the run-up and run-down it is named for.
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK54(),
        :Δt                   => 0.01,
        :tinit                => 0.0,
        :tend                 => 3.0,
        :diagnostics_at_times => (0:0.25:3.0),
        :case                 => "soliwave_island_dg",
        :lsource              => true,
        :SOL_VARS_TYPE        => TOTAL(),
        #---------------------------------------------------------------------------
        # Discretization
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        :lexact_integration   => false,
        :AD                   => DiscGal(),
        :numerical_flux       => rusanov_flux(),
        #---------------------------------------------------------------------------
        # Physical parameters / artificial viscosity.
        #
        # Inviscid: see the :lvisc note in the header -- :lvisc => true is a
        # MethodError here, not a slower run. :μ is kept at zero so that a
        # later :lvisc => true does not silently inherit the CG deck's 0.05 --
        # an inert-looking coefficient that is not inert once the flag flips.
        #
        # NOTE: :lkep must stay unset/false under DiscGal -- the KEP volume
        # path (_expansion_inviscid_KEP!) dispatches on AD and has no DiscGal
        # method.
        #---------------------------------------------------------------------------
        :lvisc                => false,
        :ivisc_equations      => [1, 2, 3],
        :μ                    => [0.0, 0.0, 0.0],
        #---------------------------------------------------------------------------
        # Filter: OFF, for the reason given in the CG deck -- the filter acts
        # on the full depth H, which is cone-shaped at rest with a kink at the
        # wet/dry ring, so re-projecting it every step perturbs the
        # lake-at-rest equilibrium the well-balanced flux/source split
        # preserves.
        #---------------------------------------------------------------------------
        :lfilter              => false,
        #---------------------------------------------------------------------------
        # Mesh: 25 x 30 structured quadrilaterals on [0,25] x [-15,15], i.e.
        # Δx = Δy = 1 m; at nop = 4 that is 750 elements and 18750 DG degrees
        # of freedom (the CG mesh carries 12221 -- DG duplicates the interface
        # nodes). Kept next to this deck so the case is self-contained.
        # Regenerate with:
        #   gmsh -2 problems/ShallowWater/SoliWaveIsland_dg/SoliWaveIsland.geo \
        #        -o problems/ShallowWater/SoliWaveIsland_dg/SoliWaveIsland.msh
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/ShallowWater/SoliWaveIsland_dg/SoliWaveIsland.msh",
        #---------------------------------------------------------------------------
        # Plotting / output. VTK rather than the CG deck's PNG: the DG solution
        # is discontinuous at element interfaces, and ParaView draws it per
        # element, whereas the nodal PNG map interpolates across the jump and
        # would hide exactly what distinguishes this run from the CG one.
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :lwrite_initial       => true,
        :output_dir           => "./output",
        :loutput_pert         => false,
        #---------------------------------------------------------------------------
        # Mesh adaptivity is refused under DiscGal (mod_inputs.jl): the mortar
        # projections the CG path applies to the mass matrix and the RHS have
        # no meaning on duplicated DOFs.
        #---------------------------------------------------------------------------
        :linitial_refine     => false,
        :init_refine_lvl     => 1,
    ) #Dict

    return inputs

end
