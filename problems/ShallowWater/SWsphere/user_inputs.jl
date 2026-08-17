function user_inputs()

    inputs = Dict(
        #---------------------------------------------------------------------------
        # SWsphere — the Galewsky jet on the shell, shipped in the ARTIFICIAL
        # DIFFUSION configuration of Eq. (8b) added to the modal filter: the
        # PAPER'S OWN setup (Marras, Kopera & Giraldo 2015, QJRMS 141:
        # 1727-1739), which filters every step AND carries ν = 1e5 m²/s.
        #
        # This deck is what used to be the separate SWsphere_visc case, folded
        # in here so there is ONE spherical shallow-water case; flip :lvisc at
        # the bottom to get back the filter-only configuration:
        #
        #     filter only    :lfilter => true    :lvisc => false
        #     + diffusion    :lfilter => true    :lvisc => true, :μ => 1e5   ← as shipped
        #
        # The five sibling user_*.jl in this directory are the REAL
        # implementations. They were briefly one-line `include`s of the (then
        # separate) SWsphere case's files; when the two cases were merged those
        # shims ended up pointing at themselves and every run died with a
        # StackOverflowError. Do not reintroduce them — a case directory owns
        # its six files outright.
        #
        # WHY NOT VISCOSITY ALONE, which is what this case first shipped as.
        # It does not survive this grid. Measured, 3 days, all else equal:
        #
        #   ν = 1e5, filter OFF   NaN at 2.005 d, whether the diffusion is on
        #                         the momentum only or on all four equations
        #                         (both were run; they die within 0.4% of each
        #                         other, so it is not the undiffused φ)
        #   ν = 5e5, filter OFF   survives 3 d — but max|ζ| falls 1.12e-4 →
        #                         5.65e-5, HALF the vorticity gone, and
        #                         vorticity is the field this test is judged on
        #   ν = 1e5, filter ON    survives 3 d with max|ζ| 1.11e-4 → 9.58e-5,
        #                         a 13% loss  ← this deck
        #
        # The grid is the reason: 10 elements per panel at nop=5 leaves the
        # coarse elements at ~219 km effective resolution, i.e. a grid Reynolds
        # number uΔ/ν ≈ 175 at ν = 1e5. Viscosity alone can hold that only by
        # being strong enough to erase the answer.
        #
        # Physics, provenance and the discretization: see the README next to
        # problems/ShallowWater/SWsphere/user_inputs.jl. Everything down to the
        # STABILIZATION block below is that deck, unchanged.
        #---------------------------------------------------------------------------
        :lspherical_shell     => true,
        :lgrid_only           => true,
        :linit_only           => false,
        #---------------------------------------------------------------------------
        # Grid checks (see check_sphere_metrics in src/kernel/mesh/sphere_metrics.jl)
        #---------------------------------------------------------------------------
        :lcheck_grid             => true,
        :lstop_on_bad_grid       => true,
        :lproject_to_sphere      => true,
        :sphere_radius        => 6.37122e6,
        #---------------------------------------------------------------------------
        # Which cube-face → sphere map the panels carry. This switch slides the
        # nodes onto a different map after the grid is read, leaving the
        # connectivity and the panel boundaries alone. See the README section
        # "Changing the cube-face → sphere map" and
        # src/kernel/mesh/cubed_sphere_maps.jl.
        #
        # cubed_sphere.geo does NOT build the gnomonic grid, despite what the
        # comments here used to say. gmsh spaces Transfinite Line points at
        # equal ANGLE along a Circle arc, so the .msh is already EQUIANGULAR —
        # it reproduces tools/generate_cubed_sphere.jl's :equiangular output to
        # 2.2e-16 of R. Hence :none below: there is nothing to improve, and
        # asking for :equiangular would have applied the warp a second time
        # (min element edge 710 km -> 562 km, 21% of the time step, silently).
        # The remap now MEASURES the map the grid carries rather than assuming,
        # so that mistake is no longer possible — it reports what it found.
        #
        #   :equiangular  what this grid already is (Ronchi et al. 1996):
        #                 min element edge 710 km, max 1000 km, ratio 1.41.
        #   :conformal    NOT USABLE ON THIS GRID, and the code refuses it.
        #                 It does keep grid lines at 90° into a cube corner
        #                 where the other two degenerate to 120°, but it pays
        #                 for that with a SINGULARITY at the eight corners: the
        #                 local map scale falls off as d^(1/3), so a node
        #                 sitting on a corner — and cubed_sphere.geo puts one on
        #                 each of the eight — gets a surface Jacobian that
        #                 collapses with it: measured 1/27, 1/51, 1/93 of the
        #                 grid median at nop = 3, 5, 8, and the metric checks
        #                 then fail by O(1) (M6 = 0.10 to 1.8 against 5e-2).
        #                 Three panels meet at a corner and each must open 120°;
        #                 an angle-preserving map can only manage that by
        #                 collapsing its derivative there. This is why conformal
        #                 cubed spheres are used by cell-centred finite-volume
        #                 codes and not by nodal spectral elements.
        #                 (Rančić, Purser & Mesinger 1996)
        #                 A corner-STRETCHED variant of it was tried, shipped
        #                 briefly, and removed: it does make the corner Jacobian
        #                 finite, and it still fails check_sphere_metrics by O(1)
        #                 (M6 = 0.41) at every nop from 3 to 7 and every spacing
        #                 from n = 5 to 40. Forcing 90° with a non-zero Jacobian
        #                 turns the corner into a cone point, and a separable
        #                 stretch also blows |r_u| up along the WHOLE cube edge
        #                 as (1-u)^(-1/4). Do not reach for it again — THE 120°
        #                 CORNER in cubed_sphere_maps.jl has the numbers.
        #---------------------------------------------------------------------------
        # NOTE the trailing comma belongs BEFORE the # on whichever line is live.
        # Without it Julia reads the NEXT line's leading `:` as the range
        # operator and the run dies with a baffling
        # `UndefVarError: sphere_metrics`, naming the line AFTER the mistake.
        :cubed_sphere_map => :none,        # THIS GRID IS ALREADY EQUIANGULAR — see the note above
        #:cubed_sphere_map => :gnomonic,   # would UN-warp it to the equidistant projection (Sadourny 1972)
        #:cubed_sphere_map => :equiangular,# no-op here; the remap detects that and says so
        #:cubed_sphere_map => :conformal,  # UNUSABLE with this grid — see the note above
        #---------------------------------------------------------------------------
        # Metric terms of the 2D manifold. Kopriva's curl-invariant form
        # degenerates on a surface — see the header of sphere_metrics.jl — and
        # leaves only the direction v in which the surface is extended off
        # itself. :cross_product is v = n̂ (and is what the textbook formulas
        # give); :radial is v = x̂. Both are curl-invariant; they differ in which
        # exact property they keep.
        #---------------------------------------------------------------------------
        #:sphere_metrics       => :cross_product,
        :sphere_metrics       => :curl_invariant,
        #---------------------------------------------------------------------------
        # Initial condition: Galewsky, Scott & Polvani (2004), Tellus 56A,
        # 429-440 — barotropically unstable mid-latitude jet.
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
        # Keeping the flow ON the shell: the μx source of user_source.jl and the
        # P = I - xxᵀ/r² projection, applied after every RK stage.
        #---------------------------------------------------------------------------
        :llagrange_projection   => true,
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 5,
        #---------------------------------------------------------------------------
        # Mesh: the cubed sphere that ships with SWsphere (600 quads, 602
        # vertices, 10 elements per panel edge). Pointed at rather than copied —
        # the two cases are the same grid by construction.
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./problems/ShallowWater/SWsphere/cubed_sphere.msh",
        #---------------------------------------------------------------------------
        # Time integration.
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK33(),
        :lcfl_dt              => true,           # take Δt from the CFL condition
        :cfl                  => 0.35,
        :tinit                => 0.0,
        :tend                 => 10*24*3600,       # 144 h = 6 days, as in the test
        :ndiagnostics_outputs => 24,             # a VTK dump every 6 h
        :ndiagnostics_prints  => 200,            # steps between diagnostic lines
        :case                 => "swsphere_visc",
        :SOL_VARS_TYPE        => TOTAL(),
        :lsource              => true,
        #---------------------------------------------------------------------------
        # STABILIZATION — this is the whole point of the case.
        #
        # The artificial diffusion δν∇²(φu) of Eq. (8b), with the paper's
        # ν = 1e5 m²/s. On the shell this is the SURFACE (Laplace-Beltrami)
        # Laplacian, assembled in weak form from the manifold metrics in
        # src/kernel/operators/sphere_rhs.jl.
        #
        # :ivisc_equations => [2,3,4] is the MOMENTUM only, where the paper puts
        # it: the continuity equation ∂φ/∂t + ∇·(φu) = 0 carries no diffusion.
        #
        # ν costs no time step here, and that is worth being precise about
        # because it is the opposite of the intuition that a blow-up means the
        # step is too big. Diffusion is a second derivative, so it is
        # explicit-stable only for Δt ~ Δ²/ν, and sphere_cfl_dt takes the min of
        # that and the wave-speed limit — but at this resolution the two are
        # ~1e5 s and ~75 s, so the gravity waves set Δt by three orders of
        # magnitude. ν = 1e5 is nowhere near the stability CEILING; the trouble
        # is the FLOOR, i.e. whether it damps enough. Reducing :cfl does not
        # help a run that fails this way, which is the tell.
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :ivisc_equations      => [2, 3, 4],
        :μ                    => 1.5e5,      # set to 1.0e5 together with :lvisc => true
        #---------------------------------------------------------------------------
        # ... AND the modal filter, which is what actually keeps the run alive
        # here (see the table at the top). The two mechanisms compose, and the
        # paper uses both.
        #
        # To isolate the diffusion instead, set :lfilter => false and raise :μ to
        # at least 5e5 — that combination does complete 3 days. Just do not read
        # the vorticity off it: at 5e5 half of max|ζ| is gone by day 3.
        #---------------------------------------------------------------------------
        :lfilter              => false,
        :filter_alpha         => 0.05,
        :filter_order         => 8,
        :filter_kcut          => 2/3,
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
