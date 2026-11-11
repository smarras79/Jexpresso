#---------------------------------------------------------------------------------
# Environment overrides, so that a Δt sweep, the other Prony set or the paper's
# load convention can be run without editing this file (same pattern as
# MHD/astroJetWuShu2018). Unset, each one takes the default shown.
#
#   JEXPRESSO_VE_SET         ve1 | ve2            Prony set (read in user_flux.jl)   [ve1]
#   JEXPRESSO_VE_DT          time step                                               [1.0]
#   JEXPRESSO_VE_LOAD        step | first_increment                                  [step]
#   JEXPRESSO_VE_INTEGRATOR  exponential | midpoint | backward_euler                [exponential]
#   JEXPRESSO_VE_NOP         polynomial order                                        [4]
#   JEXPRESSO_VE_OUTFORMAT   vtk | none   (none: solid_history.csv only)             [vtk]
#---------------------------------------------------------------------------------
_ve_env(key, default) = String(strip(get(ENV, key, default)))
_ve_dt()   = something(tryparse(Float64, _ve_env("JEXPRESSO_VE_DT", "")), 1.0)
_ve_nop()  = something(tryparse(Int,     _ve_env("JEXPRESSO_VE_NOP", "")), 4)
function _ve_initial_solve()
    load = lowercase(_ve_env("JEXPRESSO_VE_LOAD", "step"))
    load in ("step", "first_increment") ||
        error(" # ERROR viscoelastic_cylinder2d: JEXPRESSO_VE_LOAD = \"", load,
              "\" — use step or first_increment.")
    return load == "step"
end

function user_inputs()

    #---------------------------------------------------------------------------
    # Thick-walled viscoelastic cylinder under internal pressure, plane strain:
    # Artioli, Beirão da Veiga, Lovadina & Sacco, Comput. Mech. 60 (2017),
    # Sec. 4.2 — here with order-4 spectral elements.
    #
    # QUASI-STATIC: no inertia. Every time step is one linear solve; time enters
    # only through the material's memory (generalized Maxwell in the deviator).
    # See src/kernel/solvers/quasistatic_solid.jl.
    #
    # WHAT TO EXPECT (Prony set ve1): u_r at A = (2,0) creeps from the Lamé
    # value 0.0381333 at t = 0⁺ to 0.66078 at t = 20, and at B = (4,0) from
    # 0.0242667 to 0.33638. With 20 unit steps the exponential integrator is
    # within 5e-4 of that closed form at A; halve Δt and the gap quarters. The
    # spatial error at nop = 4 is orders of magnitude below that (README).
    #
    # Material, geometry and load live in user_flux.jl (cylinder_properties).
    #---------------------------------------------------------------------------
    inputs = Dict(
        #---------------------------------------------------------------------------
        # The quasi-static solid path, and its time stepping
        #---------------------------------------------------------------------------
        :lquasistatic_solid      => true,
        :tinit                   => 0.0,
        :tend                    => 20.0,
        :Δt                      => _ve_dt(),         # the paper: 20 unit steps
        #
        # "exponential": q̇ + q/λ = ė integrated exactly for a strain linear in
        # time over the step — the recursion of Zienkiewicz, Taylor & Fox (2014)
        # that FEAP, and so the paper, uses. Second order. "midpoint" is the
        # Simo & Hughes (1998) Box 10.1 variant, also second order;
        # "backward_euler" is first order and is there as a control.
        #
        :viscoelastic_integrator => _ve_env("JEXPRESSO_VE_INTEGRATOR", "exponential"),
        #
        # true:  p is a Heaviside step at t = 0 — an instantaneous (glassy)
        #        elastic solve at t = 0⁺, then the creep. This is the problem the
        #        closed form in user_analytic.jl solves, so the error against it
        #        is reported at every output time.
        # false: start unloaded and let p enter during the FIRST step, which the
        #        integrator reads as a ramp over [0, Δt]. That is how FEAP runs it
        #        and what the paper's Fig. 5.2 shows (Δt = 1).
        #
        :linitial_elastic_solve  => _ve_initial_solve(),
        #
        # The paper's control points: A on the inner, B on the outer surface,
        # both on the x axis, where u_r = ux.
        #
        :solid_probes            => [("A", 2.0, 0.0), ("B", 4.0, 0.0)],
        :diagnostics_at_times    => collect(0.0:1.0:20.0),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes     => "lgl",
        :nop                     => _ve_nop(),
        #---------------------------------------------------------------------------
        # Mesh: a quarter annulus 2 ≤ r ≤ 4 of 4 (radial) x 8 (circumferential)
        # quads; regenerate with problems/Elasticity/make_meshes.py. gmsh writes
        # it straight-sided, so the two arcs are curved back onto their circles:
        # the boundary is then the isoparametric degree-nop interpolant of the
        # circle, and the pressure is integrated on it with its exact normal.
        #---------------------------------------------------------------------------
        :lread_gmsh              => true,
        :gmsh_filename           => "./problems/Elasticity/viscoelastic_cylinder2d/quarter_annulus_4x8.msh",
        :exact_geometry          => Dict("inner" => (:circle, 0.0, 0.0, 2.0),
                                         "outer" => (:circle, 0.0, 0.0, 4.0)),
        #---------------------------------------------------------------------------
        # Output: VTK fields at every unit time, plus solid_history.csv at every
        # step (displacement at A and B, the exact one, and the error norms).
        #---------------------------------------------------------------------------
        :outformat               => _ve_env("JEXPRESSO_VE_OUTFORMAT", "vtk"),
        :loverwrite_output       => true,
        :output_dir              => "./output",
    ) #Dict

    return inputs
end
