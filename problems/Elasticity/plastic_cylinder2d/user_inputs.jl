#---------------------------------------------------------------------------------
# Environment overrides (same pattern as viscoelastic_cylinder2d):
#
#   JEXPRESSO_PC_DT          load increment in t ∈ [0, 1]                      [0.02]
#   JEXPRESSO_PC_NOP         polynomial order                                  [4]
#   JEXPRESSO_PC_OUTFORMAT   vtk | none                                        [vtk]
#   JEXPRESSO_PC_PLOT        true | false  (deformed-shape PNGs + gif)         [true]
#---------------------------------------------------------------------------------
_pc_env(key, default) = String(strip(get(ENV, key, default)))
_pc_dt()  = something(tryparse(Float64, _pc_env("JEXPRESSO_PC_DT", "")), 0.02)
_pc_nop() = something(tryparse(Int,     _pc_env("JEXPRESSO_PC_NOP", "")), 4)

function user_inputs()

    #---------------------------------------------------------------------------
    # Thick-walled cylinder, J2 elastic-perfectly plastic, plane strain, its
    # bore displaced radially from 0 to δ = 0.05 over the load parameter
    # t ∈ [0, 1]. Quasi-static; the plastic history (eᵖ, ēᵖ) is carried at
    # every quadrature point and each increment is a Newton iteration with the
    # consistent tangent. See user_flux.jl and src/kernel/solvers/quasistatic_solid.jl.
    #
    # WHAT TO EXPECT: the bore pressure p = R_inner,x / a (solid_history.csv)
    # follows Lamé up to first yield, p_e = 10.50 at δ = 5.7e-3; the plastic zone
    # (ep_bar > 0) then spreads outwards from the bore, reaches the outer
    # surface, and p levels off at the exact limit p_L = (2/√3) σy0 ln 2 = 19.45.
    # Two or three Newton iterations per increment.
    #---------------------------------------------------------------------------
    inputs = Dict(
        #---------------------------------------------------------------------------
        # Quasi-static solid path, J2 plasticity (user_solid_material), history
        # at the quadrature points (the only option for plasticity)
        #---------------------------------------------------------------------------
        :lquasistatic_solid      => true,
        :tinit                   => 0.0,
        :tend                    => 1.0,
        :Δt                      => _pc_dt(),        # 50 increments of the bore displacement
        :linitial_elastic_solve  => false,           # start unloaded; δ(0) = 0 anyway
        :solid_history           => "quadrature",
        :newton_rtol             => 1.0e-10,
        :newton_maxit            => 25,
        :solid_probes            => [("A", 2.0, 0.0), ("B", 4.0, 0.0)],
        :diagnostics_at_times    => collect(0.0:0.1:1.0),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes     => "lgl",
        :nop                     => _pc_nop(),
        #---------------------------------------------------------------------------
        # Mesh: the quarter annulus 2 ≤ r ≤ 4, 8 (radial) x 8 quads, both arcs
        # curved back onto their circles (problems/Elasticity/make_meshes.py).
        #---------------------------------------------------------------------------
        :lread_gmsh              => true,
        :gmsh_filename           => "./problems/Elasticity/plastic_cylinder2d/quarter_annulus_8x8.msh",
        :exact_geometry          => Dict("inner" => (:circle, 0.0, 0.0, 2.0),
                                         "outer" => (:circle, 0.0, 0.0, 4.0)),
        #---------------------------------------------------------------------------
        # Output. VTK: ux, uy, the stresses, svm and ep_bar, plus the vector "u"
        # for ParaView's Warp By Vector. PNG: the deformed grid, displacement
        # magnified 10x, coloured by the plastic strain — the plastic zone
        # growing from the bore.
        #---------------------------------------------------------------------------
        :outformat               => _pc_env("JEXPRESSO_PC_OUTFORMAT", "vtk"),
        :vtk_vectors             => [("u", "ux", "uy")],
        :lplot_deformed          => lowercase(_pc_env("JEXPRESSO_PC_PLOT", "true")) == "true",
        :plot_deformed_field     => "ep_bar",
        :plot_deformed_scale     => 10.0,
        :plot_deformed_clims     => (0.0, 0.03),
        :loverwrite_output       => true,
        :output_dir              => "./output",
    ) #Dict

    return inputs
end
