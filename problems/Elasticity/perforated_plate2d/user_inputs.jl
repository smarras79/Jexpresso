#---------------------------------------------------------------------------------
# Environment overrides (same pattern as the sibling solid cases):
#
#   JEXPRESSO_PP_STEPS       number of load increments                         [400]
#   JEXPRESSO_PP_NOP         polynomial order                                  [4]
#   JEXPRESSO_PP_OUTFORMAT   vtk | none                                        [vtk]
#   JEXPRESSO_PP_PLOT        true | false  (deformed-shape PNGs + gif)         [true]
#---------------------------------------------------------------------------------
_pp_env(key, default) = String(strip(get(ENV, key, default)))
_pp_steps() = something(tryparse(Int, _pp_env("JEXPRESSO_PP_STEPS", "")), 400)
_pp_nop()   = something(tryparse(Int, _pp_env("JEXPRESSO_PP_NOP", "")), 4)

function user_inputs()

    #---------------------------------------------------------------------------
    # The perforated plastic plate of Artioli et al. (2017), Sec. 4.3, with
    # order-4 spectral elements: plane strain, J2 perfectly plastic, the top
    # edge pulled up by 2 mm in 400 increments of the load parameter
    # t ∈ [0, 1]. Quasi-static; plastic history at every quadrature point,
    # Newton with the consistent tangent at every increment. See user_flux.jl
    # and src/kernel/solvers/quasistatic_solid.jl.
    #
    # WHAT TO COMPARE (paper, Table 1, end of loading): u_A = 2.741 mm (the
    # horizontal displacement of A = (R, 0)), v_B = 1.859 mm (the vertical one
    # of B = (0, R)) — A_ux and B_uy in solid_history.csv. Fig. 5.4: the
    # reaction sum on the top edge, R_top_y x 10 (the paper's plate is 10 mm
    # thick). tools/artioli2017_figures.jl draws both from the VTK output.
    #---------------------------------------------------------------------------
    inputs = Dict(
        #---------------------------------------------------------------------------
        # Quasi-static solid path, J2 plasticity, history at the quadrature points
        #---------------------------------------------------------------------------
        :lquasistatic_solid      => true,
        :tinit                   => 0.0,
        :tend                    => 1.0,
        :Δt                      => 1.0/_pp_steps(),  # δ = 2 mm in 400 equal increments
        :linitial_elastic_solve  => false,            # start unloaded; δ(0) = 0
        :solid_history           => "quadrature",
        :newton_rtol             => 1.0e-10,
        :newton_maxit            => 25,
        :solid_probes            => [("A", 50.0, 0.0), ("B", 0.0, 50.0)],
        :diagnostics_at_times    => collect(0.0:0.01:1.0),   # every 4 increments
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes     => "lgl",
        :nop                     => _pp_nop(),
        #---------------------------------------------------------------------------
        # Mesh: a quarter of the strip in 132 quads, graded towards the hole
        # (problems/Elasticity/make_meshes.py); the hole is curved back onto its
        # circle.
        #---------------------------------------------------------------------------
        :lread_gmsh              => true,
        :gmsh_filename           => "./problems/Elasticity/perforated_plate2d/quarter_plate_hole_132.msh",
        :exact_geometry          => Dict("hole" => (:circle, 0.0, 0.0, 50.0)),
        #---------------------------------------------------------------------------
        # Output: VTK (ux, uy, stresses, svm, ep_bar, the reactions rx, ry and
        # the vector u for ParaView's Warp By Vector) every 0.01 in t, and the
        # deformed grid (displacement x5) coloured by the plastic strain, as
        # PNGs and a gif. tools/artioli2017_figures.jl reads the VTK files.
        #---------------------------------------------------------------------------
        :outformat               => _pp_env("JEXPRESSO_PP_OUTFORMAT", "vtk"),
        :vtk_vectors             => [("u", "ux", "uy")],
        :lplot_deformed          => lowercase(_pp_env("JEXPRESSO_PP_PLOT", "true")) == "true",
        :plot_deformed_field     => "ep_bar",
        :plot_deformed_scale     => 5.0,
        :loverwrite_output       => true,
        :output_dir              => "./output",
    ) #Dict

    return inputs
end
