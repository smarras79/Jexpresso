#---------------------------------------------------------------------------------
# Environment overrides (same pattern as the sibling solid cases):
#
#   JEXPRESSO_SMA_STEPS        increments in each of the four load branches  [40]
#   JEXPRESSO_SMA_HEAT_STEPS   increments in the heating branch              [10]
#   JEXPRESSO_SMA_NOP          polynomial order                              [4]
#   JEXPRESSO_SMA_HYPOTHESIS   plane_stress | plane_strain                   [plane_stress]
#   JEXPRESSO_SMA_OUTFORMAT    vtk | none                                    [vtk]
#   JEXPRESSO_SMA_PLOT         true | false (deformed-shape PNGs + gif)      [true]
#   JEXPRESSO_SMA_VERIFY       elastic: the elastic arch against its exact
#                              solution, one load step (user_analytic.jl)    [off]
#---------------------------------------------------------------------------------
_sma_steps()      = something(tryparse(Int, _sma_env("JEXPRESSO_SMA_STEPS", "")), 40)
_sma_heat_steps() = something(tryparse(Int, _sma_env("JEXPRESSO_SMA_HEAT_STEPS", "")), 10)
_sma_nop()        = something(tryparse(Int, _sma_env("JEXPRESSO_SMA_NOP", "")), 4)

# Table 3: n increments in each of the branches [0,1], …, [3,4], m in [4,5]
function sma_step_times(n, m)
    t = [k/n for k = 0:4n]
    append!(t, [4 + k/m for k = 1:m])
    return t
end

function user_inputs()

    #---------------------------------------------------------------------------
    # The shape-memory-alloy arch of Artioli et al. (2017), Sec. 4.4, with
    # order-4 spectral elements, plane stress: the load–temperature history of
    # Table 3 over the load parameter t ∈ [0, 5] (see user_flux.jl).
    # Quasi-static; transformation strain at every quadrature point, Newton with
    # the consistent tangent at every increment (src/kernel/solvers/
    # quasistatic_solid.jl, quasistatic_solid_sma.jl).
    #
    # WHAT TO COMPARE: the paper's Fig. 5.6, the load q against the horizontal
    # displacement of A = (4.5, 0) — A_ux in solid_history.csv, or q from the
    # clamp reaction (−R_clamp_x) — and Fig. 5.5b, the deformed arch at t = 1
    # and t = 3. tools/artioli2017_figures.jl draws both from the VTK output.
    #---------------------------------------------------------------------------
    verify = sma_verify_mode()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # Quasi-static solid path, SMA, history at the quadrature points
        #---------------------------------------------------------------------------
        :lquasistatic_solid      => true,
        :tinit                   => 0.0,
        :tend                    => verify ? 1.0 : 5.0,
        :Δt                      => verify ? 1.0 : 1.0/_sma_steps(),
        :solid_step_times        => verify ? [0.0, 1.0] : sma_step_times(_sma_steps(), _sma_heat_steps()),
        :linitial_elastic_solve  => verify,           # the device starts unloaded: q(0) = 0
        :solid_history           => "quadrature",
        # 1e-8: near zero load the residual of points sitting on a transformation
        # kink cycles at ~1e-10 N/mm, which a tighter relative test cannot see past
        :newton_rtol             => 1.0e-8,
        :newton_maxit            => 25,
        :solid_probes            => [("A", 4.5, 0.0), ("C", 0.0, 4.5)],    # loaded end, crown
        :diagnostics_at_times    => verify ? [0.0, 1.0] : sma_step_times(_sma_steps(), _sma_heat_steps()),
        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes     => "lgl",
        :nop                     => _sma_nop(),
        #---------------------------------------------------------------------------
        # Mesh: the arch in 2 x 24 quads (problems/Elasticity/make_meshes.py), both
        # arcs curved back onto their circles.
        #---------------------------------------------------------------------------
        :lread_gmsh              => true,
        :gmsh_filename           => "./problems/Elasticity/sma_arch2d/half_annulus_2x24.msh",
        :exact_geometry          => Dict("inner" => (:circle, 0.0, 0.0, 3.5),
                                         "outer" => (:circle, 0.0, 0.0, 4.5)),
        #---------------------------------------------------------------------------
        # Output: VTK at every increment (u as a vector for ParaView's Warp By
        # Vector), and the deformed arch at true scale coloured by the
        # saturation of the transformation strain, as PNGs and a gif.
        #---------------------------------------------------------------------------
        :outformat               => _sma_env("JEXPRESSO_SMA_OUTFORMAT", "vtk"),
        :vtk_vectors             => [("u", "ux", "uy")],
        :lplot_deformed          => lowercase(_sma_env("JEXPRESSO_SMA_PLOT", "true")) == "true",
        :plot_deformed_field     => "etr_sat",
        :plot_deformed_scale     => 1.0,
        :plot_deformed_clims     => (0.0, 1.0),
        :loverwrite_output       => true,
        :output_dir              => "./output",
    ) #Dict

    return inputs
end
