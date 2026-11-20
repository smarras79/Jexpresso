#-------------------------------------------------------------------------------------------------
# 1D diffusion with insulated (zero-flux) ends: the natural condition, needing no code at all
# (docs/boundary_conditions.md, §4.1.1).
# Problem and exact solution: initialize.jl. Boundary conditions: user_bc.jl.
#
#   julia --project=. -e 'using Jexpresso; Jexpresso.run_case("AdvDiff", "1d_insulated")'
#-------------------------------------------------------------------------------------------------
function user_inputs()
    inputs = Dict(
        :ode_solver           => SSPRK54(),
        :Δt                   => 2.0e-4,     # diffusion-limited: ν Δt/Δx_min² ≈ 0.07 (10 elements, nop 4)
        :tend                 => 1.0,
        :diagnostics_at_times => (0.5, 1.0),
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        :lsource              => false,
        :lvisc                => true,
        :μ                    => [0.1],      # ν
        :bdy_fluxes           => false,      # nothing to add: zero flux is what the weak form gives by default
        :lread_gmsh           => false,      # 1D mesh from :xmin, :xmax, :nelx
        :xmin                 => 0.0,
        :xmax                 => 1.0,
        :nelx                 => 10,
        :outformat            => "png",      # plots overlay the exact solution (user_analytic.jl)
        :loverwrite_output    => true,
        :output_dir           => "./output/",
    )
    return inputs
end
