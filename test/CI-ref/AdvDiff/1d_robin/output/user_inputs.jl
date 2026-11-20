#-------------------------------------------------------------------------------------------------
# 1D diffusion with a general Robin condition  a u + b ν ∂u/∂n = g  at both ends, with different
# coefficients at each end (docs/boundary_conditions.md, §4.1.1).
# Problem and exact solution: initialize.jl. Boundary conditions: user_bc.jl.
#
#   julia --project=. -e 'using Jexpresso; Jexpresso.run_case("AdvDiff", "1d_robin")'
#-------------------------------------------------------------------------------------------------
function user_inputs()
    inputs = Dict(
        :ode_solver           => SSPRK54(),
        :Δt                   => 2.0e-4,     # diffusion-limited: ν Δt/Δx_min² ≈ 0.07 (10 elements, nop 4)
        :tend                 => 2.0,
        :diagnostics_at_times => (0.5, 1.0, 2.0),
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        :lsource              => false,
        :lvisc                => true,
        :μ                    => [0.1],      # ν
        # Robin coefficients  a u + b ν ∂u/∂n = g  (n outward), per end. b ≠ 0 (b = 0 is Dirichlet)
        # and a/b > 0 (the wall exchanges heat with an ambient; a/b < 0 would be a heat pump).
        :robin_left           => (a = 0.5, b = 1.0, g = 0.5),
        :robin_right          => (a = 1.0, b = 4.0, g = 0.0),
        :bdy_fluxes           => true,       # REQUIRED for the flux condition (user_bc_neumann!)
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
