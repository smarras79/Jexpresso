#-------------------------------------------------------------------------------------------------
# 1D diffusion with a nonlinear (radiation) flux condition on the right and a Dirichlet value on
# the left (docs/boundary_conditions.md, §4.1.1).
# Problem and exact solution: initialize.jl. Boundary conditions: user_bc.jl.
#
#   julia --project=. -e 'using Jexpresso; Jexpresso.run_case("AdvDiff", "1d_radiation")'
#-------------------------------------------------------------------------------------------------
function user_inputs()
    inputs = Dict(
        :ode_solver           => SSPRK54(),
        :Δt                   => 1.0e-3,     # diffusion-limited: ν Δt/Δx_min² ≈ 0.05 (4 elements, nop 4)
        :tend                 => 20.0,       # ~15 e-folds of the slowest mode: steady state
        :diagnostics_at_times => (1.0, 5.0, 20.0),
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        :lsource              => false,
        :lvisc                => true,
        :μ                    => [0.1],      # ν
        :rad_coeff            => 0.1,        # εσ: inward flux εσ (u∞⁴ - u⁴) at the right end
        :rad_uinf             => 1.0,        # u∞, temperature of the surroundings
        :u_left               => 2.0,        # Dirichlet value at the left end
        :bdy_fluxes           => true,       # REQUIRED for the flux condition (user_bc_neumann!)
        :lread_gmsh           => false,      # 1D mesh from :xmin, :xmax, :nelx
        :xmin                 => 0.0,
        :xmax                 => 1.0,
        :nelx                 => 4,
        :outformat            => "png",      # plots overlay the exact solution (user_analytic.jl)
        :loverwrite_output    => true,
        :output_dir           => "./output/",
    )
    return inputs
end
