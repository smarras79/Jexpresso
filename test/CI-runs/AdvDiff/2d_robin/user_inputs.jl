#-------------------------------------------------------------------------------------------------
# 2D advection-diffusion with Robin walls: a template for natural (flux) boundary conditions.
# Problem, exact solution and parameters: initialize.jl. How the Robin condition is imposed:
# user_bc.jl.
#
#   julia --project=. -e 'using Jexpresso; Jexpresso.run_case("AdvDiff", "2d_robin")'
#-------------------------------------------------------------------------------------------------
function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # Time integration
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK54(),
        :Δt                   => 1.0e-3,     # diffusion-limited (ν Δt/Δx_min² ≈ 0.08 on 10x5, nop 4); halve when refining
        :tend                 => 1.0,
        :diagnostics_at_times => (0.25, 0.5, 0.75, 1.0),
        #---------------------------------------------------------------------------
        # Discretization
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        :lsource              => false,
        #---------------------------------------------------------------------------
        # Physics: ∂u/∂t + c ∂u/∂x = ν ∇²u,  -ν ∂u/∂n = h (u - u∞) on the walls
        #---------------------------------------------------------------------------
        :lvisc                => true,
        :μ                    => [0.1],      # ν, the diffusivity
        :adv_c                => 1.0,        # c, advection speed along x (periodic direction)
        :robin_h              => 0.2,        # h, transfer coefficient (Biot number hH/ν = 2)
        :robin_uinf_bottom    => 1.0,        # u∞ at y = ymin
        :robin_uinf_top       => 0.0,        # u∞ at y = ymax
        :robin_amp            => 0.5,        # amplitude of the decaying mode in the initial condition
        :bdy_fluxes           => true,       # REQUIRED: calls user_bc_neumann! and adds ∮ ψ F dS
        #---------------------------------------------------------------------------
        # Mesh: [0, 2] x [0, 1], periodic in x, walls tagged "bottom" and "top" (robin_2d.geo)
        #---------------------------------------------------------------------------
        :lread_gmsh           => true,
        :gmsh_filename        => "./test/CI-runs/AdvDiff/2d_robin/robin_2d_10x5.msh",
        #---------------------------------------------------------------------------
        # Output
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :output_dir           => "./output/",
        :loverwrite_output    => true,
    )
    return inputs
end
