#=---------------------------------------------------------------------------------
# run_advdiff_2d.jl
#
# Stand-alone 2D example for SymbolicJE.jl — the 2D analogue of run_gaussian_1d.jl,
# set up exactly like Jexpresso's problems/AdvDiff/kopriva:
#
#       ∂q/∂t + ∇⋅(u q) = μ ∇⋅∇(q)
#
#   * a gaussian blob centred at (xc, yc) = (0, 3), σx = σy = 1, A = 1
#   * advecting velocity u = (0.5, 1.0)
#   * isotropic diffusion μ = 0.1
#   * domain x∈[-10,10], y∈[0,10], periodic
#
# The SEM backend reads the SAME gmsh grid kopriva uses
# (./meshes/gmsh_grids/kopriva_periodic.msh) through Jexpresso's sem_setup: the
# mesh, the metric terms and the dψ/ω basis come straight from Jexpresso, and the
# weak-form operator is the one in rhs.jl's _expansion_inviscid! (DSS + Minv).
#
# Run from the Jexpresso root (so the gmsh path and project resolve):
#       julia --project=. tools/SymbolicJE.jl/run_advdiff_2d.jl
#
# To run instead on a structured Cartesian grid (no gmsh / no Jexpresso needed),
# set :method => :fd (or :sem without :gmsh_filename) and give :npoinx/:npoiny
# (FD) or :nelx/:nely (SEM); see the commented block below.
#---------------------------------------------------------------------------------=#

isdefined(Main, :SymbolicJE) || include(joinpath(@__DIR__, "src", "SymbolicJE.jl"))
using .SymbolicJE

#---------------------------------------------------------------------------------
# 1. The equation, with live Julia symbols (residual form, = 0 implied):
#
#       ∂q/∂t + ∇⋅(u q) - μ ∇⋅∇(q) = 0
#
#    u is a 2-vector here (∇⋅ contracts it over both directions); μ is the scalar
#    (isotropic) diffusivity used by the ∇⋅∇ term.
#---------------------------------------------------------------------------------
@vars q u μ
equation = ∂t(q) + ∇⋅(u*q) - μ*∇⋅∇(q)

#---------------------------------------------------------------------------------
# 2. Inputs — the kopriva setup. The gmsh grid is read by Jexpresso's sem_setup.
#---------------------------------------------------------------------------------
function user_inputs()
    xc, yc, sx, sy, A = 0.0, 3.0, 1.0, 1.0, 1.0          # kopriva gaussian blob
    inputs = Dict(
        #-------------------------------------------------------------------
        # Space dimension, discretization and the EXISTING kopriva grid
        #-------------------------------------------------------------------
        :nsd                  => 2,
        :method               => :sem,         # :sem (gmsh grid) or :fd (structured)
        :nop                  => 4,            # polynomial order (LGL: nop+1 nodes/dir)
        :interpolation_nodes  => "lgl",
        :periodic             => true,
        # the SAME mesh file referenced in problems/AdvDiff/kopriva/user_inputs.jl
        :lread_gmsh           => true,
        :gmsh_filename        => "./meshes/gmsh_grids/kopriva_periodic.msh",
        #-------------------------------------------------------------------
        # Physical parameters referenced by the equation
        #-------------------------------------------------------------------
        :u                    => [0.5, 1.0],   # advecting velocity (vector)
        :μ                    => 0.1,          # isotropic diffusion coefficient
        #-------------------------------------------------------------------
        # Initial condition: gaussian centred at (xc, yc) (kopriva initialize.jl)
        #-------------------------------------------------------------------
        :q0                   => (x, y) -> A * exp(-((x - xc) / sx)^2) * exp(-((y - yc) / sy)^2),
        #-------------------------------------------------------------------
        # Time integration (kopriva: Δt = 0.005, tend = 10)
        #-------------------------------------------------------------------
        :tend                 => 10.0,
        :Δt                   => 0.005,
        #-------------------------------------------------------------------
        # Output (PNG node-map snapshots + on-the-fly window)
        #-------------------------------------------------------------------
        :output_dir           => joinpath(@__DIR__, "output_2d"),
        :outformat            => "png",
        :ndiagnostics_outputs => 20,
        :plot_live            => true,
    )
    return inputs
end

# --- structured-grid alternative (no gmsh / no Jexpresso) ------------------------
# Uncomment to run the very same equation on a Cartesian grid over the kopriva box:
#
# function user_inputs()
#     xc, yc = 0.0, 3.0
#     return Dict(:nsd => 2, :method => :fd, :periodic => true,
#                 :xmin => -10.0, :xmax => 10.0, :ymin => 0.0, :ymax => 10.0,
#                 :npoinx => 200, :npoiny => 100,
#                 :u => [0.5, 1.0], :μ => 0.1,
#                 :q0 => (x, y) -> exp(-((x-xc))^2) * exp(-((y-yc))^2),
#                 :tend => 10.0, :CFL => 0.4,
#                 :output_dir => joinpath(@__DIR__, "output_2d"),
#                 :outformat => "png", :ndiagnostics_outputs => 20, :plot_live => true)
# end

#---------------------------------------------------------------------------------
# 3. Solve.
#---------------------------------------------------------------------------------
mesh, q0, q = SymbolicJE.solve(equation, user_inputs())
