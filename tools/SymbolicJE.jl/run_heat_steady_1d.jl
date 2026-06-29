#=---------------------------------------------------------------------------------
# run_heat_steady_1d.jl
#
# Stand-alone example of a TIME-INDEPENDENT problem for SymbolicJE.jl.
#
# Because the equation has no ∂q/∂t term, it is solved as a steady (elliptic)
# problem: the discrete operator is assembled into a matrix and solved directly,
# with Dirichlet boundary data.
#
#   Steady heat equation (Poisson form) on the rod [-1, 1]:
#
#         ∇⋅∇(q) = f          with   q(-1) = q(1) = 0
#
#   and a localized heat source f(x) (so q bulges up between the cold ends).
#
# Run with:
#       julia --project=. tools/SymbolicJE.jl/run_heat_steady_1d.jl
#---------------------------------------------------------------------------------=#

# Load the module ONCE per session (see note in run_gaussian_1d.jl): re-including
# this script reuses the same module so the symbolic `Node`s and `solve` match.
isdefined(Main, :SymbolicJE) || include(joinpath(@__DIR__, "src", "SymbolicJE.jl"))
using .SymbolicJE

# Equation as live symbols, residual form (no ∂t term ⇒ steady solve):
#   ∇⋅∇(q) - f = 0   <=>   ∇²q = f.   (A string "∇⋅∇(q) = f" also works.)
@vars q f
equation = ∇⋅∇(q) - f

function user_inputs()
    inputs = Dict(
        #-------------------------------------------------------------------
        # Space dimension and grid (steady ⇒ non-periodic, with Dirichlet BCs)
        #-------------------------------------------------------------------
        :nsd       => 1,
        :xmin      => -1.0,
        :xmax      =>  1.0,
        :npoin     => 201,
        :periodic  => false,
        #-------------------------------------------------------------------
        # Right-hand-side source f(x): a localized (Gaussian) heat source.
        # ∇²q = f with f < 0 drives q upward between the cold ends.
        #-------------------------------------------------------------------
        :f         => x -> -10.0 * exp(-(x / 0.15)^2),
        #-------------------------------------------------------------------
        # Dirichlet boundary values  q(xmin), q(xmax)
        #-------------------------------------------------------------------
        :bc_left   => 0.0,
        :bc_right  => 0.0,
        #-------------------------------------------------------------------
        # Output
        #-------------------------------------------------------------------
        :output_dir => joinpath(@__DIR__, "output_steady"),
        :outformat  => "png",
        :plot_live  => true,
    )
    return inputs
end

mesh, q0, q = SymbolicJE.solve(equation, user_inputs())
