#=---------------------------------------------------------------------------------
# run_gaussian_1d.jl
#
# Stand-alone example for SymbolicJE.jl
#
#   1D advection-diffusion of a gaussian wave on the periodic grid [-1, 1]:
#
#         ∂q/∂t + ∇⋅(u q) = μ ∇⋅∇(q)        with a small μ
#
# The equation is written exactly the way the user would type it (unicode /
# LaTeX), the grid is described by a Jexpresso-style `user_inputs()` Dict, and
# the PDE is parsed, discretized with central finite differences and integrated
# in time automatically.
#
# Run with:
#       julia run_gaussian_1d.jl
#---------------------------------------------------------------------------------=#

# Load the module ONCE per Julia session. Re-`include`ing this script (e.g. after
# editing the inputs below) then reuses the same module, so the symbolic `Node`s
# in `equation` and the `solve` method always come from the same module instance.
# (Re-defining the module on every include would create two incompatible `Node`
# types and break `solve`. To pick up edits to src/SymbolicJE.jl, restart Julia
# or use Revise.jl.)
isdefined(Main, :SymbolicJE) || include(joinpath(@__DIR__, "src", "SymbolicJE.jl"))
using .SymbolicJE

#---------------------------------------------------------------------------------
# 1. The equation, written with live Julia symbols (no string).
#    Declare the symbols, then write the residual ( = 0 implied):
#
#         ∂q/∂t + ∇⋅(u q) - μ ∇⋅∇(q) = 0      <=>   ∂q/∂t + ∇⋅(u q) = μ ∇⋅∇(q)
#
#    Notation notes: use `∂t(q)` (the literal `∂q/∂t` lexes as two identifiers),
#    and a `*` after a coefficient (`μ*…`, since `μ∇` would lex as one symbol).
#    A plain string equation, e.g. "∂q/∂t + ∇⋅(\\mathbf{u}q) = \\mu∇⋅∇(q)", also works.
#---------------------------------------------------------------------------------
@vars q u μ
equation = ∂t(q) + ∇⋅(u*q) - 2*μ*∇⋅∇(q)

#---------------------------------------------------------------------------------
# 2. Inputs — same spirit as a Jexpresso `user_inputs()`.
#---------------------------------------------------------------------------------
function user_inputs()
    inputs = Dict(
        #-------------------------------------------------------------------
        # Space dimension and grid (exactly like Jexpresso's 1D path)
        #-------------------------------------------------------------------
        :nsd       => 1,
        :method    => :sem,        # discretization: :fd (this) or :sem (next step)
        :xmin      => -1.0,
        :xmax      =>  1.0,
        :npoin     => 200,        # number of grid points
        :periodic  => true,       # periodic boundary conditions
        #-------------------------------------------------------------------
        # Physical parameters referenced by the equation string
        #-------------------------------------------------------------------
        :u         => [1.0],      # advecting velocity (a vector; ∇⋅ contracts it)
        :μ         => 1.0e-3,     # small diffusion coefficient in μ∇⋅∇(q)
        #-------------------------------------------------------------------
        # Initial condition: a gaussian wave centred at x = 0
        #-------------------------------------------------------------------
        :q0        => x -> exp(-(x^2) / (2 * 0.1^2)),
        #-------------------------------------------------------------------
        # Time integration
        #-------------------------------------------------------------------
        :tend      => 2.0,        # one full period for u = 1 on a domain of length 2
        :CFL       => 0.4,        # automatic Δt = CFL * stable step
        #-------------------------------------------------------------------
        # Output (PNG + on-the-fly plotting, exactly like problems/CompEuler/sod1d)
        #-------------------------------------------------------------------
        :output_dir           => joinpath(@__DIR__, "output"),
        :outformat            => "png",   # "png" (Plots.jl) or "ascii"
        :ndiagnostics_outputs => 20,      # number of on-the-fly snapshots
        :plot_live            => true,    # update an on-screen window during the run
    )
    return inputs
end

#---------------------------------------------------------------------------------
# 3. Solve.
#---------------------------------------------------------------------------------
mesh, q0, q = SymbolicJE.solve(equation, user_inputs())
