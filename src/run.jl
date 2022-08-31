#--------------------------------------------------------
# external packages
#--------------------------------------------------------
using Crayons.Box
using DifferentialEquations
using Gridap
using GridapGmsh
using MPI
using Revise
using WriteVTK

#Plots
using Plots; gr()
plotlyjs()

#Constants
const TInt   = Int64
const TFloat = Float64

#--------------------------------------------------------
# jexpresso modules
#--------------------------------------------------------
include("./IO/mod_inputs.jl")
include("./Mesh/mod_mesh.jl")
include("./basis/basis_structs.jl")
include("./Infrastructure/Kopriva_functions.jl")
include("./Infrastructure/2D_3D_structures.jl")
include("./element_matrices.jl")
include("../tests/plot_lagrange_polynomial.jl")
#--------------------------------------------------------

struct EDGES <:At_geo_entity end
struct FACES <:At_geo_entity end
#MPI.Init()
#comm = MPI.COMM_WORLD

#if MPI.Comm_rank(comm) == 0
mod_inputs_print_welcome()

#--------------------------------------------------------
#User inputs:
#--------------------------------------------------------
inputs        = Dict{}()
inputs, nvars = mod_inputs_user_inputs()

#--------------------------------------------------------
# Create/read mesh
# return mesh::St_mesh
#--------------------------------------------------------
mesh = mod_mesh_mesh_driver(inputs)

#--------------------------------------------------------
# Problem setup
# !!!!!!
# !!!!!! WARNING: MOVE all the setup parameters to user_input.jl
# !!!!!!
#--------------------------------------------------------
exact_quadrature = false

P1  = LGL1D()
P2  = CGL1D()
T1  = NodalGalerkin()
T2  = Collocation()

#--------------------------------------------------------
# Build interpolation nodes:
#             the user decides among LGL, GL, etc. 
# Return:
# ξ = ND.ξ.ξ
# ω = ND.ξ.ω
#--------------------------------------------------------
N  = inputs[:nop]
ND = build_nodal_Storage([N],P1,T1) # --> ξ <- ND.ξ.ξ
ξ  = ND.ξ.ξ

exact_integration = inputs[:lexact_integration]
if exact_integration
    #
    # Exact quadrature:
    # Quadrature order (Q = N+1) ≠ polynomial order (N)
    #
    TP  = Exact()
    Q   = N + 1
    
    NDq = build_nodal_Storage([Q],P1,T1) # --> ξ <- ND.ξ.ξ
    ξq  = NDq.ξ.ξ
    ω   = NDq.ξ.ω
    
else  
    #
    # Inexact quadrature:
    # Quadrature and interpolation orders coincide (Q = N)
    #
    TP  = Inexact()
    Q   = N
    NDq = ND
    ξq  = ξ
    ω   = ND.ξ.ω
end



#--------------------------------------------------------
# Build Lagrange polynomials:
#
# Return:
# ψ     = basis.ψ[N+1, Q+1]
# dψ/dξ = basis.dψ[N+1, Q+1]
#--------------------------------------------------------
basis = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat)


#--------------------------------------------------------
# Build element mass matrix
#
# Return:
# el_mat.M[iel, i, j] <-- if exact (full)
# el_mat.M[iel, i]    <-- if inexact (diagonal)
# el_mat.D[iel, i, j] <-- either exact (full) OR inexact (sparse)
#--------------------------------------------------------
el_mat = build_element_matrices!(TP, basis.ψ, basis.dψ, ω, mesh.nelem, N, Q, TFloat)
