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
include("./solver/mod_solution.jl")
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
#--------------------------------------------------------
mod_mesh_mesh_driver(inputs)


#--------------------------------------------------------
# Build mass matrix
#--------------------------------------------------------
TInt=Int64
TFloat=Float64

N   = 6
Nit = 100
Tol = 0.1

P1  = LGL1D()
P2  = CGL1D()
T1  = NodalGalerkin()
T2  = Collocation()

dim  = 1
dims = [N]
ND   = build_nodal_Storage(dims,P1,T1)


#ElementMassMatrix(Int64(inputs[:nop]), Int64(inputs[:nop]), MassMatrix1D(), Float64)
        

