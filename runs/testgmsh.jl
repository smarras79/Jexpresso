using Crayons.Box
using Gridap
using GridapGmsh
using MPI

model = GmshDiscreteModel("/Users/simone/Work/Codes/jexpresso/demo/gmsh_grids/hexa_UNSTR.msh")

writevtk(model,"gmsh_grid")


order = 10
reffe = ReferenceFE(lagrangian,Float64,order)
#=
V = TestFESpace(model,reffe,dirichlet_tags=["boundary1","boundary2"])
U = TrialFESpace(V,[0,1])
Ω = Triangulation(model)
dΩ = Measure(Ω,2*order)
a(u,v) = ∫( ∇(v)⋅∇(u) )dΩ
l(v) = 0
op = AffineFEOperator(a,l,U,V)
uh = solve(op)
writevtk(Ω,"gmsh_grid"); #,cellfields=["uh"=>uh])
=#
