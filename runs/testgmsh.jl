using Crayons.Box
using Gridap
using GridapGmsh
using MPI

model = GmshDiscreteModel("/Users/simone/Work/Codes/jexpresso/demo/gmsh_grids/hexa_UNSTR.msh")

#model.grid.node_coordinates <-- coords
#model.grid.cell_node_ids    <-- conn

writevtk(model,"gmsh_grid")
