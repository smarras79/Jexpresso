# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/NodeNumbering.jl/blob/master/LICENSE

module NodeNumbering

using PyPlot: matshow

include("calc_bw.jl")
include("create_adjacency_graph.jl")
include("node_degrees.jl")
include("RCM.jl")
include("renumbering.jl")
include("create_RCM_adjacency.jl")
include("adjacency_visualization.jl")

#=
elements = Dict(
                  1 => [15, 1, 8, 4],
                  2 => [1, 3, 2, 8],
                  3 => [3, 11, 13, 2],
                  4 => [4, 8, 5, 10],
                  5 => [8, 2, 7, 5],
                  6 => [2, 13, 6, 7],
                  7 => [5, 7, 12],
                  8 => [7, 6, 14, 12],
                  9 => [12, 14, 9]);

element_types = Dict(
                       1 => :Quad4,
                       2 => :Quad4,
                       3 => :Quad4,
                       4 => :Quad4,
                       5 => :Quad4,
                       6 => :Quad4,
                       7 => :Tri3,
                       8 => :Quad4,
                       9 => :Tri3);

P = 15
=#
elements = Dict(
                  1 => [1, 2, 3, 4, 5, 6, 7, 8],
                  2 => [2, 9, 11, 3, 6, 10, 12, 7]);

element_types = Dict(
                       1 => :Hexa8,
                       2 => :Hexa8);

P = 2

adjacency = create_adjacency_graph(elements, element_types)

element_adjacencies[:Quad4]
element_types[1]

degrees = node_degrees(adjacency)
neworder = RCM(adjacency, degrees, P, P)
finalorder = renumbering(neworder)
RCM_adjacency = create_RCM_adjacency(adjacency, finalorder)
newmatrix = adjacency_visualization(RCM_adjacency)
display(heatmap(newmatrix, cgrad=(scale=:log10)))

end
