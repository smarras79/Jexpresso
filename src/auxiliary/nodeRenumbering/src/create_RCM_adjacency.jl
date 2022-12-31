# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/NodeNumbering.jl/blob/master/LICENSE

"""
    create_RCM_adjacency(adjacency::Dict{Int, Vector{Int}}, finalorder::Dict{Int, Int})

Create an adjacency Dict of the RCM ordered nodes. A part of the Reverse Cuthill-McKee Algorithm.

# Examples

```jldoctest
julia> adjacency = Dict{Int, Vector{Int}}(
                   1 => [3, 8, 9],
                   2 => [3, 8, 7, 6],
                   3 => [1, 2],
                   4 => [8, 9],
                   5 => [7, 8],
                   6 => [2, 7],
                   7 => [5, 2, 6],
                   8 => [1, 2, 4, 5],
                   9 => [1, 4]);

julia> finalorder =  Dict{Int, Int}(
                     7 => 2,
                     9 => 9,
                     4 => 8,
                     2 => 3,
                     3 => 5,
                     5 => 4,
                     8 => 6,
                     6 => 1,
                     1 => 7);

julia> create_RCM_adjacency(adjacency, finalorder)
Dict{Int64,Array{Int64,1}} with 9 entries:
  7 => [5, 6, 9]
  9 => [7, 8]
  4 => [2, 6]
  2 => [4, 3, 1]
  3 => [5, 6, 2, 1]
  8 => [6, 9]
  5 => [7, 3]
  6 => [7, 3, 8, 4]
  1 => [3, 2]

```

# References

* Wikipedia contributors. "Adjacency list". Wikipedia, The Free Encyclopedia. Wikipedia, The Free Encyclopedia, 7 Jun. 2017. Web. 17 Jul. 2017. https://en.wikipedia.org/wiki/Adjacency_list
"""
function create_RCM_adjacency(adjacency::Dict{Int,Vector{Int}}, finalorder::Dict{Int,Int})
    RCM_adjacency = Dict(finalorder[i] => [finalorder[j] for j in c] for (i,c) in adjacency)
    return RCM_adjacency
end
