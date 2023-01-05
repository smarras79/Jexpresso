# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/NodeNumbering.jl/blob/master/LICENSE

"""
    calculate_bandwidth(G::Dict{Int, Vector{Int}})

Calculate bandwidth of graph ``G``. The bandwidth of the matrix is number ``k``
such that ``A_{ij} = 0`` if ``|i-j| > k``.

# Examples

```jldoctest
julia> G = Dict(
           1 => [3, 8, 9],
           2 => [3, 8, 7],
           3 => [1, 2],
           4 => [8, 9],
           5 => [7, 8],
           6 => [2, 7],
           7 => [5, 2, 6],
           8 => [1, 2, 4, 5],
           9 => [1, 4]);

julia> calculate_bandwidth(G)
17

```

# References

* Wikipedia contributors. "Band matrix." Wikipedia, The Free Encyclopedia. Wikipedia, The Free Encyclopedia, 28 Apr. 2017. Web. 15 Jul. 2017. https://en.wikipedia.org/wiki/Band_matrix.
"""
function calculate_bandwidth(G::Dict{Int, Vector{Int}})
    return 2*maximum([maximum(j-g) for (j,g) in G]) + 1
end
