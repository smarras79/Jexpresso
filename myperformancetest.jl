function barw_array( xnodes )
    n = length(xnodes)
    barw = Array{Float64}(undef,n)
    for j ∈ 1:n
        barw[j] = 1.0
        for k ∈ 1:n
            if k != j
                barw[j] = barw[j] / (xnodes[j] - xnodes[k])
            end
        end
    end
    return sum(barw)
end

function run_barw_array()
    xnodes = LinRange(0,1,4)

    for i ∈ 1:100000
        u = barw_array( xnodes )
    end
end

function barw_workvector( barw, xnodes )
    n = length(xnodes)
    for j ∈ 1:n
        barw[j] = 1.0
        for k ∈ 1:n
            if k != j
                barw[j] = barw[j] / (xnodes[j] - xnodes[k])
            end
        end
    end
    return sum(barw)
end

function run_barw_workvector()
    xnodes = LinRange(0,1,4)

    barw = Array{Float64}(undef,length(xnodes))
    for i ∈ 1:100000
        u = barw_workvector( barw, xnodes )
    end

end

@btime run_barw_array()
@btime run_barw_workvector()
#############################

using StaticArrays

function xdiff(xnodes,i)
    xd = 1.0
    for k ∈ eachindex(xnodes)
        if k != i
            xd = xd / (xnodes[i] - xnodes[k])
        end
    end
    return xd
end

function barw_sarray_constructor( xnodes )
    n = length(xnodes)
    barw = SVector{n,Float64}( ( xdiff(xnodes,i) for i in 1:n )... )
    s = sum( barw )
    return s
end


function run_barw_sarray_constructor()
    xnodes = LinRange(0,1,4)

    for i ∈ 1:100000
        u = barw_sarray_constructor( xnodes )
    end

end

@btime run_barw_sarray_constructor()


using Revise
using BenchmarkTools
function barw_array( xnodes, ::Val{n} ) where n
    barw = zeros(MVector{n,Float64})
    for j ∈ 1:n
        barw[j] = 1.0
        for k ∈ 1:n
            if k != j
                barw[j] = barw[j] / (xnodes[j] - xnodes[k])
            end
        end
    end
    return sum(barw)
end

function run_barw_array()
    xnodes = LinRange(0,1,4)
    n = Val(length(xnodes))
    for i ∈ 1:100000
        u = barw_array( xnodes, n)
    end
end

@btime run_barw_array()
