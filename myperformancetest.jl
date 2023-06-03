using BenchmarkTools
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
#@info "NAIVE "
#@btime run_barw_array()
####

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
#@info "work vector without static arrays"
#@btime run_barw_workvector()
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

#@info "work vector WITH static arrays"
#@btime run_barw_sarray_constructor()


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
#@info "work vector WITH static arrays"
#@btime run_barw_array()

function mytest(n)
    a = zeros(Int64, n)
    b = zeros(Int64, n)
    sum = 0.0
    for i=1:n
        a[i] = i
        b[i] = i*2
    end
    return a, b
end

function mytest!(a::Array,b::Array)
    n=length(a)
    for i=1:n
        a[i] = i
        b[i] = i*2
    end
    
end


using BenchmarkTools
using StaticArrays

Base.@kwdef mutable struct St_SolutionVars{TFloat <: AbstractFloat}

    a = Array{TFloat}(undef, 0, 0, 0)
    b = Array{TFloat}(undef, 0, 0, 0)
    c = Array{TFloat}(undef, 0, 0, 0)
    d = Array{TFloat}(undef, 0, 0, 0)
    
end

function define_q(nvar, nx)
    
    q = St_SolutionVars{Float64}(a = zeros(nx,nx, nvar),
                                 b = zeros(nx,nx, nvar),
                                 c = zeros(nx,nx, nvar),
                                 d = zeros(nx,nx, nvar))
    return q
end

function user_flux!(d,a0,b0,c0, i, nvar)
    #for n=1:nvar
    #d[i,i,n] = (a0[i,i,n])#+b0[i-1,i-1,n]) # + 129.0*(c0[i,i,n]*b0[i+1,i+1,n]) + 300.0*(a0[i-1,i-1,n]*c0[i+1,i+1,n])
    d[i,i,1] = 1.0#(a0[n])#+b0[i-1,i-1,n]) # + 129.0*(c0[i,i,n]*b0[i+1,i+1,n]) + 300.0*(a0[i-1,i-1,n]*c0[i+1,i+1,n])
    d[i,i,2] = 1.0#(a0[n])#+b0[i-1,i-1,n]) # + 129.0*(c0[i,i,n]*b0[i+1,i+1,n]) + 300.0*(a0[i-1,i-1,n]*c0[i+1,i+1,n])
    d[i,i,3] = 1.0#(a0[n])#+b0[i-1,i-1,n]) # + 129.0*(c0[i,i,n]*b0[i+1,i+1,n]) + 300.0*(a0[i-1,i-1,n]*c0[i+1,i+1,n])
    d[i,i,nvar] = 1.0#(a0[n])#+b0[i-1,i-1,n]) # + 129.0*(c0[i,i,n]*b0[i+1,i+1,n]) + 300.0*(a0[i-1,i-1,n]*c0[i+1,i+1,n])
    #d[i,i,n] *= 5
    #end
end

function _build_rhs(a, b, c, d, e, nx, nvar)   
    
    for i=2:nx-1
        user_flux!(d, a, b, c, i, nvar)
        #@info "first"
        #@info d[2,2,end]
    end
    
    #=for ieq = 1:nvar
        for j=1:nx, i=1:nx
            dFdξ = 0.0
            for k = 1:nx
                dFdξ += d[k,j,ieq]
            end
        end
    end=#
    
    #for i=2:nx-1
    #@info "second"
    #@info d[i,i,end]
    #end
    #@info "dFdξ" dFdξ
end

nx = 5
nvar = 4
#q=define_q(nvar, nx)
a,b,c,d, e = [rand(nx, nx, nvar) for _=1:5]
#@btime _build_rhs($q.a,$q.b,$q.c,$q.d,$nx,$nvar)
@btime _build_rhs($a,$b,$c,$d,$e,$nx,$nvar)
