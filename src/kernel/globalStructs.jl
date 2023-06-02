#
# Space dimensions
#
abstract type AbstractSpaceDimensions end
struct NSD_1D <: AbstractSpaceDimensions end
struct NSD_2D <: AbstractSpaceDimensions end
struct NSD_3D <: AbstractSpaceDimensions end

Base.@kwdef mutable struct St_SolutionVars{TFloat <: AbstractFloat}

    qnp1 = Array{TFloat}(undef, 0, 0)       # qⁿ⁺¹
    qn   = Array{TFloat}(undef, 0, 0)       # qⁿ
    qnm1 = Array{TFloat}(undef, 0, 0)       # qⁿ⁻¹
    qnm2 = Array{TFloat}(undef, 0, 0)       # qⁿ⁻²
    qnm3 = Array{TFloat}(undef, 0, 0)       # qⁿ⁻³
    qe   = Array{TFloat}(undef, 0, 0)       # qexact    
    zb   = Array{TFloat}(undef, 0, 0)       # zb #shallow water moving bathymetry
    qnel = Array{TFloat}(undef, 0, 0, 0, 0) # Fⁿ[ngl,ngl,ngl,neqs]
    F    = Array{TFloat}(undef, 0, 0, 0, 0) # Gⁿ[ngl,ngl,ngl,neqs]
    G    = Array{TFloat}(undef, 0, 0, 0, 0) # Hⁿ[ngl,ngl,ngl,neqs]
    H    = Array{TFloat}(undef, 0, 0, 0, 0) # Sⁿ[ngl,ngl,ngl,neqs]
    S    = Array{TFloat}(undef, 0, 0, 0, 0) # qelⁿ[ngl,ngl,ngl,neqs]
    μ    = Array{TFloat}(undef, 0)          # μ (dynamic viscosity)
    neqs = UInt8(1)
    qvars= Array{String}(undef, neqs)
end

Base.@kwdef mutable struct St_PostProcessVars{TFloat <: AbstractFloat}

    qpost = Array{TFloat}(undef, 0, 0)
    
end

function allocate_post_process_vars(nelem, npoin, ngl, TFloat; neqs)

    qpost = St_SolutionVars{TFloat}(μ = zeros(npoin, neqs))    

    return qpost
end

function define_q(SD::NSD_1D, nelem, npoin, ngl, TFloat; neqs=1)

    q = St_SolutionVars{TFloat}(neqs=neqs,
                                qn   = zeros(npoin, neqs), # qn
                                qnm1 = zeros(npoin, neqs), # qⁿ
                                qnm2 = zeros(npoin, neqs), # qⁿ
                                qe   = zeros(npoin, neqs), # qexact
                                μ    = zeros(nelem)) # μ
    
    return q
end

function define_q(SD::NSD_2D, nelem, npoin, ngl, TFloat; neqs=1)

    q = St_SolutionVars{TFloat}(neqs=neqs,
                                qn   = zeros(npoin, neqs), # qn
                                qnm1 = zeros(npoin, neqs), # qⁿ
                                qnm2 = zeros(npoin, neqs), # qⁿ
                                qe   = zeros(npoin, neqs), # qexact
                                μ    = zeros(nelem)) # μ
    
    return q
end
