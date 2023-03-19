#
# Space dimensions
#
abstract type AbstractSpaceDimensions end
struct NSD_1D <: AbstractSpaceDimensions end
struct NSD_2D <: AbstractSpaceDimensions end
struct NSD_3D <: AbstractSpaceDimensions end

Base.@kwdef mutable struct St_SolutionVars{TFloat <: AbstractFloat}

    qnp1 = Array{TFloat}(undef, 0, 0)       #1 qⁿ⁺¹
    qn   = Array{TFloat}(undef, 0, 0)       #2 qⁿ
    qnm1 = Array{TFloat}(undef, 0, 0)       #3 qⁿ⁻¹
    qnm2 = Array{TFloat}(undef, 0, 0)       #4 qⁿ⁻²
    qnm3 = Array{TFloat}(undef, 0, 0)       #5 qⁿ⁻³
    qe   = Array{TFloat}(undef, 0, 0)       #6 qexact    
    qnel = Array{TFloat}(undef, 0, 0, 0, 0) #7 qelⁿ[ngl,ngl,ngl,nelem]
    F    = Array{TFloat}(undef, 0, 0, 0, 0) #8  Fⁿ
    G    = Array{TFloat}(undef, 0, 0, 0, 0) #9  Gⁿ
    H    = Array{TFloat}(undef, 0, 0, 0, 0) #10 Hⁿ
end

"""
    allocate_q(nelem, npoin, ngl, neqs)

TBW
"""
function allocate_q(nelem, npoin, ngl, neqs, TFloat)
    
    q = St_SolutionVars{TFloat}(zeros(1, 1),               # qn+1
                                zeros(1, 1),               # qn
                                zeros(1, 1),               # qn-1
                                zeros(1, 1),               # qn-2
                                zeros(1, 1),               # qn-3
                                zeros(1, 1),               # qe
                                zeros(1, 1, 1, 1),  # qelⁿ[ngl,ngl,ngl,nelem]
                                zeros(1, 1, 1, 1),  # Fⁿ
                                zeros(1, 1, 1, 1),  # Gⁿ
                                zeros(1, 1, 1, 1))  # Hⁿ
    
    return q
end

function define_q(SD::NSD_1D, nelem, npoin, ngl, neqs, TFloat)

    q = St_SolutionVars{TFloat}(qn = zeros(npoin, neqs),    # qn
                                F  = zeros(ngl, nelem))  # Fⁿ
    
    return q
end

function define_q(SD::NSD_2D, nelem, npoin, ngl, neqs, TFloat)
    
    q = St_SolutionVars{TFloat}(qn = zeros(npoin, neqs), # qⁿ
                                F  = zeros(ngl, nelem),  # Fⁿ
                                G  = zeros(ngl, nelem))  # Gⁿ
    
    return q
end


function define_q(SD::NSD_3D, nelem, npoin, ngl, neqs, TFloat)
    
    q = St_SolutionVars{TFloat}(qn = zeros(npoin, neqs), # qⁿ
                                F  = zeros(ngl, nelem, neqs),  # Fⁿ
                                G  = zeros(ngl, nelem, neqs),  # Gⁿ
                                H  = zeros(ngl, nelem, neqs))  # Hⁿ
    
    return q
end
