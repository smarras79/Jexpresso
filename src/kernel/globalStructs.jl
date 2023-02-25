Base.@kwdef mutable struct St_SolutionVars{TFloat}

    #=qnp1::Array{TFloat} #1 qⁿ⁺¹
    qn::Array{TFloat}   #2 qⁿ
    qnm1::Array{TFloat} #3 qⁿ⁻¹
    qnm2::Array{TFloat} #4 qⁿ⁻²
    qnm3::Array{TFloat} #5 qⁿ⁻³
    qe::Array{TFloat}   #6 qexact    
    qnel::Array{TFloat} #7 qelⁿ[ngl,ngl,ngl,nelem]

    F::Array{TFloat}   #8  Fⁿ
    G::Array{TFloat}   #9  Gⁿ
    H::Array{TFloat}   #10 Hⁿ
    =#

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
function allocate_q(nelem, npoin, ngl, neqs)
    
    q = St_SolutionVars{TFloat}(zeros(npoin, 3),               # qn+1
                                zeros(npoin, 3),               # qn
                                zeros(npoin, 3),               # qn-1
                                zeros(npoin, 3),               # qn-2
                                zeros(npoin, 3),               # qn-3
                                zeros(npoin, 3),               # qe
                                zeros(ngl, ngl, nelem, neqs),  # qelⁿ[ngl,ngl,ngl,nelem]
                                zeros(ngl, ngl, nelem, neqs),  # Fⁿ
                                zeros(ngl, ngl, nelem, neqs),  # Gⁿ
                                zeros(ngl, ngl, nelem, neqs))  # Hⁿ
    
    return q
end