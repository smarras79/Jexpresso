mutable struct St_SolutionVars{TFloat}

    qnp1::Array{TFloat} #1 qⁿ⁺¹
    qn::Array{TFloat}   #2 qⁿ
    qnm1::Array{TFloat} #3 qⁿ⁻¹
    qnm2::Array{TFloat} #4 qⁿ⁻²
    qnm3::Array{TFloat} #5 qⁿ⁻³
    qe::Array{TFloat}   #6 qexact    
    qnel::Array{TFloat} #7 qelⁿ[ngl,ngl,ngl,nelem]

    F::Array{TFloat}   #8  Fⁿ
    G::Array{TFloat}   #9  Gⁿ
    H::Array{TFloat}   #10 Hⁿ
end
