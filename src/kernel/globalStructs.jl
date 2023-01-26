mutable struct St_SolutionVars{TFloat}

    qnp1::Array{TFloat} #qⁿ⁺¹
    qn::Array{TFloat}   #qⁿ
    qnm1::Array{TFloat} #qⁿ⁻¹
    qnm2::Array{TFloat} #qⁿ⁻²
    qnm3::Array{TFloat} #qⁿ⁻³
    qe::Array{TFloat}   #qexact    
    qnel::Array{TFloat} #qⁿ[ngl,ngl,ngl,nelem]

    F::Array{TFloat}   #Fⁿ
    G::Array{TFloat}   #Gⁿ
    H::Array{TFloat}   #Hⁿ
end
