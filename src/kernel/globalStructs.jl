Base.@kwdef mutable struct St_SolutionVars{TFloat <: AbstractFloat}


    qnp1::Array{TFloat,2} = zeros(TFloat, 0, 0)       # qⁿ⁺¹
    qn::Array{TFloat,2}   = zeros(TFloat, 0, 0)       # qⁿ
    qq::Array{TFloat,2}   = zeros(TFloat, 0, 0)       # qⁿ
    qnm1::Array{TFloat,2} = zeros(TFloat, 0, 0)       # qⁿ⁻¹
    qnm2::Array{TFloat,2} = zeros(TFloat, 0, 0)       # qⁿ⁻²
    qnm3::Array{TFloat,2} = zeros(TFloat, 0, 0)       # qⁿ⁻³
    qe::Array{TFloat,2}   = zeros(TFloat, 0, 0)       # qexact    
    zb::Array{TFloat,2}   = zeros(TFloat, 0, 0)       # zb #shallow water moving bathymetry
    qnel::Array{TFloat,4} = zeros(TFloat, 0, 0, 0, 0)
    μ::Array{TFloat,1}    = zeros(TFloat, 0)  
    neqs = UInt8(1)
    qvars = Array{Union{Nothing, String}}(nothing, 0)
end

Base.@kwdef mutable struct St_PostProcessVars{TFloat <: AbstractFloat}

    qpost = Array{TFloat}(undef, 0, 0)
    
end

function allocate_post_process_vars(nelem, npoin, ngl, TFloat; neqs)

    qpost = St_SolutionVars{TFloat}(μ = zeros(npoin, neqs))    

    return qpost
end


function define_q(SD, nelem, npoin, ngl, qvars, TFloat; neqs=1)

    q = St_SolutionVars{TFloat}(neqs=neqs,
                                qn   = zeros(TFloat, npoin, neqs+1), # qn
                                qnm1 = zeros(TFloat, npoin, neqs+1), # qⁿ
                                qnm2 = zeros(TFloat, npoin, neqs+1), # qⁿ
                                qe   = zeros(TFloat, npoin, neqs+1), # qexact
                                qvars=qvars) # μ
    return q
end
