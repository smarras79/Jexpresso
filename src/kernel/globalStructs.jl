Base.@kwdef mutable struct St_SolutionVars{TFloat <: AbstractFloat, backend}


    qnp1  = KernelAbstractions.zeros(backend,  TFloat, 0, 0)       # qⁿ⁺¹
    qn    = KernelAbstractions.zeros(backend,  TFloat, 0, 0)       # qⁿ
    qq    = KernelAbstractions.zeros(backend,  TFloat, 0, 0)       # qⁿ
    qnm1  = KernelAbstractions.zeros(backend,  TFloat, 0, 0)       # qⁿ⁻¹
    qnm2  = KernelAbstractions.zeros(backend,  TFloat, 0, 0)       # qⁿ⁻²
    qnm3  = KernelAbstractions.zeros(backend,  TFloat, 0, 0)       # qⁿ⁻³
    qe    = KernelAbstractions.zeros(backend,  TFloat, 0, 0)       # qexact    
    zb    = KernelAbstractions.zeros(backend,  TFloat, 0, 0)       # zb #shallow water moving bathymetry
    qnel  = KernelAbstractions.zeros(backend,  TFloat, 0, 0, 0, 0)
    μ     = KernelAbstractions.zeros(backend,  TFloat, 0)  
    neqs = UInt8(1)
    qvars = Array{Union{Nothing, String}}(nothing, 0)
end

Base.@kwdef mutable struct St_PostProcessVars{TFloat <: AbstractFloat}

    qpost = Array{TFloat}(undef, 0, 0)
    
end

function allocate_post_process_vars(nelem, npoin, ngl, TFloat; neqs)

    qpost = St_SolutionVars{TFloat}(μ = KernelAbstractions.zeros(backend,  npoin, neqs))    

    return qpost
end


function define_q(SD, nelem, npoin, ngl, qvars, TFloat, backend; neqs=1)

    q = St_SolutionVars{TFloat, backend}(neqs=neqs,
                                qn   = KernelAbstractions.zeros(backend,  TFloat, Int64(npoin), neqs+1), # qn
                                qnm1 = KernelAbstractions.zeros(backend,  TFloat, Int64(npoin), neqs+1), # qⁿ
                                qnm2 = KernelAbstractions.zeros(backend,  TFloat, Int64(npoin), neqs+1), # qⁿ
                                qe   = KernelAbstractions.zeros(backend,  TFloat, Int64(npoin), neqs+1), # qexact
                                qvars=qvars) # μ
    @info typeof(q.qe)
    @info backend
    return q
end
