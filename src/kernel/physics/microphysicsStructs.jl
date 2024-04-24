Base.@kwdef mutable struct St_Microphysics{TFloat <: AbstractFloat}
    
    prod::Array{TFloat,1}  = zeros(TFloat, 0)
    rhs::Array{TFloat,1}   = zeros(TFloat, 0)
    vt::Array{TFloat,1}    = zeros(TFloat, 0)
    prodk::Array{TFloat,1} = zeros(TFloat, 0)
    vtden::Array{TFloat,1} = zeros(TFloat, 0)
    rdzw::Array{TFloat,1}  = zeros(TFloat, 0)
    rdzk::Array{TFloat,1}  = zeros(TFloat, 0)
    ρk::Array{TFloat,1}    = zeros(TFloat, 0)
    pp::Array{TFloat,1}    = zeros(TFloat, 0)  

end

function allocate_mp(SD, nelem, npoin, ngl, qvars, TFloat; neqs=1)

    mp = St_Microphysics{TFloat}(neqs=neqs,
                                 prod  = zeros(TFloat, npoin),
                                 rhs   = zeros(TFloat, npoin),
                                 vt    = zeros(TFloat, npoin),
                                 prodk = zeros(TFloat, npoin),
                                 vtden = zeros(TFloat, npoin),
                                 rdzw  = zeros(TFloat, npoin),
                                 rdzk  = zeros(TFloat, npoin),
                                 ρk    = zeros(TFloat, npoin),
                                 pp    = zeros(TFloat, npoin))
    
    return mp
end
