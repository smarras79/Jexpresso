#
# Solution arrays
#
Base.@kwdef mutable struct St_SolutionVars{TFloat <: AbstractFloat, dim1}
    
    qnp1::Array{TFloat,dim1}  # qⁿ⁺¹
    qn::Array{TFloat,dim1}         # qⁿ
    qq::Array{TFloat,dim1}      # qⁿ
    qnm1::Array{TFloat,dim1}    # qⁿ⁻¹
    qnm2::Array{TFloat,dim1}       # qⁿ⁻²
    qnm3::Array{TFloat,dim1}       # qⁿ⁻³
    qe::Array{TFloat,dim1}       # qexact    
    zb::Array{TFloat,dim1}        # zb #shallow water moving bathymetry
    qvars = Array{Union{Nothing, String}}(nothing, 0)
    neqs = UInt8(1)
end
function St_SolutionVars(T::Type, SD::AbstractSpaceDimensions, dim1::Int, dims1, qvars, neqs)
    
    St_SolutionVars{T, length(dims1)}(qnp1  = zeros(T, dims1),
                                      qn    = zeros(T, dims1),
                                      qq    = zeros(T, dims1),
                                      qnm1  = zeros(T, dims1),
                                      qnm2  = zeros(T, dims1),
                                      qnm3  = zeros(T, dims1),
                                      qe    = zeros(T, dims1),
                                      zb    = zeros(T, dims1),
                                      qvars = qvars,
                                      neqs  = neqs)

    
end
function define_q(SD, nelem, npoin, ngl, qvars, TFloat; neqs=1)

    neqs = length(qvars)
    dims1 = (npoin, neqs+1)
    
    q = St_SolutionVars(TFloat, SD, length(dims1), dims1, qvars, neqs)
    
    return q
end

#
# rhs
#
Base.@kwdef mutable struct St_rhs{TFloat <: AbstractFloat, dim1, dim2}
    
    RHS::Array{TFloat,dim1}          
    RHS_visc::Array{TFloat,dim1}     
    rhs_el::Array{TFloat,dim2}     
    rhs_diff_el::Array{TFloat, dim2} 
    rhs_diffξ_el::Array{TFloat, dim2}
    rhs_diffη_el::Array{TFloat, dim2}
    rhs_diffζ_el::Array{TFloat, dim2}
    
end
function St_rhs(T::Type, SD::AbstractSpaceDimensions, dim1::Int, dim2::Int, dims1, dims2)
    
    St_rhs{T, length(dims1), length(dims2)}(RHS          = zeros(T, dims1),
                                            RHS_visc     = zeros(T, dims1),
                                            rhs_el       = zeros(T, dims2),
                                            rhs_diff_el  = zeros(T, dims2),
                                            rhs_diffξ_el = zeros(T, dims2),
                                            rhs_diffη_el = zeros(T, dims2),
                                            rhs_diffζ_el = zeros(T, dims2))
end

function allocate_rhs(SD::NSD_1D, nelem, npoin, ngl, TFloat; neqs=1)

    dims1 = (npoin, neqs)
    dims2 = (nelem, ngl, neqs) 

    rhs = St_rhs(TFloat, SD, length(dims1), length(dims2), dims1, dims2)
                 
    return rhs
end


function allocate_rhs(SD::NSD_2D, nelem, npoin, ngl, TFloat; neqs=1)

    dims1 = (npoin, neqs)
    dims2 = (nelem, ngl, ngl, neqs) 

    rhs = St_rhs(TFloat, SD, length(dims1), length(dims2), dims1, dims2)
               
    return rhs
    
end

function allocate_rhs(SD::NSD_3D, nelem, npoin, ngl, TFloat; neqs=1)

    dims1 = (npoin, neqs)
    dims2 = (nelem, ngl, ngl, ngl, neqs) 

    rhs = St_rhs(TFloat, SD, length(dims1), length(dims2), dims1, dims2)
    
    return rhs
    
end

#
# Fluxes
#
Base.@kwdef mutable struct St_fluxes{TFloat <: AbstractFloat, dim1, dim2}
    
    F::Array{TFloat,dim1}
    G::Array{TFloat,dim1}
    H::Array{TFloat,dim1}
    S::Array{TFloat,dim1}
    uprimitive::Array{TFloat, dim2} 
    
end
function St_fluxes(T::Type, SD::AbstractSpaceDimensions, dim1::Int, dim2::Int, dims1, dims2)
    
    St_fluxes{T, length(dims1), length(dims2)}(F = zeros(T, dims1),
                                               G = zeros(T, dims1),
                                               H = zeros(T, dims1),
                                               S = zeros(T, dims1),
                                               uprimitive = zeros(T, dims2))
end

function allocate_fluxes(SD::NSD_1D, nelem, npoin, ngl, TFloat; neqs=1)

    dims1 = (ngl, neqs)
    dims2 = (ngl, neqs+1) 

    fluxes = St_fluxes(TFloat, SD, length(dims1), length(dims2), dims1, dims2)
    
    return fluxes
end

function allocate_fluxes(SD::NSD_2D, nelem, npoin, ngl, TFloat; neqs=1)

    dims1 = (ngl, ngl, neqs)
    dims2 = (ngl, ngl, neqs+1) 

    fluxes = St_fluxes(TFloat, SD, length(dims1), length(dims2), dims1, dims2)
    
    return fluxes
end

function allocate_fluxes(SD::NSD_3D, nelem, npoin, ngl, TFloat; neqs=1)

    dims1 = (ngl, ngl, ngl, neqs)
    dims2 = (ngl, ngl, ngl, neqs+1) 

    fluxes = St_fluxes(TFloat, SD, length(dims1), length(dims2), dims1, dims2)
    
    return fluxes
end

#
# rhs Laguerre
#
Base.@kwdef mutable struct St_rhs_lag{TFloat <: AbstractFloat, dim1, dim2}
    
    RHS_lag::Array{TFloat,dim1}          
    RHS_visc_lag::Array{TFloat,dim1}     
    rhs_el_lag::Array{TFloat,dim2}     
    rhs_diff_el_lag::Array{TFloat, dim2} 
    rhs_diffξ_el_lag::Array{TFloat, dim2}
    rhs_diffη_el_lag::Array{TFloat, dim2}
    rhs_diffζ_el_lag::Array{TFloat, dim2}
    
end
function St_rhs_lag(T::Type, SD::AbstractSpaceDimensions, dim1::Int, dim2::Int, dims1, dims2)
    
    St_rhs_lag{T, length(dims1), length(dims2)}(RHS_lag          = zeros(T, dims1),
                                                RHS_visc_lag     = zeros(T, dims1),
                                                rhs_el_lag       = zeros(T, dims2),
                                                rhs_diff_el_lag  = zeros(T, dims2),
                                                rhs_diffξ_el_lag = zeros(T, dims2),
                                                rhs_diffη_el_lag = zeros(T, dims2),
                                                rhs_diffζ_el_lag = zeros(T, dims2))
end

function allocate_rhs_lag(SD::NSD_1D, nelem_semi_inf, npoin, ngr, TFloat; neqs=1)

    dims1 = (npoin, neqs)
    dims2 = (nelem_semi_inf, ngr, neqs) 

    rhs_lag = St_rhs_lag(TFloat, SD, length(dims1), length(dims2), dims1, dims2)
                 
    return rhs_lag
end


function allocate_rhs_lag(SD::NSD_2D, nelem_semi_inf, npoin, ngl, ngr, TFloat; neqs=1)

    dims1 = (npoin, neqs)
    dims2 = (nelem_semi_inf, ngl, ngr, neqs) 

    rhs_lag = St_rhs_lag(TFloat, SD, length(dims1), length(dims2), dims1, dims2)
               
    return rhs_lag
    
end

function allocate_rhs_lag(SD::NSD_3D, nelem_semi_inf, npoin, ngl, ngr, TFloat; neqs=1)

     error(" src/kernel/infrastructore/params_setup.jl: 3D Laguerre arrays not coded yet!")
    
end

#
# Fluxes Laguerre
#
Base.@kwdef mutable struct St_fluxes_lag{TFloat <: AbstractFloat, dim1, dim2}
    
    F_lag::Array{TFloat,dim1}
    G_lag::Array{TFloat,dim1}
    H_lag::Array{TFloat,dim1}
    S_lag::Array{TFloat,dim1}
    uprimitive_lag::Array{TFloat, dim2} 
    
end
function St_fluxes_lag(T::Type, SD::AbstractSpaceDimensions, dim1::Int, dim2::Int, dims1, dims2)
    
    St_fluxes_lag{T, length(dims1), length(dims2)}(F_lag = zeros(T, dims1),
                                                   G_lag = zeros(T, dims1),
                                                   H_lag = zeros(T, dims1),
                                                   S_lag = zeros(T, dims1),
                                                   uprimitive_lag = zeros(T, dims2))
end

function allocate_fluxes_lag(SD::NSD_1D, nelem, npoin, ngr, TFloat; neqs=1)

    dims1 = (ngr, neqs)
    dims2 = (ngr, neqs+1) 

    fluxes_lag = St_fluxes_lag(TFloat, SD, length(dims1), length(dims2), dims1, dims2)
    
    return fluxes_lag
end

function allocate_fluxes_lag(SD::NSD_2D, nelem, npoin, ngl, ngr, TFloat; neqs=1)

    dims1 = (ngl, ngr, neqs)
    dims2 = (ngl, ngr, neqs+1) 

    fluxes_lag = St_fluxes_lag(TFloat, SD, length(dims1), length(dims2), dims1, dims2)
    
    return fluxes_lag
end

function allocate_fluxes_lag(SD::NSD_3D, nelem, npoin, ngl, TFloat; neqs=1)
    nothing
end


#
# Moist variables
#
Base.@kwdef mutable struct St_MoistVars{TFloat <: AbstractFloat, dim1}

    rainnc::Array{TFloat, dim1}
    rainncv::Array{TFloat, dim1}
    vt::Array{TFloat, dim1}
    prod::Array{TFloat, dim1}
    prodk::Array{TFloat, dim1}
    vtden::Array{TFloat, dim1}
    rdzk::Array{TFloat, dim1}
    ρk::Array{TFloat, dim1}
    
end
function St_MoistVars(T::Type, SD::AbstractSpaceDimensions, dim1::Int, dims1)
    
    St_MoistVars{T, length(dims1), length(dims2)}(zeros(T, dims1),
                                                  zeros(T, dims1),
                                                  zeros(T, dims1),
                                                  zeros(T, dims1),
                                                  zeros(T, dims1),
                                                  zeros(T, dims1),
                                                  zeros(T, dims1),
                                                  zeros(T, dims1))
end


function allocate_MoistVars(npoin, TFloat)
    
    dims1 = (npoin)

    moistvars = St_MoistVars(TFloat, SD, length(dims1), dims1)
    
    return moistvars
end
