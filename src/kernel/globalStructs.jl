#
# Solution arrays
#
Base.@kwdef mutable struct St_SolutionVars{T <: AbstractFloat, dim1}
    
    qnp1::Array{T,dim1}  # qⁿ⁺¹
    qn::Array{T,dim1}         # qⁿ
    qq::Array{T,dim1}      # qⁿ
    qnm1::Array{T,dim1}    # qⁿ⁻¹
    qnm2::Array{T,dim1}       # qⁿ⁻²
    qnm3::Array{T,dim1}       # qⁿ⁻³
    qe::Array{T,dim1}       # qexact    
    zb::Array{T,dim1}        # zb #shallow water moving bathymetry
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
function define_q(SD, nelem, npoin, ngl, qvars, T; neqs=1)

    neqs = length(qvars)
    dims1 = (npoin, neqs+1)
    
    q = St_SolutionVars(T, SD, length(dims1), dims1, qvars, neqs)
    
    return q
end

#
# rhs
#
Base.@kwdef mutable struct St_rhs{T <: AbstractFloat, dim1, dim2}
    
    RHS::Array{T,dim1}          
    RHS_visc::Array{T,dim1}     
    rhs_el::Array{T,dim2}     
    rhs_diff_el::Array{T, dim2} 
    rhs_diffξ_el::Array{T, dim2}
    rhs_diffη_el::Array{T, dim2}
    rhs_diffζ_el::Array{T, dim2}
    
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

function allocate_rhs(SD::NSD_1D, nelem, npoin, ngl, T; neqs=1)

    dims1 = (npoin, neqs)
    dims2 = (nelem, ngl, neqs) 

    rhs = St_rhs(T, SD, length(dims1), length(dims2), dims1, dims2)
                 
    return rhs
end


function allocate_rhs(SD::NSD_2D, nelem, npoin, ngl, T; neqs=1)

    dims1 = (npoin, neqs)
    dims2 = (nelem, ngl, ngl, neqs) 

    rhs = St_rhs(T, SD, length(dims1), length(dims2), dims1, dims2)
               
    return rhs
    
end

function allocate_rhs(SD::NSD_3D, nelem, npoin, ngl, T; neqs=1)

    dims1 = (npoin, neqs)
    dims2 = (nelem, ngl, ngl, ngl, neqs) 

    rhs = St_rhs(T, SD, length(dims1), length(dims2), dims1, dims2)
    
    return rhs
    
end

#
# Fluxes
#
Base.@kwdef mutable struct St_fluxes{T <: AbstractFloat, dim1, dim2}
    
    F::Array{T,dim1}
    G::Array{T,dim1}
    H::Array{T,dim1}
    S::Array{T,dim1}
    uprimitive::Array{T, dim2} 
    
end
function St_fluxes(T::Type, SD::AbstractSpaceDimensions, dim1::Int, dim2::Int, dims1, dims2)
    
    St_fluxes{T, length(dims1), length(dims2)}(F = zeros(T, dims1),
                                               G = zeros(T, dims1),
                                               H = zeros(T, dims1),
                                               S = zeros(T, dims1),
                                               uprimitive = zeros(T, dims2))
end

function allocate_fluxes(SD::NSD_1D, nelem, npoin, ngl, T; neqs=1)

    dims1 = (ngl, neqs)
    dims2 = (ngl, neqs+1) 

    fluxes = St_fluxes(T, SD, length(dims1), length(dims2), dims1, dims2)
    
    return fluxes
end

function allocate_fluxes(SD::NSD_2D, nelem, npoin, ngl, T; neqs=1)

    dims1 = (ngl, ngl, neqs)
    dims2 = (ngl, ngl, neqs+1) 

    fluxes = St_fluxes(T, SD, length(dims1), length(dims2), dims1, dims2)
    
    return fluxes
end

function allocate_fluxes(SD::NSD_3D, nelem, npoin, ngl, T; neqs=1)

    dims1 = (ngl, ngl, ngl, neqs)
    dims2 = (ngl, ngl, ngl, neqs+1) 

    fluxes = St_fluxes(T, SD, length(dims1), length(dims2), dims1, dims2)
    
    return fluxes
end

#
# rhs Laguerre
#
Base.@kwdef mutable struct St_rhs_lag{T <: AbstractFloat, dim1, dim2}
    
    RHS_lag::Array{T,dim1}          
    RHS_visc_lag::Array{T,dim1}     
    rhs_el_lag::Array{T,dim2}     
    rhs_diff_el_lag::Array{T, dim2} 
    rhs_diffξ_el_lag::Array{T, dim2}
    rhs_diffη_el_lag::Array{T, dim2}
    rhs_diffζ_el_lag::Array{T, dim2}
    
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

function allocate_rhs_lag(SD::NSD_1D, nelem_semi_inf, npoin, ngl, ngr, T; neqs=1)

    dims1 = (npoin, neqs)
    dims2 = (nelem_semi_inf, ngr, neqs) 

    rhs_lag = St_rhs_lag(T, SD, length(dims1), length(dims2), dims1, dims2)
                 
    return rhs_lag
end


function allocate_rhs_lag(SD::NSD_2D, nelem_semi_inf, npoin, ngl, ngr, T; neqs=1)

    dims1 = (npoin, neqs)
    dims2 = (nelem_semi_inf, ngl, ngr, neqs) 

    rhs_lag = St_rhs_lag(T, SD, length(dims1), length(dims2), dims1, dims2)
               
    return rhs_lag
    
end

function allocate_rhs_lag(SD::NSD_3D, nelem_semi_inf, npoin, ngl, ngr, T; neqs=1)

     error(" src/kernel/infrastructore/params_setup.jl: 3D Laguerre arrays not coded yet!")
    
end

#
# Fluxes Laguerre
#
Base.@kwdef mutable struct St_fluxes_lag{T <: AbstractFloat, dim1, dim2}
    
    F_lag::Array{T,dim1}
    G_lag::Array{T,dim1}
    H_lag::Array{T,dim1}
    S_lag::Array{T,dim1}
    uprimitive_lag::Array{T, dim2} 
    
end
function St_fluxes_lag(T::Type, SD::AbstractSpaceDimensions, dim1::Int, dim2::Int, dims1, dims2)
    
    St_fluxes_lag{T, length(dims1), length(dims2)}(F_lag = zeros(T, dims1),
                                                   G_lag = zeros(T, dims1),
                                                   H_lag = zeros(T, dims1),
                                                   S_lag = zeros(T, dims1),
                                                   uprimitive_lag = zeros(T, dims2))
end

function allocate_fluxes_lag(SD::NSD_1D, ngl, ngr, T; neqs=1)

    dims1 = (ngr, neqs)
    dims2 = (ngr, neqs+1) 

    fluxes_lag = St_fluxes_lag(T, SD, length(dims1), length(dims2), dims1, dims2)
    
    return fluxes_lag
end

function allocate_fluxes_lag(SD::NSD_2D, ngl, ngr, T; neqs=1)

    dims1 = (ngl, ngr, neqs)
    dims2 = (ngl, ngr, neqs+1) 

    fluxes_lag = St_fluxes_lag(T, SD, length(dims1), length(dims2), dims1, dims2)
    
    return fluxes_lag
end

function allocate_fluxes_lag(SD::NSD_3D, ngl, ngr, T; neqs=1)
    nothing
end

#
# Filter:
#
Base.@kwdef mutable struct St_filter{T <: AbstractFloat, dim1, dim2, dim3, dim4}
    
    q_t::Array{T,dim1}
    fqf::Array{T,dim1}
    q_ti::Array{T,dim2}
    b::Array{T,dim3}
    B::Array{T,dim4}   
    
end
function St_filter(T::Type, SD::AbstractSpaceDimensions, dim1::Int, dim2::Int, dim3::Int, dim4::Int,
                   dims1, dims2, dims3, dims4)
    
    St_filter{T,
              length(dims1), length(dims2),
              length(dims3), length(dims4)}(q_t  = zeros(T, dims1),
                                            fqf  = zeros(T, dims1),
                                            q_ti = zeros(T, dims2),
                                            b    = zeros(T, dims3),
                                            B    = zeros(T, dims4))
    
end

function allocate_filter(SD::NSD_1D, nelem, npoin, ngl, T; neqs=1, lfilter=false)

    if lfilter
        dims1 = (neqs, ngl)
        dims2 = (ngl)
        dims3 = (nelem, ngl, neqs)
        dims4 = (npoin, neqs)
    else
        dims1 = (1, 1)
        dims2 = (1)
        dims3 = (1, 1, 1)
        dims4 = (1, 1)
    end

    filter = St_filter(T, SD,
                       length(dims1), length(dims2), length(dims3), length(dims4),
                       dims1, dims2, dims3, dims4)
    
    return filter
end

function allocate_filter(SD::NSD_2D, nelem, npoin, ngl, T; neqs=1, lfilter=false)

    if lfilter
        dims1 = (neqs, ngl, ngl)
        dims2 = (ngl, ngl)
        dims3 = (nelem, ngl, ngl, neqs)
        dims4 = (npoin, neqs)
    else
        dims1 = (1, 1, 1)
        dims2 = (1, 1)
        dims3 = (1, 1, 1, 1)
        dims4 = (1, 1)
    end
        
    filter = St_filter(T, SD,
                       length(dims1), length(dims2), length(dims3), length(dims4),
                       dims1, dims2, dims3, dims4)
    return filter
        
end


function allocate_filter(SD::NSD_3D, nelem, npoin, ngl, T; neqs=1, lfilter=false)

    dims1 = (1, 1, 1, 1)
    dims2 = (1, 1, 1)
    dims3 = (1, 1, 1, 1, 1)
    dims4 = (1, 1)
    
    filter = St_filter(T, SD,
                       length(dims1), length(dims2), length(dims3), length(dims4),
                       dims1, dims2, dims3, dims4)
    if lfilter
        warning( " NOT IMPLEMENTED: 3D filter not implemented")
    end

    return filter
    
end


#
# Laguerre filter
#
Base.@kwdef mutable struct St_filter_lag{T <: AbstractFloat, dim1, dim2, dim3, dim4}
    
    q_t_lag::Array{T,dim1}
    fqf_lag::Array{T,dim1}
    q_ti_lag::Array{T,dim2}
    b_lag::Array{T,dim3}
    B_lag::Array{T,dim4}   
    
end
function St_filter_lag(T::Type, SD::AbstractSpaceDimensions, dim1::Int, dim2::Int, dim3::Int, dim4::Int,
    dims1, dims2, dims3, dims4)
    
    St_filter_lag{T,  
                  length(dims1), length(dims2),
                  length(dims3), length(dims4)}(q_t_lag  = zeros(T, dims1),
                                                fqf_lag  = zeros(T, dims1),
                                                q_ti_lag = zeros(T, dims2),
                                                b_lag    = zeros(T, dims3),
                                                B_lag    = zeros(T, dims4))
    
end

function allocate_filter_lag(SD::NSD_1D, nelem_semi_inf, npoin, ngl, ngr, T; neqs=1, lfilter=false)

    if lfilter
        dims1 = (neqs, ngr)
        dims2 = (ngr)
        dims3 = (nelem_semi_inf, ngr, neqs)
        dims4 = (npoin, neqs)
    else
        dims1 = (1, 1)
        dims2 = (1)
        dims3 = (1, 1, 1)
        dims4 = (1, 1)
    end

    filter_lag = St_filter_lag(T, SD,
                               length(dims1), length(dims2), length(dims3), length(dims4),
                               dims1, dims2, dims3, dims4)

    return filter_lag
end

function allocate_filter_lag(SD::NSD_2D, nelem_semi_inf, npoin, ngl, ngr, T; neqs=1, lfilter=false)

    if lfilter
        dims1 = (neqs, ngl, ngr)
        dims2 = (ngl, ngr)
        dims3 = (nelem_semi_inf, ngl, ngr, neqs)
        dims4 = (npoin, neqs)
    else
        dims1 = (1, 1, 1)
        dims2 = (1, 1)
        dims3 = (1, 1, 1, 1)
        dims4 = (1, 1)
    end

    filter_lag = St_filter_lag(T, SD,
                               length(dims1), length(dims2), length(dims3), length(dims4),
                               dims1, dims2, dims3, dims4)
    
    return filter_lag
end

function allocate_filter_lag(SD::NSD_3D, nelem_semi_inf, npoin, ngl, ngr, T; neqs=1, lfilter=false)
    
    dims1 = (1, 1, 1, 1)
    dims2 = (1, 1, 1)
    dims3 = (1, 1, 1, 1, 1)
    dims4 = (1, 1)
    
    filter_lag = St_filter_lag(T, SD,
                       length(dims1), length(dims2), length(dims3), length(dims4),
                       dims1, dims2, dims3, dims4)
    if lfilter
        warning( " NOT IMPLEMENTED: 3D filter not implemented")
    end

    
    if lfilter
        warning( " NOT IMPLEMENTED: 3D filter not implemented")
    end

    return filter_lag
end


#
# Moist variables
#
Base.@kwdef mutable struct St_MoistVars{T <: AbstractFloat, dim1}

    rainnc::Array{T, dim1}
    rainncv::Array{T, dim1}
    vt::Array{T, dim1}
    prod::Array{T, dim1}
    prodk::Array{T, dim1}
    vtden::Array{T, dim1}
    rdzk::Array{T, dim1}
    ρk::Array{T, dim1}
    
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


function allocate_MoistVars(npoin, T)
    
    dims1 = (npoin)

    moistvars = St_MoistVars(T, SD, length(dims1), dims1)
    
    return moistvars
end
