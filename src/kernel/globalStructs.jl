#-------------------------------------------------------------------------------------------
# Solution variables
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_SolutionVars{T <: AbstractFloat, dims1, nvars, backend}

    qnp1  = KernelAbstractions.zeros(backend,  T, dims1) # qⁿ⁺¹
    qn    = KernelAbstractions.zeros(backend,  T, dims1) # qⁿ
    qq    = KernelAbstractions.zeros(backend,  T, dims1) # qⁿ
    qnm1  = KernelAbstractions.zeros(backend,  T, dims1) # qⁿ⁻¹
    qnm2  = KernelAbstractions.zeros(backend,  T, dims1) # qⁿ⁻²
    qnm3  = KernelAbstractions.zeros(backend,  T, dims1) # qⁿ⁻³
    qe    = KernelAbstractions.zeros(backend,  T, dims1) # qexact    
    zb    = KernelAbstractions.zeros(backend,  T, dims1) # zb #shallow water moving bathymetry
    qvars = Array{Union{Nothing, String}}(nothing, nvars)
    neqs  = nvars
    
end
function define_q(SD, nelem, npoin, ngl, qvars, T, backend; neqs=1)
    
    dims1 = (npoin, neqs+1)
    
    q =  St_SolutionVars{T, dims1, neqs, backend}()
    
    return q
end

#-------------------------------------------------------------------------------------------
# rhs
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_rhs{T <: AbstractFloat, dims1, dims2, backend}
    
    RHS          = KernelAbstractions.zeros(backend,  T, dims1)         
    RHS_visc     = KernelAbstractions.zeros(backend,  T, dims1)    
    rhs_el       = KernelAbstractions.zeros(backend,  T, dims2)     
    rhs_diff_el  = KernelAbstractions.zeros(backend,  T, dims2) 
    rhs_diffξ_el = KernelAbstractions.zeros(backend,  T, dims2)
    rhs_diffη_el = KernelAbstractions.zeros(backend,  T, dims2)
    rhs_diffζ_el = KernelAbstractions.zeros(backend,  T, dims2)
    
end
function allocate_rhs(SD, nelem, npoin, ngl, T, backend; neqs=1)

    if SD == NSD_1D()
        dims1 = (npoin, neqs)
        dims2 = (nelem, ngl, neqs) 
    elseif SD == NSD_2D()
        dims1 = (npoin, neqs)
        dims2 = (nelem, ngl, ngl, neqs) 
    elseif SD == NSD_3D()
        dims1 = (npoin, neqs)
        dims2 = (nelem, ngl, ngl, ngl, neqs) 
    end
    
    rhs = St_rhs{T, dims1, dims2, backend}()
    
    return rhs
end


#-------------------------------------------------------------------------------------------
# Fluxes
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_fluxes{T <: AbstractFloat, dims1, dims2, backend}
    
    F = KernelAbstractions.zeros(backend,  T, dims1)
    G = KernelAbstractions.zeros(backend,  T, dims1)
    H = KernelAbstractions.zeros(backend,  T, dims1)
    S = KernelAbstractions.zeros(backend,  T, dims1)
    uprimitive = KernelAbstractions.zeros(backend,  T, dims2)
    
end
function allocate_fluxes(SD, npoin, ngl, T, backend; neqs=1)

    if SD == NSD_1D()
        dims1 = (ngl, neqs)
        dims2 = (ngl, neqs+1) 
    else SD == NSD_2D()
        dims1 = (ngl, ngl, neqs)
        dims2 = (ngl, ngl, neqs+1) 
    else SD == NSD_3D)(
        dims1 = (ngl, ngl, ngl, neqs)
        dims2 = (ngl, ngl, ngl, neqs+1) 
    end
        
    fluxes = St_fluxes{T, dims1, dims2, backend}()
    
    return fluxes
end

#-------------------------------------------------------------------------------------------
# rhs Laguerre
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_rhs_lag{T <: AbstractFloat, dims1, dims2, backend}
    
    RHS_lag          = KernelAbstractions.zeros(backend,  T, dims1)         
    RHS_visc_lag     = KernelAbstractions.zeros(backend,  T, dims1)    
    rhs_el_lag       = KernelAbstractions.zeros(backend,  T, dims2)     
    rhs_diff_el_lag  = KernelAbstractions.zeros(backend,  T, dims2) 
    rhs_diffξ_el_lag = KernelAbstractions.zeros(backend,  T, dims2)
    rhs_diffη_el_lag = KernelAbstractions.zeros(backend,  T, dims2)
    rhs_diffζ_el_lag = KernelAbstractions.zeros(backend,  T, dims2)
    
end
function allocate_rhs_lag(SD, nelem_semi_inf, npoin, ngl, ngr, T, backend; neqs=1)

    if SD == NSD_1D()
        dims1 = (npoin, neqs)
        dims2 = (nelem_semi_inf, ngr, neqs) 
    elseif SD == NSD_2D()
        dims1 = (npoin, neqs)
        dims2 = (nelem_semi_inf, ngl, ngr, neqs) 
    elseif SD == NSD_3D()
        error(" src/kernel/infrastructore/params_setup.jl: 3D Laguerre arrays not coded yet!")
    end
    
    rhs_lag = St_rhs_lag{T, dims1, dims2, backend}()
    
    return rhs_lag
end

#-------------------------------------------------------------------------------------------
# Fluxes Laguerre
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_fluxes_lag{T <: AbstractFloat, dims1, dims2, backend}
    
    F_lag= KernelAbstractions.zeros(backend,  T, dims1)
    G_lag= KernelAbstractions.zeros(backend,  T, dims1)
    H_lag= KernelAbstractions.zeros(backend,  T, dims1)
    S_lag= KernelAbstractions.zeros(backend,  T, dims1)
    uprimitive_lag= KernelAbstractions.zeros(backend,  T, dims2) 
    
end
function allocate_fluxes_lag(SD, ngl, ngr, T, backend; neqs=1)
    if SD == NSD_1D()
        dims1 = (ngr, neqs)
        dims2 = (ngr, neqs+1) 
    else SD == NSD_2D()
        dims1 = (ngl, ngr, neqs)
        dims2 = (ngl, ngr, neqs+1) 
    else SD == NSD_3D()
        error(" src/kernel/infrastructore/params_setup.jl: 3D Laguerre arrays not coded yet!")
    end
    
    fluxes_lag = St_fluxes_lag{T, dims1, dims2, backend}()
    
    return fluxes_lag
end

#-------------------------------------------------------------------------------------------
# Filter:
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_filter{T <: AbstractFloat, dims1, dims2, dims3, dims4, backend}   
    q_t  = KernelAbstractions.zeros(backend,  T, dims1)
    fqf  = KernelAbstractions.zeros(backend,  T, dims1)
    q_ti = KernelAbstractions.zeros(backend,  T, dims2)
    b    = KernelAbstractions.zeros(backend,  T, dims3)
    B    = KernelAbstractions.zeros(backend,  T, dims4)
end
function allocate_filter(SD, nelem, npoin, ngl, T; neqs=1, lfilter=false)

    if lfilter
        if SD == NSD_1D()
            dims1 = (neqs, ngl)
            dims2 = (ngl)
            dims3 = (nelem, ngl, neqs)
            dims4 = (npoin, neqs)
        elseif SD == NSD_2D()
            dims1 = (neqs, ngl, ngl)
            dims2 = (ngl, ngl)
            dims3 = (nelem, ngl, ngl, neqs)
            dims4 = (npoin, neqs)
        elseif SD == NSD_3D()
            # WARNING Allocate only 1 because there is no 3D filter yet
            dims1 = (1, 1, 1, 1)
            dims2 = (1, 1, 1)
            dims3 = (1, 1, 1, 1, 1)
            dims4 = (1, 1)
            warning( " 3D filter not implemented yet")
        end
    else
        if SD == NSD_1D()
            dims1 = (1, 1)
            dims2 = (1)
            dims3 = (1, 1, 1)
            dims4 = (1, 1)
        elseif SD == NSD_2D()
            dims1 = (1, 1, 1)
            dims2 = (1, 1)
            dims3 = (1, 1, 1, 1)
            dims4 = (1, 1)
        elseif SD == NSD_3D()
            dims1 = (1, 1, 1, 1)
            dims2 = (1, 1, 1)
            dims3 = (1, 1, 1, 1, 1)
            dims4 = (1, 1)
        end
    end

    filter = St_filter{T, dims1, dims2, dims3, dims4, backend}()
    
    return filter
end



#-------------------------------------------------------------------------------------------
# Laguerre filter
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_filter_lag{T <: AbstractFloat, dims1, dims2, dims3, dims4, backend}
    
    q_t_lag  = KernelAbstractions.zeros(backend,  T, dims1)
    fqf_lag  = KernelAbstractions.zeros(backend,  T, dims1)
    q_ti_lag = KernelAbstractions.zeros(backend,  T, dims2)
    b_lag    = KernelAbstractions.zeros(backend,  T, dims3)
    B_lag    = KernelAbstractions.zeros(backend,  T, dims4)   
    
end
function allocate_filter_lag(SD, nelem_semi_inf, npoin, ngl, ngr, T; neqs=1, lfilter=false)

    if lfilter
        if SD == NSD_1D()
            dims1 = (neqs, ngr)
            dims2 = (ngr)
            dims3 = (nelem_semi_inf, ngr, neqs)
            dims4 = (npoin, neqs)
        elseif SD == NSD_2D()
            dims1 = (neqs, ngl, ngr)
            dims2 = (ngl, ngr)
            dims3 = (nelem_semi_inf, ngl, ngr, neqs)
            dims4 = (npoin, neqs)
        elseif SD == NSD_3D()
            # WARNING Allocate only 1 because there is no 3D filter yet
            dims1 = (1, 1, 1, 1)
            dims2 = (1, 1, 1)
            dims3 = (1, 1, 1, 1, 1)
            dims4 = (1, 1)
            warning( " 3D laguerre filter not implemented yet")
        end
        
    else
        if SD == NSD_1D()
            dims1 = (1, 1)
            dims2 = (1)
            dims3 = (1, 1, 1)
            dims4 = (1, 1)
        elseif SD == NSD_2D()
            dims1 = (1, 1, 1)
            dims2 = (1, 1)
            dims3 = (1, 1, 1, 1)
            dims4 = (1, 1)
        elseif SD == NSD_3D()             
            dims1 = (1, 1, 1, 1)
            dims2 = (1, 1, 1)
            dims3 = (1, 1, 1, 1, 1)
            dims4 = (1, 1)            
        end
    end
    
    filter_lag = St_filter_lag{T, dims1, dims2, dims3, dims4, backend}()
    
    return filter_lag
end

#-------------------------------------------------------------------------------------------
# Moist variables WIP    
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_MoistVars{T <: AbstractFloat, dims1, backend}

    # WIP
    
    rainnc  = KernelAbstractions.zeros(backend,  T, dims1)
    rainncv = KernelAbstractions.zeros(backend,  T, dims1)
    vt      = KernelAbstractions.zeros(backend,  T, dims1)
    prod    = KernelAbstractions.zeros(backend,  T, dims1)
    prodk   = KernelAbstractions.zeros(backend,  T, dims1)
    vtden   = KernelAbstractions.zeros(backend,  T, dims1)
    rdzk    = KernelAbstractions.zeros(backend,  T, dims1)
    ρk      = KernelAbstractions.zeros(backend,  T, dims1)
    
end

function allocate_MoistVars(nelem, npoin, ngl, T; neqs=1, lfilter=false)
    
    # WIP
    
    dims1 = (npoin)
    
    moistvars = St_MoistVars{T, dims1, backend}()
    
    return moistvars
end
