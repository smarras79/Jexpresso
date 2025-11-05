#-------------------------------------------------------------------------------------------
# Solution variables
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_uODE{T <: AbstractFloat, dims1, dims2, dims3, dims4, backend}

    u    = KernelAbstractions.zeros(backend, T, dims1)
    uaux = KernelAbstractions.zeros(backend, T, dims2)
    vaux = KernelAbstractions.zeros(backend, T, dims3) #generic auxiliary array for general use
    fluxaux = KernelAbstractions.zeros(backend, T, dims4) #generic auxiliary array for general use
    
end
function allocate_uODE(SD, npoin, T, backend; neqs=1)

    dims1 = (Int64(npoin)*Int64(neqs))
    dims2 = (Int64(npoin), Int64(neqs+1))
    dims3 = (Int64(npoin))
    dims4 = (Int64(npoin), Int64(neqs+3))

    uODE = St_uODE{T, dims1, dims2, dims3, dims4, backend}()
    
    return uODE
end

Base.@kwdef mutable struct St_SolutionVars{T <: AbstractFloat, dims1, nvars, backend}

    qnp1  = KernelAbstractions.zeros(backend,  T, dims1) # qⁿ⁺¹
    qn    = KernelAbstractions.zeros(backend,  T, dims1) # qⁿ
    qout  = KernelAbstractions.zeros(backend,  T, dims1) # qⁿ
    qnm1  = KernelAbstractions.zeros(backend,  T, dims1) # qⁿ⁻¹
    qnm2  = KernelAbstractions.zeros(backend,  T, dims1) # qⁿ⁻²
    qnm3  = KernelAbstractions.zeros(backend,  T, dims1) # qⁿ⁻³
    qe    = KernelAbstractions.zeros(backend,  T, dims1) # qexact 
    press = KernelAbstractions.zeros(backend,  T, dims1) # qⁿ⁺¹   
    zb    = KernelAbstractions.zeros(backend,  T, dims1) # zb #shallow water moving bathymetry 
    
    qvars    = Array{Union{Nothing, String}}(nothing, nvars)
    qoutvars = Array{Union{Nothing, String}}(nothing, nvars)
    neqs  = nvars
    
end
function define_q(SD, nelem, npoin, ngl, qvars, T, backend; neqs=1, qoutvars=qvars)
    
    dims1 = (Int64(npoin), Int64(neqs+1))
    
    q          = St_SolutionVars{T, dims1, neqs, backend}()
    q.qvars    = qvars
    q.qoutvars = qoutvars
    
    return q
end


#-------------------------------------------------------------------------------------------
# rhs
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_rhs{T <: AbstractFloat, dims1, dims2, backend}

    #T1::Type=Float16
    
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
        dims1 = (Int64(npoin), Int64(neqs))
        dims2 = (Int64(nelem), Int64(ngl), Int64(neqs)) 
    elseif SD == NSD_2D()
        dims1 = (Int64(npoin), Int64(neqs))
        dims2 = (Int64(nelem), Int64(ngl), Int64(ngl), Int64(neqs)) 
    elseif SD == NSD_3D()
        dims1 = (Int64(npoin), Int64(neqs))
        dims2 = (Int64(nelem), Int64(ngl), Int64(ngl), Int64(ngl), Int64(neqs)) 
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
        dims1 = (Int64(ngl), Int64(neqs))
        dims2 = (Int64(ngl), Int64(neqs+1)) 
    elseif SD == NSD_2D()
        dims1 = (Int64(ngl), Int64(ngl), Int64(neqs))
        dims2 = (Int64(ngl), Int64(ngl), Int64(neqs+1)) 
    elseif SD == NSD_3D()
        dims1 = (Int64(ngl), Int64(ngl), Int64(ngl), Int64(neqs))
        dims2 = (Int64(ngl), Int64(ngl), Int64(ngl), Int64(neqs+1)) 
    end
    
    fluxes = St_fluxes{T, dims1, dims2, backend}()
    
    return fluxes
end
#-------------------------------------------------------------------------------------------
# Boundary Fluxes
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_bdy_fluxes{T <: AbstractFloat, dims1, dims2, dims3, backend}
    
    F_surf = KernelAbstractions.zeros(backend,  T, dims1)
    S_face = KernelAbstractions.zeros(backend,  T, dims2)
    S_flux = KernelAbstractions.zeros(backend,  T, dims3)

end

function allocate_bdy_fluxes(SD, nfaces, nedges, npoin, ngl, T, backend; neqs=1)

    if SD == NSD_1D()
        dims1 = (Int64(1), Int64(1))
        dims2 = (Int64(1), Int64(1))
        dims3 = (Int64(1), Int64(1))
    elseif SD == NSD_2D()
        dims1 = (Int64(ngl), Int64(neqs))
        dims2 = (Int64(nedges), Int64(ngl), Int64(neqs))
        dims3 = (Int64(npoin), Int64(neqs))
    elseif SD == NSD_3D()
        dims1 = (Int64(ngl), Int64(ngl), Int64(neqs))
        dims2 = (Int64(nfaces), Int64(ngl), Int64(ngl), Int64(neqs))
        dims3 = (Int64(npoin), Int64(neqs))
    end

    bdy_fluxes = St_bdy_fluxes{T, dims1, dims2, dims3, backend}()

    return bdy_fluxes

end
#-------------------------------------------------------------------------------------------
# Arbitrary ijk-defined quantity f:
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_fijk{T <: AbstractFloat, dims1, backend}
    fijk = KernelAbstractions.zeros(backend,  T, dims1)
end
function allocate_fijk(SD, ngl, T, backend; neqs=1)

    if SD == NSD_1D()
        dims1 = (Int64(ngl), Int64(neqs))
    elseif SD == NSD_2D()
        dims1 = (Int64(ngl), Int64(ngl), Int64(neqs))
    elseif SD == NSD_3D()
        dims1 = (Int64(ngl), Int64(ngl), Int64(ngl), Int64(neqs))
    end
    
    fijk = St_fijk{T, dims1, backend}()
    
    return fijk
end

#-------------------------------------------------------------------------------------------
# Derivative operators: e.g. gradient(f): ∇f = [∂f∂x, ∂f/∂y, ∂f/∂z]
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_∇f{T <: AbstractFloat, dims1, backend}
    ∇f_el = KernelAbstractions.zeros(backend,  T, dims1)
end
function allocate_∇f(SD, nelem, ngl, T, backend; neqs=1)
    
    if SD == NSD_1D()
        dims1 = (Int64(nelem), Int64(ngl), 1) 
    elseif SD == NSD_2D()
        dims1 = (Int64(nelem), Int64(ngl), Int64(ngl), 2) 
    elseif SD == NSD_3D()
        dims1 = (Int64(nelem), Int64(ngl), Int64(ngl), Int64(ngl), 3) 
    end
    
    ∇f = St_∇f{T, dims1, backend}()
    
    return ∇f
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
        dims1 = (Int64(npoin), Int64(neqs))
        dims2 = (Int64(nelem_semi_inf), Int64(ngr), Int64(neqs)) 
    elseif SD == NSD_2D()
        dims1 = (Int64(npoin), Int64(neqs))
        dims2 = (Int64(nelem_semi_inf), Int64(ngl), Int64(ngr), Int64(neqs)) 
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
        dims1 = (Int64(ngr), Int64(neqs))
        dims2 = (Int64(ngr), Int64(neqs+1))
    elseif SD == NSD_2D()
        dims1 = (Int64(ngl), Int64(ngr), Int64(neqs))
        dims2 = (Int64(ngl), Int64(ngr), Int64(neqs+1))
    elseif SD == NSD_3D()
        error(" src/kernel/infrastructore/params_setup.jl: 3D Laguerre arrays not coded yet!")
    end
    
    fluxes_lag = St_fluxes_lag{T, dims1, dims2, backend}()
    
    return fluxes_lag
end

#-------------------------------------------------------------------------------------------
# Filter:
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_filter{T <: AbstractFloat, dims1, dims2, dims3, dims4, backend}   
    q_t   = KernelAbstractions.zeros(backend,  T, dims1)
    fqf   = KernelAbstractions.zeros(backend,  T, dims1)
    q_ti  = KernelAbstractions.zeros(backend,  T, dims2)
    q_tij = KernelAbstractions.zeros(backend,  T, dims2)
    b     = KernelAbstractions.zeros(backend,  T, dims3)
    B     = KernelAbstractions.zeros(backend,  T, dims4)
end
function allocate_filter(SD, nelem, npoin, ngl, T, backend; neqs=1, lfilter=false)

    if lfilter
        if SD == NSD_1D()
            dims1 = (Int64(neqs), Int64(ngl))
            dims2 = (Int64(ngl))
            dims3 = (Int64(nelem), Int64(ngl), Int64(neqs))
            dims4 = (Int64(npoin), Int64(neqs))
        elseif SD == NSD_2D()
            dims1 = (Int64(neqs), Int64(ngl), Int64(ngl))
            dims2 = (Int64(ngl), Int64(ngl))
            dims3 = (Int64(nelem), Int64(ngl), Int64(ngl), Int64(neqs))
            dims4 = (Int64(npoin), Int64(neqs))
        elseif SD == NSD_3D()
            dims1 = (Int64(neqs), Int64(ngl), Int64(ngl), Int64(ngl))
            dims2 = (Int64(ngl), Int64(ngl), Int64(ngl))
            dims3 = (Int64(nelem), Int64(ngl), Int64(ngl), Int64(ngl), Int64(neqs))
            dims4 = (Int64(npoin), Int64(neqs))
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
function allocate_filter_lag(SD, nelem_semi_inf, npoin, ngl, ngr, T, backend; neqs=1, lfilter=false)

    if lfilter
        if SD == NSD_1D()
            dims1 = (Int64(neqs), Int64(ngr))
            dims2 = (Int64(ngr))
            dims3 = (Int64(nelem_semi_inf), Int64(ngr), Int64(neqs))
            dims4 = (Int64(npoin), Int64(neqs))
        elseif SD == NSD_2D()
            dims1 = (Int64(neqs), Int64(ngl), Int64(ngr))
            dims2 = (Int64(ngl), Int64(ngr))
            dims3 = (Int64(nelem_semi_inf), Int64(ngl), Int64(ngr), Int64(neqs))
            dims4 = (Int64(npoin), Int64(neqs))
        elseif SD == NSD_3D()
            # WARNING Allocate only 1 because there is no 3D filter yet
            #dims1 = (1, 1, 1, 1)
            #dims2 = (1, 1, 1)
            #dims3 = (1, 1, 1, 1, 1)
            #dims4 = (1, 1)
            #warning( " 3D laguerre filter not implemented yet")
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
# GPU auxiliary arrays
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_gpuAux{T <: AbstractFloat, dims1, dims2, dims3, backend}

    flux_gpu   = KernelAbstractions.zeros(backend, T, dims1)
    source_gpu = KernelAbstractions.zeros(backend, T, dims2)
    qbdy_gpu   = KernelAbstractions.zeros(backend, T, dims3)
    
end
function allocate_gpuAux(SD, nelem, nedges_bdy, nfaces_bdy, ngl, T, backend; neqs=1)

    if backend == CPU()
        dims1 = (1, 1, 1, 1)
        dims2 = dims1
        dims3 = (1, 1, 1)
    else
        if SD == NSD_1D()
            dims1 = (Int64(nelem),      Int64(ngl), 2*neqs)
            dims2 = (Int64(nelem),      Int64(ngl),   neqs)
            dims3 = (Int64(0))
        elseif SD == NSD_2D()
            dims1 = (Int64(nelem),      Int64(ngl), Int64(ngl), 2*neqs)
            dims2 = (Int64(nelem),      Int64(ngl), Int64(ngl),   neqs)
            dims3 = (Int64(nedges_bdy), Int64(ngl),               neqs)
        elseif SD == NSD_3D()
            dims1 = (Int64(nelem),      Int64(ngl), Int64(ngl), Int64(ngl), 3*neqs)
            dims2 = (Int64(nelem),      Int64(ngl), Int64(ngl), Int64(ngl),   neqs)
            dims3 = (Int64(nfaces_bdy), Int64(ngl), Int64(ngl),               neqs)
        end
    end
    
    gpuAux = St_gpuAux{T, dims1, dims2, dims3, backend}()
    
    return gpuAux
end
#
# GPU Laguerre
#
Base.@kwdef mutable struct St_gpuAux_lag{T <: AbstractFloat, dims1, dims2, dims3, backend}

    flux_lag_gpu   = KernelAbstractions.zeros(backend, T, dims1)
    source_lag_gpu = KernelAbstractions.zeros(backend, T, dims2)
    qbdy_lag_gpu   = KernelAbstractions.zeros(backend, T, dims3)
    
end
function allocate_gpuAux_lag(SD, nelem_semi_inf, nedges_bdy, nfaces_bdy, ngl, ngr, T, backend; neqs=1)

     if backend == CPU()
        dims1 = (1, 1, 1, 1)
        dims2 = dims1
        dims3 = (1, 1, 1)
    else
        if SD == NSD_1D()
            dims1 = (Int64(nelem_semi_inf), Int64(ngr), 2*neqs)
            dims2 = (Int64(nelem_semi_inf), Int64(ngr),   neqs)
            dims3 = (Int64(0))
        elseif SD == NSD_2D()
            dims1 = (Int64(nelem_semi_inf), Int64(ngl), Int64(ngr), 2*neqs)
            dims2 = (Int64(nelem_semi_inf), Int64(ngl), Int64(ngr),   neqs)
            dims3 = (Int64(nedges_bdy),     Int64(ngl),               neqs)
        elseif SD == NSD_3D()
            error(" globalStructs.jl: --> 3D Laguerre not implemented yet!")
        end
    end
    
    gpuAux_lag = St_gpuAux_lag{T, dims1, dims2, dims3, backend}()
        
    return gpuAux_lag
end

Base.@kwdef mutable struct St_gpuMoist{T <: AbstractFloat, dims1, dims2, dims3, dims4, backend}

    flux_micro   = KernelAbstractions.zeros(backend, T, dims1)
    source_micro = KernelAbstractions.zeros(backend, T, dims2)
    adjusted     = KernelAbstractions.zeros(backend, T, dims3)
    Pm           = KernelAbstractions.zeros(backend, T, dims4)
end

function allocate_gpuMoist(SD, npoin, nelem, ngl, T, backend, lmoist; neqs=1)

    if backend == CPU() || lmoist == false
        dims1 = (1, 1, 1, 1)
        dims2 = dims1
        dims3 = (1,1)
        dims4 = dims3
    else
        dims1 = (Int64(nelem),      Int64(ngl), Int64(ngl), Int64(ngl), 4)
        dims2 = (Int64(nelem),      Int64(ngl), Int64(ngl), Int64(ngl), 4)
        dims3 = (Int64(npoin), 9)
        dims4 = (Int64(npoin), 3)
    end

    gpuMoist = St_gpuMoist{T, dims1, dims2, dims3, dims4, backend}()

    return gpuMoist
end


