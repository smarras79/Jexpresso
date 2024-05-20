#-------------------------------------------------------------------------------------------
# Microphysics (mp) variables:
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_MicroPhysics{T <: AbstractFloat, dims1, backend}

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

function allocate_MicroPhysics(nelem, npoin, ngl, T, backend; neqs=1, lfilter=false)
    
    # WIP
    
    dims1 = (Int64(npoin))
    
    mp = St_MoistVars{T, dims1, backend}()
    
    return mp
end
