#-------------------------------------------------------------------------------------------
# Microphysics (mp) variables:
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_Microphysics{T <: AbstractFloat, dims1, backend}

    # WIP
    
    rainnc  = KernelAbstractions.zeros(backend,  T, dims1)
    rainncv = KernelAbstractions.zeros(backend,  T, dims1)
    vt      = KernelAbstractions.zeros(backend,  T, dims1)
    prod    = KernelAbstractions.zeros(backend,  T, dims1)
    prodk   = KernelAbstractions.zeros(backend,  T, dims1)
    vtden   = KernelAbstractions.zeros(backend,  T, dims1)
    rdzk    = KernelAbstractions.zeros(backend,  T, dims1)
    rdzw    = KernelAbstractions.zeros(backend,  T, dims1)
    ρk      = KernelAbstractions.zeros(backend,  T, dims1)
    temp1   = KernelAbstractions.zeros(backend,  T, dims1)
    temp2   = KernelAbstractions.zeros(backend,  T, dims1)
    
end

function allocate_Microphysics(nelem, npoin, ngl, T, backend; lmoist=false)

    if lmoist
        dims1 = (Int64(npoin))
    else
        dims1 = (Int64(1))        
    end
    
    mp = St_Microphysics{T, dims1, backend}()
    
    return mp
end
