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

Base.@kwdef mutable struct St_SamMicrophysics{T <:AbstractFloat, dims1, dims2, dims3, backend}

    Tabs    = KernelAbstractions.zeros(backend,  T, dims1) #Absolute temperature
    qn      = KernelAbstractions.zeros(backend,  T, dims1) #total cloud 
    qi      = KernelAbstractions.zeros(backend,  T, dims1) #ice cloud
    qc      = KernelAbstractions.zeros(backend,  T, dims1) #cloud water
    qr      = KernelAbstractions.zeros(backend,  T, dims1) #rain
    qs      = KernelAbstractions.zeros(backend,  T, dims1) #snow
    qg      = KernelAbstractions.zeros(backend,  T, dims1) #graupel
    Pr      = KernelAbstractions.zeros(backend,  T, dims1) #rain precipitation flux
    Ps      = KernelAbstractions.zeros(backend,  T, dims1) #snow precipitation flux
    Pg      = KernelAbstractions.zeros(backend,  T, dims1) #graupel precipitation flux
    S_micro = KernelAbstractions.zeros(backend,  T, dims1)  #microphysical source term
    dhldt   = KernelAbstractions.zeros(backend,  T, dims2, dims3, dims3, dims3) #Storage for preciptation source contributions to hl
    dqtdt   = KernelAbstractions.zeros(backend,  T, dims2, dims3, dims3, dims3) #Storage preciptation source contributions to qt
    dqpdt   = KernelAbstractions.zeros(backend,  T, dims2, dims3, dims3, dims3) #Storage preciptation source contributions to qp
end

function allocate_SamMicrophysics(nelem, npoin, ngl, T, backend; lmoist=false)

    if lmoist
        dims1 = (Int64(npoin))
        dims2 = (Int64(nelem))
        dims3 = (Int64(ngl))
    else
        dims1 = (Int64(1))
        dims2 = (Int64(1))
        dims3 = (Int64(1)) 
    end

    mp = St_SamMicrophysics{T, dims1, dims2, dims3, backend}()

    return mp
end
