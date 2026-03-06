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
        dims1 = (Int64(npoin),1)
    else
        dims1 = (Int64(1))        
    end
    
    mp = St_Microphysics{T, dims1, backend}()
    
    return mp
end

Base.@kwdef mutable struct St_SamMicrophysics{T <:AbstractFloat, dims1, dims2, dims3, backend, VT1, VT2}

    Tabs::VT1    = KernelAbstractions.zeros(backend,  T, dims1) #Absolute temperature
    qn::VT1      = KernelAbstractions.zeros(backend,  T, dims1) #total cloud
    qi::VT1      = KernelAbstractions.zeros(backend,  T, dims1) #ice cloud
    qc::VT1      = KernelAbstractions.zeros(backend,  T, dims1) #cloud water
    qr::VT1      = KernelAbstractions.zeros(backend,  T, dims1) #rain
    qs::VT1      = KernelAbstractions.zeros(backend,  T, dims1) #snow
    qg::VT1      = KernelAbstractions.zeros(backend,  T, dims1) #graupel
    Pr::VT1      = KernelAbstractions.zeros(backend,  T, dims1) #rain precipitation flux
    Ps::VT1      = KernelAbstractions.zeros(backend,  T, dims1) #snow precipitation flux
    Pg::VT1      = KernelAbstractions.zeros(backend,  T, dims1) #graupel precipitation flux
    S_micro::VT1 = KernelAbstractions.zeros(backend,  T, dims1) #microphysical source term
    qsatt::VT1   = KernelAbstractions.zeros(backend,  T, dims1) #saturation vapor fraction
    dhldt::VT2   = KernelAbstractions.zeros(backend,  T, dims3) #Storage for preciptation source contributions to hl
    dqtdt::VT2   = KernelAbstractions.zeros(backend,  T, dims3) #Storage preciptation source contributions to qt
    dqpdt::VT2   = KernelAbstractions.zeros(backend,  T, dims3) #Storage preciptation source contributions to qp
    drad_sw::VT2 = KernelAbstractions.zeros(backend,  T, dims3) #Storage longwave flux contributions
    drad_lw::VT2 = KernelAbstractions.zeros(backend,  T, dims3) #Storage shortwave flux contribution
    flux_lw::VT1 = KernelAbstractions.zeros(backend,  T, dims1) # storage for longwave flux
    flux_sw::VT1 = KernelAbstractions.zeros(backend,  T, dims1) # storage for shortwave flux
end

function allocate_SamMicrophysics(nelem, npoin, ngl, T, backend , SD; lmoist=false)

    if lmoist
        if (SD == NSD_3D())
            dims1 = (Int64(npoin))
            dims2 = (Int64(nelem))
            dims3 = (Int64(nelem), Int64(ngl), Int64(ngl), Int64(ngl))
        else
            dims1 = (Int64(npoin))
            dims2 = (Int64(nelem))
            dims3 = (Int64(nelem), Int64(ngl), Int64(ngl))
        end

    else
        dims1 = (Int64(1))
        dims2 = (Int64(1))
        dims3 = (Int64(1))
    end

    VT1 = typeof(KernelAbstractions.zeros(backend, T, dims1))
    VT2 = typeof(KernelAbstractions.zeros(backend, T, dims3))
    mp = St_SamMicrophysics{T, dims1, dims2, dims3, backend, VT1, VT2}()

    return mp
end
