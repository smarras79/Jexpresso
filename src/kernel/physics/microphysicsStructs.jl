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

Base.@kwdef mutable struct St_SamMicrophysics{T <:AbstractFloat, dims1, dims2, dims3, backend, arr1, arr3}

    Tabs::arr1    = KernelAbstractions.zeros(backend,  T, dims1) #Absolute temperature
    qn::arr1      = KernelAbstractions.zeros(backend,  T, dims1) #total cloud
    qi::arr1      = KernelAbstractions.zeros(backend,  T, dims1) #ice cloud
    qc::arr1      = KernelAbstractions.zeros(backend,  T, dims1) #cloud water
    qr::arr1      = KernelAbstractions.zeros(backend,  T, dims1) #rain
    qs::arr1      = KernelAbstractions.zeros(backend,  T, dims1) #snow
    qg::arr1      = KernelAbstractions.zeros(backend,  T, dims1) #graupel
    Pr::arr1      = KernelAbstractions.zeros(backend,  T, dims1) #rain precipitation flux
    Ps::arr1      = KernelAbstractions.zeros(backend,  T, dims1) #snow precipitation flux
    Pg::arr1      = KernelAbstractions.zeros(backend,  T, dims1) #graupel precipitation flux
    S_micro::arr1 = KernelAbstractions.zeros(backend,  T, dims1)  #microphysical source term
    qsatt::arr1   = KernelAbstractions.zeros(backend,  T, dims1)  #saturation vapor fraction
    dhldt::arr3   = KernelAbstractions.zeros(backend,  T, dims3, ) #Storage for preciptation source contributions to hl
    dqtdt::arr3   = KernelAbstractions.zeros(backend,  T, dims3) #Storage preciptation source contributions to qt
    dqpdt::arr3   = KernelAbstractions.zeros(backend,  T, dims3) #Storage preciptation source contributions to qp
    drad_sw::arr3   = KernelAbstractions.zeros(backend,  T, dims3) #Storage longwave flux contributions
    drad_lw::arr3   = KernelAbstractions.zeros(backend,  T, dims3) #Storage shortwave flux contribution
    flux_lw::arr1 = KernelAbstractions.zeros(backend,  T, dims1) # storage for longwave flux
    flux_sw::arr1 = KernelAbstractions.zeros(backend,  T, dims1) # storage for shortwave flux
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

    arr1 = typeof(KernelAbstractions.zeros(backend,  T, dims1))
    arr3 = typeof(KernelAbstractions.zeros(backend,  T, dims3))

    mp = St_SamMicrophysics{T, dims1, dims2, dims3, backend, arr1, arr3}()

    return mp
end
