#-------------------------------------------------------------------------------------------
# Fluxes
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_visc{T <: AbstractFloat, dims1, dims2, backend, VT1, VT2}
    μ::VT1     = KernelAbstractions.zeros(backend,  T, dims1)
    μ_max::VT2 = KernelAbstractions.zeros(backend,  T, dims2)
end

function allocate_visc(SD, nelem, npoin, ngl, T, backend; neqs=1)

    if SD == NSD_1D()
        dims1 = (Int64(nelem), Int64(ngl), Int64(neqs))
    elseif SD == NSD_2D()
        dims1 = (Int64(nelem), Int64(ngl), Int64(ngl), Int64(neqs))
    elseif SD == NSD_3D()
        dims1 = (Int64(nelem), Int64(ngl), Int64(ngl), Int64(ngl), Int64(neqs))
    end
    dims2 = Int64(neqs)

    VT1 = typeof(KernelAbstractions.zeros(backend, T, dims1))
    VT2 = typeof(KernelAbstractions.zeros(backend, T, dims2))
    μ = St_visc{T, dims1, dims2, backend, VT1, VT2}()

    return μ
end
