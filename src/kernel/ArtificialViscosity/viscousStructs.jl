#-------------------------------------------------------------------------------------------
# Fluxes
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_visc{T <: AbstractFloat, dims1, backend}    
    μ = KernelAbstractions.zeros(backend,  T, dims1)
end

function allocate_visc(SD, nelem, npoin, ngl, T, backend; neqs=1)

    if SD == NSD_1D()
        dims1 = (Int64(nelem), Int64(ngl), Int64(neqs))
    elseif SD == NSD_2D()
        dims1 = (Int64(nelem), Int64(ngl), Int64(ngl), Int64(neqs))
    elseif SD == NSD_3D()
        dims1 = (Int64(nelem), Int64(ngl), Int64(ngl), Int64(ngl), Int64(neqs))
    end
    
    μ = St_visc{T, dims1, backend}()
    
    return μ
end
