Base.@kwdef mutable struct St_Wall_model{T <: AbstractFloat, dims1, dims2, backend, Arr1, Arr2}

    τ_f::Arr1 = KernelAbstractions.zeros(backend,  T, dims1)
    wθ::Arr2  = KernelAbstractions.zeros(backend,  T, dims2)

end

function allocate_Wall_model(nface, ngl, T, backend; lwall_model=false)
    # Fixed: Only allocate arrays if wall model is enabled to save memory
    if lwall_model
        dims1 = (nface, ngl, ngl, 3)
        dims2 = (nface, ngl, ngl, 1)

        arr1 = typeof(KernelAbstractions.zeros(backend,  T, dims1))
        arr2 = typeof(KernelAbstractions.zeros(backend,  T, dims2))

        wm = St_Wall_model{T, dims1, dims2, backend, arr1, arr2}()
    else
        # Allocate minimal empty arrays when wall model is disabled
        dims1 = (0, 0, 0, 0)
        dims2 = (0, 0, 0, 0)

        arr1 = typeof(KernelAbstractions.zeros(backend,  T, dims1))
        arr2 = typeof(KernelAbstractions.zeros(backend,  T, dims2))

        wm = St_Wall_model{T, dims1, dims2, backend, arr1, arr2}()
    end

    return wm
end
