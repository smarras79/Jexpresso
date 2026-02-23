Base.@kwdef mutable struct St_Wall_model{T <: AbstractFloat, dims1, dims2, dims3, backend}

    τ_f = KernelAbstractions.zeros(backend,  T, dims1)
    wθ  = KernelAbstractions.zeros(backend,  T, dims2)
    wqv = KernelAbstractions.zeros(backend,  T, dims3)

end

function allocate_Wall_model(nface, ngl, T, backend; lwall_model=false, lmoist=false)
    # Fixed: Only allocate arrays if wall model is enabled to save memory
    if (lwall_model == true) && (lmoist == true)
        dims1 = (nface, ngl, ngl, 3)
        dims2 = (nface, ngl, ngl, 1)
        dims3 = (nface, ngl, ngl, 1)
    elseif (lwall_model == true) && (lmoist == false)
        dims1 = (nface, ngl, ngl, 3)
        dims2 = (nface, ngl, ngl, 1)
        dims3 = (0, 0, 0, 0)
    else
        # Allocate minimal empty arrays when wall model is disabled
        dims1 = (0, 0, 0, 0)
        dims2 = (0, 0, 0, 0)
        dims3 = (0, 0, 0, 0)
    end
    
    wm = St_Wall_model{T, dims1, dims2, dims3, backend}()
    return wm
end
