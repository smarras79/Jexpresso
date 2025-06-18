Base.@kwdef mutable struct St_Wall_model{T <: AbstractFloat, dims1, backend}

    # WIP

    τ_f  = KernelAbstractions.zeros(backend,  T, dims1)
    wθ = KernelAbstractions.zeros(backend,  T, dims1)

end

function allocate_Wall_model(nface, ngl, T, backend; lwall_model=false)

    #if lwall_model
        dims1 = (nface, ngl, ngl)
    #else
    #    dims1 = (Int64(1))
    #end

    wm = St_Wall_model{T, dims1, backend}()

    return wm
end
