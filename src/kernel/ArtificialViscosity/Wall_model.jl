Base.@kwdef mutable struct St_Wall_model{T <: AbstractFloat, dims1, backend}

    # WIP

    τ_f = KernelAbstractions.zeros(backend,  T, dims1)
    wθ  = KernelAbstractions.zeros(backend,  T, dims1)

end

function allocate_Wall_model(nface, ngl, T, backend; lwall_model=false)

    dims1 = (nface, ngl, ngl, 2)
    
    wm = St_Wall_model{T, dims1, backend}()

    return wm
end
