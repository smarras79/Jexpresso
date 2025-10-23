Base.@kwdef mutable struct St_Wall_model{T <: AbstractFloat, dims1, dims2, backend}
    
    τ_f = KernelAbstractions.zeros(backend,  T, dims1)
    wθ  = KernelAbstractions.zeros(backend,  T, dims2)

end

function allocate_Wall_model(nface, ngl, T, backend; lwall_model=false)

    dims1 = (nface, ngl, ngl, 3)
    dims2 = (nface, ngl, ngl, 1)
    
    wm = St_Wall_model{T, dims1, dims2, backend}()

    return wm
end
