#-------------------------------------------------------------------------------------------
# Element learning matrices
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_elemLearning{T <: AbstractFloat, dims1, dims2, dims3, dims4, dims5, dims6, backend}

    Avv   = KernelAbstractions.zeros(backend, T, dims1)
    Hvv   = KernelAbstractions.zeros(backend, T, dims1)
    
    A∂O∂τ = KernelAbstractions.zeros(backend, T, dims2)
    B∂O∂τ = KernelAbstractions.zeros(backend, T, dims2)
    
    A∂Ov = KernelAbstractions.zeros(backend, T, dims3)
    Av∂O = KernelAbstractions.zeros(backend, T, dims4)
    Av∂τ = KernelAbstractions.zeros(backend, T, dims5)
    
    A∂τ∂τ = KernelAbstractions.zeros(backend, T, dims6)
    
end
function allocate_elemLearning(nelem, ngl, length∂O, length∂τ, T, backend)
    
    dims1 = (Int64(ngl-2), Int(ngl-2), Int64(nelem))
    dims2 = (Int64(length∂O), Int64(length∂τ))
    dims3 = (Int64(length∂O), Int64(ngl-2)*Int64(ngl-2))
    dims4 = (Int64(ngl-2)*Int64(ngl-2), Int64(length∂O))
    dims5 = (Int64(ngl-2)*Int64(ngl-2), Int64(length∂τ))
    dims6 = (Int64(length∂τ), Int64(length∂τ))
    
    elemLearning = St_elemLearning{T, dims1, dims2, dims3, dims4, dims5, dims6, backend}()
    
    return elemLearning
end
