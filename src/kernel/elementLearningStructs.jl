#-------------------------------------------------------------------------------------------
# Element learning matrices
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_elemLearning{T <: AbstractFloat, dims1, dims2, backend}

    Avv   = KernelAbstractions.zeros(backend, T, dims1)
    Hvv   = KernelAbstractions.zeros(backend, T, dims1)
    
    A∂O∂τ = KernelAbstractions.zeros(backend, T, dims2)
    B∂O∂τ = KernelAbstractions.zeros(backend, T, dims2)
end
function allocate_elemLearning(nelem, ngl, length∂O, length∂τ, T, backend)
@info typeof(length∂τ) length∂O
    dims1 = (Int64(ngl-2), Int(ngl-2), Int64(nelem))
    dims2 = (Int64(length∂O), Int64(length∂τ))

    elemLearning = St_elemLearning{T, dims1, dims2, backend}()
    
    return elemLearning
end
