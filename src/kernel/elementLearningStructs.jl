#-------------------------------------------------------------------------------------------
# Element learning matrices
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_elemLearning{T <: AbstractFloat, dims1, dims2, dims3, backend}

    Avv   = KernelAbstractions.zeros(backend, T, dims1)
    Hvv   = KernelAbstractions.zeros(backend, T, dims1)
    
    A∂O∂τ = KernelAbstractions.zeros(backend, T, dims2)
    B∂O∂τ = KernelAbstractions.zeros(backend, T, dims2)
    
    
end
function allocate_elemLearning(SD, nelem, ngl, length∂O, length∂τ, T, backend; neqs=1)

    dims1 = (Int64(ngl), Int(ngl), Int64(nelem))    
    dims2 = (Int64(length∂O), Int64(length∂τ))

    elemLearning = St_elemLearning{T, dims1, dims2, dims3, backend}()
    
    return elemLearning
end
