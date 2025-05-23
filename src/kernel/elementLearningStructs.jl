#-------------------------------------------------------------------------------------------
# Element learning matrices
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_elemLearning{T <: AbstractFloat,
                                           dims1, dims2, dims3, dims4, dims5,
                                           dims6, dims7, dims8, dims9, dims10,
                                           dims11, dims12,
                                           backend}

    Avovo = KernelAbstractions.zeros(backend, T, dims1)
    Hvovo = KernelAbstractions.zeros(backend, T, dims1)
    Avovb = KernelAbstractions.zeros(backend, T, dims7)    
    A∂O∂τ = KernelAbstractions.zeros(backend, T, dims2)
    B∂O∂τ = KernelAbstractions.zeros(backend, T, dims2)
    B∂O∂O = KernelAbstractions.zeros(backend, T, dims8)
    B∂O∂Γ = KernelAbstractions.zeros(backend, T, dims9)
    A∂Ovo = KernelAbstractions.zeros(backend, T, dims3)
    Avo∂O = KernelAbstractions.zeros(backend, T, dims4)
    Avo∂τ = KernelAbstractions.zeros(backend, T, dims5)
    A∂τ∂τ = KernelAbstractions.zeros(backend, T, dims6)
    AIoIo = KernelAbstractions.zeros(backend, T, dims10)
    AIo∂τ = KernelAbstractions.zeros(backend, T, dims11)
    A∂OIo = KernelAbstractions.zeros(backend, T, dims12)
end

function allocate_elemLearning(nelem, ngl, length∂O, length∂τ, lengthΓ, T, backend)

    elnbdypoints = 4*Int64(ngl-2) + 4
    
    dims1 = (Int64(ngl-2)^2,  Int(ngl-2)^2, Int64(nelem))
    dims2 = (Int64(length∂O), Int64(length∂τ))
    dims3 = (Int64(length∂O), Int64(ngl-2)^2, Int64(nelem))
    dims4 = (Int64(ngl-2)^2,  Int64(length∂O), Int64(nelem))
    dims5 = (Int64(ngl-2)^2,  Int64(length∂τ), Int64(nelem))
    dims6 = (Int64(length∂τ), Int64(length∂τ))
    dims7 = (Int64(ngl-2)^2,  elnbdypoints, Int64(nelem))
    dims8 = (Int64(length∂O), Int64(length∂O))
    dims9 = (Int64(length∂O), Int64(lengthΓ))
    dims10= (Int64(ngl-2)^2*Int64(nelem), Int64(ngl-2)^2*Int64(nelem))
    dims11= (Int64(ngl-2)^2*Int64(nelem), Int64(length∂τ))
    dims12= (Int64(length∂O), Int64(ngl-2)^2*Int64(nelem))
    
    elemLearning = St_elemLearning{T,
                                   dims1, dims2, dims3, dims4, dims5,
                                   dims6, dims7, dims8, dims9, dims10,
                                   dims11, dims12,
                                   backend}()
    
    return elemLearning
end
