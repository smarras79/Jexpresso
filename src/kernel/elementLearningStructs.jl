#-------------------------------------------------------------------------------------------
# Element learning matrices
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_elemLearning{T <: AbstractFloat,
                                           dims0,  dims1,  dims2,  dims3, dims4,
                                           dims5,  dims6,  dims7,  dims8, dims9,
                                           dims10, dims11, dims12, dims13,
                                           dimsML1, dimsML2, lELTraining,
                                           backend}
    
    Avovo = KernelAbstractions.zeros(backend, T, dims1)
    Avovb = KernelAbstractions.zeros(backend, T, dims7)    
    A∂O∂τ = KernelAbstractions.zeros(backend, T, dims2)
    B∂O∂τ = KernelAbstractions.zeros(backend, T, dims2)
    B∂O∂O = KernelAbstractions.zeros(backend, T, dims8)
    B∂O∂Γ = KernelAbstractions.zeros(backend, T, dims9)
    B∂τ∂τ = KernelAbstractions.zeros(backend, T, dims6)
    Bvbvb = KernelAbstractions.zeros(backend, T, dims7)
    A∂Ovo = KernelAbstractions.zeros(backend, T, dims3)
    Avo∂O = KernelAbstractions.zeros(backend, T, dims4)
    Avo∂τ = KernelAbstractions.zeros(backend, T, dims5)
    A∂τ∂τ = KernelAbstractions.zeros(backend, T, dims6)
    AIoIo = KernelAbstractions.zeros(backend, T, dims10)
    AIo∂τ = KernelAbstractions.zeros(backend, T, dims11)
    A∂OIo = KernelAbstractions.zeros(backend, T, dims12)
    AIo∂O = KernelAbstractions.zeros(backend, T, dims13)
    lEL_Train = lELTraining

    # ML:
    input_tensor  = KernelAbstractions.zeros(backend, T, dimsML1)
    output_tensor = KernelAbstractions.zeros(backend, T, dimsML2)
    
end

function allocate_elemLearning(nelem, ngl, length∂O, length∂τ, lengthΓ, T, backend; Nsamp=1, lEL_Train=false)

    elnbdypoints = 4*Int64(ngl-2) + 4
    
    dims0 = (Int64(nelem), 2)
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
    dims13= (Int64(ngl-2)^2*Int64(nelem), Int64(length∂O))
    
    # Tensors:
    k = ngl-1
    dimsML1 = ((k+1)^2, Nsamp)     #input  tensor
    dimsML2 = (4*k*(k-1)^2, Nsamp) #output tensor
    
    elemLearning = St_elemLearning{T,
                                   dims0, dims1, dims2, dims3, dims4, dims5,
                                   dims6, dims7, dims8, dims9, dims10,
                                   dims11, dims12, dims13,
                                   dimsML1, dimsML2, lEL_Train,
                                   backend}()
    
    return elemLearning
end

function write_MLtensor(tensor_column, buffer, total_cols_written, fname)
    
    push!(buffer, tensor_column)
    data = hcat(buffer...)
    col_names = ["x$(i)" for i in (total_cols_written+1):(total_cols_written+length(buffer))]
    df = DataFrame(data, col_names)
    if total_cols_written == 0
        # First write - create file with headers
        CSV.write(fname, df, transform=(col, val) -> round(val, digits=6))
    else
        # Append columns horizontally by reading, concatenating, and writing
        existing = CSV.read(fname, DataFrame)
        combined = hcat(existing, df)
        CSV.write(fname, combined, transform=(col, val) -> round(val, digits=6))
    end
    total_cols_written += length(buffer)
    buffer = Vector{Vector{Float64}}()
    
end
