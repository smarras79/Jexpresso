include("./mesh/meshStructs.jl")

#-------------------------------------------------------------------------------------------
# Element learning matrices
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_elemLearning{T <: AbstractFloat,
                                           dims0,  dims1,  dims2,  dims3, dims4,
                                           dims5,  dims6,  dims7,  dims8, dims9,
                                           dims10, dims11, dims12, dims13, dims14,
                                           dims15, dimsML1, dimsML2, lELTraining,
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
    
    T1    = KernelAbstractions.zeros(backend, T, dims15)
    T2    = KernelAbstractions.zeros(backend, T, dims14)
    Tie   = KernelAbstractions.zeros(backend, T, dims14)
    
    lEL_Train = lELTraining

    # ML:
    input_tensor  = KernelAbstractions.zeros(backend, T, dimsML1)
    output_tensor = KernelAbstractions.zeros(backend, T, dimsML2)
    
end

function allocate_elemLearning(nelem, ngl, length∂O, length∂τ, lengthΓ, T, backend; Nsamp=1, lEL_Train=false)

    elnbdypoints = 4*Int64(ngl-2) + 4
    
    dims0  = (Int64(nelem), 2)
    dims1  = (Int64(ngl-2)^2,  Int(ngl-2)^2, Int64(nelem))
    dims2  = (Int64(length∂O), Int64(length∂τ))
    dims3  = (Int64(length∂O), Int64(ngl-2)^2, Int64(nelem))
    dims4  = (Int64(ngl-2)^2,  Int64(length∂O), Int64(nelem))
    dims5  = (Int64(ngl-2)^2,  Int64(length∂τ), Int64(nelem))
    dims6  = (Int64(length∂τ), Int64(length∂τ))
    dims7  = (Int64(ngl-2)^2,  elnbdypoints, Int64(nelem))
    dims8  = (Int64(length∂O), Int64(length∂O))
    dims9  = (Int64(length∂O), Int64(lengthΓ))
    dims10 = (Int64(ngl-2)^2*Int64(nelem), Int64(ngl-2)^2*Int64(nelem))
    dims11 = (Int64(ngl-2)^2*Int64(nelem), Int64(length∂τ))
    dims12 = (Int64(length∂O), Int64(ngl-2)^2*Int64(nelem))
    dims13 = (Int64(ngl-2)^2*Int64(nelem), Int64(length∂O))
    dims14 = (Int64(ngl-2)^2, elnbdypoints)
    dims15 = (elnbdypoints, elnbdypoints)
    
    # Tensors:
    k = ngl-1
    dimsML1 = ((k+1)^2, Nsamp)     #input  tensor
    dimsML2 = (4*k*(k-1)^2, Nsamp) #output tensor
    
    elemLearning = St_elemLearning{T,
                                   dims0, dims1, dims2, dims3, dims4, dims5,
                                   dims6, dims7, dims8, dims9, dims10,
                                   dims11, dims12, dims13, dims14, dims15,
                                   dimsML1, dimsML2, lEL_Train,
                                   backend}()
    
    return elemLearning
end

function write_MLtensor!(buffer::Vector{Vector{Float64}}, tensor_column::Vector{Float64})
    push!(buffer, tensor_column)  # accumulate in caller-owned buffer
end
function flush_MLtensor!(buffer::Vector{Vector{Float64}}, total_cols_written, fname)
    isempty(buffer) && return total_cols_written
    
    # Stack columns into a matrix without splatting
    nrows = length(buffer[1])
    ncols = length(buffer)
    data  = Matrix{Float64}(undef, nrows, ncols)
    for (j, col) in enumerate(buffer)
        data[:, j] .= col
    end

    col_names = ["x$(total_cols_written + i)" for i in 1:ncols]
    df = DataFrame(data, col_names)

    if total_cols_written == 0
        CSV.write(fname, df, transform=(col, val) -> round(val, digits=6))
    else
        existing = CSV.read(fname, DataFrame)
        combined = hcat(existing, df)
        CSV.write(fname, combined, transform=(col, val) -> round(val, digits=6))
    end

    total_cols_written += ncols
    empty!(buffer)          # ← mutates in place; caller sees the cleared buffer
    return total_cols_written
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

function elementLearning_Axb!(u, uaux, mesh::St_mesh,
                              A, ubdy, EL,
                              avisc, 
                              bufferin, bufferout,
                              ABC, BC, BOΓg, gΓ;
                              isamp=1,
                              total_cols_writtenin=0,
                              total_cols_writtenout=0) 

    mesh.lengthO =  mesh.length∂O +  mesh.lengthIo
    nelintpoints = (mesh.ngl-2)*(mesh.ngl-2)
    nelpoints    = size(mesh.conn)[2]
    elnbdypoints = nelpoints - nelintpoints

    for iel=1:mesh.nelem
        #
        # A∂oᵥₒ
        #
        ii = 1
        for i = elnbdypoints+1:nelpoints
            ipo = mesh.conn[iel, i]

            for io=1:length(mesh.∂O)
                io1 = mesh.∂O[io]
                EL.A∂Ovo[io, ii, iel] = A[io1, ipo]
            end
            
            #
            # Aᵥₒ∂τ
            #
            for jτ = 1:mesh.length∂τ
                jτ1 = mesh.∂τ[jτ]
                EL.Avo∂τ[ii, jτ, iel] = A[ipo, jτ1]
            end

            #
            # Aᵥₒᵥₒ
            #
            jj = 1
            for j = elnbdypoints+1:nelpoints
                jpo = mesh.conn[iel, j]
                EL.Avovo[ii, jj, iel] = A[ipo, jpo]
                jj += 1
            end

            #
            # Aᵥₒᵥb
            #
            for j = 1:elnbdypoints
                jpb = mesh.conn[iel, j]
                EL.Avovb[ii, j, iel] = A[ipo, jpb]
            end
            
            ii += 1
        end
    end
    
    #
    # A∂O∂τ ⊂ A∂τ∂τ
    #
    for j1=1:length(mesh.∂τ)
        jτ1 = mesh.∂τ[j1]

        for i1=1:length(mesh.∂O)
            iO1 = mesh.∂O[i1]
            EL.A∂O∂τ[i1, j1] = A[iO1, jτ1]
        end

        for j2=1:length(mesh.∂τ)
            jτ2 = mesh.∂τ[j2]
            EL.A∂τ∂τ[j1, j2] = A[jτ1, jτ2]
        end
    end

    #
    # A∂OIo
    #
    for jo=1:mesh.length∂O
        jo1 = mesh.∂O[jo]
        for io=1:mesh.lengthIo
            io1 = mesh.Io[io]
            EL.A∂OIo[jo, io] = A[jo1, io1]
        end
    end

    #
    # AIo∂O
    #
    for jo=1:mesh.length∂O
        jo1 = mesh.∂O[jo]
        for io=1:mesh.lengthIo
            io1 = mesh.Io[io]
            EL.AIo∂O[io, jo] = A[io1, jo1]
        end
    end

    #
    # AIoIo
    #
    for io = 1:mesh.lengthIo
        io1 = mesh.Io[io]
        for jo = 1:mesh.lengthIo
            jo1 = mesh.Io[jo]
            EL.AIoIo[io, jo] = A[io1, jo1]
        end
    end

    #
    # AIo∂τ
    #
    for jτ = 1:mesh.length∂τ
        jτ1 = mesh.∂τ[jτ]
        for io=1:mesh.lengthIo
            io1 = mesh.Io[io]
            EL.AIo∂τ[io, jτ] = A[io1, jτ1]
        end
    end

    # inv(AIoIo) — needed by training branch
    invAIoIo = inv(EL.AIoIo)

    # AIo,Γ — needed by training branch
    AIoΓ = similar(A, (mesh.lengthIo, mesh.lengthΓ))
    for iΓ = 1:mesh.lengthΓ
        g1 = mesh.Γ[iΓ]
        for io = 1:mesh.lengthIo
            io1 = mesh.Io[io]
            AIoΓ[io, iΓ] = A[io1, g1]
        end
    end

    # gΓ — needed by both branches
    #gΓ = zeros(mesh.lengthΓ)
    for iΓ = 1:mesh.lengthΓ
        g1 = mesh.Γ[iΓ]
        gΓ[iΓ] = ubdy[g1, 1]
    end
    
    if EL.lEL_Train

        #--------------------------------------------------------------------
        # TRAINING BRANCH — exact static condensation
        #--------------------------------------------------------------------
        # Step 4: B∂O∂τ = A∂O∂τ - Σ_iel A∂Ovo * inv(Avovo) * Avo∂τ
        for iel = 1:mesh.nelem
            LinearAlgebra.mul!(BC, inv(EL.Avovo[:,:,iel]), EL.Avo∂τ[:,:,iel])
            LinearAlgebra.mul!(@view(ABC[:,:,iel]), @view(EL.A∂Ovo[:,:,iel]), @view(BC[:,:]))

            
            
        end
        ∑el = sum(ABC, dims=3)
        EL.B∂O∂τ = EL.A∂O∂τ - ∑el

        # Step 5: Extract B∂O∂O and B∂O∂Γ from B∂O∂τ
        for i1=1:length(mesh.∂O)
            for i2=1:length(mesh.∂O)
                j2 = findall(x->x==mesh.∂O[i2], mesh.∂τ)[1]
                EL.B∂O∂O[i1, i2] = EL.B∂O∂τ[i1, j2]
            end
        end
        for iΓ = 1:mesh.lengthΓ
            jτ = findall(x->x==mesh.Γ[iΓ], mesh.∂τ)[1]
            EL.B∂O∂Γ[:, iΓ] .= EL.B∂O∂τ[:, jτ]
        end
        
        # Step 6: u∂O = -inv(B∂O∂O) * B∂O∂Γ * gΓ
        LinearAlgebra.mul!(BOΓg, EL.B∂O∂Γ, gΓ)
        u∂O      = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(mesh.length∂O))
        invB∂O∂O = inv(EL.B∂O∂O)
        LinearAlgebra.mul!(u∂O, -invB∂O∂O, BOΓg)
        
        # Step 7: uIo = -inv(AIoIo) * (AIo∂O * u∂O + AIoΓ * gΓ)
        AIoΓg  = similar(AIoΓ, (mesh.lengthIo,))
        AIou∂O = similar(AIoΓg)
        LinearAlgebra.mul!(AIoΓg,  AIoΓ,     gΓ)
        LinearAlgebra.mul!(AIou∂O, EL.AIo∂O, u∂O)
        uIo = similar(u∂O, mesh.lengthIo)
        LinearAlgebra.mul!(uIo, -invAIoIo, AIou∂O + AIoΓg)

        # Fill u
        for io = 1:mesh.lengthIo
            u[mesh.Io[io]] = uIo[io]
        end
        for io = 1:mesh.length∂O
            u[mesh.∂O[io]] = u∂O[io]
        end
        for io = 1:mesh.lengthΓ
            u[mesh.Γ[io]] = gΓ[io]
        end

        # Record input/output tensors for training.
        # Input  = local viscosity coefficient vector (a_j)_{j=1}^{(k+1)²}
        # Output = flattened T^{ie} = inv(Avovo)*Avovb for element 1
        #          (column-major, matching Julia's vec() convention)
        EL.input_tensor[:, isamp] .= avisc[:]
        
        for iel = 1:1
            Avbvo = transpose(EL.Avovb[:,:,iel])
            LinearAlgebra.mul!(EL.Tie, -inv(EL.Avovo[:,:,iel]), EL.Avovb[:,:,iel])
            #EL.Tie .= -EL.T2
            LinearAlgebra.mul!(@view(EL.T1[:,:]), @view(Avbvo[:,:]), -@view(EL.Tie[:,:]))
            EL.output_tensor[:, isamp] .= -vec(EL.Tie)
        end

        buffer = Vector{Vector{Float64}}()
        total_cols_written = 0
        
    #    write_MLtensor!(bufferin,  EL.input_tensor[:, isamp])
    #    write_MLtensor!(bufferout, EL.output_tensor[:, isamp])
        
    else

        #--------------------------------------------------------------------
        # INFERENCE BRANCH — NN-predicted T^{ie,nn} replaces exact T^{ie}
        #
        # NOTE: the NN must be trained to sufficient accuracy.
        # With cond(B∂O∂O) ≈ 69, the NN relative error on T^{ie} must be
        # well below 1/69 ≈ 1.4% to avoid corrupting the skeleton solve.
        # Use Nsamp >= 10000 and num_epochs >= 10000 in training.
        #--------------------------------------------------------------------
        @info " # RUN INFERENCE ......................................."

        # Load ONNX model
        sess        = ONNXRunTime.load_inference("./JX_NN_model.onnx")
        input_name  = first(sess.input_names)
        output_name = first(sess.output_names)

        # avisc is [1, ngl²] — same coefficient field for all elements.
        # Flatten to 1D once; reuse for every element.
        avisc_local = Float32.(vec(avisc))

        # Storage for per-element NN predictions (avoid calling NN twice)
        Tie_nn_all = zeros(Float64, nelintpoints, elnbdypoints, mesh.nelem)
        Tie_nn     = zeros(Float64, nelintpoints, elnbdypoints)
        M          = zeros(Float64, elnbdypoints, elnbdypoints)

        #--------------------------------------------------------------------
        # Steps 3+4: for each element get T^{ie,nn} from NN, update B∂τ∂τ
        #   B_{v^{ie,b}, v^{ie,b}} ← B_{v^{ie,b}, v^{ie,b}} - A_{v^{ie,b},v^{ie,o}} * T^{ie,nn}
        #--------------------------------------------------------------------
        EL.B∂τ∂τ .= EL.A∂τ∂τ   # initialise: B∂τ∂τ := A∂τ∂τ

        for iel = 1:mesh.nelem

            # Step 3: NN inference → flat prediction of T^{ie,nn}
            # Output is a column-major flattened (nelintpoints × elnbdypoints) matrix,
            # consistent with how vec(Tie) was written during training.
            y  = sess(Dict(input_name => Float32.(avisc)))
            ŷ  = y[output_name]

            Tie_nn .= reshape(vec(Float64.(ŷ)), nelintpoints, elnbdypoints)
            Tie_nn_all[:, :, iel] .= Tie_nn

            # Step 4: M = A_{v^{ie,b}, v^{ie,o}} * T^{ie,nn}  (elnbdypoints × elnbdypoints)
            Avbvo = transpose(EL.Avovb[:, :, iel])   # elnbdypoints × nelintpoints
            LinearAlgebra.mul!(M, Avbvo, Tie_nn)      # elnbdypoints × elnbdypoints

            # Update B∂τ∂τ[i', j'] -= M[i, j]
            # i', j' are the positions of v^{ie,b}(i), v^{ie,b}(j) in mesh.∂τ
            for i = 1:elnbdypoints
                vbi     = mesh.conn[iel, i]
                i_prime = findall(x -> x == vbi, mesh.∂τ)[1]
                for j = 1:elnbdypoints
                    vbj     = mesh.conn[iel, j]
                    j_prime = findall(x -> x == vbj, mesh.∂τ)[1]
                    EL.B∂τ∂τ[i_prime, j_prime] -= M[i, j]
                end
            end
        end

        #--------------------------------------------------------------------
        # Step 5: Extract B∂O∂O and B∂O∂Γ from the NN-assembled B∂τ∂τ
        #   (valid because ∂τ = ∂O ∪ Γ)
        #--------------------------------------------------------------------
        for i1 = 1:mesh.length∂O
            i_prime = findall(x -> x == mesh.∂O[i1], mesh.∂τ)[1]
            for i2 = 1:mesh.length∂O
                j_prime = findall(x -> x == mesh.∂O[i2], mesh.∂τ)[1]
                EL.B∂O∂O[i1, i2] = EL.B∂τ∂τ[i_prime, j_prime]
            end
            for iΓ = 1:mesh.lengthΓ
                j_prime = findall(x -> x == mesh.Γ[iΓ], mesh.∂τ)[1]
                EL.B∂O∂Γ[i1, iΓ] = EL.B∂τ∂τ[i_prime, j_prime]
            end
        end

        #--------------------------------------------------------------------
        # Step 6: u∂O = -inv(B∂O∂O) * B∂O∂Γ * gΓ
        #--------------------------------------------------------------------
        BOΓg_nn = zeros(mesh.length∂O)
        LinearAlgebra.mul!(BOΓg_nn, EL.B∂O∂Γ, gΓ)
        u∂O_nn  = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(mesh.length∂O))
        LinearAlgebra.mul!(u∂O_nn, -inv(EL.B∂O∂O), BOΓg_nn)

        #--------------------------------------------------------------------
        # Step 7: fill u_∂O and u_Γ so that u_∂τ is complete.
        # Must happen BEFORE uvb gather (Step 8a) since that reads from u.
        #--------------------------------------------------------------------
        for io = 1:mesh.length∂O
            u[mesh.∂O[io]] = u∂O_nn[io]
        end
        for io = 1:mesh.lengthΓ
            u[mesh.Γ[io]] = gΓ[io]
        end

        #--------------------------------------------------------------------
        # Step 8a: Gather u_{v^{ie,b}} for each element from the complete u_∂τ
        #--------------------------------------------------------------------
        uvb_nn = zeros(Float64, mesh.nelem, elnbdypoints)
        for iel = 1:mesh.nelem
            for ibdy = 1:elnbdypoints
                uvb_nn[iel, ibdy] = u[mesh.conn[iel, ibdy]]
            end
        end

        #--------------------------------------------------------------------
        # Step 8b: Local recovery using T^{ie,nn}
        #   u_{v^{ie,o}} = -T^{ie,nn} * u_{v^{ie,b}}   ∀ ie = 1:Nel
        #   This fills u_Io since Io = ⊕_{ie} v^{ie,o}
        #--------------------------------------------------------------------
        uvo_nn = zeros(Float64, nelintpoints)
        for iel = 1:mesh.nelem
            LinearAlgebra.mul!(uvo_nn, -Tie_nn_all[:, :, iel], uvb_nn[iel, :])
            for i = 1:nelintpoints
                u[mesh.conn[iel, elnbdypoints + i]] = uvo_nn[i]
            end
        end
        
        # u_∂τ (= u_∂O ∪ u_Γ) and u_Io are now both filled.
        # Since I = ∂τ ∪ Io, the full solution vector u is complete.
        
        @info " # RUN INFERENCE ....................................... END"
    end
end
