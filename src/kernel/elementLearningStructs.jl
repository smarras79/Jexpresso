include("./mesh/meshStructs.jl")

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

function elementLearning_Axb!(u, uaux, mesh::St_mesh,
                              A, ubdy, EL,
                              avisc, 
                              bufferin, bufferout;
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

    # inv(AIoIo) — needed by training branch only, but cheap to keep here
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
    gΓ = zeros(mesh.lengthΓ)
    for iΓ = 1:mesh.lengthΓ
        g1 = mesh.Γ[iΓ]
        gΓ[iΓ] = ubdy[g1, 1]
    end

    #
    # ML: input/output tensors to use in training / inference
    #
    T2  = zeros(size(EL.Avovo)[1], size(EL.Avovb)[2])
    T1  = zeros(size(EL.Avovb)[2], size(EL.Avovb)[2])
    Tie = similar(T2)
    
    if EL.lEL_Train

        #--------------------------------------------------------------------
        # TRAINING BRANCH — exact static condensation (sec:static_alg)
        #--------------------------------------------------------------------

        # Step 4: B∂O∂τ = A∂O∂τ - Σ_iel A∂Ovo * inv(Avovo) * Avo∂τ  (Eq. 13)
        ABC = zeros(mesh.length∂O, mesh.length∂τ, mesh.nelem)
        BC  = zeros(size(EL.Avo∂τ)[1], size(EL.Avo∂τ)[2])
        for iel = 1:mesh.nelem
            LinearAlgebra.mul!(BC, inv(EL.Avovo[:,:,iel]), EL.Avo∂τ[:,:,iel])
            LinearAlgebra.mul!(@view(ABC[:,:,iel]), @view(EL.A∂Ovo[:,:,iel]), @view(BC[:,:]))
        end
        ∑el = sum(ABC, dims=3)
        EL.B∂O∂τ = EL.A∂O∂τ - ∑el  # (13)

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

        # Step 6: u∂O = -inv(B∂O∂O) * B∂O∂Γ * gΓ  (Eq. 11)
        BOΓg = zeros(mesh.length∂O)
        LinearAlgebra.mul!(BOΓg, EL.B∂O∂Γ, gΓ)
        u∂O      = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(mesh.length∂O))
        invB∂O∂O = inv(EL.B∂O∂O)
        LinearAlgebra.mul!(u∂O, -invB∂O∂O, BOΓg)

        # Step 7 (via Eq. 12): uIo = -inv(AIoIo) * (AIo∂O * u∂O + AIoΓ * gΓ)
        AIoΓg  = similar(AIoΓ, (mesh.lengthIo,))
        AIou∂O = similar(AIoΓg)
        LinearAlgebra.mul!(AIoΓg,  AIoΓ,     gΓ)
        LinearAlgebra.mul!(AIou∂O, EL.AIo∂O, u∂O)
        uIo = similar(u∂O, mesh.lengthIo)
        LinearAlgebra.mul!(uIo, -invAIoIo, AIou∂O + AIoΓg)

        # Fill u (Steps 2, 6, 7)
        for io = 1:mesh.lengthIo
            u[mesh.Io[io]] = uIo[io]
        end
        for io = 1:mesh.length∂O
            u[mesh.∂O[io]] = u∂O[io]
        end
        for io = 1:mesh.lengthΓ
            u[mesh.Γ[io]] = gΓ[io]
        end

        # Record input/output tensors for training
        EL.input_tensor[:, isamp] .= avisc[:]
        for iel = 1:1
            Avbvo = transpose(EL.Avovb[:,:,iel])
            LinearAlgebra.mul!(T2, -inv(EL.Avovo[:,:,iel]), EL.Avovb[:,:,iel]) 
            Tie .= -T2
            LinearAlgebra.mul!(@view(T1[:,:]), @view(Avbvo[:,:]), @view(T2[:,:]))
            EL.output_tensor[:, isamp] .= vec(Tie)
        end

        write_MLtensor(@view(EL.input_tensor[:, isamp]), bufferin,  total_cols_writtenin,  "input_tensor.csv")
        write_MLtensor(@view(EL.output_tensor[:, isamp]), bufferout, total_cols_writtenout, "output_tensor.csv")

    else

        #--------------------------------------------------------------------
        # INFERENCE BRANCH — NN-predicted T^{ie,nn} replaces exact T^{ie}
        # Follows sec:static_alg Steps 2,4,5,6,7,8 with Step 3 replaced by NN
        #--------------------------------------------------------------------
        @info "RUN INFERENCE"

        # Load ONNX model
        sess        = ONNXRunTime.load_inference("./JX_NN_model.onnx")
        input_name  = first(sess.input_names)
        output_name = first(sess.output_names)

        # avisc is [1, ngl²] — same coefficient field for all elements.
        # Flatten to 1D once here; reuse for every element.
        # (Per the algorithm, step 5.a retrieves (a_j)_{j=1}^{(k+1)²} for element ie.
        #  Here avisc already contains those values, uniform across elements.)
        avisc_local = Float32.(vec(avisc))

        # Storage for all per-element NN predictions (avoids running NN twice)
        Tie_nn_all = zeros(Float64, nelintpoints, elnbdypoints, mesh.nelem)
        Tie_nn     = zeros(Float64, nelintpoints, elnbdypoints)
        M          = zeros(Float64, elnbdypoints, elnbdypoints)

        #--------------------------------------------------------------------
        # Step 3 (NN) + Step 4: for each element get T^{ie,nn}, update B∂τ∂τ
        #   B_{v^{ie,b}, v^{ie,b}} ←  B_{v^{ie,b}, v^{ie,b}} - A_{v^{ie,b},v^{ie,o}} * T^{ie,nn}
        #--------------------------------------------------------------------
        EL.B∂τ∂τ .= EL.A∂τ∂τ   # initialise B∂τ∂τ := A∂τ∂τ

        for iel = 1:mesh.nelem

            # Step 5.b: run NN → flat prediction of T^{ie,nn}
            #y  = sess(Dict(input_name => avisc_local))
            y  = sess(Dict(input_name => Float32.(avisc)))
            ŷ  = y[output_name]

            # Step 5.c: reshape to nelintpoints × elnbdypoints and cache
            Tie_nn .= reshape(Float64.(ŷ), nelintpoints, elnbdypoints)
            Tie_nn_all[:, :, iel] .= Tie_nn

            # Step 4 (element contribution):
            #   M = A_{v^{ie,b}, v^{ie,o}} * T^{ie,nn}   (elnbdypoints × elnbdypoints)
            Avbvo = transpose(EL.Avovb[:, :, iel])   # elnbdypoints × nelintpoints
            LinearAlgebra.mul!(M, Avbvo, Tie_nn)      # elnbdypoints × elnbdypoints

            #   B∂τ∂τ[i', j'] -= M[i, j]  where i',j' are positions of v^{ie,b}(i,j) in ∂τ
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
        #   (possible because ∂τ = ∂O ∪ Γ)
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
        # Step 6: u∂O = -inv(B∂O∂O) * B∂O∂Γ * gΓ  (Eq. glb_red_sol)
        #--------------------------------------------------------------------
        BOΓg_nn = zeros(mesh.length∂O)
        LinearAlgebra.mul!(BOΓg_nn, EL.B∂O∂Γ, gΓ)
        u∂O_nn  = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(mesh.length∂O))
        LinearAlgebra.mul!(u∂O_nn, -inv(EL.B∂O∂O), BOΓg_nn)

        #--------------------------------------------------------------------
        # Step 2 & 7: fill u_Γ and u_∂O so that u_∂τ = u_∂O ∪ u_Γ is complete.
        # This must happen BEFORE uvb gather and local recovery (Step 8)
        # because uvb reads from u.
        #--------------------------------------------------------------------
        for io = 1:mesh.length∂O
            u[mesh.∂O[io]] = u∂O_nn[io]
        end
        for io = 1:mesh.lengthΓ
            u[mesh.Γ[io]] = gΓ[io]
        end

        #--------------------------------------------------------------------
        # Step 8a: Gather u_{v^{ie,b}} for each element from the now-complete u_∂τ
        #--------------------------------------------------------------------
        uvb_nn = zeros(Float64, mesh.nelem, elnbdypoints)
        for iel = 1:mesh.nelem
            for ibdy = 1:elnbdypoints
                uvb_nn[iel, ibdy] = u[mesh.conn[iel, ibdy]]
            end
        end

        #--------------------------------------------------------------------
        # Step 8b: Local recovery using T^{ie,nn}  (Eq. sol_rec_elem)
        #   u_{v^{ie,o}} = -T^{ie,nn} * u_{v^{ie,b}}   ∀ ie = 1:Nel
        # This updates u_I^o since I^o = ⊕_{ie} v^{ie,o}
        #--------------------------------------------------------------------
        uvo_nn = zeros(Float64, nelintpoints)
        for iel = 1:mesh.nelem
            LinearAlgebra.mul!(uvo_nn, -Tie_nn_all[:, :, iel], uvb_nn[iel, :])
            for i = 1:nelintpoints
                u[mesh.conn[iel, elnbdypoints + i]] = uvo_nn[i]
            end
        end

        # After Steps 7 and 8, u_∂τ and u_I^o are both filled.
        # Since I = ∂τ ∪ I^o, the entire solution vector u is now complete.
        @info "INFERENCE COMPLETE — solution stored in u"
    end
end

# Point evaluation: interpolate at a single point (ξ, η)
