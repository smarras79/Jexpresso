include("./mesh/meshStructs.jl")

#-------------------------------------------------------------------------------------------
# Element learning matrices
#-------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_elemLearning{T <: AbstractFloat,
                                           dims0,  dims1,  dims2,  dims3, dims4,
                                           dims5,  dims6,  dims7,  dims8, dims9,
                                           dims10, dims11, dims12, dims13,
                                           dimsML1, dimsML2,
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

    # ML:
    input_tensor  = KernelAbstractions.zeros(backend, T, dimsML1)
    output_tensor = KernelAbstractions.zeros(backend, T, dimsML2)
    
end

function allocate_elemLearning(nelem, ngl, length∂O, length∂τ, lengthΓ, T, backend; Nsamp=1)

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
                                   dimsML1, dimsML2,
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
    
    #=EL = allocate_elemLearning(mesh.nelem, mesh.ngl,
    mesh.length∂O,
    mesh.length∂τ,
    mesh.lengthΓ,
    TFloat, inputs[:backend];
    Nsamp=inputs[:Nsamp])
    =#
    
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
    # inv(AiIoIo)
    invAIoIo = similar(EL.AIoIo)
    invAIoIo = inv(EL.AIoIo)
    
    dims = (mesh.lengthIo, mesh.lengthΓ)
    AIoΓ = similar(EL.AIoIo, dims);
    
    #------------------------------------------------------------------------
    # Eq. (13)
    #------------------------------------------------------------------------    
    #  B∂O∂τ[:,:] = A∂O∂τ - Sum_{iel} A∂Oᵥₒ[:,:,iel]*A⁻¹ᵥₒᵥₒ[:,:,iel]*Aᵥₒ∂τ[:,:,iel] -> A∂O∂τ - Sum_{iel}A⋅B⋅C
    #
    
    #
    # LOCAL VERSION (eq 13)
    #
    ABC = zeros(mesh.length∂O, mesh.length∂τ, mesh.nelem)
    BC  = zeros(size(EL.Avo∂τ)[1], size(EL.Avo∂τ)[2])
    for iel = 1:mesh.nelem
        
        # BC = A⁻¹ᵥₒᵥₒ[:,:,iel]⋅Aᵥₒ∂τ[:,:,iel]
        LinearAlgebra.mul!(BC, inv(EL.Avovo[:,:,iel]), EL.Avo∂τ[:,:,iel])
        
        # ABC = A∂Oᵥₒ[:,:,iel]⋅BC
        LinearAlgebra.mul!(@view(ABC[:,:,iel]), @view(EL.A∂Ovo[:,:,iel]), @view(BC[:,:]))
    end
    ∑el = similar(EL.A∂O∂τ)
    ∑el = sum(ABC, dims=3)
    EL.B∂O∂τ = EL.A∂O∂τ - ∑el # (13)
    
    #
    # WARNING: for large grids this double loop may be a bottleneck
    #
    for i1=1:length(mesh.∂O)      #row    B[i1][i2]        
        for i2=1:length(mesh.∂O)  #column B[i1][i2]
            
            j2 = findall(x->x==mesh.∂O[i2], mesh.∂τ)[1]
            EL.B∂O∂O[i1, i2] = EL.B∂O∂τ[i1, j2]
        end        
    end
    
    gΓ = zeros(mesh.lengthΓ)
    for iΓ = 1:mesh.lengthΓ
        g1=mesh.Γ[iΓ]
        
        jτ = findall(x->x==mesh.Γ[iΓ], mesh.∂τ)[1]
        EL.B∂O∂Γ[:, iΓ] .= EL.B∂O∂τ[:, jτ]

        gΓ[iΓ] = ubdy[g1, 1]
    end

    #------------------------------------------------------------------------
    # Eq. (11)
    #------------------------------------------------------------------------    
    #
    # B∂O∂Γ⋅gΓ
    #
    BOΓg = zeros(mesh.length∂O)
    LinearAlgebra.mul!(BOΓg, EL.B∂O∂Γ, gΓ)

    u∂O      = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(mesh.length∂O))
    invB∂O∂O = similar(EL.B∂O∂O)
    invB∂O∂O = inv(EL.B∂O∂O)
    #
    #  u∂O = -B⁻¹∂O∂O⋅(B∂OΓ⋅gΓ)
    #
    LinearAlgebra.mul!(u∂O, -invB∂O∂O, BOΓg)

    #
    #  AIo,Γ
    #
    AIoΓ = similar(A, (mesh.lengthIo, mesh.lengthΓ))
    for iΓ = 1:mesh.lengthΓ
        g1=mesh.Γ[iΓ]        
        for io = 1:mesh.lengthIo
            io1 = mesh.Io[io]
            
            AIoΓ[io, iΓ] = A[io1, g1]
        end
    end
    
    #
    # Eq (12)
    #
    AIoΓg = similar(AIoΓ, (mesh.lengthIo))
    LinearAlgebra.mul!(AIoΓg, AIoΓ, gΓ)

    AIou∂O = similar(AIoΓg)
    LinearAlgebra.mul!(AIou∂O, EL.AIo∂O, u∂O)

    dims = (mesh.lengthIo)
    uIo = similar(u∂O, dims)
    LinearAlgebra.mul!(uIo, -invAIoIo, (AIou∂O + AIoΓg))

    for io = 1:mesh.lengthIo # all internal without edges
        io1 = mesh.Io[io]
        u[io1] = uIo[io]
    end
    for io = 1:mesh.length∂O # all skeleton no boundaries
        io1 = mesh.∂O[io]
        u[io1] = u∂O[io]
    end
    for io = 1:mesh.lengthΓ # # all boundaries
        io1 = mesh.Γ[io]
        u[io1] = gΓ[io]
    end
    #@info u
    #@mystop
    # ∂O U Γ
    skeletonAndbdy                                              = zeros(Int64, length(mesh.∂τ))
    skeletonAndbdy[1:mesh.length∂O]                            .= mesh.∂O[:]
    skeletonAndbdy[mesh.length∂O+1:mesh.length∂O+mesh.lengthΓ] .= mesh.Γ[:]

    
    #    #for iel=1:mesh.nelem
    #        for iedge = 1:4
    #            @info iel, mesh.edge2pedge[iedge]
    #        end
    #    end
    #    @mystop
    #    bdy_edge_in_elem
    #    edge2pedge

    #
    # uvb ⊂ u∂τ  SHUKAI meeting
    #
    #=
    @info "mesh.cell_face_ids "
    @info mesh.cell_face_ids, size(mesh.cell_face_ids)
    @info "mesh.facet_cell_ids"
    @info mesh.facet_cell_ids, size(mesh.facet_cell_ids)

    for iedge=1:mesh.nedges
    @info " iedge ", iedge, "belongs to element ",  mesh.facet_cell_ids[iedge]
    end
    for iedge=1:mesh.nedges
    @info " edge ", iedge, " belong to elem ", mesh.edge_in_elem[iedge]
    end

    @info "-----"
    =#
    uvb = zeros(Float64, mesh.nelem, elnbdypoints)
    u∂τ = zeros(Float64, length(mesh.∂τ))
    for iskel = 1:length(mesh.∂τ)
        is = skeletonAndbdy[iskel]
        #        @info iskel, is
        u∂τ[iskel] = u[is]
        #x   @info iskel, is, u∂τ[iskel]
    end

    
    for iel=1:mesh.nelem
        #
        # 
        #
        ii = 1
        #
        # Aᵥₒᵥb
        #
        for isk = 1:length(mesh.∂τ)
            ipsk = skeletonAndbdy[isk]
            
            for ibdy = 1:elnbdypoints
                #jpb = mesh.conn[iel, j]
                uvb[iel, ibdy] = u∂τ[ipsk]
                # @info iel, ibdy, uvb[iel, ibdy]
            end
            ii += 1
        end
    end

    
    #=for iel = 1:mesh.nelem
    for ibdyel = 1:
    for j1=1:length(mesh.∂τ)
    jτ1 = skeletonAndbdy[j1]
    
    uvb[iel, ibdyel] = u∂τ[jτ1]
    end
    end=#

    
    
    #
    # ML: input/outpute tensors to use in training (?):
    #
    # 1. Set B∂τ∂τ := A∂τ∂τ
    #    
    T2  = zeros(size(EL.Avovo)[1], size(EL.Avovb)[2])
    T1  = zeros(size(EL.Avovb)[2], size(EL.Avovb)[2])
    Bie = similar(T2)
    
    # 2.c
    EL.input_tensor[:, isamp] .= avisc[:]

    # 2.d        
    for iel = 1:1 #mesh.nelem
        
        Avbvo = transpose(EL.Avovb[:,:,iel])
        
        # T2 = -A⁻¹ᵥₒᵥₒ[:,:,iel]⋅Avovb[:,:,iel]
        LinearAlgebra.mul!(T2, -inv(EL.Avovo[:,:,iel]), EL.Avovb[:,:,iel])
        Bie .= -T2
        
        # T1 = Avbvo[:,:,iel]⋅T2 = - Avbvo⋅A⁻¹ᵥₒᵥₒ⋅Avovb
        LinearAlgebra.mul!(@view(T1[:,:]), @view(Avbvo[:,:]), @view(T2[:,:]))
        
        # 2.e
        # Output tensor:
        EL.output_tensor[:, isamp] .= vec(Bie)  # Bie = -T2ie
    end

    #------------------------------------------------------------------------
    # Write input/output_bufferin.csv
    #------------------------------------------------------------------------
    write_MLtensor(@view(EL.input_tensor[:, isamp]), bufferin, total_cols_writtenin, "input_tensor.csv")
    write_MLtensor(@view(EL.output_tensor[:, isamp]), bufferout, total_cols_writtenout, "output_tensor.csv")
    #------------------------------------------------------------------------
    # END write input/output_buffer.csv
    #------------------------------------------------------------------------
    
    #!!!!!  ML code GOES HERE !!!!!
    
end

