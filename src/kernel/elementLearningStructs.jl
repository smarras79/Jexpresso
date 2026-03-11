include("./mesh/meshStructs.jl")
using SparseArrays
using LinearAlgebra

# =============================================================================
#  DIAGNOSTIC — call this before allocating to verify sizes
# =============================================================================
function diagnose_elemLearning(nelem, ngl, length∂O, length∂τ, lengthΓ, T; Nsamp=1)
    elnbdypoints = 4*(ngl-2) + 4
    nvo          = (ngl-2)^2
    k            = ngl - 1
    bytes        = sizeof(T)
    MB(dims)     = prod(Int64.(dims)) * bytes / 1024^2

    println("\n========== ElemLearning allocation diagnostic ==========")
    println("  nelem=$(nelem)  ngl=$(ngl)  nvo=$(nvo)  elnbdypoints=$(elnbdypoints)")
    println("  length∂O=$(length∂O)  length∂τ=$(length∂τ)  lengthΓ=$(lengthΓ)\n")

    println("  --- STRUCT (persistent) ---")
    struct_arrays = [
        ("Avovo/AIoIo ×2",   (nvo, nvo,          nelem), 2),
        ("Avovb/Avo∂τ/AIo∂τ",(nvo, elnbdypoints, nelem), 3),
        ("Avo∂O/AIo∂O",      (nvo, elnbdypoints, nelem), 2),
        ("A∂Ovo/A∂OIo",      (elnbdypoints, nvo, nelem), 2),
        ("T1",               (elnbdypoints, elnbdypoints), 1),
        ("T2/Tie",           (nvo, elnbdypoints),          2),
        ("input_tensor",     ((k+1)^2, Nsamp),             1),
        ("output_tensor",    (4*k*(k-1)^2, Nsamp),         1),
    ]
    struct_total = 0.0
    for (name, dims, count) in struct_arrays
        m = MB(dims) * count
        struct_total += m
        @printf("  %-28s  %s ×%d  →  %8.1f MB\n", name, string(dims), count, m)
    end
    @printf("  %-28s  %8.1f MB\n", "STRUCT TOTAL", struct_total)

    println("\n  --- LOCAL (sparse, per call) ---")
    nnz = nelem * elnbdypoints^2
    local_arrays = [
        ("A∂τ∂τ / B∂τ∂τ (sparse)",   nnz, 2),
        ("A∂O∂τ / B∂O∂τ (sparse)",   nnz, 2),
        ("B∂O∂O (sparse)",            nnz, 1),
        ("B∂O∂Γ (dense)",             (length∂O, lengthΓ), 1),
    ]
    local_total = 0.0
    for (name, dims, count) in local_arrays
        m = (isa(dims, Int) ? dims * bytes / 1024^2 : MB(dims)) * count
        local_total += m
        @printf("  %-28s  NNZ≤%s ×%d  →  %8.1f MB\n",
                name, isa(dims,Int) ? string(dims) : string(dims), count, m)
    end
    @printf("  %-28s  %8.1f MB\n", "LOCAL TOTAL", local_total)

    println("\n  --- WHAT WOULD HAVE CRASHED (dense) ---")
    crash_arrays = [
        ("A∂τ∂τ dense",  (length∂τ, length∂τ)),
        ("B∂τ∂τ dense",  (length∂τ, length∂τ)),
        ("A∂O∂τ dense",  (length∂O, length∂τ)),
        ("B∂O∂τ dense",  (length∂O, length∂τ)),
        ("B∂O∂O dense",  (length∂O, length∂O)),
    ]
    for (name, dims) in crash_arrays
        @printf("  %-28s  %s  →  %8.1f MB  ◄◄◄\n", name, string(dims), MB(dims))
    end
    println("========================================================\n")
end


# =============================================================================
#  ROOT CAUSE (complete picture)
#
#  FOUR categories of arrays were wrong:
#
#  1. Per-element with global ∂τ/∂O dimension — FIXED in previous iteration:
#       Avo∂τ, AIo∂τ, Avo∂O, AIo∂O, A∂Ovo, A∂OIo
#       Fix: use elnbdypoints instead of length∂τ / length∂O
#
#  2. Global skeleton matrices stored in struct — FIXED HERE:
#       A∂τ∂τ (6 GB), B∂τ∂τ (6 GB), A∂O∂τ (5.8 GB), B∂O∂τ (5.8 GB),
#       B∂O∂O (5.5 GB), B∂O∂Γ (246 MB)
#       Fix: remove from struct entirely; allocate LOCALLY as sparse
#       submatrices of A inside elementLearning_Axb! each call.
#       Since A is a SparseMatrixCSC, A[mesh.∂τ, mesh.∂τ] is sparse (11 MB).
#       NNZ ≤ nelem × elnbdypoints² = 2500 × 576 = 1,440,000.
#
#  Key insight: these matrices have NO persistent state between calls.
#  Section 2 of elementLearning_Axb! fills them from A every single call.
#  There is no reason to store them in the struct at all.
# =============================================================================


# =============================================================================
#  St_elemLearning — struct holds ONLY per-element blocks and ML tensors.
#  All global skeleton matrices (A∂τ∂τ, B∂τ∂τ, A∂O∂τ, B∂O∂τ, B∂O∂O, B∂O∂Γ)
#  are removed.  They are local sparse variables inside elementLearning_Axb!.
# =============================================================================
Base.@kwdef mutable struct St_elemLearning{T <: AbstractFloat,
                                            dims0,
                                            dims_vovo,
                                            dims_∂Ovo,
                                            dims_vovb,
                                            dims_T2,
                                            dims_T1,
                                            dimsML1,
                                            dimsML2,
                                            lELSample,
                                            backend}

    # ── Per-element: interior × interior  (nvo × nvo × nelem) ────────────────
    Avovo   = KernelAbstractions.zeros(backend, T, dims_vovo)
    AIoIo   = KernelAbstractions.zeros(backend, T, dims_vovo)

    # ── Per-element: interior × local-boundary  (nvo × elnbdy × nelem) ───────
    # All use elnbdypoints as boundary dim — NOT length∂τ / length∂O
    Avovb   = KernelAbstractions.zeros(backend, T, dims_vovb)
    Avo∂O   = KernelAbstractions.zeros(backend, T, dims_vovb)   # + ∂O mask at runtime
    Avo∂τ   = KernelAbstractions.zeros(backend, T, dims_vovb)   # ≡ Avovb
    AIo∂τ   = KernelAbstractions.zeros(backend, T, dims_vovb)
    AIo∂O   = KernelAbstractions.zeros(backend, T, dims_vovb)   # + ∂O mask at runtime

    # ── Per-element: local-boundary × interior  (elnbdy × nvo × nelem) ───────
    A∂Ovo   = KernelAbstractions.zeros(backend, T, dims_∂Ovo)
    A∂OIo   = KernelAbstractions.zeros(backend, T, dims_∂Ovo)

    # ── Local temporaries  (no nelem, no global DOF dimensions) ──────────────
    T1      = KernelAbstractions.zeros(backend, T, dims_T1)   # elnbdy × elnbdy
    T2      = KernelAbstractions.zeros(backend, T, dims_T2)   # nvo × elnbdy
    Tie     = KernelAbstractions.zeros(backend, T, dims_T2)   # nvo × elnbdy

    lEL_Sample = lELSample

    # ── ML tensors ────────────────────────────────────────────────────────────
    input_tensor  = KernelAbstractions.zeros(backend, T, dimsML1)
    output_tensor = KernelAbstractions.zeros(backend, T, dimsML2)
end


# =============================================================================
#  allocate_elemLearning
# =============================================================================
function allocate_elemLearning(nelem, ngl, length∂O, length∂τ, lengthΓ,
                                T, backend;
                                Nsamp=1, lEL_Sample=false)
    elnbdypoints = 4*(ngl-2) + 4
    nvo          = (ngl-2)^2
    k            = ngl - 1

    dims_vovo  = (nvo,          nvo,          nelem)
    dims_vovb  = (nvo,          elnbdypoints, nelem)
    dims_∂Ovo  = (elnbdypoints, nvo,          nelem)
    dims_T1    = (elnbdypoints, elnbdypoints)
    dims_T2    = (nvo,          elnbdypoints)
    dimsML1    = ((k+1)^2,        Nsamp)
    dimsML2    = (4*k*(k-1)^2,    Nsamp)
    dims0      = (nelem, 2)

    return St_elemLearning{T,
                            dims0,
                            dims_vovo,
                            dims_∂Ovo,
                            dims_vovb,
                            dims_T2,
                            dims_T1,
                            dimsML1,
                            dimsML2,
                            lEL_Sample,
                            backend}()
end


# =============================================================================
#  elementLearning_Axb!
#
#  Global skeleton matrices are now LOCAL sparse variables:
#
#    A_∂τ∂τ  = A[mesh.∂τ, mesh.∂τ]         SparseMatrixCSC  ~11 MB
#    B_∂τ∂τ  = copy(A_∂τ∂τ)                mutable working copy
#    A_∂O∂τ  = A[mesh.∂O, mesh.∂τ]         SparseMatrixCSC  ~11 MB
#    B_∂O∂τ  computed via scatter
#    B_∂O∂O  = B_∂τ∂τ[∂O_in_∂τ, ∂O_in_∂τ]  sparse submatrix
#    B_∂O∂Γ  = B_∂τ∂τ[∂O_in_∂τ, Γ_in_∂τ]   sparse submatrix
#
#  All of the above are freed when the function returns.
#  A must be a SparseMatrixCSC (from build_laplace_matrix / assemble_and_DSS).
# =============================================================================
function elementLearning_Axb!(u, uaux, mesh::St_mesh,
                               A::SparseMatrixCSC,
                               ubdy, EL,
                               avisc,
                               bufferin, bufferout,
                               BOΓg, gΓ;
                               isamp=1,
                               total_cols_writtenin=0,
                               total_cols_writtenout=0)

    mesh.lengthO  = mesh.length∂O + mesh.lengthIo
    nelintpoints  = (mesh.ngl - 2)^2
    nelpoints     = size(mesh.conn)[2]
    elnbdypoints  = nelpoints - nelintpoints

    # ── DOF → position lookup tables  (built once, O(length∂τ)) ─────────────
    ∂O_pos = Dict{Int,Int}(mesh.∂O[i] => i for i in 1:mesh.length∂O)
    ∂τ_pos = Dict{Int,Int}(mesh.∂τ[j] => j for j in 1:mesh.length∂τ)

    # Per-element scratch masks (small, reused each element)
    conn_∂O_idx = zeros(Int, elnbdypoints)
    conn_∂τ_idx = zeros(Int, elnbdypoints)

    # =========================================================================
    # SECTION 1: Extract global skeleton matrices as SPARSE submatrices of A.
    #
    # Julia sparse fancy indexing: A[idx, idx] returns SparseMatrixCSC.
    # For 50×50 N=6: each ~11 MB instead of ~6 GB.
    # These are LOCAL variables — freed on function return.
    # =========================================================================
    A_∂τ∂τ = A[mesh.∂τ, mesh.∂τ]     # SparseMatrixCSC (length∂τ × length∂τ)
    A_∂O∂τ = A[mesh.∂O, mesh.∂τ]     # SparseMatrixCSC (length∂O × length∂τ)
    # Note: A[∂O, ∂τ] is already the right submatrix — no manual loop needed.
    # The skeleton portion of A is exactly what was assembled by DSS.

    # =========================================================================
    # SECTION 2: Fill per-element 3D blocks from A
    # All boundary dimensions run over [1:elnbdypoints].
    # =========================================================================
    for iel = 1:mesh.nelem
        for j = 1:elnbdypoints
            gnode          = mesh.conn[iel, j]
            conn_∂O_idx[j] = get(∂O_pos, gnode, 0)
            conn_∂τ_idx[j] = get(∂τ_pos, gnode, 0)
        end

        ii = 1
        for i = elnbdypoints+1:nelpoints
            ipo = mesh.conn[iel, i]

            for j = 1:elnbdypoints
                gj = mesh.conn[iel, j]
                # Avovb / Avo∂τ / AIo∂τ  (all local bdy nodes are in ∂τ)
                val = A[ipo, gj]
                EL.Avovb[ii, j, iel]  = val
                EL.Avo∂τ[ii, j, iel]  = val
                EL.AIo∂τ[ii, j, iel]  = val

                # Avo∂O / AIo∂O  (only if gj ∈ ∂O)
                if conn_∂O_idx[j] != 0
                    val2 = A[ipo, gj]
                    EL.Avo∂O[ii, j, iel]  = val2
                    EL.AIo∂O[ii, j, iel]  = val2
                end

                # A∂Ovo / A∂OIo  (only if gj ∈ ∂O; row = gj, col = ipo)
                if conn_∂O_idx[j] != 0
                    val3 = A[gj, ipo]
                    EL.A∂Ovo[j, ii, iel]  = val3
                    EL.A∂OIo[j, ii, iel]  = val3
                end
            end

            # Avovo / AIoIo  (interior × interior)
            jj = 1
            for j = elnbdypoints+1:nelpoints
                val = A[ipo, mesh.conn[iel, j]]
                EL.Avovo[ii, jj, iel] = val
                EL.AIoIo[ii, jj, iel] = val
                jj += 1
            end

            ii += 1
        end
    end

    # =========================================================================
    # SECTION 3: Gather Dirichlet data gΓ
    # =========================================================================
    for iΓ = 1:mesh.lengthΓ
        gΓ[iΓ] = ubdy[mesh.Γ[iΓ], 1]
    end

    # =========================================================================
    # TRAINING BRANCH
    # =========================================================================
    if EL.lEL_Sample
        # Build B_∂O∂τ = A_∂O∂τ - Σ_iel A∂Ovo * inv(Avovo) * Avo∂τ  (sparse ← sparse - dense_scatter)
        # We accumulate the element corrections into a dense (length∂O × length∂τ) delta,
        # but only the nonzero ∂O×∂τ entries touched by mesh.conn — still sparse in practice.
        ΔB = spzeros(eltype(A), mesh.length∂O, mesh.length∂τ)

        for iel = 1:mesh.nelem
            for j = 1:elnbdypoints
                conn_∂τ_idx[j] = get(∂τ_pos, mesh.conn[iel, j], 0)
                conn_∂O_idx[j] = get(∂O_pos, mesh.conn[iel, j], 0)
            end

            invAvovo = inv(EL.Avovo[:, :, iel])
            BC_local = invAvovo * EL.Avo∂τ[:, :, iel]   # nvo × elnbdy

            for j_loc = 1:elnbdypoints
                jτ = conn_∂τ_idx[j_loc];  jτ == 0 && continue
                for i_loc = 1:elnbdypoints
                    io = conn_∂O_idx[i_loc];  io == 0 && continue
                    s = zero(eltype(A))
                    for ii = 1:nelintpoints
                        s += EL.A∂Ovo[i_loc, ii, iel] * BC_local[ii, j_loc]
                    end
                    ΔB[io, jτ] += s
                end
            end
        end

        B_∂O∂τ = A_∂O∂τ - ΔB   # SparseMatrixCSC (length∂O × length∂τ)

        # Extract ∂O and Γ position indices within ∂τ (built once)
        ∂O_in_∂τ = [∂τ_pos[mesh.∂O[i]]  for i in 1:mesh.length∂O]
        Γ_in_∂τ  = [∂τ_pos[mesh.Γ[iΓ]]  for iΓ in 1:mesh.lengthΓ]

        B_∂O∂O = B_∂O∂τ[:, ∂O_in_∂τ]   # SparseMatrixCSC (length∂O × length∂O)
        B_∂O∂Γ = B_∂O∂τ[:, Γ_in_∂τ]    # SparseMatrixCSC (length∂O × lengthΓ)

        # Step 6: u∂O = -inv(B∂O∂O) * B∂O∂Γ * gΓ
        # B_∂O∂Γ and B_∂O∂O are SparseMatrixCSC — use sparse * and sparse \.
        # NEVER call Matrix() on these: B_∂O∂O is length∂O × length∂O = 5.5 GB dense.
        BOΓg  = B_∂O∂Γ * gΓ        # sparse × dense vector → dense vector  (cheap)
        u∂O   = -(B_∂O∂O \ BOΓg)   # sparse direct solve via UMFPACK  (no densification)
print
        for io = 1:mesh.length∂O;  u[mesh.∂O[io]] = u∂O[io];  end
        for io = 1:mesh.lengthΓ;   u[mesh.Γ[io]]  = gΓ[io];   end

        # Step 7: per-element interior recovery
        AIoΓ_ie   = zeros(eltype(A), nelintpoints, mesh.lengthΓ)
        AIou∂O_ie = zeros(eltype(A), nelintpoints)
        AIoΓg_ie  = zeros(eltype(A), nelintpoints)
        uvo_ie    = zeros(eltype(A), nelintpoints)

        for iel = 1:mesh.nelem
            for j = 1:elnbdypoints
                conn_∂O_idx[j] = get(∂O_pos, mesh.conn[iel, j], 0)
            end

            for iΓ = 1:mesh.lengthΓ
                g1 = mesh.Γ[iΓ]
                for ii = 1:nelintpoints
                    AIoΓ_ie[ii, iΓ] = A[mesh.conn[iel, elnbdypoints+ii], g1]
                end
            end

            fill!(AIou∂O_ie, zero(eltype(A)))
            for j_loc = 1:elnbdypoints
                io = conn_∂O_idx[j_loc];  io == 0 && continue
                for ii = 1:nelintpoints
                    AIou∂O_ie[ii] += EL.AIo∂O[ii, j_loc, iel] * u∂O[io]
                end
            end

            LinearAlgebra.mul!(AIoΓg_ie, AIoΓ_ie, gΓ)
            rhs_ie = AIou∂O_ie .+ AIoΓg_ie
            LinearAlgebra.mul!(uvo_ie, -inv(EL.AIoIo[:, :, iel]), rhs_ie)

            for ii = 1:nelintpoints
                u[mesh.conn[iel, elnbdypoints+ii]] = uvo_ie[ii]
            end
        end

        # ML tensor recording (element 1 only)
        EL.input_tensor[:, isamp] .= vec(avisc)
        let iel = 1
            LinearAlgebra.mul!(EL.Tie, -inv(EL.Avovo[:, :, iel]), EL.Avovb[:, :, iel])
            LinearAlgebra.mul!(EL.T1, transpose(EL.Avovb[:, :, iel]), -EL.Tie)
            EL.output_tensor[:, isamp] .= -vec(EL.Tie)
        end
        write_MLtensor!(bufferin,  EL.input_tensor[:,  isamp])
        write_MLtensor!(bufferout, EL.output_tensor[:, isamp])
        
    # =========================================================================
    # INFERENCE BRANCH
    # =========================================================================
    else
        println(GREEN_FG(string(" # INFERENCE — solution stored in u ..........")))
                
        sess        = ONNXRunTime.load_inference("./JX_NN_model.onnx")
        input_name  = first(sess.input_names)
        output_name = first(sess.output_names)

        Tie_nn_all = zeros(Float64, nelintpoints, elnbdypoints, mesh.nelem)
        Tie_nn     = zeros(Float64, nelintpoints, elnbdypoints)
        M          = zeros(Float64, elnbdypoints, elnbdypoints)

        # B∂τ∂τ starts as the sparse skeleton submatrix of A, then gets
        # element corrections subtracted in-place.  ~11 MB, not 6 GB.
        B_∂τ∂τ = copy(A_∂τ∂τ)   # SparseMatrixCSC — mutable working copy

        for iel = 1:mesh.nelem
            for j = 1:elnbdypoints
                conn_∂τ_idx[j] = get(∂τ_pos, mesh.conn[iel, j], 0)
            end

            y      = sess(Dict(input_name => Float32.(avisc)))
            ŷ      = y[output_name]
            Tie_nn .= reshape(vec(Float64.(ŷ)), nelintpoints, elnbdypoints)
            Tie_nn_all[:, :, iel] .= Tie_nn

            LinearAlgebra.mul!(M, transpose(EL.Avovb[:, :, iel]), Tie_nn)

            for i = 1:elnbdypoints
                i_prime = conn_∂τ_idx[i];  i_prime == 0 && continue
                for j = 1:elnbdypoints
                    j_prime = conn_∂τ_idx[j];  j_prime == 0 && continue
                    B_∂τ∂τ[i_prime, j_prime] -= M[i, j]
                end
            end
        end

        # Extract B∂O∂O and B∂O∂Γ as sparse submatrices
        ∂O_in_∂τ = [∂τ_pos[mesh.∂O[i]]  for i in 1:mesh.length∂O]
        Γ_in_∂τ  = [∂τ_pos[mesh.Γ[iΓ]]  for iΓ in 1:mesh.lengthΓ]

        B_∂O∂O = B_∂τ∂τ[∂O_in_∂τ, ∂O_in_∂τ]   # SparseMatrixCSC
        B_∂O∂Γ = B_∂τ∂τ[∂O_in_∂τ, Γ_in_∂τ]    # SparseMatrixCSC

     
        # Step 6: u∂O = -inv(B∂O∂O) * B∂O∂Γ * gΓ
        # NEVER call Matrix() on B_∂O∂O — that densifies 26901×26901 = 5.5 GB.
        # sparse * dense vector → dense vector, sparse \ dense vector → dense vector.
        BOΓg_nn = B_∂O∂Γ * gΓ           # SparseMatrixCSC × Vector → Vector  (cheap)
        u∂O_nn  = -(B_∂O∂O \ BOΓg_nn)   # UMFPACK sparse direct solve
        
        # Step 7
        for io = 1:mesh.length∂O;  u[mesh.∂O[io]] = u∂O_nn[io];  end
        for io = 1:mesh.lengthΓ;   u[mesh.Γ[io]]  = gΓ[io];      end

        # Step 8a: gather boundary solution for each element
        uvb_nn = zeros(Float64, mesh.nelem, elnbdypoints)
        for iel = 1:mesh.nelem
            for ibdy = 1:elnbdypoints
                uvb_nn[iel, ibdy] = u[mesh.conn[iel, ibdy]]
            end
        end

        # Step 8b: local recovery  u_vo = -T^{ie,nn} * u_vb
        uvo_nn = zeros(Float64, nelintpoints)
        for iel = 1:mesh.nelem
            LinearAlgebra.mul!(uvo_nn, -Tie_nn_all[:, :, iel], uvb_nn[iel, :])
            for i = 1:nelintpoints
                u[mesh.conn[iel, elnbdypoints+i]] = uvo_nn[i]
            end
        end

        println(GREEN_FG(string(" # INFERENCE — solution stored in u .......... DONE")))
        
    end
    # A_∂τ∂τ, A_∂O∂τ, B_∂τ∂τ, B_∂O∂τ, B_∂O∂O, B_∂O∂Γ all freed here
end


# =============================================================================
#  write_MLtensor helpers (unchanged)
# =============================================================================
function write_MLtensor!(buffer::Vector{Vector{Float64}}, tensor_column::Vector{Float64})
    push!(buffer, copy(tensor_column))
end

function flush_MLtensor!(buffer::Vector{Vector{Float64}}, total_cols_written, fname)
    isempty(buffer) && return total_cols_written
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
    empty!(buffer)
    return total_cols_written + ncols
end
