include("./mesh/meshStructs.jl")
using SparseArrays
using BenchmarkTools
using Profile, PProf   # PProf gives a flame graph in the browser

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

struct EL_InferBuffers
    # ONNX staging
    avisc_f32    :: Matrix{Float32}    # (1,        nfeatures)  — Case A (shared)
    avisc_batch  :: Matrix{Float32}    # (nelem,    nfeatures)  — Case B (per-element)
    ŷ_f64_buf    :: Vector{Float64}    # (nout,)
    # Element assembly
    Tie_nn_all   :: Array{Float64, 3}  # (nelintpoints, elnbdypoints, nelem)
    conn_∂τ_idx  :: Vector{Int}        # (elnbdypoints,)
    M            :: Matrix{Float64}    # (elnbdypoints, elnbdypoints)
    B_∂τ∂τ       :: SparseMatrixCSC{Float64, Int32}   # mutable copy of A_∂τ∂τ skeleton
    # Index maps
    ∂O_in_∂τ     :: Vector{Int}        # (length∂O,)
    Γ_in_∂τ      :: Vector{Int}        # (lengthΓ,)
    # Gather / recovery
    uvb_nn       :: Matrix{Float64}    # (nelem, elnbdypoints)
    uvo_nn       :: Vector{Float64}    # (nelintpoints,)
end

"""
        EL_InferBuffers(mesh, EL, A_∂τ∂τ, nfeatures, nelintpoints, elnbdypoints)

    Allocate all working arrays for `elementLearning_infer!` once.
    Pass the returned struct to every subsequent call — nothing is allocated inside
    the hot path.
    """
function EL_InferBuffers(mesh, A_∂τ∂τ::SparseMatrixCSC,
                         nfeatures::Int, nelintpoints::Int, elnbdypoints::Int)
    nout = nelintpoints * elnbdypoints
    return EL_InferBuffers(
        Matrix{Float32}(undef, 1,          nfeatures),           # avisc_f32
        Matrix{Float32}(undef, mesh.nelem, nfeatures),           # avisc_batch
        Vector{Float64}(undef, nout),                            # ŷ_f64_buf
        Array{Float64, 3}(undef, nelintpoints, elnbdypoints, mesh.nelem),  # Tie_nn_all
        Vector{Int}(undef, elnbdypoints),                        # conn_∂τ_idx
        Matrix{Float64}(undef, elnbdypoints, elnbdypoints),      # M
        copy(A_∂τ∂τ),                                            # B_∂τ∂τ
        Vector{Int}(undef, mesh.length∂O),                       # ∂O_in_∂τ
        Vector{Int}(undef, mesh.lengthΓ),                        # Γ_in_∂τ
        Matrix{Float64}(undef, mesh.nelem, elnbdypoints),        # uvb_nn
        Vector{Float64}(undef, nelintpoints),                    # uvo_nn
    )
end

struct EL_WorkBuffers
    # Section 2 scratch
    conn_∂O_idx  :: Vector{Int}
    conn_∂τ_idx  :: Vector{Int}

    # Sampling block scratch
    ΔB           :: SparseMatrixCSC{Float64, Int32}   # (length∂O × length∂τ)
    invAvovo_buf :: Matrix{Float64}                    # (nelintpoints × nelintpoints)
    BC_local     :: Matrix{Float64}                    # (nelintpoints × elnbdypoints)
    ∂O_in_∂τ     :: Vector{Int}
    Γ_in_∂τ      :: Vector{Int}
    u∂O          :: Vector{Float64}

    # Recovery block scratch
    AIoΓ_ie      :: Matrix{Float64}                    # (nelintpoints × lengthΓ)
    AIou∂O_ie    :: Vector{Float64}
    AIoΓg_ie     :: Vector{Float64}
    rhs_ie       :: Vector{Float64}                    # replaces .+ allocation
    uvo_ie       :: Vector{Float64}
    invAIoIo_buf :: Matrix{Float64}                    # (nelintpoints × nelintpoints)

    # Inference buffers
    infer        :: EL_InferBuffers

    # ONNX session (loaded once)
    sess
    input_name   :: String
    output_name  :: String
end

"""
        EL_WorkBuffers(mesh, A, A_∂τ∂τ, nfeatures, nelintpoints, elnbdypoints)

    Allocate every working array for `elementLearning_Axb!` once at setup.
    """
function EL_WorkBuffers(mesh, A::SparseMatrixCSC, A_∂τ∂τ::SparseMatrixCSC,
                        nfeatures::Int, nelintpoints::Int, elnbdypoints::Int,
                        onnx_path::String)
    T = eltype(A)

    if isfile(onnx_path)
        
        # Must be set before the ORT session is created
        ENV["OMP_NUM_THREADS"]        = "1"
        ENV["ORT_NUM_THREADS"]        = "1"
        ENV["OPENBLAS_NUM_THREADS"]   = "1"   # avoids contention with Julia's BLAS
        #ENV["OMP_NUM_THREADS"] = string(Sys.CPU_THREADS)
        
        sess        = ONNXRunTime.load_inference(onnx_path)
        input_name  = first(sess.input_names)
        output_name = first(sess.output_names)
        
    else
        sess = nothing
        input_name  = "nothing"
        output_name = "nothing"
    end
   
    return EL_WorkBuffers(
        zeros(Int, elnbdypoints),                                   # conn_∂O_idx
        zeros(Int, elnbdypoints),                                   # conn_∂τ_idx
        spzeros(T, mesh.length∂O, mesh.length∂τ),                  # ΔB
        Matrix{T}(undef, nelintpoints, nelintpoints),               # invAvovo_buf
        Matrix{T}(undef, nelintpoints, elnbdypoints),               # BC_local
        Vector{Int}(undef, mesh.length∂O),                         # ∂O_in_∂τ
        Vector{Int}(undef, mesh.lengthΓ),                          # Γ_in_∂τ
        Vector{T}(undef, mesh.length∂O),                           # u∂O
        Matrix{T}(undef, nelintpoints, mesh.lengthΓ),              # AIoΓ_ie
        Vector{T}(undef, nelintpoints),                            # AIou∂O_ie
        Vector{T}(undef, nelintpoints),                            # AIoΓg_ie
        Vector{T}(undef, nelintpoints),                            # rhs_ie
        Vector{T}(undef, nelintpoints),                            # uvo_ie
        Matrix{T}(undef, nelintpoints, nelintpoints),               # invAIoIo_buf
        EL_InferBuffers(mesh, A_∂τ∂τ, nfeatures,
                        nelintpoints, elnbdypoints),                # infer
        sess, input_name, output_name,
    )
end

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  Main function                                                               ║
# ╚══════════════════════════════════════════════════════════════════════════════╝
function elementLearning_Axb!(u, uaux, mesh::St_mesh,
                              A::SparseMatrixCSC,
                              ubdy, EL,
                              avisc,
                              bufferin, bufferout,
                              BOΓg, gΓ,
                              wbuf::EL_WorkBuffers;       # ← replaces all internal allocs
                              isamp=1,
                              total_cols_writtenin=0,
                              total_cols_writtenout=0)

    mesh.lengthO  = mesh.length∂O + mesh.lengthIo
    nelintpoints  = (mesh.ngl - 2)^2
    nelpoints     = size(mesh.conn, 2)
    elnbdypoints  = nelpoints - nelintpoints

    # ── DOF → position lookup tables ─────────────────────────────────────────
    # These are O(length∂τ) Dicts.  If this function is called many times,
    # promote ∂O_pos and ∂τ_pos into EL_WorkBuffers and build them once.
    ∂O_pos = Dict{Int,Int}(mesh.∂O[i] => i for i in 1:mesh.length∂O)
    ∂τ_pos = Dict{Int,Int}(mesh.∂τ[j] => j for j in 1:mesh.length∂τ)

    # =========================================================================
    # SECTION 1: Sparse skeleton submatrices
    # =========================================================================
    A_∂τ∂τ = A[mesh.∂τ, mesh.∂τ]
    A_∂O∂τ = A[mesh.∂O, mesh.∂τ]

    # =========================================================================
    # SECTION 2: Fill per-element 3D blocks from A (no changes needed here)
    # =========================================================================
    nelem = mesh.nelem
    @inbounds for iel = 1:nelem
        for j = 1:elnbdypoints
            gnode                  = mesh.conn[iel, j]
            wbuf.conn_∂O_idx[j]   = get(∂O_pos, gnode, 0)
            wbuf.conn_∂τ_idx[j]   = get(∂τ_pos, gnode, 0)
        end

        ii = 1
        for i = elnbdypoints+1:nelpoints
            ipo = mesh.conn[iel, i]
            for j = 1:elnbdypoints
                gj  = mesh.conn[iel, j]
                val = A[ipo, gj]
                EL.Avovb[ii, j, iel] = val
                EL.Avo∂τ[ii, j, iel] = val
                EL.AIo∂τ[ii, j, iel] = val
                if wbuf.conn_∂O_idx[j] != 0
                    EL.Avo∂O[ii, j, iel] = val
                    EL.AIo∂O[ii, j, iel] = val
                    val3 = A[gj, ipo]
                    EL.A∂Ovo[j, ii, iel] = val3
                    EL.A∂OIo[j, ii, iel] = val3
                end
            end
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
    # SECTION 3: Gather Dirichlet data
    # =========================================================================
    lengthΓ = mesh.lengthΓ
    @inbounds for iΓ = 1:lengthΓ
        gΓ[iΓ] = ubdy[mesh.Γ[iΓ], 1]
    end

    if EL.lEL_Sample

        # ── Build ΔB: reset sparsity pattern in-place then accumulate ─────────
        # spzeros every call was the hidden allocator here — reuse wbuf.ΔB
        fill!(nonzeros(wbuf.ΔB), zero(eltype(A)))

        @inbounds for iel = 1:nelem
            for j = 1:elnbdypoints
                gnode               = mesh.conn[iel, j]
                wbuf.conn_∂τ_idx[j] = get(∂τ_pos, gnode, 0)
                wbuf.conn_∂O_idx[j] = get(∂O_pos, gnode, 0)
            end

            # inv() into pre-allocated buffer — no allocation
            # Use copyto! + lu factorization to avoid inv() heap alloc
            copyto!(wbuf.invAvovo_buf, @view(EL.Avovo[:, :, iel]))
            invAvovo = inv(wbuf.invAvovo_buf)   # in-place LU, result in wbuf.invAvovo_buf

            # BC_local = invAvovo * Avo∂τ — write into pre-allocated buffer
            LinearAlgebra.mul!(wbuf.BC_local, invAvovo,
                               @view(EL.Avo∂τ[:, :, iel]))

            for j_loc = 1:elnbdypoints
                jτ = wbuf.conn_∂τ_idx[j_loc];  jτ == 0 && continue
                for i_loc = 1:elnbdypoints
                    io = wbuf.conn_∂O_idx[i_loc];  io == 0 && continue
                    s  = zero(eltype(A))
                    for ii = 1:nelintpoints
                        s += EL.A∂Ovo[i_loc, ii, iel] * wbuf.BC_local[ii, j_loc]
                    end
                    wbuf.ΔB[io, jτ] += s
                end
            end
        end

        B_∂O∂τ = A_∂O∂τ - wbuf.ΔB   # one alloc per call — unavoidable

        # ── Index maps in-place ───────────────────────────────────────────────
        @inbounds for i  = 1:mesh.length∂O
            wbuf.∂O_in_∂τ[i]  = ∂τ_pos[mesh.∂O[i]]
        end
        @inbounds for iΓ = 1:mesh.lengthΓ
            wbuf.Γ_in_∂τ[iΓ]  = ∂τ_pos[mesh.Γ[iΓ]]
        end

        B_∂O∂O = B_∂O∂τ[:, wbuf.∂O_in_∂τ]
        B_∂O∂Γ = B_∂O∂τ[:, wbuf.Γ_in_∂τ]

        BOΓg_tmp          = B_∂O∂Γ * gΓ
        wbuf.u∂O         .= -(B_∂O∂O \ BOΓg_tmp)  # eq. (1.4)

        @inbounds for io = 1:mesh.length∂O;  u[mesh.∂O[io]] = wbuf.u∂O[io];  end
        @inbounds for io = 1:mesh.lengthΓ;   u[mesh.Γ[io]]  = gΓ[io];        end

        # ── Per-element interior recovery ─────────────────────────────────────
        @inbounds for iel = 1:mesh.nelem
            for j = 1:elnbdypoints
                wbuf.conn_∂O_idx[j] = get(∂O_pos, mesh.conn[iel, j], 0)
            end

            # Fill AIoΓ_ie in-place — no zeros() allocation
            for iΓ = 1:mesh.lengthΓ
                g1 = mesh.Γ[iΓ]
                for ii = 1:nelintpoints
                    wbuf.AIoΓ_ie[ii, iΓ] = A[mesh.conn[iel, elnbdypoints+ii], g1]
                end
            end

            fill!(wbuf.AIou∂O_ie, zero(eltype(A)))
            for j_loc = 1:elnbdypoints
                io = wbuf.conn_∂O_idx[j_loc];  io == 0 && continue
                for ii = 1:nelintpoints
                    wbuf.AIou∂O_ie[ii] += EL.AIo∂O[ii, j_loc, iel] * wbuf.u∂O[io]
                end
            end

            LinearAlgebra.mul!(wbuf.AIoΓg_ie, wbuf.AIoΓ_ie, gΓ)

            # rhs_ie = AIou∂O_ie + AIoΓg_ie — in-place, no .+ allocation
            wbuf.rhs_ie .= wbuf.AIou∂O_ie .+ wbuf.AIoΓg_ie

            # -inv(AIoIo) * rhs_ie — inv into pre-allocated buffer
            copyto!(wbuf.invAIoIo_buf, @view(EL.AIoIo[:, :, iel]))
            invAIoIo = inv(wbuf.invAIoIo_buf)
            LinearAlgebra.mul!(wbuf.uvo_ie, invAIoIo, wbuf.rhs_ie, -1.0, 0.0)

            for ii = 1:nelintpoints
                u[mesh.conn[iel, elnbdypoints+ii]] = wbuf.uvo_ie[ii]
            end
        end

        # ── ML tensor recording ───────────────────────────────────────────────
        EL.input_tensor[:, isamp] .= vec(avisc)
        let iel = 1
            # inv into pre-allocated buffer, mul! with @view — no temp allocs
            copyto!(wbuf.invAvovo_buf, @view(EL.Avovo[:, :, iel]))
            invAvovo = inv(wbuf.invAvovo_buf)
            LinearAlgebra.mul!(EL.Tie,  invAvovo, @view(EL.Avovb[:, :, iel]), -1.0, 0.0)
            LinearAlgebra.mul!(EL.T1,   transpose(@view(EL.Avovb[:, :, iel])), EL.Tie, -1.0, 0.0)
            EL.output_tensor[:, isamp] .= -vec(EL.Tie)
        end
        write_MLtensor!(bufferin,  EL.input_tensor[:,  isamp])
        write_MLtensor!(bufferout, EL.output_tensor[:, isamp])
        
    else
        # ── Inference (sess already loaded in wbuf — no load_inference here) ──
        println(YELLOW_FG(string(" # --- INFERENCE — solution stored in u ..........")))
        elementLearning_infer!(u, mesh,
                               wbuf.sess, wbuf.input_name, wbuf.output_name,
                               avisc, EL, A_∂τ∂τ, ∂τ_pos, gΓ, wbuf.infer,
                               nelintpoints, elnbdypoints)

        @btime elementLearning_infer!($u, $mesh,
                                      $wbuf.sess, $wbuf.input_name, $wbuf.output_name,
                                      $avisc, $EL, $A_∂τ∂τ, $∂τ_pos, $gΓ, $wbuf.infer,
                                      $nelintpoints, $elnbdypoints)
        println(YELLOW_FG(string(" # --- INFERENCE — solution stored in u .......... DONE")))
    end

    return nothing
end

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  Main inference function                                                    ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

"""
        elementLearning_infer!(u, mesh, sess, input_name, output_name,
                               avisc, EL, A_∂τ∂τ, ∂τ_pos, gΓ, buf,
                               nelintpoints, elnbdypoints)

    Element-learning inference step. Computes the NN-based transformation matrices
    T^{ie} for all elements, assembles the skeleton system, solves for interior and
    boundary DOFs, and writes the result into `u`.

    # Arguments
    - `u`             : global solution vector/matrix — written in place
    - `mesh`          : mesh struct (.nelem, .conn, .∂O, .Γ, .length∂O, .lengthΓ)
    - `sess`          : ONNXRunTime.InferenceSession
    - `input_name`    : ONNX input tensor name  (String)
    - `output_name`   : ONNX output tensor name (String)
    - `avisc`         : NN input features — (1, nfeatures) shared OR (nelem, nfeatures)
    - `EL`            : element-learning struct (.Avovb[:,:,iel])
    - `A_∂τ∂τ`        : sparse skeleton submatrix of A (read-only reference)
    - `∂τ_pos`        : Dict mapping global node → index in ∂τ numbering
    - `gΓ`            : Dirichlet values on Γ  (Vector{Float64})
    - `buf`           : EL_InferBuffers — all pre-allocated working arrays
    - `nelintpoints`  : number of interior quadrature points per element
    - `elnbdypoints`  : number of boundary nodes per element
    """
function elementLearning_infer!(
    u            :: AbstractMatrix{Float64},
    mesh,
    sess,
    input_name   :: String,
    output_name  :: String,
    avisc        :: Matrix{Float64},
    EL,
    A_∂τ∂τ       :: SparseMatrixCSC,
    ∂τ_pos       :: Dict,
    gΓ           :: Vector{Float64},
    buf          :: EL_InferBuffers,
    nelintpoints :: Int,
    elnbdypoints :: Int,
    )
    nfeatures    = size(avisc, 2)
    nout         = nelintpoints * elnbdypoints
    nelem_avisc  = size(avisc, 1)   # 1 → shared input; mesh.nelem → per-element

    # ── Reset B_∂τ∂τ to the skeleton matrix (in-place, no allocation) ─────────
    copyto!(buf.B_∂τ∂τ, A_∂τ∂τ)

    # ══════════════════════════════════════════════════════════════════════════
    # INFERENCE
    # ══════════════════════════════════════════════════════════════════════════

    if nelem_avisc == 1
        # ── CASE A: single shared avisc row — one inference call, broadcast ───
        # Float64 → Float32 in-place, no allocation
        @inbounds for k = 1:nfeatures
            buf.avisc_f32[1, k] = Float32(avisc[1, k])
        end

        # One ONNX call: (1, nfeatures) → (1, nout)
        y_single = sess(Dict(input_name => buf.avisc_f32))
        ŷ_single = y_single[output_name]

        # Cast output to Float64 in-place
        @inbounds for k = 1:nout
            buf.ŷ_f64_buf[k] = Float64(ŷ_single[1, k])
        end

        # Reshape is a free view — broadcast same Tie_nn to every element
        Tie_nn = reshape(buf.ŷ_f64_buf, nelintpoints, elnbdypoints)
        @inbounds for iel = 1:mesh.nelem
            buf.Tie_nn_all[:, :, iel] .= Tie_nn
        end

    else
        # ── CASE B: per-element avisc — one batched inference call ─────────────
        # Pack avisc → Float32 batch in-place; avisc is (nelem, nfeatures)
        @inbounds for iel = 1:mesh.nelem
            for k = 1:nfeatures
                buf.avisc_batch[iel, k] = Float32(avisc[iel, k])
            end
        end

        # One ONNX call: (nelem, nfeatures) → (nelem, nout)
        y_batch = sess(Dict(input_name => buf.avisc_batch))
        ŷ_batch = y_batch[output_name]

        # Unpack into Tie_nn_all in-place; ŷ_f64_buf used as row staging buffer
        @inbounds for iel = 1:mesh.nelem
            for k = 1:nout
                buf.ŷ_f64_buf[k] = Float64(ŷ_batch[iel, k])
            end
            # reshape of a contiguous column view — allocation-free
            buf.Tie_nn_all[:, :, iel] .= reshape(buf.ŷ_f64_buf,
                                                 nelintpoints, elnbdypoints)
        end
    end

    # ══════════════════════════════════════════════════════════════════════════
    # STEP 4 — Assemble B_∂τ∂τ
    # transpose(@view(...)) is a lazy no-copy wrapper
    # ══════════════════════════════════════════════════════════════════════════
    @inbounds for iel = 1:mesh.nelem
        for j = 1:elnbdypoints
            buf.conn_∂τ_idx[j] = get(∂τ_pos, mesh.conn[iel, j], 0)
        end
        LinearAlgebra.mul!(buf.M,
                           transpose(@view(EL.Avovb[:, :, iel])),
                           @view(buf.Tie_nn_all[:, :, iel]))
        for i = 1:elnbdypoints
            i_prime = buf.conn_∂τ_idx[i];  i_prime == 0 && continue
            for j = 1:elnbdypoints
                j_prime = buf.conn_∂τ_idx[j];  j_prime == 0 && continue
                buf.B_∂τ∂τ[i_prime, j_prime] -= buf.M[i, j]
            end
        end
    end

    # ══════════════════════════════════════════════════════════════════════════
    # STEP 5 — Build ∂O and Γ index maps in-place
    # ══════════════════════════════════════════════════════════════════════════
    @inbounds for i  = 1:mesh.length∂O
        buf.∂O_in_∂τ[i]  = ∂τ_pos[mesh.∂O[i]]
    end
    @inbounds for iΓ = 1:mesh.lengthΓ
        buf.Γ_in_∂τ[iΓ]  = ∂τ_pos[mesh.Γ[iΓ]]
    end

    # ══════════════════════════════════════════════════════════════════════════
    # STEP 6 — Extract sparse submatrices and solve
    # Sparse indexing allocates O(nnz) — unavoidable, happens once per call
    # NEVER call Matrix() on B_∂O∂O — that densifies to potentially GBs
    # ══════════════════════════════════════════════════════════════════════════
    B_∂O∂O  = buf.B_∂τ∂τ[buf.∂O_in_∂τ, buf.∂O_in_∂τ]   # SparseMatrixCSC
    B_∂O∂Γ  = buf.B_∂τ∂τ[buf.∂O_in_∂τ, buf.Γ_in_∂τ]    # SparseMatrixCSC

    BOΓg_nn = B_∂O∂Γ * gΓ             # sparse × dense → dense  (one alloc)
    u∂O_nn  = -(B_∂O∂O \ BOΓg_nn)     # UMFPACK sparse direct solve (one alloc) eq. (1.4)

    # ══════════════════════════════════════════════════════════════════════════
    # STEP 7 — Scatter solution into u
    # ══════════════════════════════════════════════════════════════════════════
    @inbounds for io = 1:mesh.length∂O
        u[mesh.∂O[io]] = u∂O_nn[io]
    end
    @inbounds for io = 1:mesh.lengthΓ
        u[mesh.Γ[io]]  = gΓ[io]
    end

    # ══════════════════════════════════════════════════════════════════════════
    # STEP 8a — Gather boundary solution for each element
    # ══════════════════════════════════════════════════════════════════════════
    @inbounds for iel = 1:mesh.nelem
        for ibdy = 1:elnbdypoints
            buf.uvb_nn[iel, ibdy] = u[mesh.conn[iel, ibdy]]
        end
    end

    # ══════════════════════════════════════════════════════════════════════════
    # STEP 8b — Local interior recovery  u_vo = -T^{ie,nn} * u_vb
    # mul!(C, A, B, α, β) → C = α*A*B + β*C  (folds negation, no extra alloc)
    # @view(uvb_nn[iel,:]) is a contiguous row — no copy
    # ══════════════════════════════════════════════════════════════════════════
    @inbounds for iel = 1:mesh.nelem
        LinearAlgebra.mul!(buf.uvo_nn,
                           @view(buf.Tie_nn_all[:, :, iel]),
                           @view(buf.uvb_nn[iel, :]),
                           -1.0, 0.0)
        for i = 1:nelintpoints
            u[mesh.conn[iel, elnbdypoints + i]] = buf.uvo_nn[i]
        end
    end

    return nothing
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
