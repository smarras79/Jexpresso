#include("../../EL_Jexpresso/NN_RFRC.jl")
using SparseArrays
using JLD2
using ONNXRunTime
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
    Avovb   = KernelAbstractions.zeros(backend, T, dims_vovb)
    Avo∂O   = KernelAbstractions.zeros(backend, T, dims_vovb)
    Avo∂τ   = KernelAbstractions.zeros(backend, T, dims_vovb)
    AIo∂τ   = KernelAbstractions.zeros(backend, T, dims_vovb)
    AIo∂O   = KernelAbstractions.zeros(backend, T, dims_vovb)

    # ── Per-element: local-boundary × interior  (elnbdy × nvo × nelem) ───────
    A∂Ovo   = KernelAbstractions.zeros(backend, T, dims_∂Ovo)
    A∂OIo   = KernelAbstractions.zeros(backend, T, dims_∂Ovo)

    # ── Local temporaries ─────────────────────────────────────────────────────
    T1      = KernelAbstractions.zeros(backend, T, dims_T1)
    T2      = KernelAbstractions.zeros(backend, T, dims_T2)
    Tie     = KernelAbstractions.zeros(backend, T, dims_T2)

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
#  EL_InferBuffers — pre-allocated working arrays for inference
# =============================================================================
struct EL_InferBuffers
    # ── ONNX staging (row-major: samples × features) ─────────────────────────
    avisc_f32    :: Matrix{Float32}    # (1,          nfeatures)  — Case A shared
    avisc_batch  :: Matrix{Float32}    # (nelem,      nfeatures)  — Case B per-element

    # ── JLD2/RFRC staging (column-major: features × samples) ─────────────────
    avisc_f32_T  :: Matrix{Float32}    # (nfeatures,  nelem)      — transposed for NNRFRC
    ŷ_f32_batch  :: Matrix{Float32}    # (nout,       nelem)      — JLD2 output buffer

    # ── Shared staging ────────────────────────────────────────────────────────
    ŷ_f64_buf    :: Vector{Float64}    # (nout,)                  — per-element cast buffer

    # ── Element assembly ──────────────────────────────────────────────────────
    Tie_nn_all   :: Array{Float64, 3}  # (nelintpoints, elnbdypoints, nelem)
    conn_∂τ_idx  :: Vector{Int}        # (elnbdypoints,)
    M            :: Matrix{Float64}    # (elnbdypoints, elnbdypoints)
    B_∂τ∂τ       :: SparseMatrixCSC{Float64, Int32}

    # ── Index maps ────────────────────────────────────────────────────────────
    ∂O_in_∂τ     :: Vector{Int}
    Γ_in_∂τ      :: Vector{Int}

    # ── Gather / recovery ─────────────────────────────────────────────────────
    uvb_nn       :: Matrix{Float64}    # (nelem, elnbdypoints)
    uvo_nn       :: Vector{Float64}    # (nelintpoints,)
end

"""
    EL_InferBuffers(mesh, A_∂τ∂τ, nfeatures, nelintpoints, elnbdypoints)

Allocate all working arrays for `elementLearning_infer!` once.
"""
function EL_InferBuffers(mesh, A_∂τ∂τ::SparseMatrixCSC,
                         nfeatures::Int, nelintpoints::Int, elnbdypoints::Int)
    nout = nelintpoints * elnbdypoints
    return EL_InferBuffers(
        # ONNX staging
        Matrix{Float32}(undef, 1,          nfeatures),
        Matrix{Float32}(undef, mesh.nelem, nfeatures),
        # JLD2/RFRC staging
        Matrix{Float32}(undef, nfeatures,  mesh.nelem),
        Matrix{Float32}(undef, nout,       mesh.nelem),
        # Shared staging
        Vector{Float64}(undef, nout),
        # Element assembly
        Array{Float64, 3}(undef, nelintpoints, elnbdypoints, mesh.nelem),
        Vector{Int}(undef, elnbdypoints),
        Matrix{Float64}(undef, elnbdypoints, elnbdypoints),
        copy(A_∂τ∂τ),
        # Index maps
        Vector{Int}(undef, mesh.length∂O),
        Vector{Int}(undef, mesh.lengthΓ),
        # Gather / recovery
        Matrix{Float64}(undef, mesh.nelem, elnbdypoints),
        Vector{Float64}(undef, nelintpoints),
    )
end


# =============================================================================
#  EL_WorkBuffers — all scratch arrays + model session
# =============================================================================
struct EL_WorkBuffers
    # Section 2 scratch
    conn_∂O_idx  :: Vector{Int}
    conn_∂τ_idx  :: Vector{Int}

    # Sampling block scratch
    ΔB           :: SparseMatrixCSC{Float64, Int32}
    invAvovo_buf :: Matrix{Float64}
    BC_local     :: Matrix{Float64}
    ∂O_in_∂τ     :: Vector{Int}
    Γ_in_∂τ      :: Vector{Int}
    u∂O          :: Vector{Float64}

    # Recovery block scratch
    AIoΓ_ie      :: Matrix{Float64}
    AIou∂O_ie    :: Vector{Float64}
    AIoΓg_ie     :: Vector{Float64}
    rhs_ie       :: Vector{Float64}
    uvo_ie       :: Vector{Float64}
    invAIoIo_buf :: Matrix{Float64}

    # Inference buffers
    infer        :: EL_InferBuffers

    # ── Model (ONNX or JLD2) ─────────────────────────────────────────────────
    model                              # ONNXRunTime.InferenceSession OR NNRFRC.RFRC
    model_type   :: Symbol             # :ONNX or :JLD2
    input_name   :: String             # ONNX input tensor name   (empty for JLD2)
    output_name  :: String             # ONNX output tensor name  (empty for JLD2)
end

"""
    EL_WorkBuffers(mesh, A, A_∂τ∂τ, nfeatures, nelintpoints, elnbdypoints, NNfile)

Allocate every working array for `elementLearning_Axb!` once at setup.
Loads the model from `NNfile` — supports both `.onnx` and `.jld2` extensions.
"""
function EL_WorkBuffers(mesh, A::SparseMatrixCSC, A_∂τ∂τ::SparseMatrixCSC,
                        nfeatures::Int, nelintpoints::Int, elnbdypoints::Int,
                        NNfile::String)
    T = eltype(A)

    local model
    local m_type    = :NONE
    local m_inname  = ""
    local m_outname = ""

    if !isfile(NNfile)
        error("Model file not found: $NNfile")
    end

    ext = lowercase(splitext(NNfile)[2])

    if ext == ".onnx"
        # ── ONNX: load session, cache tensor names ────────────────────────────
        ENV["OMP_NUM_THREADS"]      = "1"
        ENV["ORT_NUM_THREADS"]      = "1"
        ENV["OPENBLAS_NUM_THREADS"] = "1"

        model      = ONNXRunTime.load_inference(NNfile)
        m_type     = :ONNX
        m_inname   = first(model.input_names)
        m_outname  = first(model.output_names)

    elseif ext == ".jld2"
        # ── JLD2: load and reconstruct into the local NNRFRC.RFRC type ────────
        raw = jldopen(NNfile, "r") do file; file["model"]; end
        if raw isa NNRFRC.RFRC
            model = raw
        else
            # JLD2 returns ReconstructedMutable when module path differs;
            # rebuild field-by-field matching the RFRC struct definition:
            #   struct RFRC{T}
            #       W_in, W_res, W_out, b_out, ridge_alpha, activation
            #   end
            model = NNRFRC.RFRC{Float32}(
                Matrix{Float32}(raw.W_in),       # 1st: W_in
                Matrix{Float32}(raw.W_res),      # 2nd: W_res
                Matrix{Float32}(raw.W_out),      # 3rd: W_out
                Vector{Float32}(raw.b_out),      # 4th: b_out
                Float32(raw.ridge_alpha),         # 5th: ridge_alpha
                raw.activation                    # 6th: activation (Function)
            )
        end
        m_type = :JLD2

    else
        error("Unsupported model extension: \"$ext\" — use .onnx or .jld2")
    end

    return EL_WorkBuffers(
        zeros(Int, elnbdypoints),                                   # conn_∂O_idx
        zeros(Int, elnbdypoints),                                   # conn_∂τ_idx
        spzeros(T, mesh.length∂O, mesh.length∂τ),                   # ΔB
        Matrix{T}(undef, nelintpoints, nelintpoints),               # invAvovo_buf
        Matrix{T}(undef, nelintpoints, elnbdypoints),               # BC_local
        Vector{Int}(undef, mesh.length∂O),                          # ∂O_in_∂τ
        Vector{Int}(undef, mesh.lengthΓ),                           # Γ_in_∂τ
        Vector{T}(undef, mesh.length∂O),                            # u∂O
        Matrix{T}(undef, nelintpoints, mesh.lengthΓ),               # AIoΓ_ie
        Vector{T}(undef, nelintpoints),                             # AIou∂O_ie
        Vector{T}(undef, nelintpoints),                             # AIoΓg_ie
        Vector{T}(undef, nelintpoints),                             # rhs_ie
        Vector{T}(undef, nelintpoints),                             # uvo_ie
        Matrix{T}(undef, nelintpoints, nelintpoints),               # invAIoIo_buf
        EL_InferBuffers(mesh, A_∂τ∂τ, nfeatures,
                        nelintpoints, elnbdypoints),                # infer
        model, m_type, m_inname, m_outname,
    )
end


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  Main function                                                              ║
# ╚══════════════════════════════════════════════════════════════════════════════╝
function elementLearning_Axb!(u, uaux, mesh::St_mesh,
                              A::SparseMatrixCSC,
                              ubdy, EL,
                              avisc,
                              bufferin, bufferout,
                              BOΓg, gΓ,
                              wbuf::EL_WorkBuffers;
                              isamp=1,
                              total_cols_writtenin=0,
                              total_cols_writtenout=0)

    mesh.lengthO  = mesh.length∂O + mesh.lengthIo
    nelintpoints  = (mesh.ngl - 2)^2
    nelpoints     = size(mesh.conn, 2)
    elnbdypoints  = nelpoints - nelintpoints

    # ── DOF → position lookup tables ─────────────────────────────────────────
    ∂O_pos = Dict{Int,Int}(mesh.∂O[i] => i for i in 1:mesh.length∂O)
    ∂τ_pos = Dict{Int,Int}(mesh.∂τ[j] => j for j in 1:mesh.length∂τ)

    # =========================================================================
    # SECTION 1: Sparse skeleton submatrices
    # =========================================================================
    A_∂τ∂τ = A[mesh.∂τ, mesh.∂τ]
    A_∂O∂τ = A[mesh.∂O, mesh.∂τ]

    # =========================================================================
    # SECTION 2: Fill per-element 3D blocks from A
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

        # ── Build ΔB ──────────────────────────────────────────────────────────
        fill!(nonzeros(wbuf.ΔB), zero(eltype(A)))

        @inbounds for iel = 1:nelem
            for j = 1:elnbdypoints
                gnode               = mesh.conn[iel, j]
                wbuf.conn_∂τ_idx[j] = get(∂τ_pos, gnode, 0)
                wbuf.conn_∂O_idx[j] = get(∂O_pos, gnode, 0)
            end

            copyto!(wbuf.invAvovo_buf, @view(EL.Avovo[:, :, iel]))
            invAvovo = inv(wbuf.invAvovo_buf)

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

        B_∂O∂τ = A_∂O∂τ - wbuf.ΔB

        @inbounds for i  = 1:mesh.length∂O
            wbuf.∂O_in_∂τ[i]  = ∂τ_pos[mesh.∂O[i]]
        end
        @inbounds for iΓ = 1:mesh.lengthΓ
            wbuf.Γ_in_∂τ[iΓ]  = ∂τ_pos[mesh.Γ[iΓ]]
        end

        B_∂O∂O = B_∂O∂τ[:, wbuf.∂O_in_∂τ]
        B_∂O∂Γ = B_∂O∂τ[:, wbuf.Γ_in_∂τ]

        BOΓg_tmp          = B_∂O∂Γ * gΓ
        wbuf.u∂O         .= -(B_∂O∂O \ BOΓg_tmp)

        @inbounds for io = 1:mesh.length∂O;  u[mesh.∂O[io]] = wbuf.u∂O[io];  end
        @inbounds for io = 1:mesh.lengthΓ;   u[mesh.Γ[io]]  = gΓ[io];        end

        # ── Per-element interior recovery ─────────────────────────────────────
        @inbounds for iel = 1:mesh.nelem
            for j = 1:elnbdypoints
                wbuf.conn_∂O_idx[j] = get(∂O_pos, mesh.conn[iel, j], 0)
            end

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
            wbuf.rhs_ie .= wbuf.AIou∂O_ie .+ wbuf.AIoΓg_ie

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
            copyto!(wbuf.invAvovo_buf, @view(EL.Avovo[:, :, iel]))
            invAvovo = inv(wbuf.invAvovo_buf)
            LinearAlgebra.mul!(EL.Tie,  invAvovo, @view(EL.Avovb[:, :, iel]), -1.0, 0.0)
            LinearAlgebra.mul!(EL.T1,   transpose(@view(EL.Avovb[:, :, iel])), EL.Tie, -1.0, 0.0)
            EL.output_tensor[:, isamp] .= -vec(EL.Tie)
        end
        write_MLtensor!(bufferin,  EL.input_tensor[:,  isamp])
        write_MLtensor!(bufferout, EL.output_tensor[:, isamp])
        
    else
        # ── Inference ─────────────────────────────────────────────────────────
        println(YELLOW_FG(string(" # --- INFERENCE — solution stored in u ..........")))
        elementLearning_infer!(u, mesh,
                               wbuf.model, wbuf.model_type,
                               wbuf.input_name, wbuf.output_name,
                               avisc, EL, A_∂τ∂τ, ∂τ_pos, gΓ, wbuf.infer,
                               nelintpoints, elnbdypoints)

        # PERF: ad-hoc @btime instrumentation removed. The line below was a
        # *duplicate* call of elementLearning_infer! for timing during
        # development; it forced `using BenchmarkTools` at the top of
        # src/Jexpresso.jl just so the @btime macro could be expanded at
        # module load. BenchmarkTools adds ~50 MB to the per-rank load
        # baseline. The real inference call above is unchanged.
        # @btime elementLearning_infer!($u, $mesh,
        #                               $wbuf.model, $wbuf.model_type,
        #                               $wbuf.input_name, $wbuf.output_name,
        #                               $avisc, $EL, $A_∂τ∂τ, $∂τ_pos, $gΓ, $wbuf.infer,
        #                               $nelintpoints, $elnbdypoints)
        println(YELLOW_FG(string(" # --- INFERENCE — solution stored in u .......... DONE")))
    end

    return nothing
end


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  Inference engine — supports both ONNX and JLD2/RFRC models                ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

"""
    elementLearning_infer!(u, mesh, model, model_type, input_name, output_name,
                           avisc, EL, A_∂τ∂τ, ∂τ_pos, gΓ, buf,
                           nelintpoints, elnbdypoints)

Element-learning inference step.  Dispatches on `model_type`:
  - `:ONNX` — calls the ONNXRunTime InferenceSession
  - `:JLD2`  — calls NNRFRC native Julia model

# Arguments
- `u`             : global solution vector/matrix — written in place
- `mesh`          : mesh struct
- `model`         : ONNXRunTime.InferenceSession  OR  NNRFRC.RFRC model
- `model_type`    : `:ONNX` or `:JLD2`
- `input_name`    : ONNX input tensor name   (unused for JLD2)
- `output_name`   : ONNX output tensor name  (unused for JLD2)
- `avisc`         : NN input features — (1, nfeatures) shared OR (nelem, nfeatures)
- `EL`            : element-learning struct (.Avovb[:,:,iel])
- `A_∂τ∂τ`        : sparse skeleton submatrix of A (read-only)
- `∂τ_pos`        : Dict mapping global node → index in ∂τ numbering
- `gΓ`            : Dirichlet values on Γ
- `buf`           : EL_InferBuffers — all pre-allocated working arrays
- `nelintpoints`  : number of interior points per element
- `elnbdypoints`  : number of boundary nodes per element
"""
function elementLearning_infer!(
    u            :: AbstractMatrix{Float64},
    mesh,
    model,
    model_type   :: Symbol,
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
    nelem        = mesh.nelem
    nelem_avisc  = size(avisc, 1)   # 1 → shared input; nelem → per-element

    # ── Reset B_∂τ∂τ to the skeleton matrix ───────────────────────────────────
    copyto!(buf.B_∂τ∂τ, A_∂τ∂τ)

    # ══════════════════════════════════════════════════════════════════════════
    # INFERENCE — dispatch on model_type
    # ══════════════════════════════════════════════════════════════════════════

    if model_type == :ONNX
        _infer_onnx!(buf, model, input_name, output_name,
                     avisc, nfeatures, nout, nelem, nelem_avisc,
                     nelintpoints, elnbdypoints)

    elseif model_type == :JLD2
        _infer_jld2!(buf, model,
                     avisc, nfeatures, nout, nelem, nelem_avisc,
                     nelintpoints, elnbdypoints)
    else
        error("Unknown model_type: $model_type — expected :ONNX or :JLD2")
    end

    # ══════════════════════════════════════════════════════════════════════════
    # STEP 4 — Assemble B_∂τ∂τ (Schur complement)
    # ══════════════════════════════════════════════════════════════════════════
    @inbounds for iel = 1:nelem
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
    # ══════════════════════════════════════════════════════════════════════════
    B_∂O∂O  = buf.B_∂τ∂τ[buf.∂O_in_∂τ, buf.∂O_in_∂τ]
    B_∂O∂Γ  = buf.B_∂τ∂τ[buf.∂O_in_∂τ, buf.Γ_in_∂τ]

    BOΓg_nn = B_∂O∂Γ * gΓ
    u∂O_nn  = -(B_∂O∂O \ BOΓg_nn)

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
    @inbounds for iel = 1:nelem
        for ibdy = 1:elnbdypoints
            buf.uvb_nn[iel, ibdy] = u[mesh.conn[iel, ibdy]]
        end
    end

    # ══════════════════════════════════════════════════════════════════════════
    # STEP 8b — Local interior recovery  u_vo = -T^{ie,nn} * u_vb
    # ══════════════════════════════════════════════════════════════════════════
    @inbounds for iel = 1:nelem
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
#  _infer_onnx!  — ONNX inference sub-routine
#
#  ONNXRunTime.jl may return outputs transposed relative to the Python
#  convention, so we auto-detect the layout in every code path.
# =============================================================================
function _infer_onnx!(buf, model, input_name, output_name,
                      avisc, nfeatures, nout, nelem, nelem_avisc,
                      nelintpoints, elnbdypoints)

    if nelem_avisc == 1
        # ── CASE A: single shared avisc row — one call, broadcast result ──────
        @inbounds for k = 1:nfeatures
            buf.avisc_f32[1, k] = Float32(avisc[1, k])
        end

        y_single = model(Dict(input_name => buf.avisc_f32))
        ŷ_single = y_single[output_name]

        # Auto-detect output layout
        if size(ŷ_single, 1) == 1 && size(ŷ_single, 2) == nout
            # Shape (1, nout) — row-per-sample
            @inbounds for k = 1:nout
                buf.ŷ_f64_buf[k] = Float64(ŷ_single[1, k])
            end
        elseif size(ŷ_single, 1) == nout && size(ŷ_single, 2) == 1
            # Shape (nout, 1) — column-per-sample
            @inbounds for k = 1:nout
                buf.ŷ_f64_buf[k] = Float64(ŷ_single[k, 1])
            end
        elseif ndims(ŷ_single) == 1 && length(ŷ_single) == nout
            # Shape (nout,) — flat vector
            @inbounds for k = 1:nout
                buf.ŷ_f64_buf[k] = Float64(ŷ_single[k])
            end
        else
            error("ONNX output shape $(size(ŷ_single)) incompatible — " *
                  "expected (1,$nout), ($nout,1), or ($nout,)")
        end

        Tie_nn = reshape(buf.ŷ_f64_buf, nelintpoints, elnbdypoints)
        @inbounds for iel = 1:nelem
            buf.Tie_nn_all[:, :, iel] .= Tie_nn
        end

    else
        # ── CASE B: per-element avisc — one batched call ──────────────────────
        @inbounds for iel = 1:nelem, k = 1:nfeatures
            buf.avisc_batch[iel, k] = Float32(avisc[iel, k])
        end

        y_batch = model(Dict(input_name => buf.avisc_batch))
        ŷ_batch = y_batch[output_name]

        # Auto-detect output layout
        if size(ŷ_batch, 1) == nelem && size(ŷ_batch, 2) == nout
            # Shape (nelem, nout) — row-per-sample
            @inbounds for iel = 1:nelem
                for k = 1:nout
                    buf.ŷ_f64_buf[k] = Float64(ŷ_batch[iel, k])
                end
                buf.Tie_nn_all[:, :, iel] .= reshape(buf.ŷ_f64_buf,
                                                     nelintpoints, elnbdypoints)
            end
        elseif size(ŷ_batch, 1) == nout && size(ŷ_batch, 2) == nelem
            # Shape (nout, nelem) — transposed by ONNXRunTime
            @inbounds for iel = 1:nelem
                for k = 1:nout
                    buf.ŷ_f64_buf[k] = Float64(ŷ_batch[k, iel])
                end
                buf.Tie_nn_all[:, :, iel] .= reshape(buf.ŷ_f64_buf,
                                                     nelintpoints, elnbdypoints)
            end
        else
            error("ONNX output shape $(size(ŷ_batch)) incompatible — " *
                  "expected ($nelem,$nout) or ($nout,$nelem)")
        end
    end

    return nothing
end


# =============================================================================
#  _infer_jld2!  — JLD2/RFRC inference sub-routine
#
#  Forward pass (from NN_RFRC.jl / train_RFRC.jl):
#    H = get_features(model, X)       # H = activation(W_in * X)
#    Y = W_out * H + b_out
#
#  get_features expects X as Matrix{Float32} of shape (nfeatures, nbatch).
#  We pass raw features; get_features handles the W_in projection internally.
# =============================================================================
function _infer_jld2!(buf, model,
                      avisc, nfeatures, nout, nelem, nelem_avisc,
                      nelintpoints, elnbdypoints)

    if nelem_avisc == 1
        # ── CASE A: single shared avisc row — evaluate once, broadcast ────────
        # Transpose (1, nfeatures) → (nfeatures, 1) into buffer
        @inbounds for k = 1:nfeatures
            buf.avisc_f32_T[k, 1] = Float32(avisc[1, k])
        end
        # get_features wants Matrix{T}, not SubArray — use slicing (allocates a small copy)
        avisc_in = buf.avisc_f32_T[:, 1:1]   # (nfeatures, 1) Matrix{Float32}

        # Forward pass: H = activation(W_in * X), then Y = W_out * H + b_out
        H_res = NNRFRC.get_features(model, avisc_in)
        ŷ_col = model.W_out * H_res .+ model.b_out   # (nout, 1)

        @inbounds for k = 1:nout
            buf.ŷ_f64_buf[k] = Float64(ŷ_col[k, 1])
        end

        Tie_nn = reshape(buf.ŷ_f64_buf, nelintpoints, elnbdypoints)
        @inbounds for iel = 1:nelem
            buf.Tie_nn_all[:, :, iel] .= Tie_nn
        end

    else
        # ── CASE B: per-element — batch inference in one shot ─────────────────
        # Transpose avisc (nelem, nfeatures) → (nfeatures, nelem) into buffer
        @inbounds for iel = 1:nelem, k = 1:nfeatures
            buf.avisc_f32_T[k, iel] = Float32(avisc[iel, k])
        end

        # Forward pass: H = activation(W_in * X)
        H_res = NNRFRC.get_features(model, buf.avisc_f32_T)

        # Output: W_out * H + b_out → (nout, nelem)
        LinearAlgebra.mul!(buf.ŷ_f32_batch, model.W_out, H_res)
        buf.ŷ_f32_batch .+= model.b_out

        # Unpack columns into Tie_nn_all
        @inbounds for iel = 1:nelem
            for k = 1:nout
                buf.ŷ_f64_buf[k] = Float64(buf.ŷ_f32_batch[k, iel])
            end
            buf.Tie_nn_all[:, :, iel] .= reshape(buf.ŷ_f64_buf,
                                                 nelintpoints, elnbdypoints)
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

#---------------------------------------------------------------------------------------
# Element learning (standalone Lx = RHS path).
#
# Consolidates the full element-learning workflow — allocation, the sampling
# branch (writes input/output tensors) and the inference branch (calls the
# trained NN and writes the solution) — that used to live inline inside the
# `if inputs[:lelementLearning]` block of driver(). Called from there.
#---------------------------------------------------------------------------------------
function element_learning_linsolve!(sem, params, qp, inputs, OUTPUT_DIR, TFloat, rank)

    if rank == 0 println(BLUE_FG(string(" # ALLOCATE FOR ELEMENT LEARNING ......."))) end

    ngl            = sem.mesh.ngl
    ngr            = sem.mesh.ngr
    npoin          = sem.mesh.npoin
    nelem          = sem.mesh.nelem
    nelem_semi_inf = sem.mesh.nelem_semi_inf
    nelintpoints   = (ngl - 2)^2
    nelpoints      = ngl^2
    elnbdypoints   = nelpoints - nelintpoints
    length∂O       = sem.mesh.length∂O
    length∂τ       = sem.mesh.length∂τ
    lengthΓ        = sem.mesh.lengthΓ
    
    RHS   = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(npoin))
    Mdiag = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(npoin))
    
    EL = @time allocate_elemLearning(nelem, ngl,
                                     length∂O,
                                     length∂τ,
                                     lengthΓ,
                                     TFloat, inputs[:backend];
                                     Nsamp=inputs[:Nsamp],
                                     lEL_Sample=inputs[:lEL_Sample])

    if rank == 0 println(BLUE_FG(string(" # ALLOCATE FOR ELEMENT LEARNING ....... DONE"))) end

    BOΓg        = zeros(length∂O)
    gΓ          = zeros(lengthΓ)
    lvtk_sample = false

    if EL.lEL_Sample
        #-----------------------------------------------------
        # 1. Sampling
        #-----------------------------------------------------
        bufferin  = Vector{Vector{Float64}}()
        bufferout = Vector{Vector{Float64}}()
        total_cols_writtenin  = 0
        total_cols_writtenout = 0

        if isfile("input_tensor.csv");  rm("input_tensor.csv");  end
        if isfile("output_tensor.csv"); rm("output_tensor.csv"); end

        # ── Allocate ONCE outside the loop ────────────────────────────────────────
        A       = sem.matrix.L
        A_∂τ∂τ  = A[sem.mesh.∂τ, sem.mesh.∂τ]
        avisc   = zeros(TFloat, 1, ngl^2)          # shape fixed, values change each iter
        nfeatures = size(avisc, 2)

        wbuf = EL_WorkBuffers(params.mesh, A, A_∂τ∂τ, nfeatures,
                              nelintpoints, elnbdypoints,
                              inputs[:NNfile])  # load_inference called ONCE here

        Nsamples = inputs[:Nsamp]
        for isamp = 1:Nsamples
            println(" # --- sample = $isamp")

            # avisc changes each sample — update values in-place, no reallocation
            ranvisc      = 0.5 + rand()
            avisc[1, :] .= ranvisc

            for ip = 1:npoin
                user_source!(RHS[ip], params.qp.qn[ip], params.qp.qe[ip],
                             npoin, inputs[:CL], inputs[:SOL_VARS_TYPE];
                             neqs=1, x=sem.mesh.x[ip], y=sem.mesh.y[ip],
                             xmax=sem.mesh.xmax, xmin=sem.mesh.xmin,
                             ymax=sem.mesh.ymax, ymin=sem.mesh.ymin)
            end
            RHS = sem.matrix.M .* RHS

            apply_boundary_conditions_lin_solve!(sem.matrix.L,
                                                 0.0, params.qp.qe,
                                                 params.mesh.coords,
                                                 params.metrics.nx,
                                                 params.metrics.ny,
                                                 params.metrics.nz,
                                                 npoin,
                                                 params.mesh.npoin_linear,
                                                 params.mesh.poin_in_bdy_edge,
                                                 params.mesh.poin_in_bdy_face,
                                                 params.mesh.nedges_bdy,
                                                 params.mesh.nfaces_bdy,
                                                 ngl, ngr,
                                                 nelem_semi_inf,
                                                 params.basis.ψ, params.basis.dψ,
                                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                 RHS, 0.0, params.ubdy,
                                                 params.mesh.connijk_lag,
                                                 params.mesh.bdy_edge_in_elem,
                                                 params.mesh.bdy_edge_type,
                                                 params.ω, qp.neqs,
                                                 params.inputs, params.AD, sem.mesh.SD)

            # wbuf reused — no new allocations, no new ONNX sessions
            elementLearning_Axb!(params.qp.qn, params.uaux, sem.mesh,
                                 A, RHS, EL,
                                 avisc,
                                 bufferin, bufferout,
                                 BOΓg, gΓ, wbuf;
                                 isamp=isamp,
                                 total_cols_writtenin=total_cols_writtenin,
                                 total_cols_writtenout=total_cols_writtenout)
        end # isamp loop

        total_cols_writtenin  = flush_MLtensor!(bufferin,  total_cols_writtenin,  "input_tensor.csv")
        total_cols_writtenout = flush_MLtensor!(bufferout, total_cols_writtenout, "output_tensor.csv")

        if rank == 0 println(BLUE_FG(" # EL SAMPLING .......... DONE")) end

    else
        #-----------------------------------------------------
        # 2. Inference:
        #-----------------------------------------------------
        #
        # L*q = M*RHS   See algo 12.18 of Giraldo's book
        #
        # 2.a/b
        μ        = 1
        #â        = zeros(TFloat, ngl, ngl)
        avisc      = zeros(TFloat, 1, ngl^2)
        avisc[1,:].= 0.5 + rand() #Uniform distribution between 0.5 and 1.5
        nfeatures  = size(avisc, 2)
        #ψ        = sem.basis.ψ
        #expansion_2d!(â, ψ)

        for ip =1:npoin
            RHS[ip] = user_source!(RHS[ip],
                                   params.qp.qn[ip],
                                   params.qp.qe[ip],
                                   npoin,
                                   inputs[:CL], inputs[:SOL_VARS_TYPE];
                                   neqs=1, x=sem.mesh.x[ip], y=sem.mesh.y[ip],
                                   xmax=sem.mesh.xmax, xmin=sem.mesh.xmin,
                                   ymax=sem.mesh.ymax, ymin=sem.mesh.ymin)
        end
        RHS = sem.matrix.M.*RHS

        apply_boundary_conditions_lin_solve!(sem.matrix.L,
                                             0.0, params.qp.qe,
                                             params.mesh.coords,
                                             params.metrics.nx,
                                             params.metrics.ny,
                                             params.metrics.nz,
                                             npoin,
                                             params.mesh.npoin_linear,
                                             params.mesh.poin_in_bdy_edge,
                                             params.mesh.poin_in_bdy_face,
                                             params.mesh.nedges_bdy,
                                             params.mesh.nfaces_bdy,
                                             ngl, ngr,
                                             nelem_semi_inf,
                                             params.basis.ψ, params.basis.dψ,
                                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                             RHS, 0.0, params.ubdy,
                                             params.mesh.connijk_lag,
                                             params.mesh.bdy_edge_in_elem,
                                             params.mesh.bdy_edge_type,
                                             params.ω, qp.neqs,
                                             params.inputs, params.AD, sem.mesh.SD)

        #-----------------------------------------------------
        # Element-learning infrastructure
        #-----------------------------------------------------
        nfeatures    = size(avisc, 2)
        A            = sem.matrix.L
        A_∂τ∂τ       = A[sem.mesh.∂τ, sem.mesh.∂τ]   # needed by EL_WorkBuffers constructor

        wbuf = EL_WorkBuffers(params.mesh, A, A_∂τ∂τ, nfeatures,
                              nelintpoints, elnbdypoints,
                              inputs[:NNfile])

        total_cols_writtenin  = 0
        total_cols_writtenout = 0

        println(GREEN_FG(string(" # INFERENCE: call to elementLearning_Axb! .......... ")))
        elementLearning_Axb!(params.qp.qn, params.uaux, sem.mesh,
                             A, RHS, EL,
                             avisc,
                             [0.0], [0.0],
                             BOΓg, gΓ, wbuf;
                             isamp=1,
                             total_cols_writtenin=total_cols_writtenin,
                             total_cols_writtenout=total_cols_writtenout)


        println(GREEN_FG(string(" # INFERENCE: call to elementLearning_Axb! .......... DONE")))
        usol = params.qp.qn
        neqs = params.qp.neqs
        args = (params.SD, usol, params.uaux, 1, 1,
                sem.mesh, nothing,
                nothing, nothing,
                0.0, 0.0, 0.0,
                OUTPUT_DIR, inputs,
                params.qp.qvars,
                params.qp.qoutvars,
                inputs[:outformat])

        write_output(args...; nvar=neqs, qexact=params.qp.qe)
        #-----------------------------------------------------
        # END Element-learning infrastructure
        #-----------------------------------------------------
    end

    return nothing
end


# Point evaluation: interpolate at a single point (ξ, η)
function expansion_2d!(a::Matrix, ψ::Matrix)

    # Tensor product form: ψᵀ * A * ψ
    return dot(ψ, a * ψ)

end
