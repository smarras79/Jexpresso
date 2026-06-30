# =============================================================================
#  Classical FFT (Fourier spectral) solver for the Laplace / Poisson equation
# =============================================================================
#
#  Solves, on a PERIODIC rectangular box (2D or 3D),
#
#       -∇²u = f,          u periodic,
#
#  by diagonalizing the Laplacian with the Fourier transform. On a uniform grid
#  the Fourier modes exp(i k·x) are eigenfunctions of -∇² with eigenvalue |k|²,
#  so the solve is purely algebraic in Fourier space:
#
#       f̂ = FFT(f)                         (forward FFT, Fourier coeffs)
#       û_k = f̂_k / |k|²                   (and û_0 = 0 — the constant null space)
#       u = real( FFT⁻¹(û) )               (inverse FFT)
#
#  This is the "classical FFT Poisson solver". It is an ALTERNATIVE to the
#  standard SEM direct solve (standard_linsolve!) and the element-learning
#  solve (element_learning_linsolve!); the three are selected in
#  problems/drivers.jl. It is spectrally exact: for a band-limited manufactured
#  solution the error is at machine precision.
#
#  Infrastructure: the building-block transforms are Kopriva's FFT routines
#  (InitializeFFT, Radix2FFT, Forward2DFFT, Backward2DFFT) in
#  src/kernel/infrastructure/Kopriva_functions.jl, adjusted there to a clean,
#  self-consistent (and unit-tested) convention. The N-D core here (fft_nd,
#  fft_poisson_solve_nd) applies the 1-D Radix2FFT along each axis, so the same
#  code path serves 2-D and 3-D.
#
#  TWO grid sources, selected by :fft_use_mesh in user_inputs.jl:
#
#   • :fft_use_mesh => false (default) — SYNTHETIC grid (2D, mesh-independent):
#       :fft_N   => N           points per direction (power of 2)
#       :fft_Lx, :fft_Ly        domain lengths            (default 2π)
#       :fft_x0, :fft_y0        lower-left corner         (default 0)
#
#   • :fft_use_mesh => true — solve ON the structured periodic GMSH mesh that
#     Jexpresso reads (2D or 3D, by mesh.nsd; requires :nop => 1 so the global
#     nodes are equispaced). Period per axis from :fft_Lx/:fft_Ly[/:fft_Lz]. The
#     solution is scattered back to the mesh nodes and written by the normal
#     Jexpresso output.
#
#  The problem supplies the data through plain functions in its user_source.jl
#  (give the arity matching the mesh dimension; z may default for 2D/3D reuse):
#       user_fft_rhs(x,y[,z])   the right-hand side f          (REQUIRED)
#       user_fft_exact(x,y[,z]) the exact solution u_ex        (optional; enables
#                               the automatic L2-error verification)
# =============================================================================

# ── 1-D signed Fourier wavenumbers for an N-point periodic grid of length L ──
# FFT index m = 0:N-1 maps to the signed mode m̃ = m (m ≤ N/2) or m-N (m > N/2);
# the physical wavenumber is k = 2π m̃ / L. The Nyquist mode (m = N/2) is kept at
# +N/2; its k² is the same either way, so the Laplacian eigenvalue is unaffected.
function fft_wavenumbers(N::Int, L::Real)
    k = Vector{Float64}(undef, N)
    @inbounds for m = 0:N-1
        m̃ = m <= N ÷ 2 ? m : m - N
        k[m+1] = 2π * m̃ / L
    end
    return k
end

# ── Core: spectral solve of -∇²u = f on a periodic N×M grid ──────────────────
# `F` is the sampled RHS (N×M). `kx`,`ky` are the signed wavenumbers from
# fft_wavenumbers. Returns the real solution with zero mean (the constant null
# space of the periodic Laplacian is fixed by setting the (0,0) mode to zero).
function fft_poisson_solve(F::AbstractMatrix, kx::AbstractVector, ky::AbstractVector)
    N = size(F, 1); M = size(F, 2)
    wf1 = InitializeFFT(N,  1); wf2 = InitializeFFT(M,  1)
    wb1 = InitializeFFT(N, -1); wb2 = InitializeFFT(M, -1)

    F̂ = Forward2DFFT(F, wf1, wf2)
    Û = zeros(ComplexF64, N, M)
    @inbounds for j = 1:M, i = 1:N
        λ = kx[i]*kx[i] + ky[j]*ky[j]      # eigenvalue of -∇² for this mode
        Û[i,j] = λ == 0.0 ? 0.0 : F̂[i,j] / λ
    end
    return real.(Backward2DFFT(Û, wb1, wb2))
end

# ── Legacy-VTK STRUCTURED_POINTS writer for the uniform FFT grid ──────────────
# Writes u (and, when available, the exact field and the point-wise error) as
# point data on the regular grid so the result opens directly in ParaView,
# matching the :outformat => "vtk" convention used elsewhere.
function write_fft_vtk(path, x, y, u, uex, err)
    N = length(x); M = length(y)
    dx = N > 1 ? x[2]-x[1] : 1.0
    dy = M > 1 ? y[2]-y[1] : 1.0
    open(path, "w") do io
        println(io, "# vtk DataFile Version 3.0")
        println(io, "Jexpresso FFT Laplace/Poisson solution")
        println(io, "ASCII")
        println(io, "DATASET STRUCTURED_POINTS")
        println(io, "DIMENSIONS $N $M 1")
        println(io, "ORIGIN $(x[1]) $(y[1]) 0.0")
        println(io, "SPACING $dx $dy 1.0")
        println(io, "POINT_DATA $(N*M)")
        println(io, "SCALARS u double 1")
        println(io, "LOOKUP_TABLE default")
        @inbounds for j = 1:M, i = 1:N
            println(io, u[i,j])
        end
        if uex !== nothing
            println(io, "SCALARS u_exact double 1")
            println(io, "LOOKUP_TABLE default")
            @inbounds for j = 1:M, i = 1:N
                println(io, uex[i,j])
            end
            println(io, "SCALARS error double 1")
            println(io, "LOOKUP_TABLE default")
            @inbounds for j = 1:M, i = 1:N
                println(io, err[i,j])
            end
        end
    end
    return path
end

# ── Recover the uniform tensor grid from a structured PERIODIC mesh ──────────
# Classical FFT needs an equispaced grid. The global Jexpresso SEM grid is
# equispaced only at :nop => 1 (LGL with 2 nodes/element ⇒ the element corners);
# at nop>1 the in-element LGL nodes are non-uniform and an FFT cannot use them.
#
# This maps every mesh node ip to a tensor index on the periodic lattice by
# folding each of its coordinates into the period: φ_d = mod(x_d - min_d, L_d)/L_d.
# Folding makes a "closed" mesh (a node at both x=min and x=min+L) and a
# Jexpresso periodic-merged mesh (only x=min kept) collapse to the SAME lines —
# the seam node lands on line 1 either way. Dimension-agnostic: ND=2 uses
# (x,y); ND=3 uses (x,y,z). Returns
#   (dims::NTuple{ND,Int}, lines::NTuple{ND,Vector}, idxof::Vector{NTuple{ND,Int}})
# where lines[d][i]=min_d+(i-1)L_d/dims[d] and idxof[ip] is the node's lattice
# index. Errors if the folded lines are not uniformly spaced (⇒ nop>1 or a
# non-uniform mesh) or any axis count is not a power of two.
function fft_grid_from_mesh(mesh, Ls::NTuple{ND,Float64}) where {ND}
    npoin  = Int(mesh.npoin)
    coords = ND == 3 ? (mesh.x, mesh.y, mesh.z) : (mesh.x, mesh.y)
    mins   = ntuple(d -> minimum(@view coords[d][1:npoin]), ND)
    dirs   = ("x", "y", "z")
    snap   = 5e-8                          # fold the period seam (φ≈1) back to 0

    ph = ntuple(d -> Vector{Float64}(undef, npoin), ND)
    @inbounds for d = 1:ND, ip = 1:npoin
        p = mod(coords[d][ip] - mins[d], Ls[d]) / Ls[d]
        ph[d][ip] = (p > 1 - snap || p < snap) ? 0.0 : p
    end

    dims = ntuple(ND) do d
        u = sort(unique(round.(ph[d], digits = 7)))   # distinct periodic lines in [0,1)
        N = length(u)
        _fft_assert_uniform(u, N, dirs[d])
        (N > 1 && (N & (N-1)) == 0) ||
            error(" # fft_linsolve!: mesh has N=$N points along $(dirs[d]); an FFT needs a power of 2.")
        N
    end

    idxof = Vector{NTuple{ND,Int}}(undef, npoin)
    @inbounds for ip = 1:npoin
        idxof[ip] = ntuple(d -> mod(round(Int, ph[d][ip]*dims[d]), dims[d]) + 1, ND)
    end

    # every lattice cell must be reached by a node (i.e. a genuine tensor grid)
    seen = falses(dims...)
    @inbounds for ip = 1:npoin
        seen[idxof[ip]...] = true
    end
    all(seen) || error(" # fft_linsolve!: the mesh nodes do not fill a $(join(dims, "×")) "*
                       "tensor grid (is it structured and periodic?).")

    lines = ntuple(d -> [ mins[d] + (i-1)*Ls[d]/dims[d] for i = 1:dims[d] ], ND)
    return dims, lines, idxof
end

# Assert the folded grid lines `u` are at the uniform phases (k-1)/N.
function _fft_assert_uniform(u, N, dir)
    @inbounds for k = 1:N
        if abs(u[k] - (k-1)/N) > 1e-4
            error(" # fft_linsolve!: the mesh is NOT equispaced along $dir (line $k at "*
                  "phase $(u[k]) ≠ $((k-1)/N)). Classical FFT needs a uniform grid — use "*
                  ":nop => 1 on a uniform structured periodic mesh, or :fft_use_mesh => false.")
        end
    end
end

# Project the RHS onto the zero-mean subspace (the periodic Poisson solvability /
# compatibility condition ∫f = 0); warn if a non-negligible mean was removed.
function fft_enforce_zero_mean!(F)
    fbar = sum(F) / length(F)
    if abs(fbar) > 1e-10
        println(string(" # fft_linsolve!: RHS mean = ", fbar,
                       " ≠ 0; removing it (periodic compatibility condition)."))
        F .-= fbar
    end
    return F
end

# Evaluate the user RHS / exact field at an nd-tuple of coordinates (2D or 3D).
# The case's user_fft_rhs / user_fft_exact must accept the matching arity.
_eval_fft_rhs(c::NTuple{2}) = user_fft_rhs(c[1], c[2])
_eval_fft_rhs(c::NTuple{3}) = user_fft_rhs(c[1], c[2], c[3])
_eval_fft_exact(c::NTuple{2}) = user_fft_exact(c[1], c[2])
_eval_fft_exact(c::NTuple{3}) = user_fft_exact(c[1], c[2], c[3])

# Separable N-D FFT: apply the 1-D Radix2FFT along each axis in turn. `ws[d]` are
# the twiddles for axis d (forward s=+1 or backward s=-1, from InitializeFFT).
# `normalize` divides by the total point count — used on the forward transform so
# it returns the Fourier coefficients (matching the 2-D Forward2DFFT convention).
function fft_nd(A, ws; normalize::Bool)
    B = ComplexF64.(A)
    for d = 1:ndims(B)
        B = mapslices(v -> Radix2FFT(v, ws[d]), B; dims = d)
    end
    normalize && (B ./= length(B))
    return B
end

# Solve -∇²u = f on an N-D periodic grid by spectral diagonalization. `F` is the
# sampled RHS (its size sets the per-axis point counts); `Ls` the per-axis period.
# The constant null mode (|k|=0) is set to zero ⇒ the returned u is zero-mean.
function fft_poisson_solve_nd(F, Ls)
    dims = size(F)
    ND   = ndims(F)
    ks   = ntuple(d -> fft_wavenumbers(dims[d], Ls[d]), ND)
    wf   = ntuple(d -> InitializeFFT(dims[d],  1), ND)
    wb   = ntuple(d -> InitializeFFT(dims[d], -1), ND)
    F̂    = fft_nd(F, wf; normalize = true)
    Û    = Array{ComplexF64}(undef, dims)
    @inbounds for I in CartesianIndices(F̂)
        λ = 0.0
        for d = 1:ND
            λ += ks[d][I[d]]^2
        end
        Û[I] = λ == 0.0 ? zero(ComplexF64) : F̂[I] / λ
    end
    return real.(fft_nd(Û, wb; normalize = false))
end

# L2 / L∞ error of the grid solution `u` vs the exact field sampled on the same
# lines. `lines` = per-axis coordinate vectors, `Ls` = per-axis periods (2D or
# 3D). The periodic solution is unique only up to a constant, so the constant is
# pinned to the exact field's mean before comparing.
function fft_report_grid_error(u, lines, Ls)
    ND  = ndims(u)
    uex = Array{Float64}(undef, size(u))
    @inbounds for I in CartesianIndices(u)
        coord  = ntuple(d -> lines[d][I[d]], ND)
        uex[I] = Float64(_eval_fft_exact(coord))
    end
    u .+= (sum(uex) - sum(u)) / length(u)
    err  = u .- uex
    dV   = prod(Ls) / length(u)
    l2   = sqrt(sum(abs2, err) * dV)
    ref2 = sqrt(sum(abs2, uex) * dV)
    linf = maximum(abs, err)
    relstr = ref2 > 0 ? string(" , relative ‖e‖_L2 = ", l2/ref2) : ""
    println(GREEN_FG(string(" # MMS verification: FFT solve vs exact  →  ‖e‖_L2 = ", l2,
                            relstr, " , ‖e‖_∞ = ", linf)))
    return uex, err
end

# ── Driver: classical FFT solve of the Laplace/Poisson equation ──────────────
# Same call signature as standard_linsolve! / element_learning_linsolve! so the
# driver dispatch in problems/drivers.jl can swap solvers transparently.
#
# Two grid sources, selected by :fft_use_mesh:
#   false (default) : build a synthetic uniform periodic grid from :fft_N,
#                     :fft_Lx/:fft_Ly, :fft_x0/:fft_y0 (mesh unused).
#   true            : solve ON the structured periodic mesh read by Jexpresso
#                     (requires :nop => 1; period set by :fft_Lx/:fft_Ly), then
#                     scatter the solution back to the mesh nodes and write it
#                     through the standard Jexpresso output.
function fft_linsolve!(sem, params, qp, inputs, OUTPUT_DIR)

    if inputs[:backend] != CPU()
        error(" # fft_linsolve!: the FFT Laplace solver is CPU-only.")
    end
    isdefined(@__MODULE__, :user_fft_rhs) ||
        error(" # fft_linsolve!: define user_fft_rhs(x,y) in the case's user_source.jl")
    has_exact = isdefined(@__MODULE__, :user_fft_exact)

    if get(inputs, :fft_use_mesh, false)
        return fft_linsolve_on_mesh!(sem, params, inputs, OUTPUT_DIR, has_exact)
    end

    #=====================================================================
      Synthetic-grid mode (mesh-independent)
    =====================================================================#
    N  = Int(get(inputs, :fft_N, 64))
    M  = Int(get(inputs, :fft_M, N))
    (N > 0 && (N & (N-1)) == 0) || error(" # fft_linsolve!: :fft_N=$N must be a power of 2")
    (M > 0 && (M & (M-1)) == 0) || error(" # fft_linsolve!: :fft_M=$M must be a power of 2")

    Lx = Float64(get(inputs, :fft_Lx, 2π)); Ly = Float64(get(inputs, :fft_Ly, 2π))
    x0 = Float64(get(inputs, :fft_x0, 0.0)); y0 = Float64(get(inputs, :fft_y0, 0.0))

    x = [ x0 + (i-1)*Lx/N for i = 1:N ]     # periodic: node N+1 ≡ node 1 (omitted)
    y = [ y0 + (j-1)*Ly/M for j = 1:M ]

    println(YELLOW_FG(string(" # Solve -∇²u = f by classical FFT (Fourier spectral): ",
                             N, "×", M, " synthetic periodic grid ..............")))

    F = Matrix{Float64}(undef, N, M)
    @inbounds for j = 1:M, i = 1:N
        F[i,j] = Float64(user_fft_rhs(x[i], y[j]))
    end
    fft_enforce_zero_mean!(F)

    u = fft_poisson_solve(F, fft_wavenumbers(N, Lx), fft_wavenumbers(M, Ly))

    println(YELLOW_FG(string(" # Solve -∇²u = f by classical FFT ............................... DONE")))

    uex = nothing; err = nothing
    has_exact && ((uex, err) = fft_report_grid_error(u, (x, y), (Lx, Ly)))

    vtkpath = joinpath(OUTPUT_DIR, "fft_laplace.vtk")
    write_fft_vtk(vtkpath, x, y, u, uex, err)
    println(string(" # FFT solution written to ", vtkpath))
    return u
end

# ── Mesh-grid mode: solve on the actual structured periodic mesh (2D or 3D) ──
function fft_linsolve_on_mesh!(sem, params, inputs, OUTPUT_DIR, has_exact)
    mesh = sem.mesh
    ND   = Int(mesh.nsd)
    (ND == 2 || ND == 3) ||
        error(" # fft_linsolve!: mesh nsd=$ND not supported by the FFT solver (need 2 or 3).")

    # Per-axis period: from inputs if given, else the mesh extent.
    Ls = ND == 3 ?
        (Float64(get(inputs, :fft_Lx, mesh.xmax - mesh.xmin)),
         Float64(get(inputs, :fft_Ly, mesh.ymax - mesh.ymin)),
         Float64(get(inputs, :fft_Lz, mesh.zmax - mesh.zmin))) :
        (Float64(get(inputs, :fft_Lx, mesh.xmax - mesh.xmin)),
         Float64(get(inputs, :fft_Ly, mesh.ymax - mesh.ymin)))

    dims, lines, idxof = fft_grid_from_mesh(mesh, Ls)

    println(YELLOW_FG(string(" # Solve -∇²u = f by classical FFT (Fourier spectral) on the ",
                             join(dims, "×"), " periodic mesh grid ..............")))

    # RHS sampled on the canonical grid lines (periodic ⇒ seam value is unique)
    F = Array{Float64}(undef, dims)
    @inbounds for I in CartesianIndices(F)
        coord = ntuple(d -> lines[d][I[d]], ND)
        F[I]  = Float64(_eval_fft_rhs(coord))
    end
    fft_enforce_zero_mean!(F)

    ugrid = fft_poisson_solve_nd(F, Ls)

    println(YELLOW_FG(string(" # Solve -∇²u = f by classical FFT ............................... DONE")))

    has_exact && fft_report_grid_error(ugrid, lines, Ls)

    # Scatter the grid solution back onto every mesh node (seam duplicates get the
    # same periodic value), then write through the standard Jexpresso output so the
    # field is visualized on the real mesh.
    npoin = Int(mesh.npoin)
    sol   = Vector{TFloat}(undef, npoin)
    @inbounds for ip = 1:npoin
        sol[ip] = ugrid[idxof[ip]...]
    end

    args = (params.SD, sol, params.uaux, 1, 1,
            mesh, nothing, nothing, nothing,
            0.0, 0.0, 0.0, OUTPUT_DIR, inputs,
            params.qp.qvars, params.qp.qoutvars, inputs[:outformat])
    write_output(args...; nvar=params.qp.neqs, qexact=params.qp.qe, metrics=params.metrics)
    println(string(" # FFT solution written to ", OUTPUT_DIR, " (on the mesh grid)"))
    return sol
end
