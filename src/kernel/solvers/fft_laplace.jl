# =============================================================================
#  Classical FFT (Fourier spectral) solver for the Laplace / Poisson equation
# =============================================================================
#
#  Solves, on a rectangular PERIODIC domain [x0, x0+Lx] × [y0, y0+Ly],
#
#       -∇²u(x,y) = f(x,y),          u periodic,
#
#  by diagonalizing the Laplacian with the Fourier transform. On a uniform
#  N×M grid the Fourier modes  exp(i(kx x + ky y))  are eigenfunctions of -∇²
#  with eigenvalue (kx² + ky²), so the solve is purely algebraic in Fourier
#  space:
#
#       f̂ = F2D(f)                         (forward 2-D FFT, Fourier coeffs)
#       û_{kx,ky} = f̂_{kx,ky} / (kx²+ky²)  (and û_{0,0}=0 — the constant null space)
#       u = real( F2D⁻¹(û) )               (inverse 2-D FFT)
#
#  This is the "classical FFT Poisson solver". It is an ALTERNATIVE to the
#  standard SEM direct solve (standard_linsolve!) and the element-learning
#  solve (element_learning_linsolve!); the three are selected in
#  problems/drivers.jl. It is spectrally exact: for a band-limited manufactured
#  solution the error is at machine precision.
#
#  Infrastructure: the forward/backward transforms are Kopriva's FFT routines
#  (InitializeFFT, Radix2FFT, Forward2DFFT, Backward2DFFT) in
#  src/kernel/infrastructure/Kopriva_functions.jl, adjusted there to a clean,
#  self-consistent (and unit-tested) convention:
#       Forward2DFFT  : normalized forward  →  Fourier coefficients
#       Backward2DFFT : un-normalized inverse (exact inverse of the forward)
#
#  Because an FFT needs a UNIFORM, periodic, power-of-two grid (which the
#  unstructured GMSH / LGL mesh used by the SEM path is not), this solver builds
#  its OWN tensor-product grid. Geometry and discretization come from
#  user_inputs.jl:
#       :fft_N   => N           number of points per direction (power of 2)
#       :fft_Lx, :fft_Ly        domain lengths            (default 2π)
#       :fft_x0, :fft_y0        lower-left corner         (default 0)
#  and the problem supplies the data through two plain functions in its
#  user_source.jl:
#       user_fft_rhs(x,y)       the right-hand side f          (REQUIRED)
#       user_fft_exact(x,y)     the exact solution u_ex        (optional; enables
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

# ── Driver: classical FFT solve of the Laplace/Poisson equation ──────────────
# Same call signature as standard_linsolve! / element_learning_linsolve! so the
# driver dispatch in problems/drivers.jl can swap solvers transparently.
function fft_linsolve!(sem, params, qp, inputs, OUTPUT_DIR)

    if inputs[:backend] != CPU()
        error(" # fft_linsolve!: the FFT Laplace solver is CPU-only.")
    end

    isdefined(@__MODULE__, :user_fft_rhs) ||
        error(" # fft_linsolve!: define user_fft_rhs(x,y) in the case's user_source.jl")
    has_exact = isdefined(@__MODULE__, :user_fft_exact)

    #-----------------------------------------------------
    # Uniform periodic grid (FFT requires a structured grid)
    #-----------------------------------------------------
    N  = Int(get(inputs, :fft_N, 64))
    M  = Int(get(inputs, :fft_M, N))
    (N > 0 && (N & (N-1)) == 0) || error(" # fft_linsolve!: :fft_N=$N must be a power of 2")
    (M > 0 && (M & (M-1)) == 0) || error(" # fft_linsolve!: :fft_M=$M must be a power of 2")

    Lx = Float64(get(inputs, :fft_Lx, 2π)); Ly = Float64(get(inputs, :fft_Ly, 2π))
    x0 = Float64(get(inputs, :fft_x0, 0.0)); y0 = Float64(get(inputs, :fft_y0, 0.0))

    x = [ x0 + (i-1)*Lx/N for i = 1:N ]     # periodic: node N+1 ≡ node 1 (omitted)
    y = [ y0 + (j-1)*Ly/M for j = 1:M ]

    println(YELLOW_FG(string(" # Solve -∇²u = f by classical FFT (Fourier spectral): ",
                             N, "×", M, " periodic grid ..............")))

    #-----------------------------------------------------
    # Sample the RHS (and exact field) on the grid
    #-----------------------------------------------------
    F = Matrix{Float64}(undef, N, M)
    @inbounds for j = 1:M, i = 1:N
        F[i,j] = Float64(user_fft_rhs(x[i], y[j]))
    end

    # Periodic Poisson is solvable only for a zero-mean RHS (∫f = 0). Project the
    # mean out and warn if it was non-negligible (the (0,0) mode is dropped anyway).
    fbar = sum(F) / (N*M)
    if abs(fbar) > 1e-10
        println(string(" # fft_linsolve!: RHS mean = ", fbar,
                       " ≠ 0; removing it (periodic compatibility condition)."))
        F .-= fbar
    end

    #-----------------------------------------------------
    # Spectral solve
    #-----------------------------------------------------
    kx = fft_wavenumbers(N, Lx)
    ky = fft_wavenumbers(M, Ly)
    u  = fft_poisson_solve(F, kx, ky)       # zero-mean solution

    println(YELLOW_FG(string(" # Solve -∇²u = f by classical FFT ............................... DONE")))

    #-----------------------------------------------------
    # Verification against the manufactured solution (if provided)
    #-----------------------------------------------------
    uex = nothing
    err = nothing
    if has_exact
        uex = Matrix{Float64}(undef, N, M)
        @inbounds for j = 1:M, i = 1:N
            uex[i,j] = Float64(user_fft_exact(x[i], y[j]))
        end
        # The periodic solution is defined up to a constant; pin it to the exact
        # field's mean so the comparison measures the genuine (mode-by-mode) error.
        u .+= (sum(uex) - sum(u)) / (N*M)

        err   = u .- uex
        dA    = (Lx*Ly) / (N*M)             # uniform-grid quadrature weight
        l2    = sqrt(sum(abs2, err) * dA)
        ref2  = sqrt(sum(abs2, uex) * dA)
        linf  = maximum(abs, err)
        relstr = ref2 > 0 ? string(" , relative ‖e‖_L2 = ", l2/ref2) : ""
        println(GREEN_FG(string(" # MMS verification: FFT solve vs exact  →  ‖e‖_L2 = ", l2,
                                relstr, " , ‖e‖_∞ = ", linf)))
    end

    #-----------------------------------------------------
    # Output
    #-----------------------------------------------------
    vtkpath = joinpath(OUTPUT_DIR, "fft_laplace.vtk")
    write_fft_vtk(vtkpath, x, y, u, uex, err)
    println(string(" # FFT solution written to ", vtkpath))

    return u
end
