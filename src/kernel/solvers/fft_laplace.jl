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

# ── Recover the uniform tensor grid from a structured PERIODIC mesh ──────────
# Classical FFT needs an equispaced grid. The global Jexpresso SEM grid is
# equispaced only at :nop => 1 (LGL with 2 nodes/element ⇒ the element corners);
# at nop>1 the in-element LGL nodes are non-uniform and an FFT cannot use them.
#
# This maps every mesh node ip to a tensor index (ix,iy) on the Nx×Ny periodic
# lattice by folding its coordinate into the period: φ = mod(x-xmin, Lx)/Lx.
# Folding makes a "closed" mesh (a node at both x=xmin and x=xmin+Lx) and a
# Jexpresso periodic-merged mesh (only x=xmin kept) collapse to the SAME Nx
# lines — the seam node lands on line 1 either way. Returns
#   (Nx, Ny, xline, yline, ixof, iyof)
# with xline[i]=xmin+(i-1)Lx/Nx, ixof[ip]∈1:Nx. Errors if the folded lines are
# not uniformly spaced (⇒ nop>1 or a non-uniform mesh) or not a power of two.
function fft_grid_from_mesh(mesh, Lx::Float64, Ly::Float64)
    npoin = Int(mesh.npoin)
    xmin  = minimum(@view mesh.x[1:npoin])
    ymin  = minimum(@view mesh.y[1:npoin])
    snap  = 5e-8                          # fold the period seam (φ≈1) back to 0

    phx = Vector{Float64}(undef, npoin)
    phy = Vector{Float64}(undef, npoin)
    @inbounds for ip = 1:npoin
        px = mod(mesh.x[ip]-xmin, Lx)/Lx
        py = mod(mesh.y[ip]-ymin, Ly)/Ly
        phx[ip] = (px > 1-snap || px < snap) ? 0.0 : px
        phy[ip] = (py > 1-snap || py < snap) ? 0.0 : py
    end

    ux = sort(unique(round.(phx, digits=7)))   # distinct periodic grid lines in [0,1)
    uy = sort(unique(round.(phy, digits=7)))
    Nx = length(ux); Ny = length(uy)

    _fft_assert_uniform(ux, Nx, "x")
    _fft_assert_uniform(uy, Ny, "y")
    (Nx > 1 && (Nx & (Nx-1)) == 0) ||
        error(" # fft_linsolve!: mesh has Nx=$Nx points along x; an FFT needs a power of 2.")
    (Ny > 1 && (Ny & (Ny-1)) == 0) ||
        error(" # fft_linsolve!: mesh has Ny=$Ny points along y; an FFT needs a power of 2.")

    ixof = Vector{Int}(undef, npoin); iyof = Vector{Int}(undef, npoin)
    @inbounds for ip = 1:npoin
        ixof[ip] = mod(round(Int, phx[ip]*Nx), Nx) + 1
        iyof[ip] = mod(round(Int, phy[ip]*Ny), Ny) + 1
    end

    # every lattice cell must be reached by a node (i.e. a genuine tensor grid)
    seen = falses(Nx, Ny)
    @inbounds for ip = 1:npoin
        seen[ixof[ip], iyof[ip]] = true
    end
    all(seen) || error(" # fft_linsolve!: the mesh nodes do not fill an $Nx×$Ny tensor "*
                       "grid (is it structured and periodic?).")

    xline = [ xmin + (i-1)*Lx/Nx for i = 1:Nx ]
    yline = [ ymin + (j-1)*Ly/Ny for j = 1:Ny ]
    return Nx, Ny, xline, yline, ixof, iyof
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

# L2 / L∞ error of the grid solution `u` (Nx×Ny) vs the exact field sampled on
# the same lines. The periodic solution is unique only up to a constant, so the
# constant is pinned to the exact field's mean before comparing.
function fft_report_grid_error(u, xline, yline, Lx, Ly)
    Nx = length(xline); Ny = length(yline)
    uex = Matrix{Float64}(undef, Nx, Ny)
    @inbounds for j = 1:Ny, i = 1:Nx
        uex[i,j] = Float64(user_fft_exact(xline[i], yline[j]))
    end
    u .+= (sum(uex) - sum(u)) / (Nx*Ny)
    err  = u .- uex
    dA   = (Lx*Ly) / (Nx*Ny)
    l2   = sqrt(sum(abs2, err) * dA)
    ref2 = sqrt(sum(abs2, uex) * dA)
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
    has_exact && ((uex, err) = fft_report_grid_error(u, x, y, Lx, Ly))

    vtkpath = joinpath(OUTPUT_DIR, "fft_laplace.vtk")
    write_fft_vtk(vtkpath, x, y, u, uex, err)
    println(string(" # FFT solution written to ", vtkpath))
    return u
end

# ── Mesh-grid mode: solve on the actual structured periodic mesh ─────────────
function fft_linsolve_on_mesh!(sem, params, inputs, OUTPUT_DIR, has_exact)
    mesh = sem.mesh
    Lx = Float64(get(inputs, :fft_Lx, mesh.xmax - mesh.xmin))
    Ly = Float64(get(inputs, :fft_Ly, mesh.ymax - mesh.ymin))

    Nx, Ny, xline, yline, ixof, iyof = fft_grid_from_mesh(mesh, Lx, Ly)

    println(YELLOW_FG(string(" # Solve -∇²u = f by classical FFT (Fourier spectral) on the ",
                             Nx, "×", Ny, " periodic mesh grid ..............")))

    # RHS sampled on the canonical grid lines (periodic ⇒ seam value is unique)
    F = Matrix{Float64}(undef, Nx, Ny)
    @inbounds for j = 1:Ny, i = 1:Nx
        F[i,j] = Float64(user_fft_rhs(xline[i], yline[j]))
    end
    fft_enforce_zero_mean!(F)

    ugrid = fft_poisson_solve(F, fft_wavenumbers(Nx, Lx), fft_wavenumbers(Ny, Ly))

    println(YELLOW_FG(string(" # Solve -∇²u = f by classical FFT ............................... DONE")))

    has_exact && fft_report_grid_error(ugrid, xline, yline, Lx, Ly)

    # Scatter the grid solution back onto every mesh node (seam duplicates get the
    # same periodic value), then write through the standard Jexpresso output so the
    # field is visualized on the real mesh.
    npoin = Int(mesh.npoin)
    sol   = Vector{TFloat}(undef, npoin)
    @inbounds for ip = 1:npoin
        sol[ip] = ugrid[ixof[ip], iyof[ip]]
    end

    args = (params.SD, sol, params.uaux, 1, 1,
            mesh, nothing, nothing, nothing,
            0.0, 0.0, 0.0, OUTPUT_DIR, inputs,
            params.qp.qvars, params.qp.qoutvars, inputs[:outformat])
    write_output(args...; nvar=params.qp.neqs, qexact=params.qp.qe, metrics=params.metrics)
    println(string(" # FFT solution written to ", OUTPUT_DIR, " (on the mesh grid)"))
    return sol
end
