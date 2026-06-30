# =============================================================================
#  Chebyshev spectral (collocation) solver for the Laplace / Poisson equation
# =============================================================================
#
#  Solves the NON-PERIODIC problem
#
#       -∇²u(x,y) = f(x,y)   on a box [ax,bx] × [ay,by],
#        u = g(x,y)          on the boundary (Dirichlet),
#
#  with a global Chebyshev spectral method. The Fourier solver (fft_laplace.jl)
#  is the right tool for PERIODIC BCs (uniform grid); for NON-periodic BCs the
#  analogous spectral basis is Chebyshev, whose Chebyshev–Gauss–Lobatto (CGL)
#  nodes  x_j = cos(jπ/N)  cluster toward the boundaries (that clustering is what
#  defeats the Runge phenomenon and resolves boundary layers — see Kopriva).
#
#  Method: Chebyshev COLLOCATION. The 1-D CGL differentiation matrix comes from
#  Jexpresso's existing infrastructure — Kopriva's `CGLDerivativeMatrix`
#  (src/kernel/infrastructure/Kopriva_functions.jl). The 2-D Laplacian on the
#  tensor grid is  ∇² = D²ₓ ⊗ I + I ⊗ D²ᵧ ; Dirichlet rows are replaced by the
#  identity, and the resulting linear system is solved directly. It is spectrally
#  accurate (exact to round-off for polynomial data; exponential convergence for
#  smooth data).
#
#  This is selected in problems/drivers.jl by :lcheb (alongside :lfft and the SEM
#  / element-learning solves). Discretization & geometry come from user_inputs.jl:
#       :cheb_N                      nodes per direction (or :cheb_Nx / :cheb_Ny)
#       :cheb_xmin,:cheb_xmax        x-extent of the box  (default [-1, 1])
#       :cheb_ymin,:cheb_ymax        y-extent of the box  (default [-1, 1])
#  and the problem supplies, in its user_source.jl:
#       user_cheb_rhs(x,y)     the right-hand side f               (REQUIRED)
#       user_cheb_bc(x,y)      the Dirichlet boundary value g      (optional; if
#                              absent, falls back to user_cheb_exact, else 0)
#       user_cheb_exact(x,y)   the exact solution u_ex             (optional;
#                              enables the automatic L∞-error verification)
# =============================================================================

using LinearAlgebra

# ── CGL nodes on the reference interval [-1,1]:  ξ_j = cos(jπ/N), j = 0…N ─────
cheb_cgl_nodes(N::Int) = [ cospi(j/N) for j = 0:N ]

# ── 1-D second-derivative collocation matrix on [a,b] ────────────────────────
# Built from Kopriva's CGL first-derivative matrix D (reference [-1,1]) squared,
# then scaled by (dξ/dx)² = (2/(b-a))² for the physical interval [a,b].
function cheb_second_derivative_matrix(ξ::AbstractVector, N::Int, a::Real, b::Real)
    D = CGLDerivativeMatrix(ξ, N)             # Kopriva infrastructure
    s = 2.0 / (b - a)
    return (s*s) .* (D * D)
end

# ── Core: 2-D Chebyshev-collocation solve of -∇²u = f, Dirichlet g ───────────
# Returns (x, y, U) with x (length Nx+1), y (length Ny+1) the physical CGL node
# coordinates and U[i,j] = u(x_i, y_j).
function cheb_poisson_solve_2d(Nx::Int, Ny::Int,
                               ax::Real, bx::Real, ay::Real, by::Real,
                               frhs, gbc)
    ξ = cheb_cgl_nodes(Nx);  η = cheb_cgl_nodes(Ny)
    D2x = cheb_second_derivative_matrix(ξ, Nx, ax, bx)
    D2y = cheb_second_derivative_matrix(η, Ny, ay, by)

    # physical node coordinates (ξ=+1 ↦ b, ξ=-1 ↦ a)
    x = ax .+ (bx - ax) .* (ξ .+ 1) ./ 2
    y = ay .+ (by - ay) .* (η .+ 1) ./ 2

    nx1 = Nx + 1;  ny1 = Ny + 1;  npts = nx1*ny1
    # 2-D Laplacian on the tensor grid, node order P = (i-1)*ny1 + j
    Lap = kron(D2x, Matrix(I, ny1, ny1)) .+ kron(Matrix(I, nx1, nx1), D2y)
    M   = -Lap
    rhs = zeros(npts)

    @inbounds for i = 1:nx1, j = 1:ny1
        p = (i-1)*ny1 + j
        if i == 1 || i == nx1 || j == 1 || j == ny1
            M[p, :] .= 0.0;  M[p, p] = 1.0       # Dirichlet row
            rhs[p]   = Float64(gbc(x[i], y[j]))
        else
            rhs[p]   = Float64(frhs(x[i], y[j]))
        end
    end

    u = M \ rhs
    U = Matrix(reshape(u, ny1, nx1)')            # U[i,j] = u(x_i, y_j)
    return x, y, U
end

# ── STRUCTURED_GRID VTK writer (the CGL grid is non-uniform) ─────────────────
function write_cheb_vtk(path, x, y, U, Uex, Err)
    nx1 = length(x); ny1 = length(y)
    open(path, "w") do io
        println(io, "# vtk DataFile Version 3.0")
        println(io, "Jexpresso Chebyshev Laplace/Poisson solution")
        println(io, "ASCII")
        println(io, "DATASET STRUCTURED_GRID")
        println(io, "DIMENSIONS $ny1 $nx1 1")
        println(io, "POINTS $(nx1*ny1) double")
        @inbounds for i = 1:nx1, j = 1:ny1
            println(io, x[i], " ", y[j], " 0.0")
        end
        println(io, "POINT_DATA $(nx1*ny1)")
        println(io, "SCALARS u double 1")
        println(io, "LOOKUP_TABLE default")
        @inbounds for i = 1:nx1, j = 1:ny1
            println(io, U[i,j])
        end
        if Uex !== nothing
            println(io, "SCALARS u_exact double 1")
            println(io, "LOOKUP_TABLE default")
            @inbounds for i = 1:nx1, j = 1:ny1; println(io, Uex[i,j]); end
            println(io, "SCALARS error double 1")
            println(io, "LOOKUP_TABLE default")
            @inbounds for i = 1:nx1, j = 1:ny1; println(io, Err[i,j]); end
        end
    end
    return path
end

# ── Driver: Chebyshev collocation solve of the Laplace/Poisson equation ──────
# Same call signature as the other linear solvers so problems/drivers.jl can
# dispatch transparently.
function cheb_linsolve!(sem, params, qp, inputs, OUTPUT_DIR)

    if inputs[:backend] != CPU()
        error(" # cheb_linsolve!: the Chebyshev Laplace solver is CPU-only.")
    end
    isdefined(@__MODULE__, :user_cheb_rhs) ||
        error(" # cheb_linsolve!: define user_cheb_rhs(x,y) in the case's user_source.jl")

    has_exact = isdefined(@__MODULE__, :user_cheb_exact)
    # Dirichlet data: explicit user_cheb_bc, else the exact field, else 0.
    gbc = isdefined(@__MODULE__, :user_cheb_bc) ? user_cheb_bc :
          (has_exact ? user_cheb_exact : (x, y) -> 0.0)

    Nx = Int(get(inputs, :cheb_Nx, get(inputs, :cheb_N, 24)))
    Ny = Int(get(inputs, :cheb_Ny, get(inputs, :cheb_N, Nx)))
    ax = Float64(get(inputs, :cheb_xmin, -1.0)); bx = Float64(get(inputs, :cheb_xmax, 1.0))
    ay = Float64(get(inputs, :cheb_ymin, -1.0)); by = Float64(get(inputs, :cheb_ymax, 1.0))

    println(YELLOW_FG(string(" # Solve -∇²u = f by Chebyshev collocation (spectral): ",
                             Nx, "×", Ny, " CGL grid on [", ax, ",", bx, "]×[",
                             ay, ",", by, "] ..............")))

    x, y, U = jx_time_solve("Chebyshev collocation solve", () ->
                  cheb_poisson_solve_2d(Nx, Ny, ax, bx, ay, by, user_cheb_rhs, gbc))

    println(YELLOW_FG(string(" # Solve -∇²u = f by Chebyshev collocation ...................... DONE")))

    Uex = nothing; Err = nothing
    if has_exact
        Uex = [ Float64(user_cheb_exact(x[i], y[j])) for i = 1:Nx+1, j = 1:Ny+1 ]
        Err = U .- Uex
        linf  = maximum(abs, Err)
        ref   = maximum(abs, Uex)
        l2rel = sqrt(sum(abs2, Err) / sum(abs2, Uex))   # relative discrete L2
        relstr = ref > 0 ? string(" , relative L∞ = ", linf/ref) : ""
        println(GREEN_FG(string(" # MMS verification: Chebyshev solve vs exact  →  ‖e‖_∞ = ",
                                linf, relstr)))
        jx_record_solve_error(; linf = linf, l2rel = l2rel, npts = (Nx+1)*(Ny+1))
    end

    vtkpath = joinpath(OUTPUT_DIR, "cheb_laplace.vtk")
    write_cheb_vtk(vtkpath, x, y, U, Uex, Err)
    println(string(" # Chebyshev solution written to ", vtkpath))
    return U
end
