# =============================================================================
#  FFT (Fourier spectral) solver for the periodic Laplace / Poisson equation —
#  the Jexpresso driver
# =============================================================================
#
#  Solves, on a PERIODIC rectangular box (2D or 3D),
#
#       -∇²u = f,          u periodic, zero mean,
#
#  with the FFTW-based solver in fft_poisson_core.jl (FFTPoissonSolver /
#  fft_poisson_solve!): forward rfft, division by the eigenvalue of -∇², inverse
#  rfft. The mathematics, the null-space handling and the two discretisations of
#  -∇² are documented there; this file only feeds that core from a case deck and
#  writes the result. It is an ALTERNATIVE to the SEM direct solve
#  (standard_linsolve!) and the element-learning solve
#  (element_learning_linsolve!); the three are selected in problems/drivers.jl.
#
#  Case-deck inputs (user_inputs.jl):
#
#   :fft_laplacian  => "spectral" (default) | "fd2"
#                      spectral : Fourier-exact -∇² (round-off error for a
#                                 band-limited solution)
#                      fd2      : 2nd-order finite-difference -∇²; the FFT then
#                                 returns the solution of the sparse FD system
#                                 (the system an AMG solver is compared against)
#   :fft_plan       => "measure" (default) | "estimate"   FFTW planner effort
#
#  and ONE of two grid sources, selected by :fft_use_mesh:
#
#   • :fft_use_mesh => false (default) — SYNTHETIC grid (2D, mesh-independent):
#       :fft_N, :fft_M          points per direction (any size; :fft_M defaults
#                               to :fft_N; products of 2,3,5,7 are fastest)
#       :fft_Lx, :fft_Ly        domain lengths            (default 2π)
#       :fft_x0, :fft_y0        lower-left corner         (default 0)
#
#   • :fft_use_mesh => true — solve ON the structured periodic GMSH mesh that
#     Jexpresso reads (2D or 3D, by mesh.nsd; requires :nop => 1 so the global
#     nodes are equispaced). Period per axis from :fft_Lx/:fft_Ly[/:fft_Lz]
#     (default: the mesh extent). The solution is scattered back to the mesh
#     nodes and written by the normal Jexpresso output.
#
#  The problem supplies the data through plain functions in its user_source.jl
#  (give the arity matching the mesh dimension):
#       user_fft_rhs(x,y[,z])   the right-hand side f          (REQUIRED)
#       user_fft_exact(x,y[,z]) the exact solution u_ex        (optional; enables
#                               the automatic error verification)
# =============================================================================

# :fft_laplacian / :fft_plan from the deck → the core's Symbol / FFTW flag.
function _fft_inputs(inputs)
    lap = Symbol(get(inputs, :fft_laplacian, "spectral"))
    lap in FFT_POISSON_LAPLACIANS ||
        error(" # fft_linsolve!: :fft_laplacian => \"$lap\"; expected one of $(FFT_POISSON_LAPLACIANS).")
    plan  = lowercase(string(get(inputs, :fft_plan, "measure")))
    flags = plan == "measure"  ? FFTW.MEASURE :
            plan == "estimate" ? FFTW.ESTIMATE :
            error(" # fft_linsolve!: :fft_plan => \"$plan\"; expected \"measure\" or \"estimate\".")
    return lap, flags
end

# Plan, solve and time the solve (planning excluded: it is a one-time setup
# cost, like a factorisation). Reports a removed RHS mean that is more than
# round-off — a periodic problem with ∫f ≠ 0 has no solution, so the solver
# answers the projected problem -∇²u = f - mean(f).
function _fft_solve_timed(F::Array{Float64}, Ls, inputs)
    lap, flags = _fft_inputs(inputs)
    S = FFTPoissonSolver(size(F), Ls; laplacian = lap, flags = flags)
    u = similar(F)
    jx_robust_solve(string("FFT (FFTW, ", lap, ") solve"), () -> fft_poisson_solve!(u, S, F);
                    robust  = get(inputs, :lbenchmark_solve, true),
                    seconds = Float64(get(inputs, :EL_timing_seconds, 2.0)))
    fmean = S.fmean[]
    if abs(fmean) > 1e-10 * max(1.0, maximum(abs, F))
        println(string(" # fft_linsolve!: RHS mean = ", fmean,
                       " ≠ 0; solved the projected problem -∇²u = f - mean(f) ",
                       "(periodic compatibility condition)."))
    end
    return u, lap
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
# The FFT needs an equispaced grid. The global Jexpresso SEM grid is equispaced
# only at :nop => 1 (LGL with 2 nodes/element ⇒ the element corners); at nop>1
# the in-element LGL nodes are non-uniform and an FFT cannot use them.
#
# This maps every mesh node ip to a tensor index on the periodic lattice by
# folding each of its coordinates into the period: φ_d = mod(x_d - min_d, L_d)/L_d.
# Folding makes a "closed" mesh (a node at both x=min and x=min+L) and a
# Jexpresso periodic-merged mesh (only x=min kept) collapse to the SAME lines —
# the seam node lands on line 1 either way. Dimension-agnostic (ND = 2 or 3).
# Returns
#   (dims::NTuple{ND,Int}, lines::NTuple{ND,Vector}, idxof::Vector{NTuple{ND,Int}})
# where lines[d][i]=min_d+(i-1)L_d/dims[d] and idxof[ip] is the node's lattice
# index. Errors if the folded lines are not uniformly spaced (⇒ nop>1 or a
# non-uniform mesh) or the nodes do not fill a tensor grid.
function fft_grid_from_mesh(mesh, Ls::NTuple{ND,Float64}) where {ND}
    npoin  = Int(mesh.npoin)
    coords = ntuple(d -> @view(mesh.coords[d,:]), ND)
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
        N > 1 || error(" # fft_linsolve!: the mesh has a single node line along $(dirs[d]).")
        _fft_assert_uniform(u, N, dirs[d])
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
                  "phase $(u[k]) ≠ $((k-1)/N)). The FFT needs a uniform grid — use "*
                  ":nop => 1 on a uniform structured periodic mesh, or :fft_use_mesh => false.")
        end
    end
end

# The user RHS / exact field at one grid point (arity = grid dimension).
_fft_user_rhs(c::NTuple{2})   = Float64(user_fft_rhs(c[1], c[2]))
_fft_user_rhs(c::NTuple{3})   = Float64(user_fft_rhs(c[1], c[2], c[3]))
_fft_user_exact(c::NTuple{2}) = Float64(user_fft_exact(c[1], c[2]))
_fft_user_exact(c::NTuple{3}) = Float64(user_fft_exact(c[1], c[2], c[3]))

function _fft_sample_rhs(lines::NTuple{ND}) where {ND}
    F = Array{Float64}(undef, map(length, lines))
    @inbounds for I in CartesianIndices(F)
        F[I] = _fft_user_rhs(ntuple(d -> lines[d][I[d]], ND))
    end
    return F
end

# L2 / L∞ error of the grid solution `u` vs the exact field sampled on the same
# lines. `lines` = per-axis coordinate vectors, `Ls` = per-axis periods (2D or
# 3D). The periodic solution is unique only up to a constant, so the constant is
# pinned to the exact field's mean before comparing (u is shifted in place).
function fft_report_grid_error(u, lines, Ls)
    ND  = ndims(u)
    uex = Array{Float64}(undef, size(u))
    @inbounds for I in CartesianIndices(u)
        uex[I] = _fft_user_exact(ntuple(d -> lines[d][I[d]], ND))
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
    jx_record_solve_error(; linf = linf, l2rel = (ref2 > 0 ? l2/ref2 : NaN), npts = length(u))
    return uex, err
end

# ── Driver: FFT solve of the Laplace/Poisson equation ────────────────────────
# Same call signature as standard_linsolve! / element_learning_linsolve! so the
# driver dispatch in problems/drivers.jl can swap solvers transparently.
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
      Synthetic-grid mode (mesh-independent, 2D)
    =====================================================================#
    N  = Int(get(inputs, :fft_N, 64))
    M  = Int(get(inputs, :fft_M, N))
    (N > 0 && M > 0) || error(" # fft_linsolve!: :fft_N=$N, :fft_M=$M must be positive")

    Lx = Float64(get(inputs, :fft_Lx, 2π)); Ly = Float64(get(inputs, :fft_Ly, 2π))
    x0 = Float64(get(inputs, :fft_x0, 0.0)); y0 = Float64(get(inputs, :fft_y0, 0.0))
    x, y = periodic_grid_lines((N, M), (Lx, Ly), (x0, y0))

    println(YELLOW_FG(string(" # Solve -∇²u = f by FFT (FFTW): ",
                             N, "×", M, " synthetic periodic grid ..............")))

    F = _fft_sample_rhs((x, y))
    u, lap = _fft_solve_timed(F, (Lx, Ly), inputs)

    println(YELLOW_FG(string(" # Solve -∇²u = f by FFT (", lap, ") ............................ DONE")))

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

    println(YELLOW_FG(string(" # Solve -∇²u = f by FFT (FFTW) on the ",
                             join(dims, "×"), " periodic mesh grid ..............")))

    # RHS sampled on the canonical grid lines (periodic ⇒ seam value is unique)
    F = _fft_sample_rhs(lines)
    ugrid, lap = _fft_solve_timed(F, Ls, inputs)

    println(YELLOW_FG(string(" # Solve -∇²u = f by FFT (", lap, ") ............................ DONE")))

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
