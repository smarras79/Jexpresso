#=---------------------------------------------------------------------------------
# SymbolicFD.jl
#
# Write a PDE with Julia's unicode characters (a few LaTeX spellings are accepted
# too) and have it discretized and solved automatically. For example
#
#       ∂q/∂t + ∇⋅(\mathbf{u}q) = \mu∇⋅∇(q)
#
# is read, each differential operator is replaced by its own numerical
# discretization (finite differences for now), the right-hand side ∂q/∂t = …
# is assembled by composing those operators, and the solution is integrated in
# time. PNG/on-the-fly output mirrors Jexpresso's problems/CompEuler/sod1d.
#
# This is NOT a string-to-equation lookup and uses NO AI: ∇, ∇⋅ and ∇² are
# *standalone, composable operators* acting on discrete `Field`s, exactly the
# spirit of Gridap.jl's differential operators (here in strong form, by finite
# differences, rather than weak form). The equation is parsed into an expression
# tree whose nodes are operators and fields; evaluating that tree on the current
# state field produces ∂q/∂t. Because the operators are generic, ANY composition
# works -- ∇⋅(u q), ∇⋅(D∇q), u⋅∇q, ∇²q, … -- not just advection-diffusion.
#
# Only dependency is Plots.jl (already in Jexpresso's Project.toml). Designed so
# the operator layer extends to 2D/3D: fields carry a tensor rank and operators
# loop over `nsd` directions through a single directional-derivative primitive;
# only a multi-dimensional mesh + that primitive need to be added.
#---------------------------------------------------------------------------------=#
module SymbolicFD

using Printf
using Krylov               # steady (time-independent) solve: matrix-free GMRES
using LinearOperators      # wrap the discrete operator as a matrix-free LinearOperator
using Plots                # PNG output + on-the-fly window, exactly as in Jexpresso
import LinearAlgebra: ⋅    # overloaded for the symbolic DSL: ∇⋅F = divergence

export solve, parse_equation, FDMesh1D, FDMesh2D, Field, ascii_plot
export ∇, Δ, ∂t, ⋅, @vars      # symbolic equation DSL

#---------------------------------------------------------------------------------
# 1. Normalization: turn LaTeX shortcuts into the unicode the parser expects
#---------------------------------------------------------------------------------
const LATEX2UNICODE = Dict(
    "\\nabla"   => "∇",
    "\\partial" => "∂",
    "\\cdot"    => "⋅",
    "\\Delta"   => "Δ",
    "\\mu"      => "μ",
    "\\nu"      => "ν",
    "\\rho"     => "ρ",
    "\\kappa"   => "κ",
    "\\alpha"   => "α",
    "\\beta"    => "β",
    "\\gamma"   => "γ",
    "\\lambda"  => "λ",
    "\\sigma"   => "σ",
)

"""
    normalize_equation(eqn) -> String

Strip LaTeX decorations (`\\mathbf{u}` → `u`, `\\vec{u}` → `u`, …), translate
LaTeX greek/operator commands to unicode, collapse the contiguous Laplacian
spellings `∇⋅∇`, `∇^2`, `Δ` to the single token `∇²`, and remove whitespace.

Note: only the *contiguous* `∇⋅∇` becomes the compact Laplacian operator. An
explicit composition such as `∇⋅(∇q)` is left as `divergence(gradient(...))`
and discretized as the genuine composition of the two first-difference
operators (a deliberately honest, if wider, stencil).
"""
function normalize_equation(eqn::AbstractString)
    s = String(eqn)
    for cmd in ("mathbf", "vec", "mathrm", "boldsymbol", "mathbb", "mathcal")
        s = replace(s, Regex("\\\\$cmd\\{([^}]*)\\}") => s"\1")
    end
    for (k, v) in LATEX2UNICODE
        s = replace(s, k => v)
    end
    for junk in ("\\left", "\\right", "\\,", "\\;", "\\!", "\\ ", " ", "\t", "\n")
        s = replace(s, junk => "")
    end
    s = replace(s, "∇⋅∇" => "∇²")
    s = replace(s, "∇^2" => "∇²")
    s = replace(s, "Δ"   => "∇²")
    return s
end

#---------------------------------------------------------------------------------
# 2. Discretization backend (user-selectable) and the 1D grid
#
# The numerical method is pluggable. Everything above the directional-derivative
# primitive (the operators ∇/∇⋅/∇², the parser, the RHS/steady assembly, the
# time loop) is method-agnostic; only `deriv1`/`deriv2` change with the backend,
# dispatched on the mesh's `disc::Discretization`:
#
#   :method => :fd   ->  FDMethod   : 2nd-order central finite differences (this file)
#   :method => :sem  ->  SEMMethod  : spectral element, REUSING Jexpresso's basis
#
# The SEM backend does not reimplement the spectral basis: it calls Jexpresso's
# basis_structs_ξ_ω!/build_Interpolation_basis! (src/kernel/bases/basis_structs.jl)
# for the LGL nodes/weights and the dψ differentiation matrix, and applies dψ with
# the same kernel as rhs.jl's _expansion_inviscid! (DSS-averaged to C0). DSGS
# shock capturing (src/kernel/physics/SGS.jl), as in sod1d, remains a follow-up.
#---------------------------------------------------------------------------------
abstract type Discretization end
struct FDMethod  <: Discretization end
Base.@kwdef struct SEMMethod <: Discretization
    nop::Int = 4                       # polynomial order per element (LGL: nop+1 nodes)
end

method_name(::FDMethod)  = "finite differences (2nd-order central)"
method_name(d::SEMMethod) = "spectral element (LGL, nop = $(d.nop); reuses Jexpresso basis)"

function discretization(inputs::Dict)
    m = get(inputs, :method, :fd)
    m === :fd  && return FDMethod()
    m === :sem && return SEMMethod(nop = Int(get(inputs, :nop, 4)))
    error("SymbolicFD: unknown :method `$m`; choose :fd (finite differences) or " *
          ":sem (spectral element).")
end

# Supertype for every mesh. The operator algebra (§3), the parser, the RHS/time
# loop and the output are written against this abstract type and `nsd(m)`, so a
# new dimension only needs its mesh struct + the `deriv1/deriv2` primitives.
abstract type AbstractFDMesh end

"""
    SEMData

Per-mesh spectral-element data for the SEM backend: `nelem` elements of order
`nop` (`ngl = nop+1` LGL nodes each), reference LGL nodes `ξ`, the nodal
differentiation matrix `D` (Jexpresso's `basis.dψ`, indexed `[basis, point]`),
the LGL quadrature weights `ω`, element→global connectivity `conn`, and the
assembled inverse lumped mass `Minv` (`Minv_i = 1 / Σ_{e∋i} ω_i J_e`). These are
exactly the ingredients of Jexpresso's weak-form RHS: `Minv · DSS(ω ⊙ Dᵀ f)`.
"""
struct SEMData
    nelem::Int
    ngl::Int
    ξ::Vector{Float64}
    D::Matrix{Float64}
    ω::Vector{Float64}
    conn::Matrix{Int}
    Minv::Vector{Float64}
end

"""
    FDMesh1D

1D grid container with `npoin` nodes on `[xmin, xmax]` and the chosen
discretization `disc`. For `FDMethod` the nodes are uniform; for `SEMMethod`
they are element-wise LGL nodes and `sem` holds the spectral-element data.
"""
struct FDMesh1D <: AbstractFDMesh
    nsd::Int
    npoin::Int
    xmin::Float64
    xmax::Float64
    Δx::Float64
    x::Vector{Float64}
    periodic::Bool
    disc::Discretization
    sem::Union{Nothing,SEMData}
end

function FDMesh1D(inputs::Dict)
    xmin     = Float64(get(inputs, :xmin, -1.0))
    xmax     = Float64(get(inputs, :xmax,  1.0))
    periodic = Bool(get(inputs, :periodic, true))
    disc     = discretization(inputs)

    disc isa SEMMethod && return build_sem_mesh(inputs, disc, xmin, xmax, periodic)

    npoin = Int(get(inputs, :npoin, get(inputs, :nelx, 100)))
    Δx = periodic ? (xmax - xmin) / npoin : (xmax - xmin) / (npoin - 1)
    x  = [xmin + i * Δx for i in 0:npoin-1]
    return FDMesh1D(1, npoin, xmin, xmax, Δx, x, periodic, disc, nothing)
end

nsd(m::FDMesh1D) = m.nsd

@inline _ip(i, m::FDMesh1D) = i == m.npoin ? (m.periodic ? 1 : m.npoin) : i + 1
@inline _im(i, m::FDMesh1D) = i == 1       ? (m.periodic ? m.npoin : 1) : i - 1

#---------------------------------------------------------------------------------
# 2b. Two-dimensional meshes.
#
# Two flavours, both <: AbstractFDMesh and both implementing only `nsd` and the
# `deriv1/deriv2(f, m, dir)` primitives (the operator algebra is unchanged):
#
#   FDMesh2D   : a structured Cartesian grid (tensor product of two 1D grids).
#                Supports FDMethod (central differences) and SEMMethod
#                (tensor-product LGL, reusing the SAME 1D Jexpresso basis).
#   JexSEMMesh : wraps the spectral-element bundle Jexpresso's `sem_setup`
#                builds from a *gmsh* mesh (connijk, x/y, metrics dξdx…, dψ, ω).
#                This is the path that runs on the EXISTING kopriva grid.
#
# Nodes are numbered ip = (jy-1)*nx + ix (x fastest) so a structured field
# reshapes to an nx×ny image for heat-map output.
#---------------------------------------------------------------------------------
"per-mesh tensor-product spectral-element data for a structured `FDMesh2D`."
struct SEMData2D
    nelx::Int
    nely::Int
    ngl::Int
    ξ::Vector{Float64}
    D::Matrix{Float64}
    ω::Vector{Float64}
    conn::Array{Int,3}        # conn[e,i,j] -> global node; e = (ey-1)*nelx + ex
    Jx::Float64
    Jy::Float64
    Minv::Vector{Float64}     # 1 / DSS(ω_i ω_j Jx Jy)
end

"""
    FDMesh2D

Structured Cartesian 2D mesh on `[xmin,xmax]×[ymin,ymax]`. `nx`,`ny` are the
unique node counts per direction (`npoin = nx*ny`), `xu`,`yu` the 1D coordinate
lines, and `x`,`y` the per-node coordinates (`ip = (jy-1)*nx + ix`). `dx`,`dy`
are the FD spacings; `Δx = min(dx,dy)` feeds the CFL estimate.
"""
struct FDMesh2D <: AbstractFDMesh
    nsd::Int
    npoin::Int
    nx::Int
    ny::Int
    xmin::Float64; xmax::Float64
    ymin::Float64; ymax::Float64
    dx::Float64;   dy::Float64
    Δx::Float64
    x::Vector{Float64}
    y::Vector{Float64}
    xu::Vector{Float64}
    yu::Vector{Float64}
    periodic::Bool
    disc::Discretization
    sem::Union{Nothing,SEMData2D}
end

"""
    JexSEMMesh

Spectral-element mesh read from a *gmsh* file through Jexpresso's `sem_setup`.
Holds references to the pieces of the returned `sem` bundle that the weak-form
operator needs: `connijk[e,i,j]`, node coords `x`,`y`, the `dψ` matrix, weights
`ω`, the metric terms `dξdx,dξdy,dηdx,dηdy` and Jacobian `Je`, and the assembled
inverse lumped mass `Minv = 1/DSS(ω_iω_j Je)`. `nx`,`ny`,`xu`,`yu` are filled
when the grid happens to be a structured quad mesh (used only for plotting).
"""
struct JexSEMMesh <: AbstractFDMesh
    nsd::Int
    npoin::Int
    nelem::Int
    ngl::Int
    x::Vector{Float64}
    y::Vector{Float64}
    Δx::Float64
    periodic::Bool
    disc::Discretization
    connijk::Array{Int,3}            # [iel, i, j] -> global node
    D::Matrix{Float64}               # dψ[k,i]
    ω::Vector{Float64}
    Je::Array{Float64,3}             # [iel,i,j]
    dξdx::Array{Float64,3}; dξdy::Array{Float64,3}
    dηdx::Array{Float64,3}; dηdy::Array{Float64,3}
    Minv::Vector{Float64}
    nx::Int; ny::Int                 # structured plotting grid (0 if unstructured)
    xu::Vector{Float64}; yu::Vector{Float64}
end

nsd(::FDMesh2D)   = 2
nsd(::JexSEMMesh) = 2

# directional-derivative forwarders: the operator algebra calls these on any mesh
deriv1(f::Vector{Float64}, m::AbstractFDMesh, dir::Int) = deriv1(m.disc, f, m, dir)
deriv2(f::Vector{Float64}, m::AbstractFDMesh, dir::Int) = deriv2(m.disc, f, m, dir)

#---------------------------------------------------------------------------------
# SEM backend: REUSE Jexpresso's spectral basis (no reimplementation).
#
# We do NOT recompute LGL nodes or the differentiation matrix: they come straight
# from Jexpresso's src/kernel/bases/basis_structs.jl —
#     basis_structs_ξ_ω!(LGL(), nop, CPU())            -> LGL nodes ξ and weights ω
#     build_Interpolation_basis!(LagrangeBasis(),ξ,ξ,…) -> the dψ differentiation matrix
# The package is loaded lazily (only when :method => :sem) so the FD path keeps
# its light dependency footprint.
#---------------------------------------------------------------------------------
const _JEX_UUID = Base.UUID("b6c1962e-2b89-11f0-398d-4f31b209be71")   # Jexpresso
const _KA_UUID  = Base.UUID("63c18a36-062a-441e-b654-da1e3ab1ce7c")   # KernelAbstractions
const _JEX = Ref{Module}()
const _KA  = Ref{Module}()
_jexpresso() = isassigned(_JEX) ? _JEX[] : (_JEX[] = Base.require(Base.PkgId(_JEX_UUID, "Jexpresso")))
_kabstr()    = isassigned(_KA)  ? _KA[]  : (_KA[]  = Base.require(Base.PkgId(_KA_UUID,  "KernelAbstractions")))

# The actual basis build. Calls methods from packages that were just loaded with
# `Base.require`, so it must run in the latest world age (see `sem_basis`).
function _sem_basis_impl(jex::Module, ka::Module, nop::Int)
    # Jexpresso's basis/IO helpers query the global communicator, so MPI must be
    # initialised first — normally src/run.jl does this, but the SEM path loads
    # the basis directly. Singleton (no-arg) init is fine for this serial use.
    if !jex.MPI.Initialized()
        jex.MPI.Init()
    end
    cpu = ka.CPU()
    lgl = jex.basis_structs_ξ_ω!(jex.LGL(), nop, cpu)         # St_lgl: ξ, ω
    ξ   = Float64.(collect(lgl.ξ))
    ω   = Float64.(collect(lgl.ω))
    basis = jex.build_Interpolation_basis!(jex.LagrangeBasis(), ξ, ξ, Float64, cpu)
    D   = Matrix{Float64}(basis.dψ)                           # dψ[k,i] = dψ_k(ξ_i)
    return ξ, ω, D
end

"reuse Jexpresso's LGL nodes/weights and the dψ differentiation matrix."
function sem_basis(nop::Int)
    jex = _jexpresso()
    ka  = _kabstr()
    # `_jexpresso()/_kabstr()` may have just loaded the packages, advancing the
    # world age past this method; `invokelatest` runs the build in the new world.
    return Base.invokelatest(_sem_basis_impl, jex, ka, nop)::Tuple{Vector{Float64},Vector{Float64},Matrix{Float64}}
end

function build_sem_mesh(inputs::Dict, disc::SEMMethod, xmin, xmax, periodic)
    nop = disc.nop; ngl = nop + 1
    ξ, ω, D = sem_basis(nop)
    nelem = haskey(inputs, :nelx) ? Int(inputs[:nelx]) :
            max(1, round(Int, (Int(get(inputs, :npoin, nop + 1)) - (periodic ? 0 : 1)) / nop))
    h  = (xmax - xmin) / nelem
    Je = h / 2                                              # affine metric dx/dξ
    npoin = periodic ? nelem * nop : nelem * nop + 1
    conn  = zeros(Int, nelem, ngl)
    x     = zeros(npoin)
    for e in 1:nelem
        xL = xmin + (e - 1) * h
        for i in 1:ngl
            gf  = (e - 1) * nop + i
            gid = periodic ? mod1(gf, npoin) : gf
            conn[e, i] = gid
            (periodic && gf == npoin + 1) || (x[gid] = xL + (ξ[i] + 1) / 2 * h)
        end
    end
    M = zeros(npoin)                                        # assembled lumped mass
    for e in 1:nelem, i in 1:ngl
        M[conn[e, i]] += ω[i] * Je
    end
    Minv = 1.0 ./ M
    Δx = minimum(diff(sort(x)))                             # min node spacing (CFL)
    sem = SEMData(nelem, ngl, ξ, D, ω, conn, Minv)
    return FDMesh1D(1, npoin, xmin, xmax, Δx, x, periodic, disc, sem)
end

#---------------------------------------------------------------------------------
# 2D mesh construction
#---------------------------------------------------------------------------------
# number of nodes / spacing for one structured direction
function _axis_1d(np::Int, lo::Float64, hi::Float64, periodic::Bool)
    d  = periodic ? (hi - lo) / np : (hi - lo) / (np - 1)
    xu = [lo + i * d for i in 0:np-1]
    return d, xu
end

function FDMesh2D(inputs::Dict)
    disc     = discretization(inputs)
    periodic = Bool(get(inputs, :periodic, true))
    xmin = Float64(get(inputs, :xmin, -1.0)); xmax = Float64(get(inputs, :xmax, 1.0))
    ymin = Float64(get(inputs, :ymin, -1.0)); ymax = Float64(get(inputs, :ymax, 1.0))

    disc isa SEMMethod && return build_sem_mesh_2d(inputs, disc, xmin, xmax, ymin, ymax, periodic)

    nx = Int(get(inputs, :npoinx, get(inputs, :nelx, 64)))
    ny = Int(get(inputs, :npoiny, get(inputs, :nely, 64)))
    dx, xu = _axis_1d(nx, xmin, xmax, periodic)
    dy, yu = _axis_1d(ny, ymin, ymax, periodic)
    npoin  = nx * ny
    x = zeros(npoin); y = zeros(npoin)
    for jy in 1:ny, ix in 1:nx
        ip = (jy - 1) * nx + ix
        x[ip] = xu[ix]; y[ip] = yu[jy]
    end
    return FDMesh2D(2, npoin, nx, ny, xmin, xmax, ymin, ymax, dx, dy,
                    min(dx, dy), x, y, xu, yu, periodic, disc, nothing)
end

# structured tensor-product SEM: reuse the SAME 1D LGL basis in x and y
function build_sem_mesh_2d(inputs, disc::SEMMethod, xmin, xmax, ymin, ymax, periodic)
    nop = disc.nop; ngl = nop + 1
    ξ, ω, D = sem_basis(nop)
    nelx = Int(get(inputs, :nelx, 8)); nely = Int(get(inputs, :nely, 8))
    hx = (xmax - xmin) / nelx; hy = (ymax - ymin) / nely
    Jx = hx / 2; Jy = hy / 2
    npx = periodic ? nelx * nop : nelx * nop + 1
    npy = periodic ? nely * nop : nely * nop + 1
    npoin = npx * npy
    gx(ex, i) = periodic ? mod1((ex - 1) * nop + i, npx) : (ex - 1) * nop + i
    gy(ey, j) = periodic ? mod1((ey - 1) * nop + j, npy) : (ey - 1) * nop + j
    # 1D coordinate lines
    xu = zeros(npx); yu = zeros(npy)
    for ex in 1:nelx, i in 1:ngl
        (periodic && (ex - 1) * nop + i == npx + 1) ||
            (xu[gx(ex, i)] = xmin + (ex - 1) * hx + (ξ[i] + 1) / 2 * hx)
    end
    for ey in 1:nely, j in 1:ngl
        (periodic && (ey - 1) * nop + j == npy + 1) ||
            (yu[gy(ey, j)] = ymin + (ey - 1) * hy + (ξ[j] + 1) / 2 * hy)
    end
    conn = zeros(Int, nelx * nely, ngl, ngl)
    x = zeros(npoin); y = zeros(npoin)
    for ey in 1:nely, ex in 1:nelx
        e = (ey - 1) * nelx + ex
        for j in 1:ngl, i in 1:ngl
            ix = gx(ex, i); iy = gy(ey, j); ip = (iy - 1) * npx + ix
            conn[e, i, j] = ip; x[ip] = xu[ix]; y[ip] = yu[iy]
        end
    end
    M = zeros(npoin)
    for e in 1:nelx*nely, j in 1:ngl, i in 1:ngl
        M[conn[e, i, j]] += ω[i] * ω[j] * Jx * Jy
    end
    Minv = 1.0 ./ M
    Δx = min(minimum(diff(sort(xu))), minimum(diff(sort(yu))))
    sem = SEMData2D(nelx, nely, ngl, ξ, D, ω, conn, Jx, Jy, Minv)
    return FDMesh2D(2, npoin, npx, npy, xmin, xmax, ymin, ymax, hx, hy,
                    Δx, x, y, xu, yu, periodic, disc, sem)
end

#---------------------------------------------------------------------------------
# Jexpresso-mesh SEM (gmsh): read the EXISTING grid through `sem_setup`.
#
# We do not reimplement mesh reading or metric terms: Jexpresso's `sem_setup`
# (src/kernel/infrastructure/sem_setup.jl) reads the *.msh file named in the
# inputs, builds the high-order connectivity `connijk`, the metric terms
# (dξdx,…,Je) and the basis dψ/ω. We pull those arrays out of the returned `sem`
# bundle and apply the same weak-form operator as rhs.jl's _expansion_inviscid!.
#---------------------------------------------------------------------------------
const _PA_UUID = Base.UUID("5a9dfac6-5c52-46f7-8278-5e2210713be9")   # PartitionedArrays
const _PA = Ref{Module}()
_partarrays() = isassigned(_PA) ? _PA[] : (_PA[] = Base.require(Base.PkgId(_PA_UUID, "PartitionedArrays")))

# Build the sem bundle by replaying Jexpresso's own preprocessing: fill input
# defaults (mod_inputs_user_inputs!), then run sem_setup inside a with_mpi block
# (PartitionedArrays creates the serial `distribute`). Runs in the latest world.
function _jex_sem_setup_impl(jex::Module, ka::Module, pa::Module, inputs::Dict)
    if !jex.MPI.Initialized()
        jex.MPI.Init()
    end
    rank = jex.MPI.Comm_rank(jex.MPI.COMM_WORLD)
    jin = copy(inputs)
    get!(jin, :backend, ka.CPU())
    jex.mod_inputs_user_inputs!(jin, rank)        # fill all the defaults sem_setup needs
    semref = Ref{Any}(nothing)
    pa.with_mpi(; comm = jex.MPI.COMM_WORLD) do distribute
        semref[] = jex.sem_setup(jin, 1, distribute)   # (sem_namedtuple, partitioned_model[, …])
    end
    return semref[]
end

function build_jex_sem_mesh(inputs::Dict, disc::SEMMethod)
    jex = _jexpresso(); ka = _kabstr(); pa = _partarrays()
    res  = Base.invokelatest(_jex_sem_setup_impl, jex, ka, pa, inputs)
    sem  = res isa Tuple ? res[1] : res
    return Base.invokelatest(_assemble_jex_mesh, sem, disc, inputs)
end

# Pull the arrays we need out of the `sem` bundle and assemble a JexSEMMesh.
function _assemble_jex_mesh(sem, disc::SEMMethod, inputs::Dict)
    mesh    = sem.mesh
    metrics = sem.metrics
    npoin   = Int(mesh.npoin)
    nelem   = Int(mesh.nelem)
    ngl     = Int(mesh.ngl)
    x = Float64.(collect(mesh.x)); y = Float64.(collect(mesh.y))
    connijk = Int.(collect(mesh.connijk))[1:nelem, 1:ngl, 1:ngl]
    D  = Matrix{Float64}(collect(sem.basis.dψ))
    ω  = Float64.(collect(sem.ω))
    Je   = Float64.(collect(metrics.Je))
    dξdx = Float64.(collect(metrics.dξdx)); dξdy = Float64.(collect(metrics.dξdy))
    dηdx = Float64.(collect(metrics.dηdx)); dηdy = Float64.(collect(metrics.dηdy))
    # assembled inverse lumped mass  Minv_ip = 1 / Σ_{e∋ip} ω_i ω_j Je_ij
    M = zeros(npoin)
    @inbounds for e in 1:nelem, j in 1:ngl, i in 1:ngl
        M[connijk[e, i, j]] += ω[i] * ω[j] * Je[e, i, j]
    end
    Minv = 1.0 ./ M
    Δx = _min_node_distance(connijk, x, y, nelem, ngl)
    periodic = Bool(get(inputs, :periodic, true))
    nx, ny, xu, yu = _structured_axes(x, y)        # for plotting, if the grid is a quad lattice
    return JexSEMMesh(2, npoin, nelem, ngl, x, y, Δx, periodic, disc,
                      connijk, D, ω, Je, dξdx, dξdy, dηdx, dηdy, Minv, nx, ny, xu, yu)
end

# smallest distance between the first two nodes of an element edge (CFL scale)
function _min_node_distance(connijk, x, y, nelem, ngl)
    dmin = Inf
    @inbounds for e in 1:nelem
        a = connijk[e, 1, 1]; b = connijk[e, 2, 1]; c = connijk[e, 1, 2]
        dmin = min(dmin, hypot(x[a]-x[b], y[a]-y[b]), hypot(x[a]-x[c], y[a]-y[c]))
    end
    return isfinite(dmin) ? dmin : 1.0
end

# recover the (nx,ny) lattice + axes if the node cloud is a structured grid,
# else return zeros (heat-map output then falls back to a scatter).
function _structured_axes(x, y)
    xu = sort(unique(round.(x; digits = 10)))
    yu = sort(unique(round.(y; digits = 10)))
    (length(xu) * length(yu) == length(x)) || return 0, 0, Float64[], Float64[]
    return length(xu), length(yu), xu, yu
end

#---------------------------------------------------------------------------------
# Directional-derivative primitives — the ONLY method-specific numerics.
# The operators call deriv1/deriv2(f, m, dir) (forwarders above route to the
# backend, dispatching on the mesh's `disc` and the mesh type).
#---------------------------------------------------------------------------------
# --- spectral-element backend (SEMMethod): weak-form derivative, as in rhs.jl ----
# Galerkin weak first derivative with lumped LGL mass, exactly the structure of
# rhs.jl's _expansion_inviscid! + DSS_rhs! + Minv:
#     (∂f/∂x)_i = Minv_i · DSS( ω_i · Σ_k dψ[k,i] f_k )
# (the affine metric dξ/dx cancels between the ω-weighted weak derivative and the
# lumped mass, so it does not appear explicitly). Σ_k dψ[k,i] f_k is the same
# dF/dξ|_i kernel rhs.jl uses. deriv2 composes two weak first derivatives.
function deriv1(::SEMMethod, f::Vector{Float64}, m::FDMesh1D, dir::Int)
    dir == 1 || error("SymbolicFD: direction $dir is unavailable on a 1D mesh.")
    s = m.sem::SEMData
    out = zeros(m.npoin)
    @inbounds for e in 1:s.nelem
        for i in 1:s.ngl
            dFdξ = 0.0
            for k in 1:s.ngl
                dFdξ += s.D[k, i] * f[s.conn[e, k]]      # dF/dξ at node i (rhs.jl kernel)
            end
            out[s.conn[e, i]] += s.ω[i] * dFdξ            # ω-weighted weak derivative, DSS-summed
        end
    end
    @inbounds for ip in 1:m.npoin
        out[ip] *= s.Minv[ip]                             # inverse lumped mass
    end
    return out
end

deriv2(d::SEMMethod, f::Vector{Float64}, m::FDMesh1D, dir::Int) =
    deriv1(d, deriv1(d, f, m, dir), m, dir)

# --- finite-difference backend (FDMethod) --------------------------------------
"first derivative ∂f/∂x_dir at every node (2nd-order central; periodic wrap)."
function deriv1(::FDMethod, f::Vector{Float64}, m::FDMesh1D, dir::Int)
    dir == 1 || error("SymbolicFD: direction $dir is unavailable on a 1D mesh.")
    n, h = m.npoin, m.Δx
    out = similar(f)
    @inbounds for i in 1:n
        out[i] = (f[_ip(i,m)] - f[_im(i,m)]) / (2h)
    end
    return out
end

"second derivative ∂²f/∂x_dir² at every node (compact 3-point central)."
function deriv2(::FDMethod, f::Vector{Float64}, m::FDMesh1D, dir::Int)
    dir == 1 || error("SymbolicFD: direction $dir is unavailable on a 1D mesh.")
    n, h = m.npoin, m.Δx
    out = similar(f)
    @inbounds for i in 1:n
        out[i] = (f[_ip(i,m)] - 2f[i] + f[_im(i,m)]) / (h*h)
    end
    return out
end

# --- finite-difference backend on the structured 2D grid -----------------------
# ip = (jy-1)*nx + ix; x-neighbours step ±1 (wrap in ix), y-neighbours step ±nx.
@inline function _nb2d(ip, nx, ny, dir, off, periodic)
    ix = (ip - 1) % nx + 1; jy = (ip - 1) ÷ nx + 1
    if dir == 1
        ix2 = ix + off
        ix2 = ix2 < 1 ? (periodic ? ix2 + nx : 1) : ix2 > nx ? (periodic ? ix2 - nx : nx) : ix2
        return (jy - 1) * nx + ix2
    else
        jy2 = jy + off
        jy2 = jy2 < 1 ? (periodic ? jy2 + ny : 1) : jy2 > ny ? (periodic ? jy2 - ny : ny) : jy2
        return (jy2 - 1) * nx + ix
    end
end

function deriv1(::FDMethod, f::Vector{Float64}, m::FDMesh2D, dir::Int)
    (dir == 1 || dir == 2) || error("SymbolicFD: direction $dir is invalid on a 2D mesh.")
    h = dir == 1 ? m.dx : m.dy
    out = similar(f)
    @inbounds for ip in 1:m.npoin
        out[ip] = (f[_nb2d(ip, m.nx, m.ny, dir,  1, m.periodic)] -
                   f[_nb2d(ip, m.nx, m.ny, dir, -1, m.periodic)]) / (2h)
    end
    return out
end

function deriv2(::FDMethod, f::Vector{Float64}, m::FDMesh2D, dir::Int)
    (dir == 1 || dir == 2) || error("SymbolicFD: direction $dir is invalid on a 2D mesh.")
    h = dir == 1 ? m.dx : m.dy
    out = similar(f)
    @inbounds for ip in 1:m.npoin
        out[ip] = (f[_nb2d(ip, m.nx, m.ny, dir,  1, m.periodic)] - 2f[ip] +
                   f[_nb2d(ip, m.nx, m.ny, dir, -1, m.periodic)]) / (h * h)
    end
    return out
end

# --- tensor-product SEM backend on the structured 2D grid ----------------------
# Weak directional derivative, the 2D analogue of the 1D form: the metric factors
# of the affine tensor-product element fold into the lumped mass, so
#   (∂f/∂x)_ip = Minv_ip · DSS( ω_i ω_j · Jy · Σ_k dψ[k,i] f[k,j] )      (dir = 1)
#   (∂f/∂y)_ip = Minv_ip · DSS( ω_i ω_j · Jx · Σ_k dψ[k,j] f[i,k] )      (dir = 2)
# with Minv = 1/DSS(ω_i ω_j Jx Jy).  divergence = ∂F/∂x + ∂G/∂y by linearity.
function deriv1(::SEMMethod, f::Vector{Float64}, m::FDMesh2D, dir::Int)
    s = m.sem::SEMData2D; ngl = s.ngl
    geo = dir == 1 ? s.Jy : s.Jx          # transverse Jacobian (the other direction)
    out = zeros(m.npoin)
    @inbounds for e in 1:s.nelx*s.nely, j in 1:ngl, i in 1:ngl
        dfdξ = 0.0
        if dir == 1
            for k in 1:ngl; dfdξ += s.D[k, i] * f[s.conn[e, k, j]]; end
        else
            for k in 1:ngl; dfdξ += s.D[k, j] * f[s.conn[e, i, k]]; end
        end
        out[s.conn[e, i, j]] += s.ω[i] * s.ω[j] * geo * dfdξ
    end
    @inbounds for ip in 1:m.npoin; out[ip] *= s.Minv[ip]; end
    return out
end

deriv2(d::SEMMethod, f::Vector{Float64}, m::FDMesh2D, dir::Int) =
    deriv1(d, deriv1(d, f, m, dir), m, dir)

# --- SEM backend on the Jexpresso (gmsh) mesh: weak deriv with metric terms -----
# Exactly rhs.jl's _expansion_inviscid! per-direction (CL, Inexact, ContGal, 2D):
#   ∂f/∂ξ|ij = Σ_k dψ[k,i] f[k,j] ,  ∂f/∂η|ij = Σ_k dψ[k,j] f[i,k]
#   ∂f/∂x    = ∂f/∂ξ·dξdx + ∂f/∂η·dηdx   (dir = 1; dξdy/dηdy for dir = 2)
#   (∂f/∂x)_ip = Minv_ip · DSS( ω_i ω_j Je_ij · ∂f/∂x|ij )
# divergence ∇⋅(F,G) = ∂F/∂x + ∂G/∂y reproduces the fused kernel by linearity.
function deriv1(::SEMMethod, f::Vector{Float64}, m::JexSEMMesh, dir::Int)
    ngl = m.ngl
    mx  = dir == 1 ? m.dξdx : m.dξdy
    mη  = dir == 1 ? m.dηdx : m.dηdy
    out = zeros(m.npoin)
    @inbounds for e in 1:m.nelem, j in 1:ngl, i in 1:ngl
        dfdξ = 0.0; dfdη = 0.0
        for k in 1:ngl
            dfdξ += m.D[k, i] * f[m.connijk[e, k, j]]
            dfdη += m.D[k, j] * f[m.connijk[e, i, k]]
        end
        dfdx = dfdξ * mx[e, i, j] + dfdη * mη[e, i, j]
        out[m.connijk[e, i, j]] += m.ω[i] * m.ω[j] * m.Je[e, i, j] * dfdx
    end
    @inbounds for ip in 1:m.npoin; out[ip] *= m.Minv[ip]; end
    return out
end

deriv2(d::SEMMethod, f::Vector{Float64}, m::JexSEMMesh, dir::Int) =
    deriv1(d, deriv1(d, f, m, dir), m, dir)

#---------------------------------------------------------------------------------
# 3. Discrete fields and the differential operators that act on them
#
# A `Field` of tensor `rank` r stores nsd^r nodal-value components. Operators
# raise/lower the rank and are written generically over `nsd`, so they already
# read as the textbook vector-calculus identities:
#
#     ∇f        (gradient)   : rank r -> r+1
#     ∇⋅F       (divergence) : rank r -> r-1
#     ∇²f = Δf  (Laplacian)  : rank r -> r   (Σ_d ∂²/∂x_d², compact stencil)
#
# This is the "each symbol has its own numerical counterpart" layer: the
# discretization lives here, in the operator, independent of any equation.
#---------------------------------------------------------------------------------
"""
    Field(mesh, rank, comp)

Discrete tensor field. `comp` holds `nsd^rank` component vectors (each of length
`npoin`); component `c-1` decodes to a base-`nsd` multi-index of length `rank`.
"""
struct Field
    mesh::AbstractFDMesh
    rank::Int
    comp::Vector{Vector{Float64}}
end

scalar_field(m::AbstractFDMesh, v::Vector{Float64}) = Field(m, 0, [v])
const_field(m::AbstractFDMesh, c::Real)             = Field(m, 0, [fill(Float64(c), m.npoin)])
function vec_const_field(m::AbstractFDMesh, v::AbstractVector)
    length(v) == nsd(m) ||
        error("SymbolicFD: a vector constant needs nsd = $(nsd(m)) components, got $(length(v)).")
    return Field(m, 1, [fill(Float64(v[d]), m.npoin) for d in 1:nsd(m)])
end

# --- gradient: ∇ ---------------------------------------------------------------
function gradient(f::Field)
    m, d, nc = f.mesh, nsd(f.mesh), length(f.comp)
    out = Vector{Vector{Float64}}(undef, d * nc)
    for dir in 1:d, c in 1:nc           # new (leading) index = derivative direction
        out[(dir-1)*nc + c] = deriv1(f.comp[c], m, dir)
    end
    return Field(m, f.rank + 1, out)
end

# --- divergence: ∇⋅  (contracts the leading index with its derivative dir) -----
function divergence(f::Field)
    m, d = f.mesh, nsd(f.mesh)
    if f.rank == 0
        # 1D convenience: ∇⋅ of a scalar flux is just ∂/∂x (the single direction).
        d == 1 || error("SymbolicFD: ∇⋅ needs a vector (rank ≥ 1) field when nsd > 1.")
        return Field(m, 0, [deriv1(f.comp[1], m, 1)])
    end
    nc_out = length(f.comp) ÷ d
    out = [zeros(m.npoin) for _ in 1:nc_out]
    for c in 1:nc_out, dir in 1:d
        out[c] .+= deriv1(f.comp[(dir-1)*nc_out + c], m, dir)
    end
    return Field(m, f.rank - 1, out)
end

# --- Laplacian: ∇² = Δ  (compact, per component, summed over directions) -------
function laplacian(f::Field)
    m, d = f.mesh, nsd(f.mesh)
    out = [zeros(m.npoin) for _ in 1:length(f.comp)]
    for c in 1:length(f.comp), dir in 1:d
        out[c] .+= deriv2(f.comp[c], m, dir)
    end
    return Field(m, f.rank, out)
end

# --- field algebra -------------------------------------------------------------
"scalar × field (or field × scalar): scales every component."
function fmul(a::Field, b::Field)
    if a.rank == 0
        return Field(b.mesh, b.rank, [a.comp[1] .* bc for bc in b.comp])
    elseif b.rank == 0
        return Field(a.mesh, a.rank, [b.comp[1] .* ac for ac in a.comp])
    end
    error("SymbolicFD: `*` is scalar×field; use `⋅` for a vector·vector contraction.")
end

"inner product `⋅` : equal-rank fields contracted to a scalar (rank 0 = product)."
function fdot(a::Field, b::Field)
    a.rank == b.rank || error("SymbolicFD: `⋅` needs equal-rank operands.")
    a.rank == 0 && return Field(a.mesh, 0, [a.comp[1] .* b.comp[1]])
    out = zeros(a.mesh.npoin)
    for c in 1:length(a.comp)
        out .+= a.comp[c] .* b.comp[c]
    end
    return Field(a.mesh, 0, [out])
end

function fadd(a::Field, b::Field, s::Float64)
    a.rank == b.rank ||
        error("SymbolicFD: cannot add a rank-$(a.rank) and a rank-$(b.rank) field " *
              "(the equation is dimensionally inconsistent).")
    return Field(a.mesh, a.rank, [a.comp[c] .+ s .* b.comp[c] for c in 1:length(a.comp)])
end

fneg(a::Field) = Field(a.mesh, a.rank, [-ac for ac in a.comp])

#---------------------------------------------------------------------------------
# 4. The expression tree (AST) and the recursive-descent parser
#
# Grammar (on the normalized, whitespace-free string; identifiers are single
# characters so juxtaposition `uq` means u*q, as in handwritten maths):
#
#   expr   := term (('+'|'-') term)*
#   term   := factor ( ('*'|'⋅'|<juxtaposition>) factor )*
#   factor := '∇⋅' factor | '∇²' factor | '∇' factor | '∂' id '/∂t' | atom
#   atom   := '(' expr ')' | number | id
#---------------------------------------------------------------------------------
struct Node
    op::Symbol           # :const :func :vecconst :var :grad :div :lap :dt :mul :dot :add :sub :neg
    val::Float64
    vec::Vector{Float64}
    name::String
    fun::Union{Nothing,Function}   # for :func -- a spatially varying parameter/source x->value
    args::Vector{Node}
end
Node(op::Symbol; val = 0.0, vec = Float64[], name = "", fun = nothing, args = Node[]) =
    Node(op, val, vec, name, fun, args)

mutable struct PState
    chars::Vector{Char}
    pos::Int
    var::String
    inputs::Dict
end

_peek(p::PState)  = p.pos <= length(p.chars) ? p.chars[p.pos]     : '\0'
_peek2(p::PState) = p.pos+1 <= length(p.chars) ? p.chars[p.pos+1] : '\0'
_adv!(p::PState)  = (c = _peek(p); p.pos += 1; c)
_expect!(p::PState, c::Char) = _adv!(p) == c ||
    error("SymbolicFD: expected `$c` near position $(p.pos) in the equation.")

_isidchar(c::Char) = isletter(c) || c == '_'
function _is_factor_start(c::Char)
    c == '\0' && return false
    c in ('+', '-', '*', '⋅', ')', '=', '/', '²') && return false
    return true
end

function _read_number(p::PState)
    a = p.pos
    while isdigit(_peek(p)) || _peek(p) == '.'; _adv!(p); end
    if _peek(p) == 'e' || _peek(p) == 'E'
        _adv!(p)
        (_peek(p) == '+' || _peek(p) == '-') && _adv!(p)
        while isdigit(_peek(p)); _adv!(p); end
    end
    return parse(Float64, String(p.chars[a:p.pos-1]))
end

# resolve a single-character identifier to a leaf node
function _leaf(p::PState, name::String)
    if name == p.var
        return Node(:var; name = name)
    elseif haskey(p.inputs, Symbol(name))
        v = p.inputs[Symbol(name)]
        v isa Function       && return Node(:func;     name = name, fun = v)
        v isa AbstractVector && return Node(:vecconst; vec = Float64.(v))
        return Node(:const; val = Float64(v))
    end
    error("SymbolicFD: unknown symbol `$name` in the equation. " *
          "Add `:$name => value` to your inputs Dict (a vector for a vector quantity, " *
          "a function `x->…` for a spatially varying coefficient/source).")
end

function parse_atom(p::PState)
    c = _peek(p)
    if c == '('
        _adv!(p)
        node = parse_expr(p)
        _expect!(p, ')')
        return node
    elseif isdigit(c) || c == '.'
        return Node(:const; val = _read_number(p))
    elseif _isidchar(c)
        return _leaf(p, string(_adv!(p)))     # single-character identifier
    end
    error("SymbolicFD: unexpected character `$c` near position $(p.pos).")
end

function parse_timederiv(p::PState)
    _expect!(p, '∂')
    _isidchar(_peek(p)) || error("SymbolicFD: expected a variable name after ∂.")
    name = string(_adv!(p))
    _expect!(p, '/'); _expect!(p, '∂'); _expect!(p, 't')
    return Node(:dt; name = name)
end

function parse_factor(p::PState)
    c = _peek(p)
    if c == '∇'
        _adv!(p)
        n = _peek(p)
        if n == '⋅'
            _adv!(p); return Node(:div; args = [parse_factor(p)])
        elseif n == '²'
            _adv!(p); return Node(:lap; args = [parse_factor(p)])
        else
            return Node(:grad; args = [parse_factor(p)])
        end
    elseif c == '∂'
        return parse_timederiv(p)
    else
        return parse_atom(p)
    end
end

function parse_term(p::PState)
    node = parse_factor(p)
    while true
        c = _peek(p)
        if c == '*'
            _adv!(p); node = Node(:mul; args = [node, parse_factor(p)])
        elseif c == '⋅'
            _adv!(p); node = Node(:dot; args = [node, parse_factor(p)])
        elseif _is_factor_start(c)
            node = Node(:mul; args = [node, parse_factor(p)])   # juxtaposition
        else
            break
        end
    end
    return node
end

function parse_expr(p::PState)
    node = parse_term(p)
    while _peek(p) == '+' || _peek(p) == '-'
        opc = _adv!(p)
        node = Node(opc == '+' ? :add : :sub; args = [node, parse_term(p)])
    end
    return node
end

parse_side(s::AbstractString, var, inputs) =
    parse_expr(PState(collect(s), 1, var, inputs))

# flatten an expression into signed additive leaves
function flatten!(acc::Vector{Tuple{Float64,Node}}, n::Node, s::Float64)
    if n.op == :add
        flatten!(acc, n.args[1], s); flatten!(acc, n.args[2], s)
    elseif n.op == :sub
        flatten!(acc, n.args[1], s); flatten!(acc, n.args[2], -s)
    elseif n.op == :neg
        flatten!(acc, n.args[1], -s)
    else
        push!(acc, (s, n))
    end
    return acc
end

# Infer the unknown of a *steady* (no ∂/∂t) equation: the single-character
# symbol that is not a provided parameter. Overridable with `:unknown => "q"`.
function find_unknown(s::AbstractString, inputs::Dict)
    haskey(inputs, :unknown) && return String(inputs[:unknown])
    cands = String[]
    for c in s
        if isletter(c)
            str = string(c)
            (!haskey(inputs, Symbol(str)) && str ∉ cands) && push!(cands, str)
        end
    end
    length(cands) == 1 ||
        error("SymbolicFD: could not infer the unknown (candidates: $cands). " *
              "Provide `:unknown => \"q\"`, or make sure every parameter is in the inputs Dict.")
    return cands[1]
end

"""
    parse_equation(eqn, inputs) -> (var::String, mode::Symbol, node::Node)

Parse a unicode/LaTeX PDE into the unknown variable name, a `mode`, and an
expression tree, with NO operator tied to a specific equation:

- if a `∂var/∂t` term is present, `mode = :transient` and `node` is the
  right-hand side of `∂var/∂t = node` (the time term is removed and every other
  term moved across `=` with the proper sign);
- otherwise `mode = :steady` and `node` is the residual `lhs - rhs`, to be
  driven to zero by a direct solve.
"""
function parse_equation(eqn::AbstractString, inputs::Dict)
    s = normalize_equation(eqn)
    occursin("=", s) ||
        error("SymbolicFD: the equation must contain `=` (e.g. `∂q/∂t + ∇⋅(uq) = μ∇²q`).")
    lhs, rhs = split(s, "=", limit = 2)
    tm = match(r"∂([^/]+)/∂t", s)

    if tm === nothing
        # ---- steady / time-independent (e.g. the steady heat equation) ----
        var = find_unknown(s, inputs)
        length(var) == 1 ||
            error("SymbolicFD: use a single-character unknown name (got `$var`).")
        resid = Node(:sub; args = [parse_side(lhs, var, inputs),
                                   parse_side(rhs, var, inputs)])
        return var, :steady, resid
    end

    # ---- transient ----
    var = String(tm.captures[1])
    length(var) == 1 ||
        error("SymbolicFD: use single-character variable/parameter names (got `$var`); " *
              "juxtaposition like `uq` is read as u*q.")
    lhs_terms = flatten!(Tuple{Float64,Node}[], parse_side(lhs, var, inputs), 1.0)
    rhs_terms = flatten!(Tuple{Float64,Node}[], parse_side(rhs, var, inputs), 1.0)

    signed = Tuple{Float64,Node}[]          # dq/dt terms (rhs keep sign, lhs flip)
    dt_found = false
    for (sg, nd) in rhs_terms
        nd.op == :dt ? error("SymbolicFD: the ∂q/∂t term must be on the left of `=`.") :
                       push!(signed, (sg, nd))
    end
    for (sg, nd) in lhs_terms
        if nd.op == :dt
            sg == 1.0 || error("SymbolicFD: the ∂q/∂t term must appear with a + sign.")
            nd.name == var || error("SymbolicFD: time-derivative variable mismatch.")
            dt_found = true
        else
            push!(signed, (-sg, nd))
        end
    end
    dt_found || error("SymbolicFD: no `∂$var/∂t` term found.")

    isempty(signed) && return var, :transient, Node(:const; val = 0.0)
    s1, n1 = signed[1]
    dqdt = s1 == 1.0 ? n1 : Node(:neg; args = [n1])
    for (sg, nd) in signed[2:end]
        dqdt = Node(sg == 1.0 ? :add : :sub; args = [dqdt, nd])
    end
    return var, :transient, dqdt
end

#---------------------------------------------------------------------------------
# 4b. Symbolic front-end: write the equation as live Julia symbols, no string.
#
#   @vars q u μ
#   eq = ∂t(q) + ∇⋅(u*q) - μ*∇⋅∇(q)          # residual form,  = 0 implied
#   solve(eq, inputs)
#
# `∇`, `Δ`, `∂t` are real objects and `+ - * ⋅` are overloaded on `Node`, so the
# expression *evaluates* to the same AST the string parser produces. Notes on
# Julia lexing: write `∂t(q)` (the literal `∂q/∂t` lexes as two identifiers
# `∂q`,`∂t`), and put a `*` between a coefficient and an operator (`μ*∇⋅∇(q)`,
# since `μ∇` would lex as one identifier). Parameter symbols are resolved
# (scalar / vector / function of x) from the inputs Dict at solve time; the
# unknown is the symbol carrying `∂t(·)`, or the one symbol absent from inputs.
#---------------------------------------------------------------------------------
struct Nabla     end                       # ∇
struct Laplacian end                       # Δ  (and ∇⋅∇)
struct TimeDeriv end                       # ∂t
struct ScaledNabla; coef::Node; end          # c*∇  (see note on precedence below)

const ∇  = Nabla()
const Δ  = Laplacian()
const ∂t = TimeDeriv()

tonode(x::Node)           = x
tonode(x::Real)           = Node(:const;    val = Float64(x))
tonode(x::AbstractVector) = Node(:vecconst; vec = Float64.(x))
tonode(x) = error("SymbolicFD: cannot use $(typeof(x)) in a symbolic equation; " *
                  "use field/parameter symbols (see @vars), numbers, or vectors.")

# gradient ∇(f) ; divergence ∇⋅(F) ; Laplacian Δ(f) and the composed ∇⋅∇(f)
(::Nabla)(x)     = Node(:grad; args = [tonode(x)])
(::Laplacian)(x) = Node(:lap;  args = [tonode(x)])
⋅(::Nabla, n::Node) = n.op == :grad ? Node(:lap; args = n.args) :   # ∇⋅∇ → compact Δ
                                      Node(:div; args = [n])
⋅(a::Node, b::Node) = Node(:dot; args = [a, b])

# Julia lexes/associates `μ*∇⋅∇(q)` as `(μ*∇)⋅∇(q)` (`*` and `⋅` share precedence,
# left-assoc), so `μ*∇` must build a "scaled ∇" that carries the coefficient into
# the divergence/Laplacian. This is what lets the user write `μ*∇⋅∇(q)` verbatim.
Base.:*(c::Node, ::Nabla) = ScaledNabla(c)
Base.:*(c,       ::Nabla) = ScaledNabla(tonode(c))
⋅(s::ScaledNabla, n::Node) =
    n.op == :grad ? Node(:mul; args = [s.coef, Node(:lap; args = n.args)]) :
                    Node(:mul; args = [s.coef, Node(:div; args = [n])])

# time derivative ∂t(q)
function (::TimeDeriv)(n::Node)
    n.op == :sym || error("SymbolicFD: ∂t(...) expects the unknown field symbol, e.g. ∂t(q).")
    return Node(:dt; name = n.name)
end

# field algebra on the AST
Base.:+(a::Node, b::Node) = Node(:add; args = [a, b])
Base.:+(a::Node, b)       = Node(:add; args = [a, tonode(b)])
Base.:+(a,       b::Node) = Node(:add; args = [tonode(a), b])
Base.:-(a::Node, b::Node) = Node(:sub; args = [a, b])
Base.:-(a::Node, b)       = Node(:sub; args = [a, tonode(b)])
Base.:-(a,       b::Node) = Node(:sub; args = [tonode(a), b])
Base.:-(a::Node)          = Node(:neg; args = [a])
Base.:*(a::Node, b::Node) = Node(:mul; args = [a, b])
Base.:*(a::Node, b)       = Node(:mul; args = [a, tonode(b)])
Base.:*(a,       b::Node) = Node(:mul; args = [tonode(a), b])
Base.:/(a::Node, b::Real) = Node(:mul; args = [a, tonode(1 / b)])

"""
    @vars q u μ ...

Declare each name as a symbolic field/parameter for building an equation, e.g.
`@vars q u μ` then `eq = ∂t(q) + ∇⋅(u*q) - μ*∇⋅∇(q)`.
"""
macro vars(names...)
    blk = Expr(:block)
    for nm in names
        push!(blk.args, :($(esc(nm)) = $(Node)(:sym; name = $(string(nm)))))
    end
    push!(blk.args, :nothing)
    return blk
end

# walk a symbolic tree --------------------------------------------------------
function _find_dt_name(n::Node)
    n.op == :dt && return n.name
    for a in n.args
        r = _find_dt_name(a)
        r != "" && return r
    end
    return ""
end

function _collect_syms!(n::Node, acc::Set{String})
    n.op == :sym && push!(acc, n.name)
    for a in n.args; _collect_syms!(a, acc); end
    return acc
end

# resolve a parameter symbol from the inputs Dict (scalar / vector / function)
function _resolve_param(name::String, inputs::Dict)
    haskey(inputs, Symbol(name)) ||
        error("SymbolicFD: unknown symbol `$name`. Add `:$name => value` to inputs " *
              "(a vector for a vector quantity, a function `x->…` for a field).")
    v = inputs[Symbol(name)]
    v isa Function       && return Node(:func;     name = name, fun = v)
    v isa AbstractVector && return Node(:vecconst; vec  = Float64.(v))
    return Node(:const; val = Float64(v))
end

# replace :sym leaves by :var (the unknown) or resolved parameter nodes
function resolve_syms(n::Node, var::String, inputs::Dict)
    n.op == :sym && return n.name == var ? Node(:var; name = n.name) :
                                           _resolve_param(n.name, inputs)
    isempty(n.args) && return n
    return Node(n.op; val = n.val, vec = n.vec, name = n.name, fun = n.fun,
                args = [resolve_syms(a, var, inputs) for a in n.args])
end

"""
    build_from_residual(residual::Node, inputs) -> (var, mode, node, repr)

Turn a symbolic residual `R` (meaning `R = 0`) into the same `(var, mode, node)`
the string parser yields: transient `∂q/∂t = …` when an `∂t(q)` term is present,
else a steady residual to drive to zero.
"""
function build_from_residual(residual::Node, inputs::Dict)
    dtname = _find_dt_name(residual)
    if dtname != ""
        var = dtname
    else
        cands = collect(setdiff(_collect_syms!(residual, Set{String}()),
                                String[String(k) for k in keys(inputs)]))
        var = haskey(inputs, :unknown) ? String(inputs[:unknown]) :
              (length(cands) == 1 ? cands[1] :
               error("SymbolicFD: could not infer the unknown (candidates: $cands); " *
                     "pass `:unknown => \"q\"`."))
    end
    length(var) == 1 ||
        error("SymbolicFD: use a single-character unknown name (got `$var`).")

    resolved = resolve_syms(residual, var, inputs)
    repr = string(ast_str(resolved), " = 0")

    terms = flatten!(Tuple{Float64,Node}[], resolved, 1.0)
    if dtname == ""
        return var, :steady, resolved, repr
    end
    signed = Tuple{Float64,Node}[]
    for (sg, nd) in terms
        if nd.op == :dt
            sg == 1.0 || error("SymbolicFD: the ∂t(q) term must enter the residual with + sign.")
            nd.name == var || error("SymbolicFD: time-derivative variable mismatch.")
        else
            push!(signed, (-sg, nd))            # move to RHS of ∂q/∂t = …
        end
    end
    isempty(signed) && return var, :transient, Node(:const; val = 0.0), repr
    s1, n1 = signed[1]
    dqdt = s1 == 1.0 ? n1 : Node(:neg; args = [n1])
    for (sg, nd) in signed[2:end]
        dqdt = Node(sg == 1.0 ? :add : :sub; args = [dqdt, nd])
    end
    return var, :transient, dqdt, repr
end

#---------------------------------------------------------------------------------
# 5. Evaluate the expression tree on the current state field -> ∂q/∂t
#---------------------------------------------------------------------------------
# evaluate a spatially-varying coefficient/source: f(x) in 1D, f(x,y) in 2D
_eval_func(fun, m::AbstractFDMesh) =
    nsd(m) == 1 ? Float64[fun(xi) for xi in m.x] :
                  Float64[fun(m.x[ip], m.y[ip]) for ip in 1:m.npoin]

function eval_node(n::Node, qf::Field)
    op = n.op
    op == :const    && return const_field(qf.mesh, n.val)
    op == :func     && return Field(qf.mesh, 0, [_eval_func(n.fun, qf.mesh)])
    op == :vecconst && return vec_const_field(qf.mesh, n.vec)
    op == :var      && return qf
    op == :grad     && return gradient(eval_node(n.args[1], qf))
    op == :div      && return divergence(eval_node(n.args[1], qf))
    op == :lap      && return laplacian(eval_node(n.args[1], qf))
    op == :mul      && return fmul(eval_node(n.args[1], qf), eval_node(n.args[2], qf))
    op == :dot      && return fdot(eval_node(n.args[1], qf), eval_node(n.args[2], qf))
    op == :add      && return fadd(eval_node(n.args[1], qf), eval_node(n.args[2], qf),  1.0)
    op == :sub      && return fadd(eval_node(n.args[1], qf), eval_node(n.args[2], qf), -1.0)
    op == :neg      && return fneg(eval_node(n.args[1], qf))
    error("SymbolicFD: cannot evaluate node $(op).")
end

"pretty-print an AST node back to unicode, for the run banner."
function ast_str(n::Node)
    op = n.op
    op == :const    && return @sprintf("%.4g", n.val)
    op == :func     && return string(n.name, "(x)")
    op == :vecconst && return string("[", join(map(x -> @sprintf("%.4g", x), n.vec), ","), "]")
    op == :var      && return n.name
    op == :sym      && return n.name
    op == :dt       && return string("∂", n.name, "/∂t")
    op == :grad     && return string("∇(", ast_str(n.args[1]), ")")
    op == :div      && return string("∇⋅(", ast_str(n.args[1]), ")")
    op == :lap      && return string("∇²(", ast_str(n.args[1]), ")")
    op == :mul      && return string(ast_str(n.args[1]), "·", ast_str(n.args[2]))
    op == :dot      && return string(ast_str(n.args[1]), "⋅", ast_str(n.args[2]))
    op == :add      && return string(ast_str(n.args[1]), " + ", ast_str(n.args[2]))
    op == :sub      && return string(ast_str(n.args[1]), " - ", ast_str(n.args[2]))
    op == :neg      && return string("-", ast_str(n.args[1]))
    return string(op)
end

#---------------------------------------------------------------------------------
# 6. Right-hand side closure: wrap q in a Field, evaluate the tree, return ∂q/∂t
#---------------------------------------------------------------------------------
function build_rhs(node::Node, mesh::AbstractFDMesh)
    function rhs!(dq, q)
        qf = scalar_field(mesh, q)        # alias q; operators allocate, never mutate q
        r  = eval_node(node, qf)
        r.rank == 0 || error("SymbolicFD: the equation evaluated to a rank-$(r.rank) field; " *
                             "it is not a scalar balance.")
        copyto!(dq, r.comp[1])
        return dq
    end
    return rhs!
end

#---------------------------------------------------------------------------------
# 6b. Steady / time-independent solve (no ∂/∂t): MATRIX-FREE.
#
# The discrete operator is NEVER stored (not even sparse). The same operator
# evaluator is reused: for a linear PDE the residual is affine,
# residual(q) = A q + c, so a matrix-vector product A·x is just
#     A·x = residual(x) - residual(0).
# We wrap that as a `LinearOperator` (with identity rows at the two boundary
# nodes for the Dirichlet conditions) and hand it to Krylov.jl's GMRES. Only
# operator *applications* ever happen -- nothing equation-specific lives here.
#---------------------------------------------------------------------------------
# Build the matrix-free Dirichlet operator L and right-hand side b such that
# L x = b is the discrete steady BVP. Returns (L::LinearOperator, b, resid!).
function steady_operator(resid::Node, mesh::FDMesh1D, inputs::Dict)
    mesh.periodic &&
        error("SymbolicFD: a steady (no ∂/∂t) problem needs Dirichlet boundary data; " *
              "set `:periodic => false` and provide `:bc_left` / `:bc_right`.")
    n = mesh.npoin
    resid! = build_rhs(resid, mesh)

    c = zeros(n);  resid!(c, zeros(n))             # c = residual(0)

    # matrix-free product  y = L x :  interior rows = A x = residual(x) - c,
    # boundary rows = x (identity ⇒ Dirichlet). No matrix is materialized.
    function matvec!(y, x)
        resid!(y, x)                               # y = residual(x) = A x + c
        @inbounds @. y -= c                        # y = A x
        @inbounds y[1] = x[1];  y[n] = x[n]        # Dirichlet identity rows
        return y
    end
    L = LinearOperator(Float64, n, n, false, false, matvec!)

    bcl = Float64(get(inputs, :bc_left,  0.0))
    bcr = Float64(get(inputs, :bc_right, 0.0))
    b = -c;  b[1] = bcl;  b[n] = bcr
    return L, b, resid!
end

function solve_steady(resid::Node, mesh::FDMesh1D, inputs::Dict)
    L, b, resid! = steady_operator(resid, mesh, inputs)
    n = length(b)
    atol = Float64(get(inputs, :ksp_atol, 1e-12))
    rtol = Float64(get(inputs, :ksp_rtol, 1e-10))
    # full GMRES (memory = n) converges in ≤ n iterations even on the
    # ill-conditioned Laplacian; bump :ksp_memory down for a restarted variant.
    mem  = Int(get(inputs, :ksp_memory, n))
    q, stats = Krylov.gmres(L, b; atol = atol, rtol = rtol, memory = mem, itmax = 4n)
    stats.solved || @warn "SymbolicFD: GMRES did not fully converge" stats.status

    # report the interior residual of the original (unconstrained) operator
    r = zeros(n); resid!(r, q)
    rmax = n > 2 ? maximum(abs, @view r[2:n-1]) : 0.0
    return q, rmax
end

#---------------------------------------------------------------------------------
# 7. Explicit RK4 time integration with an automatic stable Δt
#---------------------------------------------------------------------------------
# Auto Δt from the advective (:u) and diffusive (:μ) scales when present; the
# user can always override with :Δt. (RK4 stays stable while dt·max|λ| ≲ 2.8.)
function stable_dt(inputs::Dict, m::AbstractFDMesh, cfl::Float64)
    h, dt = m.Δx, Inf
    if haskey(inputs, :u)
        uv = inputs[:u]
        umax = uv isa AbstractVector ? maximum(abs, uv) : abs(Float64(uv))
        umax > 0 && (dt = min(dt, cfl * h / umax))
    end
    if haskey(inputs, :μ)
        μv = inputs[:μ]
        μmax = μv isa AbstractVector ? maximum(abs, μv) : abs(Float64(μv))
        μmax > 0 && (dt = min(dt, cfl * 0.5 * h * h / μmax))
    end
    isfinite(dt) || (dt = cfl * h)
    return dt
end

struct RK4Work
    k1::Vector{Float64}; k2::Vector{Float64}
    k3::Vector{Float64}; k4::Vector{Float64}; qt::Vector{Float64}
end
RK4Work(n::Int) = RK4Work(zeros(n), zeros(n), zeros(n), zeros(n), zeros(n))

function rk4_step!(q, rhs!, dt, w::RK4Work)
    rhs!(w.k1, q)
    @. w.qt = q + 0.5dt*w.k1 ; rhs!(w.k2, w.qt)
    @. w.qt = q + 0.5dt*w.k2 ; rhs!(w.k3, w.qt)
    @. w.qt = q +     dt*w.k3 ; rhs!(w.k4, w.qt)
    @. q += (dt/6) * (w.k1 + 2w.k2 + 2w.k3 + w.k4)
    return q
end

#---------------------------------------------------------------------------------
# 8. Output: PNG + on-the-fly plotting, exactly as in Jexpresso's
#    src/io/plotting/jeplots.jl (used by problems/CompEuler/sod1d). CSV is always
#    written; a terminal ASCII plot is a display-free fallback (:outformat="ascii").
#---------------------------------------------------------------------------------
function _savefig_silent(plt, fout_name)
    plt[:overwrite_figure] = false
    Plots.savefig(plt, string(fout_name))
    return nothing
end

function plot_initial_png(x, q, var, OUTPUT_DIR::String)
    npoin = length(q)
    plt = Plots.scatter(x[1:npoin], q[1:npoin];
                        markersize = 5, color = :blue,
                        xlabel = "x", ylabel = "$var(x)", title = "$var",
                        titlefontsize = 24, guidefontsize = 18,
                        legendfontsize = 14, tickfontsize = 14,
                        legend = false, size = (800, 600))
    _savefig_silent(plt, string(OUTPUT_DIR, "/INIT-", var, ".png"))
    return plt
end

function plot_results_png(x, q, var, ttl::String, OUTPUT_DIR::String, iout::Int;
                          live::Bool = true, wfig = 600, hfig = 400)
    sort_idx = sortperm(x)
    plt = Plots.plot(x[sort_idx], q[sort_idx];
                     line = (:blue, 2), marker = (:circle, 5, :blue),
                     title = string(var, "  ", ttl), xlabel = "x",
                     titlefontsize = 22, guidefontsize = 18,
                     legendfontsize = 14, tickfontsize = 14,
                     legend = false, show = false, size = (wfig, hfig))
    fout = string(OUTPUT_DIR, "/fields-it", iout, ".png")
    if live
        try; Plots.savefig(plt, fout); catch; end   # GR: repaints live window AND writes file
        try; display(plt);             catch; end
    else
        _savefig_silent(plt, fout)
    end
    return plt
end

# 2D field as a coloured node map (works for any node cloud: structured FD,
# tensor-product SEM, or the unstructured gmsh SEM grid).
function plot_field_png_2d(m::AbstractFDMesh, q, var, ttl::String, OUTPUT_DIR::String,
                           iout::Int; live::Bool = true, init::Bool = false)
    plt = Plots.scatter(m.x, m.y; marker_z = q, markershape = :rect, markersize = 3,
                        markerstrokewidth = 0, c = :viridis, aspect_ratio = :equal,
                        xlabel = "x", ylabel = "y", title = string(var, "  ", ttl),
                        colorbar = true, legend = false, show = false, size = (700, 600))
    fout = init ? string(OUTPUT_DIR, "/INIT-", var, ".png") :
                  string(OUTPUT_DIR, "/fields-it", iout, ".png")
    if live && !init
        try; Plots.savefig(plt, fout); catch; end
        try; display(plt);             catch; end
    else
        _savefig_silent(plt, fout)
    end
    return plt
end

function save_csv(path, m::FDMesh1D, q0, q)
    open(path, "w") do io
        println(io, "x,q_initial,q_final")
        for i in 1:m.npoin
            @printf(io, "%.10e,%.10e,%.10e\n", m.x[i], q0[i], q[i])
        end
    end
    return path
end

function save_csv(path, m::AbstractFDMesh, q0, q)            # 2D: x,y,q_initial,q_final
    open(path, "w") do io
        println(io, "x,y,q_initial,q_final")
        for ip in 1:m.npoin
            @printf(io, "%.10e,%.10e,%.10e,%.10e\n", m.x[ip], m.y[ip], q0[ip], q[ip])
        end
    end
    return path
end

"""
    ascii_plot(x, q; rows=15, cols=70, label="q")

Render a quick line plot of `q(x)` to the terminal (display-free fallback).
"""
function ascii_plot(x, q; rows = 15, cols = 70, label = "q")
    n = length(q)
    qmin, qmax = minimum(q), maximum(q)
    span = qmax - qmin; span == 0 && (span = 1.0)
    grid = fill(' ', rows, cols)
    for c in 1:cols
        idx = round(Int, 1 + (c-1) * (n-1) / (cols-1))
        r   = clamp(rows - round(Int, (q[idx] - qmin) / span * (rows-1)), 1, rows)
        grid[r, c] = '*'
    end
    println()
    @printf("  %s   [min = %.4g, max = %.4g]\n", label, qmin, qmax)
    for r in 1:rows; println("  |", String(grid[r, :])); end
    println("  +", repeat("-", cols))
    @printf("   x = %.3g %s x = %.3g\n", x[1], repeat(" ", cols-12), x[end])
    return nothing
end

#---------------------------------------------------------------------------------
# 9. Top-level driver
#---------------------------------------------------------------------------------
"""
    solve(eqn, inputs) -> (mesh, q0, q)

Solve a PDE given either as a **symbolic expression** built from `@vars` and the
operators `∇`, `Δ`, `∂t`, `⋅` (residual form, e.g. `∂t(q) + ∇⋅(u*q) - μ*∇⋅∇(q)`),
or as a **string** (`"∂q/∂t + ∇⋅(u q) = μ∇²q"`, LaTeX accepted). The equation is
turned into composed finite-difference operators; with a `∂t` term it is
integrated in time (return = initial & final fields), otherwise it is solved as a
steady problem. Returns the mesh and the (initial, final) solution vectors.

Recognised `inputs` keys
- `:method`                   `:fd` finite differences (default) or `:sem`
                              spectral element (Jexpresso; not yet implemented)
- `:nop`                      polynomial order per element when `:method => :sem`
- `:xmin`, `:xmax`            domain (default [-1, 1])
- `:npoin` (or `:nelx`)       number of grid points (default 100)
- `:periodic`                 true/false (default true)
- `:nsd`                      spatial dimension: 1, or 2 (structured FD/SEM, or
                              a gmsh grid read through Jexpresso when
                              `:method => :sem` and `:gmsh_filename` is set)
- `:ymin`, `:ymax`            y-extent of the 2D domain (default [-1, 1])
- `:nelx`, `:nely`            elements (SEM) / points (FD) per direction in 2D
- `:u`, `:μ`, …               parameter symbols in the equation (vector for a
                              vector quantity `:u => [1.0]`; a function `x->…`
                              for a spatially varying coefficient/source)
- `:q0`                       initial condition `x -> value` (transient only)
- `:tend`, `:Δt`, `:CFL`      time-integration controls (transient only)
- `:output_dir`               where to write figures / `solution.csv` (default ".")
- `:outformat`                "png" (Plots.jl, like sod1d) or "ascii" (default "png")
- `:ndiagnostics_outputs`     number of on-the-fly plot snapshots (default 10)
- `:plot_live`                update an on-screen window on the fly (default true)

If the equation has **no `∂q/∂t` term** the problem is solved as a steady
(time-independent) elliptic problem, e.g. the steady heat equation `μ∇²q = f`.
In that case the discrete operator is assembled into a matrix and solved
directly; extra keys:
- `:bc_left`, `:bc_right`     Dirichlet boundary values (default 0.0)
- `:unknown`                  name of the unknown if it cannot be inferred
- (the grid defaults to non-periodic; provide `:periodic => false` explicitly
  if you also set other periodic options)
"""
solve(eqn::AbstractString, inputs::Dict) =      # string input: parse to (var,mode,node)
    _run(parse_equation(eqn, inputs)..., String(eqn), inputs)
function solve(residual::Node, inputs::Dict)    # symbolic Node input: build_from_residual
    var, mode, node, repr = build_from_residual(residual, inputs)
    return _run(var, mode, node, repr, inputs)
end

# mesh factory: 1D, structured 2D, or the gmsh-read Jexpresso 2D SEM mesh
function build_mesh(inputs::Dict)
    nsd_in = Int(get(inputs, :nsd, 1))
    nsd_in == 1 && return FDMesh1D(inputs)
    nsd_in == 2 || error("SymbolicFD: nsd must be 1 or 2 (got $nsd_in).")
    disc = discretization(inputs)
    if disc isa SEMMethod && (Bool(get(inputs, :lread_gmsh, false)) || haskey(inputs, :gmsh_filename))
        return build_jex_sem_mesh(inputs, disc)     # read the EXISTING gmsh grid
    end
    return FDMesh2D(inputs)
end

# discrete integral ∫q dV (mass): cell volume for FD, lumped mass for SEM
integrate(m::FDMesh1D, q)  = m.sem === nothing ? sum(q) * m.Δx        : sum(q ./ m.sem.Minv)
integrate(m::FDMesh2D, q)  = m.sem === nothing ? sum(q) * m.dx * m.dy : sum(q ./ m.sem.Minv)
integrate(m::JexSEMMesh, q) = sum(q ./ m.Minv)

_domain_str(m::FDMesh1D) = "[$(m.xmin), $(m.xmax)]"
_domain_str(m::FDMesh2D) = "[$(m.xmin), $(m.xmax)]×[$(m.ymin), $(m.ymax)]"
_domain_str(m::JexSEMMesh) =
    @sprintf("[%.3g, %.3g]×[%.3g, %.3g] (gmsh)", minimum(m.x), maximum(m.x), minimum(m.y), maximum(m.y))

function _run(var, mode, node, eqn_repr::AbstractString, inputs::Dict)
    nsd_in = Int(get(inputs, :nsd, 1))

    # mesh: default to periodic for transient, non-periodic for steady (unless set)
    inputs_mesh = haskey(inputs, :periodic) ? inputs :
                  merge(inputs, Dict(:periodic => (mode == :transient)))
    mesh = build_mesh(inputs_mesh)

    println("="^72)
    println(" SymbolicFD : composing operators for")
    println("   ", eqn_repr)
    println("   method          : ", method_name(mesh.disc))
    if mode == :transient
        println("   discretized RHS : ∂$var/∂t = ", ast_str(node))
    else
        println("   steady residual : ", ast_str(node), " = 0   (solve for $var)")
    end
    println("   nsd = $nsd_in, npoin = $(mesh.npoin), domain = ", _domain_str(mesh), ", ",
            mesh.periodic ? "periodic" : "non-periodic")
    println("="^72)

    outdir = String(get(inputs, :output_dir, "."))
    isdir(outdir) || mkpath(outdir)
    usepng = lowercase(String(get(inputs, :outformat, "png"))) == "png"
    live   = Bool(get(inputs, :plot_live, true))

    # ============================ steady / elliptic ============================
    if mode == :steady
        nsd_in == 1 || error("SymbolicFD: steady (no ∂/∂t) solves are 1D-only so far.")
        q, rmax = solve_steady(node, mesh, inputs)
        @printf("   max interior residual = %.3e\n", rmax)
        @printf("   q: min = %.6g, max = %.6g\n", minimum(q), maximum(q))
        usepng && plot_results_png(mesh.x, q, var, "steady", outdir, 0; live = live)
        open(joinpath(outdir, "solution.csv"), "w") do io
            println(io, "x,$var")
            for i in 1:mesh.npoin
                @printf(io, "%.10e,%.10e\n", mesh.x[i], q[i])
            end
        end
        println("   wrote ", joinpath(outdir, "solution.csv"))
        usepng ? println("   wrote ", outdir, "/fields-it0.png") :
                 ascii_plot(mesh.x, q; label = "$var(x)")
        println("="^72)
        return mesh, q, q
    end

    # =============================== transient =================================
    dqdt = node
    is2d = nsd_in == 2

    # initial condition: q0(x) in 1D, q0(x,y) in 2D
    default_q0 = is2d ? ((x, y) -> exp(-(x^2 + y^2) / (2 * 0.1^2))) :
                        (x -> exp(-(x^2) / (2 * 0.1^2)))
    q0fun = get(inputs, :q0, default_q0)
    q0 = _eval_func(q0fun, mesh)
    q  = copy(q0)

    # time stepping
    cfl  = Float64(get(inputs, :CFL, 0.5))
    tend = Float64(get(inputs, :tend, 1.0))
    dt   = haskey(inputs, :Δt) ? Float64(inputs[:Δt]) : stable_dt(inputs, mesh, cfl)
    nsteps = max(1, ceil(Int, tend / dt))
    dt     = tend / nsteps
    rhs! = build_rhs(dqdt, mesh)
    @printf("   Δt = %.3e, nsteps = %d, t_end = %.4g\n", dt, nsteps, tend)

    # output helpers (1D line plot vs 2D node map)
    plot_init(q)        = is2d ? plot_field_png_2d(mesh, q, var, "t=0", outdir, 0; init = true) :
                                 plot_initial_png(mesh.x, q, var, outdir)
    plot_snap(q, t, it) = is2d ? plot_field_png_2d(mesh, q, var, @sprintf("t=%.4g", t), outdir, it; live = live) :
                                 plot_results_png(mesh.x, q, var, @sprintf("t=%.4g", t), outdir, it; live = live)

    nout      = max(1, Int(get(inputs, :ndiagnostics_outputs, 10)))
    out_every = max(1, cld(nsteps, nout))
    if usepng
        plot_init(q0)
        plot_snap(q0, 0.0, 0)
    end

    # time loop with on-the-fly diagnostics
    w, iout = RK4Work(mesh.npoin), 0
    for sstep in 1:nsteps
        rk4_step!(q, rhs!, dt, w)
        if sstep % out_every == 0 || sstep == nsteps
            iout += 1
            t = sstep * dt
            usepng && plot_snap(q, t, iout)
            @printf("   output %3d : t = %.4g, peak = %.6g\n", iout, t, maximum(q))
        end
    end

    # diagnostics & data dump
    mass0, mass1 = integrate(mesh, q0), integrate(mesh, q)
    @printf("   mass: initial = %.6g, final = %.6g  (Δ = %.3e)\n", mass0, mass1, mass1-mass0)
    @printf("   peak: initial = %.6g, final = %.6g\n", maximum(q0), maximum(q))
    csv = save_csv(joinpath(outdir, "solution.csv"), mesh, q0, q)
    println("   wrote ", csv)
    usepng && println("   wrote ", outdir, "/INIT-", var, ".png and ",
                      outdir, "/fields-it{0..", iout, "}.png")
    if !usepng && !is2d
        ascii_plot(mesh.x, q0; label = "$var(x, t=0)")
        ascii_plot(mesh.x, q;  label = "$var(x, t=$(round(tend, digits=3)))")
    end
    println("="^72)

    return mesh, q0, q
end

end # module
