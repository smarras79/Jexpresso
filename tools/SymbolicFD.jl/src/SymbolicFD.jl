#=---------------------------------------------------------------------------------
# SymbolicFD.jl
#
# A small, self-contained "write-the-equation-and-solve-it" engine.
#
# The user writes a PDE using Julia's unicode characters (and a few common
# LaTeX shortcuts), e.g.
#
#       ∂q/∂t + ∇⋅(\mathbf{u}q) = \mu∇⋅∇(q)
#
# defines a grid (number of points and the spatial dimension, exactly like
# in Jexpresso through a `user_inputs`-style Dict), and the equation is parsed,
# discretized with a simple finite-difference scheme and integrated in time
# in the background.
#
# This module has NO external dependencies (only LinearAlgebra/Printf from the
# Julia standard library) so it can be dropped into any environment and later
# folded into Jexpresso. The public entry point is `solve(eqn, inputs)`.
#
# Currently supported (1D): conservative advection-diffusion-reaction of a
# single scalar on a periodic grid:
#
#       ∂q/∂t + ∇⋅(u q) = μ ∇⋅∇(q) ( + reaction/source terms )
#
# The design (a list of typed `PDETerm`s assembled into a right-hand-side
# operator) is meant to be extended to 2D/3D and to the Jexpresso
# CL()/PERT() conventions later on.
#---------------------------------------------------------------------------------=#
module SymbolicFD

using Printf
using Plots          # PNG output + on-the-fly window, exactly as in Jexpresso

export solve, parse_equation, FDMesh1D, PDETerm, ascii_plot

#---------------------------------------------------------------------------------
# 1. Symbolic representation of a parsed PDE
#---------------------------------------------------------------------------------
"""
    PDETerm

One additive term of a PDE, already moved to the right-hand side of

    ∂q/∂t = Σ_i  coeff_i * Op_i(q).

Fields
- `kind`   : `:advection`, `:diffusion`, `:reaction` or `:source`.
- `coeff`  : signed scalar multiplier (includes the sign coming from moving the
             term across the `=`, and any numeric/parameter factor such as μ).
- `params` : extra per-term data (e.g. the advecting velocity `u`).
"""
struct PDETerm
    kind::Symbol
    coeff::Float64
    params::Dict{Symbol,Float64}
end

#---------------------------------------------------------------------------------
# 2. Normalization: turn LaTeX shortcuts into the unicode the parser expects
#---------------------------------------------------------------------------------
# We accept a mix of unicode (∂ ∇ ⋅ ∇² Δ μ ν …) and the most common LaTeX
# spellings so the user can copy/paste straight from a paper.
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
LaTeX greek/operator commands to unicode, collapse `∇⋅∇`, `Δ`, `∇^2` to a single
canonical Laplacian token `∇²`, and remove all whitespace.
"""
function normalize_equation(eqn::AbstractString)
    s = String(eqn)

    # \mathbf{u}, \vec{u}, \mathrm{u}, \boldsymbol{u}  ->  u
    for cmd in ("mathbf", "vec", "mathrm", "boldsymbol", "mathbb", "mathcal")
        s = replace(s, Regex("\\\\$cmd\\{([^}]*)\\}") => s"\1")
    end

    # LaTeX greek/operators -> unicode
    for (k, v) in LATEX2UNICODE
        s = replace(s, k => v)
    end

    # spacing / sizing helpers that carry no meaning
    for junk in ("\\left", "\\right", "\\,", "\\;", "\\!", "\\ ", " ", "\t", "\n")
        s = replace(s, junk => "")
    end

    # Canonicalize every spelling of the Laplacian to "∇²"
    s = replace(s, "∇⋅∇" => "∇²")
    s = replace(s, "∇^2" => "∇²")
    s = replace(s, "∇2"  => "∇²")   # tolerate a missing caret
    s = replace(s, "Δ"   => "∇²")

    return s
end

#---------------------------------------------------------------------------------
# 3. Tiny scalar-coefficient evaluator
#---------------------------------------------------------------------------------
# Resolves a leading factor such as "", "-", "μ", "0.1", "2μ", "u" into a number,
# looking parameters up in the `inputs` Dict (e.g. inputs[:μ], inputs[:u]).
const TOKEN_RE = r"[0-9]+\.?[0-9]*|[A-Za-zμνρκαβγλσ_][A-Za-zμνρκαβγλσ_0-9]*"

function eval_scalar(str::AbstractString, inputs::Dict)
    s = strip(String(str))
    s = replace(s, "*" => "")
    isempty(s) && return 1.0

    sign = 1.0
    while !isempty(s) && (s[1] == '+' || s[1] == '-')
        s[1] == '-' && (sign = -sign)
        s = s[2:end]
    end
    isempty(s) && return sign

    val = sign
    for m in eachmatch(TOKEN_RE, s)
        tok = m.match
        if occursin(r"^[0-9]", tok)
            val *= parse(Float64, tok)
        else
            key = Symbol(tok)
            haskey(inputs, key) ||
                error("SymbolicFD: unknown parameter `$tok` in the equation. " *
                      "Add `:$tok => value` to your inputs Dict.")
            val *= Float64(inputs[key])
        end
    end
    return val
end

#---------------------------------------------------------------------------------
# 4. The parser:  equation string  ->  (variable name, Vector{PDETerm})
#---------------------------------------------------------------------------------
# Split a side of the equation into signed top-level additive terms, honouring
# parentheses so that "∇⋅(uq)" is not chopped at an inner sign.
function split_terms(side::AbstractString)
    terms = Tuple{Float64,String}[]
    depth = 0
    sign  = 1.0
    buf   = IOBuffer()
    flush_term!() = begin
        t = String(take!(buf))
        isempty(t) || push!(terms, (sign, t))
    end
    for c in side
        if c == '(' ; depth += 1; print(buf, c)
        elseif c == ')' ; depth -= 1; print(buf, c)
        elseif depth == 0 && (c == '+' || c == '-')
            flush_term!()
            sign = (c == '-') ? -1.0 : 1.0
        else
            print(buf, c)
        end
    end
    flush_term!()
    return terms
end

# Regexes evaluated on the normalized, whitespace-free string.
const TIME_RE = r"^∂(\w+)/∂t$"                 # ∂q/∂t
const DIFF_RE = r"^(.*?)∇²\(?(\w+)\)?$"        # μ∇²(q)  or  μ∇²q
const ADV_RE  = r"^(.*?)∇⋅\((.+)\)$"           # ∇⋅(uq)

"""
    parse_equation(eqn, inputs) -> (var::String, terms::Vector{PDETerm})

Parse a unicode/LaTeX PDE into the unknown variable name and the list of
right-hand-side terms of `∂var/∂t = Σ coeff_i Op_i(var)`.
"""
function parse_equation(eqn::AbstractString, inputs::Dict)
    s = normalize_equation(eqn)
    occursin("=", s) ||
        error("SymbolicFD: the equation must contain `=` (e.g. `∂q/∂t + ∇⋅(uq) = μ∇²q`).")
    lhs, rhs = split(s, "=", limit = 2)

    var   = ""
    terms = PDETerm[]

    # A term is parsed once; `to_rhs` carries the sign needed to move it to the
    # right-hand side of ∂q/∂t = …  (LHS terms flip, RHS terms keep their sign).
    function classify!(sign::Float64, body::String, to_rhs::Float64)
        # --- time derivative: defines the unknown, not added to the RHS ---
        m = match(TIME_RE, body)
        if m !== nothing
            var == "" || var == m.captures[1] ||
                error("SymbolicFD: more than one time-derivative variable found.")
            var = m.captures[1]
            return
        end

        # --- diffusion:  c ∇²q ---
        m = match(DIFF_RE, body)
        if m !== nothing
            c = eval_scalar(m.captures[1], inputs)
            push!(terms, PDETerm(:diffusion, to_rhs * sign * c, Dict{Symbol,Float64}()))
            return
        end

        # --- advection (flux divergence):  c ∇⋅(u q) ---
        m = match(ADV_RE, body)
        if m !== nothing
            lead  = eval_scalar(m.captures[1], inputs)
            inner = m.captures[2]
            # the flux is velocity * var; peel the var out to read the velocity
            velstr = replace(inner, string(var) => "")
            u = eval_scalar(velstr, inputs)
            push!(terms, PDETerm(:advection, to_rhs * sign * lead, Dict(:u => u)))
            return
        end

        # --- linear reaction (c·var) vs. spatially constant source (c) ---
        if endswith(body, var)
            coeffstr = chop(body; tail = length(var))
            c = eval_scalar(coeffstr, inputs)
            push!(terms, PDETerm(:reaction, to_rhs * sign * c, Dict{Symbol,Float64}()))
            return
        end
        c = eval_scalar(body, inputs)
        push!(terms, PDETerm(:source, to_rhs * sign * c, Dict{Symbol,Float64}()))
        return
    end

    for (sgn, body) in split_terms(lhs); classify!(sgn, body, -1.0); end   # LHS flips
    for (sgn, body) in split_terms(rhs); classify!(sgn, body, +1.0); end   # RHS keeps

    var == "" && error("SymbolicFD: no `∂q/∂t`-type time derivative found in the equation.")
    return var, terms
end

#---------------------------------------------------------------------------------
# 5. Finite-difference mesh (1D, Jexpresso-style inputs)
#---------------------------------------------------------------------------------
"""
    FDMesh1D

Uniform 1D finite-difference grid with `npoin` nodes on `[xmin, xmax]`.
For a periodic grid the right boundary node is identified with the left one, so
the nodes are `x_i = xmin + i*Δx`, `i = 0 … npoin-1`, with `Δx = (xmax-xmin)/npoin`.
"""
struct FDMesh1D
    npoin::Int
    xmin::Float64
    xmax::Float64
    Δx::Float64
    x::Vector{Float64}
    periodic::Bool
end

function FDMesh1D(inputs::Dict)
    xmin     = Float64(get(inputs, :xmin, -1.0))
    xmax     = Float64(get(inputs, :xmax,  1.0))
    npoin    = Int(get(inputs, :npoin, get(inputs, :nelx, 100)))
    periodic = Bool(get(inputs, :periodic, true))
    if periodic
        Δx = (xmax - xmin) / npoin
        x  = [xmin + i * Δx for i in 0:npoin-1]
    else
        Δx = (xmax - xmin) / (npoin - 1)
        x  = [xmin + i * Δx for i in 0:npoin-1]
    end
    return FDMesh1D(npoin, xmin, xmax, Δx, x, periodic)
end

#---------------------------------------------------------------------------------
# 6. Finite-difference spatial operators (2nd order central, periodic)
#---------------------------------------------------------------------------------
@inline ip1(i, n) = i == n ? 1 : i + 1
@inline im1(i, n) = i == 1 ? n : i - 1

"central first derivative  d(f)/dx"
function ddx!(out, f, m::FDMesh1D)
    n, h = m.npoin, m.Δx
    @inbounds for i in 1:n
        out[i] = (f[ip1(i,n)] - f[im1(i,n)]) / (2h)
    end
    return out
end

"central second derivative  d²(f)/dx²"
function d2dx2!(out, f, m::FDMesh1D)
    n, h = m.npoin, m.Δx
    @inbounds for i in 1:n
        out[i] = (f[ip1(i,n)] - 2f[i] + f[im1(i,n)]) / (h*h)
    end
    return out
end

#---------------------------------------------------------------------------------
# 7. Assemble the right-hand side  ∂q/∂t = rhs(q)
#---------------------------------------------------------------------------------
function build_rhs(terms::Vector{PDETerm}, m::FDMesh1D)
    n    = m.npoin
    flux = zeros(n)
    tmp  = zeros(n)
    function rhs!(dq, q)
        fill!(dq, 0.0)
        for t in terms
            if t.kind == :advection
                u = t.params[:u]
                @inbounds @. flux = u * q          # conservative flux F = u q
                ddx!(tmp, flux, m)                  # ∇⋅(u q)
            elseif t.kind == :diffusion
                d2dx2!(tmp, q, m)                   # ∇²q
            elseif t.kind == :reaction
                @inbounds @. tmp = q
            else # :source (spatially constant)
                fill!(tmp, 1.0)
            end
            @inbounds @. dq += t.coeff * tmp
        end
        return dq
    end
    return rhs!
end

#---------------------------------------------------------------------------------
# 8. Explicit RK4 time integration with an automatic stable Δt
#---------------------------------------------------------------------------------
function stable_dt(terms::Vector{PDETerm}, m::FDMesh1D, cfl::Float64)
    h = m.Δx
    dt = Inf
    for t in terms
        if t.kind == :advection
            a = abs(t.coeff * t.params[:u])
            a > 0 && (dt = min(dt, cfl * h / a))
        elseif t.kind == :diffusion
            μ = abs(t.coeff)
            μ > 0 && (dt = min(dt, cfl * 0.5 * h*h / μ))
        end
    end
    isfinite(dt) || (dt = cfl * h)   # fallback (pure reaction/source)
    return dt
end

# RK4 work arrays, allocated once and reused across steps.
struct RK4Work
    k1::Vector{Float64}; k2::Vector{Float64}
    k3::Vector{Float64}; k4::Vector{Float64}; qt::Vector{Float64}
end
RK4Work(n::Int) = RK4Work(zeros(n), zeros(n), zeros(n), zeros(n), zeros(n))

"advance `q` by one explicit RK4 step of size `dt`."
function rk4_step!(q, rhs!, dt, w::RK4Work)
    rhs!(w.k1, q)
    @. w.qt = q + 0.5dt*w.k1 ; rhs!(w.k2, w.qt)
    @. w.qt = q + 0.5dt*w.k2 ; rhs!(w.k3, w.qt)
    @. w.qt = q +     dt*w.k3 ; rhs!(w.k4, w.qt)
    @. q += (dt/6) * (w.k1 + 2w.k2 + 2w.k3 + w.k4)
    return q
end

#---------------------------------------------------------------------------------
# 9. Output helpers
#
# PNG figures and the on-the-fly window are produced with Plots.jl in exactly
# the same way as Jexpresso's src/io/plotting/jeplots.jl (used e.g. by
# problems/CompEuler/sod1d): an `INIT-<var>.png` at t = 0 and a `fields-it<iout>.png`
# at every diagnostic output. A CSV is always written and a terminal ASCII plot
# is kept as a display-free fallback (:outformat => "ascii").
#---------------------------------------------------------------------------------
# Write a figure to file without flashing a window (mirrors jeplots._savefig_silent):
# with the GR backend, :overwrite_figure = false routes the export through a
# dedicated file workstation instead of the active (on-screen) one.
function _savefig_silent(plt, fout_name)
    plt[:overwrite_figure] = false
    Plots.savefig(plt, string(fout_name))
    return nothing
end

# Initial condition figure, same look as jeplots.plot_initial(NSD_1D, …).
function plot_initial_png(x, q, var, OUTPUT_DIR::String)
    npoin = length(q)
    plt = Plots.scatter(x[1:npoin], q[1:npoin];
                        markersize = 5,
                        color = :blue,
                        xlabel = "x",
                        ylabel = "$var(x)",
                        title = "$var",
                        titlefontsize = 24,
                        guidefontsize = 18,
                        legendfontsize = 14,
                        tickfontsize = 14,
                        legend = false,
                        size = (800, 600))
    _savefig_silent(plt, string(OUTPUT_DIR, "/INIT-", var, ".png"))
    return plt
end

# One diagnostic output: build the per-variable curve (same style as
# jeplots.plot_results), save it as `fields-it<iout>.png` and, when `live` is
# true, paint the active GR workstation so the window updates on the fly --
# this is the very mechanism jeplots.render_plot_matrix relies on.
function plot_results_png(x, q, var, ttl::String, OUTPUT_DIR::String, iout::Int;
                          live::Bool = true, wfig = 600, hfig = 400)
    sort_idx = sortperm(x)
    plt = Plots.plot(x[sort_idx], q[sort_idx];
                     line = (:blue, 2),
                     marker = (:circle, 5, :blue),
                     title = string(var, "  ", ttl),
                     xlabel = "x",
                     titlefontsize = 22,
                     guidefontsize = 18,
                     legendfontsize = 14,
                     tickfontsize = 14,
                     legend = false,
                     show = false,
                     size = (wfig, hfig))
    fout = string(OUTPUT_DIR, "/fields-it", iout, ".png")
    if live
        # plain savefig on the GR backend repaints the on-screen workstation
        # (live window) AND prints the same canvas to file -- no flicker.
        try; Plots.savefig(plt, fout); catch; end
        try; display(plt); catch; end
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

"""
    ascii_plot(x, q; rows=15, cols=70, label="q")

Render a quick line plot of `q(x)` to the terminal, no plotting package needed.
"""
function ascii_plot(x, q; rows = 15, cols = 70, label = "q")
    n  = length(q)
    qmin, qmax = minimum(q), maximum(q)
    span = qmax - qmin
    span == 0 && (span = 1.0)
    grid = fill(' ', rows, cols)
    for c in 1:cols
        # nearest sample for this column
        idx = round(Int, 1 + (c-1) * (n-1) / (cols-1))
        r   = rows - round(Int, (q[idx] - qmin) / span * (rows-1))
        r   = clamp(r, 1, rows)
        grid[r, c] = '*'
    end
    println()
    @printf("  %s   [min = %.4g, max = %.4g]\n", label, qmin, qmax)
    for r in 1:rows
        println("  |", String(grid[r, :]))
    end
    println("  +", repeat("-", cols))
    @printf("   x = %.3g %s x = %.3g\n", x[1], repeat(" ", cols-12), x[end])
    return nothing
end

#---------------------------------------------------------------------------------
# 10. Top-level driver
#---------------------------------------------------------------------------------
"""
    solve(eqn, inputs) -> (mesh, q0, q)

Parse the unicode/LaTeX PDE string `eqn`, build a 1D finite-difference
discretization from the `user_inputs`-style `inputs` Dict, integrate in time
and return the mesh together with the initial and final solution vectors.

Required / recognised `inputs` keys
- `:xmin`, `:xmax`            domain (default [-1, 1])
- `:npoin` (or `:nelx`)       number of grid points (default 100)
- `:periodic`                 true/false (default true)
- `:nsd`                      spatial dimension (only 1 supported for now)
- `:u`, `:μ`, …               any parameter symbol appearing in the equation
- `:q0`                       function `x -> value` for the initial condition
- `:tend`                     final time
- `:Δt`                       time step (optional; otherwise CFL-derived)
- `:CFL`                      CFL number used for the automatic Δt (default 0.5)
- `:output_dir`              where to write the figures / `solution.csv` (default ".")
- `:outformat`               "png" (Plots.jl, like sod1d) or "ascii" (default "png")
- `:ndiagnostics_outputs`     number of on-the-fly plot snapshots (default 10)
- `:plot_live`                update an on-screen window on the fly (default true)
"""
function solve(eqn::AbstractString, inputs::Dict)
    nsd = Int(get(inputs, :nsd, 1))
    nsd == 1 || error("SymbolicFD: only nsd = 1 is implemented so far (got nsd = $nsd).")

    var, terms = parse_equation(eqn, inputs)
    mesh = FDMesh1D(inputs)

    println("="^72)
    println(" SymbolicFD : solving the user-written equation")
    println("   ", eqn)
    println("   normalized form : ", normalize_equation(eqn))
    println("   unknown         : ", var)
    println("   nsd = $nsd, npoin = $(mesh.npoin), domain = [$(mesh.xmin), $(mesh.xmax)], ",
            mesh.periodic ? "periodic" : "non-periodic")
    print("   detected terms  : ∂$var/∂t =")
    for (k, t) in enumerate(terms)
        op = t.kind == :advection ? "∇⋅(u q)" :
             t.kind == :diffusion ? "∇²q"     :
             t.kind == :reaction  ? "q"       : "1"
        @printf(" %s%.4g·%s", t.coeff ≥ 0 ? (k==1 ? "" : "+ ") : "- ", abs(t.coeff), op)
    end
    println()
    println("="^72)

    # initial condition
    q0fun = get(inputs, :q0, x -> exp(-(x^2) / (2 * 0.1^2)))   # default gaussian
    q0 = Float64[q0fun(xi) for xi in mesh.x]
    q  = copy(q0)

    # time stepping
    cfl  = Float64(get(inputs, :CFL, 0.5))
    tend = Float64(get(inputs, :tend, 1.0))
    dt   = haskey(inputs, :Δt) ? Float64(inputs[:Δt]) : stable_dt(terms, mesh, cfl)
    nsteps = max(1, ceil(Int, tend / dt))
    dt     = tend / nsteps   # land exactly on tend

    rhs! = build_rhs(terms, mesh)
    @printf("   Δt = %.3e, nsteps = %d, t_end = %.4g\n", dt, nsteps, tend)

    # ---- output set-up -------------------------------------------------------
    outdir = String(get(inputs, :output_dir, "."))
    isdir(outdir) || mkpath(outdir)
    outformat = String(get(inputs, :outformat, "png"))
    usepng    = lowercase(outformat) == "png"
    live      = Bool(get(inputs, :plot_live, true))
    nout      = max(1, Int(get(inputs, :ndiagnostics_outputs, 10)))
    out_every = max(1, cld(nsteps, nout))

    # initial snapshot
    if usepng
        plot_initial_png(mesh.x, q0, var, outdir)
        plot_results_png(mesh.x, q0, var, @sprintf("t=%.4g", 0.0), outdir, 0; live = live)
    end

    # ---- time loop with on-the-fly diagnostics ------------------------------
    w    = RK4Work(mesh.npoin)
    iout = 0
    for s in 1:nsteps
        rk4_step!(q, rhs!, dt, w)
        if s % out_every == 0 || s == nsteps
            iout += 1
            t = s * dt
            if usepng
                plot_results_png(mesh.x, q, var, @sprintf("t=%.4g", t), outdir, iout;
                                 live = live)
            end
            @printf("   output %3d : t = %.4g, peak = %.6g\n", iout, t, maximum(q))
        end
    end

    # ---- diagnostics & data dump --------------------------------------------
    mass0 = sum(q0) * mesh.Δx
    mass1 = sum(q)  * mesh.Δx
    @printf("   mass: initial = %.6g, final = %.6g  (Δ = %.3e)\n",
            mass0, mass1, mass1 - mass0)
    @printf("   peak: initial = %.6g, final = %.6g\n", maximum(q0), maximum(q))

    csv = save_csv(joinpath(outdir, "solution.csv"), mesh, q0, q)
    println("   wrote ", csv)
    usepng && println("   wrote ", outdir, "/INIT-", var, ".png and ",
                      outdir, "/fields-it{0..", iout, "}.png")

    # terminal ASCII fallback (kept for headless / dependency-free use)
    if !usepng
        ascii_plot(mesh.x, q0; label = "$var(x, t=0)")
        ascii_plot(mesh.x, q;  label = "$var(x, t=$(round(tend, digits=3)))")
    end
    println("="^72)

    return mesh, q0, q
end

end # module
