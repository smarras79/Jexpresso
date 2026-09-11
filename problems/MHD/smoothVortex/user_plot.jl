#---------------------------------------------------------------------------------
# Convergence history of the smooth MHD vortex, in the layout of
# Dao & Nazarov, J. Sci. Comput. 92:77 (2022), Fig. 1: the error of the
# velocity against the EXACT solution, against 1/sqrt(#DOFs) on log-log axes,
# one line per polynomial order, with slope guides.
#
# At the final time every run stores its error in
# `errors/nop<N>_nelx<M>_<visc>.dat` in this case directory, and the figure is
# drawn from EVERY error stored there — so a sweep over orders and meshes
# builds the whole figure and each run replaces only its own point:
#
#     tools/smooth_vortex_mesh.sh     # the meshes, once
#     tools/smooth_vortex_scan.sh     # the sweep
#
# `rm -r problems/MHD/smoothVortex/errors` starts a fresh comparison. Only
# errors from the same final time are drawn together.
#
# The hook is src/io/plotting/jeplots.jl (plot_triangulation, NSD_2D):
#   mesh     this rank's mesh (coordinates, connectivity, extents)
#   q        flat npoin*nvar vector of the OUTPUT variables (`outvar`)
#   t        simulation time
#   Minv     the solver's assembled inverse lumped mass — the quadrature
#            weights the norms below are taken with
#---------------------------------------------------------------------------------
const SV_ERR_DIR = joinpath(@__DIR__, "errors")

# One (colour, marker) per order, fixed so an order looks the same from one
# figure to the next.
const SV_STYLE = Dict(
    2 => (:goldenrod, :utriangle),
    3 => (:purple,    :rect),
    4 => (:seagreen,  :diamond),
    5 => (:royalblue, :circle),
    6 => (:crimson,   :star5),
    7 => (:black,     :xcross),
)
_sv_style(nop) = get(SV_STYLE, nop, (:gray, :cross))

_sv_tag(inputs) = (get(inputs, :lvisc, true) ? "dsgs" : "galerkin")

#---------------------------------------------------------------------------------
# Error of the velocity (u, v) against the exact solution: the vortex of
# initialize.jl translated by v₀t and wrapped periodically. Integral norms are
# taken with the nodal quadrature weights of the mesh, which for LGL nodes is
# the mass-matrix lumping the solver itself uses.
#---------------------------------------------------------------------------------
function _sv_velocity_error(mesh, q, t, outvar, inputs, Minv)
    names = string.(outvar)
    iu = findfirst(==("u"), names)
    iv = findfirst(==("v"), names)
    (iu === nothing || iv === nothing) && return nothing
    npoin = mesh.npoin

    γ  = γ_mhd
    Lx = mesh.xmax - mesh.xmin
    Ly = mesh.ymax - mesh.ymin

    # Nodal quadrature weights: the SOLVER'S OWN assembled lumped mass,
    # w_i = M_ii = 1/Minv_i, handed down by the plotting hook. That is the
    # same quadrature the DSS and the mass-matrix inversion of the run use, so
    # the norms below are the solver's own integrals and not a rule invented
    # here. On this mesh they sum to the domain area to 14 digits.
    (Minv !== nothing && length(Minv) >= npoin) ||
        error("smoothVortex: the error norms need the solver's lumped mass (Minv) from the plotting hook")
    # One weight per DEGREE OF FREEDOM. On a doubly periodic mesh two local
    # nodes can be the same unknown — the right column of the box is the left
    # one — and they carry the same assembled mass, so summing over all npoin
    # counts the periodic edges twice (the weights integrate 105 instead of
    # the box's 100 at 4x4). The solver's own local-to-global map,
    # mesh.ip2gip, is what says which nodes are the same unknown; use it
    # rather than any geometric test of ours.
    seen = Set{eltype(mesh.ip2gip)}()
    w    = zeros(npoin)
    for ip = 1:npoin
        g = mesh.ip2gip[ip]
        g in seen && continue
        push!(seen, g)
        w[ip] = 1.0/Minv[ip]
    end

    s1 = 0.0; s2 = 0.0; si = 0.0
    r1 = 0.0; r2 = 0.0; ri = 0.0
    for ip = 1:npoin
        xr = sv_wrap(mesh.x[ip] - SV_XC - SV_U0*t, Lx)
        yr = sv_wrap(mesh.y[ip] - SV_YC - SV_V0*t, Ly)
        se = sv_state(xr, yr, γ)
        ue = se[2]/se[1]; ve = se[3]/se[1]
        uh = q[(iu - 1)*npoin + ip]
        vh = q[(iv - 1)*npoin + ip]
        e  = sqrt((uh - ue)^2 + (vh - ve)^2)
        m  = sqrt(ue*ue + ve*ve)
        s1 += w[ip]*e;   s2 += w[ip]*e*e;   si = max(si, e)
        r1 += w[ip]*m;   r2 += w[ip]*m*m;   ri = max(ri, m)
    end

    # The weights must integrate the domain: a cheap check that the lumped
    # mass handed down is this rank's and complete.
    if get(ENV, "JEXPRESSO_SV_DEBUG", "") == "1"
        @info "smoothVortex quadrature check" sum_w = sum(w) area = Lx*Ly npoin = npoin nunique = length(seen)
    end
    comm = get_mpi_comm()
    ndofs = npoin
    if MPI.Comm_size(comm) > 1
        sums  = MPI.Allreduce([s1, s2, r1, r2, Float64(npoin)], MPI.SUM, comm)
        maxs  = MPI.Allreduce([si, ri], MPI.MAX, comm)
        s1, s2, r1, r2 = sums[1], sums[2], sums[3], sums[4]
        si, ri = maxs[1], maxs[2]
        ndofs  = round(Int, sums[5])
    end
    return (ndofs = ndofs,
            l1 = s1/max(r1, eps()),
            l2 = sqrt(s2)/max(sqrt(r2), eps()),
            linf = si/max(ri, eps()))
end

function _sv_save_error(e, inputs, t)
    nop = Int(get(inputs, :nop, 0))
    # :nelx carries mod_inputs' placeholder for a gmsh case, so take the
    # element count from the mesh file name.
    m    = match(r"vortex_(\d+)x", string(get(inputs, :gmsh_filename, "")))
    nelx = m === nothing ? Int(get(inputs, :nelx, 0)) : parse(Int, m.captures[1])
    try
        mkpath(SV_ERR_DIR)
        f = joinpath(SV_ERR_DIR, string("nop", nop, "_nelx", nelx, "_", _sv_tag(inputs), ".dat"))
        open(f, "w") do io
            println(io, "# smooth MHD vortex: relative velocity error against the exact solution")
            println(io, "# nop=", nop, " nelx=", nelx, " ndofs=", e.ndofs,
                        " t=", t, " visc=", _sv_tag(inputs),
                        " Cmin=", Float64(get(inputs, :dsgs_Cmin, 0.0)),
                        " dt=", Float64(get(inputs, :Δt, 0.0)))
            println(io, "# L1 L2 Linf")
            println(io, e.l1, " ", e.l2, " ", e.linf)
        end
    catch err
        @warn "smoothVortex: could not store the error" exception=err
    end
    return nothing
end

function _sv_load_errors(t)
    rows = NamedTuple[]
    isdir(SV_ERR_DIR) || return rows
    for fname in readdir(SV_ERR_DIR)
        (endswith(fname, ".dat") && startswith(fname, "nop")) || continue
        meta = Dict{String,String}(); vals = Float64[]
        try
            for line in eachline(joinpath(SV_ERR_DIR, fname))
                if startswith(line, "#")
                    for tok in split(line)
                        occursin('=', tok) || continue
                        k, v = split(tok, '=', limit = 2); meta[k] = v
                    end
                    continue
                end
                isempty(strip(line)) && continue
                append!(vals, parse.(Float64, split(line)))
            end
        catch err
            @warn "smoothVortex: skipping an unreadable error file" file=fname exception=err
            continue
        end
        length(vals) >= 3 || continue
        nop   = tryparse(Int, get(meta, "nop", ""));   nop   === nothing && continue
        ndofs = tryparse(Int, get(meta, "ndofs", "")); ndofs === nothing && continue
        tc    = tryparse(Float64, get(meta, "t", ""))
        (tc === nothing || abs(tc - t) > 1.0e-8*max(1.0, abs(t))) && continue
        push!(rows, (nop = nop, ndofs = ndofs, visc = get(meta, "visc", "dsgs"),
                     nelx = something(tryparse(Int, get(meta, "nelx", "")), 0),
                     l1 = vals[1], l2 = vals[2], linf = vals[3]))
    end
    return sort!(rows, by = r -> (r.visc, r.nop, r.ndofs))
end

_sv_rate(xs, ys) = (length(xs) < 2 || ys[end-1] <= 0 || ys[end] <= 0) ? NaN :
                   log(ys[end-1]/ys[end])/log(xs[end-1]/xs[end])

#---------------------------------------------------------------------------------
# The figure: one panel per norm, one line per order, one figure per
# stabilization (DynSGS / plain Galerkin) — the two panels of the paper's
# Fig. 1. The abscissa is 1/sqrt(#DOFs) ∝ h, as in the paper.
#---------------------------------------------------------------------------------
function _sv_plot(rows, OUTPUT_DIR, iout)
    for tag in unique(r.visc for r in rows)
        sub  = filter(r -> r.visc == tag, rows)
        nops = sort(unique(r.nop for r in sub))
        any(n -> count(r -> r.nop == n, sub) >= 2, nops) || continue

        panels = Plots.Plot[]
        for (fld, nm) in ((:l1, "L^1"), (:l2, "L^2"), (:linf, "L^\\infty"))
            allx = Float64[]; ally = Float64[]
            pl = Plots.plot(; xscale = :log10, yscale = :log10,
                            xlabel = LaTeXStrings.L"1/\sqrt{\#\mathrm{DOFs}}",
                            ylabel = LaTeXStrings.latexstring(string(
                                "\\|\\mathbf{u}_h-\\mathbf{u}_{exact}\\|_{", nm, "}\\ /\\ \\|\\mathbf{u}_{exact}\\|_{", nm, "}")),
                            framestyle = :box, grid = true,
                            legend = :bottomright, legendfontsize = 8,
                            titlefontsize = 13, guidefontsize = 11, tickfontsize = 10,
                            title = LaTeXStrings.latexstring(string(nm, "\\mathrm{-error},\\ \\mathrm{",
                                     tag == "dsgs" ? "RV\\ (DynSGS)" : "Galerkin", "}")),
                            show = false)
            for nop in nops
                g  = sort(filter(r -> r.nop == nop, sub), by = r -> r.ndofs)
                xs = [1.0/sqrt(r.ndofs) for r in g]
                ys = [getfield(r, fld) for r in g]
                keep = isfinite.(ys) .& (ys .> 0)
                any(keep) || continue
                xs = xs[keep]; ys = ys[keep]
                append!(allx, xs); append!(ally, ys)
                col, mk = _sv_style(nop)
                p = _sv_rate(xs, ys)
                lab = isfinite(p) ?
                      LaTeXStrings.latexstring(string("\\mathrm{nop}\\ ", nop, "\\ (p=", round(p; digits = 2), ")")) :
                      LaTeXStrings.latexstring(string("\\mathrm{nop}\\ ", nop))
                Plots.plot!(pl, xs, ys; line = (col, 1.8, :solid), marker = (mk, 5),
                            markerstrokewidth = 0.8, color = col, label = lab)
            end
            if !isempty(allx)
                x2 = maximum(allx); ymax = maximum(ally); ymin = minimum(ally)
                for (sl, col, anchor) in ((2, :gray40, 2.0*ymax), (5, :gray70, 0.5*ymin))
                    xg = [minimum(allx), x2]
                    Plots.plot!(pl, xg, [anchor*(xi/x2)^sl for xi in xg];
                                line = (col, 1.4, :dash),
                                label = LaTeXStrings.latexstring(string("\\mathrm{slope}\\ ", sl)))
                end
            end
            push!(panels, pl)
        end
        plt = Plots.plot(panels...; layout = (1, 3), size = (1500, 450),
                         left_margin = 9Plots.mm, bottom_margin = 8Plots.mm, show = false)
        _savefig_silent(plt, string(OUTPUT_DIR, "/convergence_", tag, "-it", iout, ".png"))
    end
    return nothing
end

function _sv_report(rows)
    println(" # smoothVortex: relative velocity error against the exact solution")
    println(" #   visc      nop  nelx   DOFs           L1           L2         Linf")
    for r in rows
        println(@sprintf(" #   %-8s  %3d  %4d %6d   %10.3e   %10.3e   %10.3e",
                         r.visc, r.nop, r.nelx, r.ndofs, r.l1, r.l2, r.linf))
    end
    return nothing
end

function user_plot_2d(mesh, q, t, outvar, inputs, OUTPUT_DIR, iout; Minv = nothing)
    isfinite(t) && t > 0.0 || return nothing         # the error figure is a final-time product
    abs(t - Float64(get(inputs, :tend, t))) < 1.0e-8*max(1.0, abs(t)) || return nothing

    e = _sv_velocity_error(mesh, q, t, outvar, inputs, Minv)
    e === nothing && return nothing
    _sv_save_error(e, inputs, t)

    if MPI.Comm_rank(get_mpi_comm()) == 0
        rows = _sv_load_errors(t)
        isempty(rows) || (_sv_report(rows); _sv_plot(rows, OUTPUT_DIR, iout))
    end
    return nothing
end
