#!/usr/bin/env julia
#
# LES Statistics: compute resolved vertical velocity variance  w'w'
# from simulation.pvd output (CompEuler / Jexpresso)
#
# Usage:  julia les_statistics.jl [n_last_steps]
#   n_last_steps : number of final time‑steps to average over (default 10)
#

using LightXML
using CairoMakie
using Statistics
using Printf

include("../src/kernel/physics/globalConstantsPhysics.jl")
include("../src/kernel/physics/constitutiveLaw.jl")


# ──────────────────────────────────────────────────────────────────
# 1.  Low‑level VTU reader  (appended‑binary, raw encoding, UInt64 header)
# ──────────────────────────────────────────────────────────────────

"""
Read a single VTU file and return a Dict with:
  :points   => (3, npts) Float64 array
  :pointdata => Dict{String, Vector{Float64}}  (scalar fields)
"""
function read_vtu(filepath::String)
    raw = read(filepath)                       # full file as bytes

    # ---- locate the appended‑data binary blob ------------------
    marker = Vector{UInt8}(codeunits("<AppendedData encoding=\"raw\">\n_"))
    idx = nothing
    for i in 1:(length(raw) - length(marker) + 1)
        if @view(raw[i:i+length(marker)-1]) == marker
            idx = i + length(marker)           # first byte after the '_'
            break
        end
    end
    if idx === nothing
        # try without the newline
        marker2 = Vector{UInt8}(codeunits("<AppendedData encoding=\"raw\">_"))
        for i in 1:(length(raw) - length(marker2) + 1)
            if @view(raw[i:i+length(marker2)-1]) == marker2
                idx = i + length(marker2)
                break
            end
        end
    end
    idx === nothing && error("Could not find AppendedData section in $filepath")
    blob_start = idx   # 1‑based index into `raw`

    # ---- parse only the XML header (before the binary blob) ----
    # Find the start of <AppendedData ...> tag and truncate there,
    # then close the XML so LightXML can parse it.
    appended_tag = Vector{UInt8}(codeunits("<AppendedData"))
    xml_end = nothing
    for i in 1:(length(raw) - length(appended_tag) + 1)
        if @view(raw[i:i+length(appended_tag)-1]) == appended_tag
            xml_end = i - 1
            break
        end
    end
    xml_end === nothing && error("Could not find <AppendedData> tag in $filepath")
    xml_header = String(raw[1:xml_end]) * "\n</VTKFile>\n"
    doc  = parse_string(xml_header)
    root = LightXML.root(doc)
    ugrid = find_element(root, "UnstructuredGrid")
    piece = find_element(ugrid, "Piece")
    npts  = parse(Int, LightXML.attribute(piece, "NumberOfPoints"))

    # helper: read one DataArray from the blob
    function read_array(da_elem)
        T = LightXML.attribute(da_elem, "type")
        jltype = T == "Float64" ? Float64 :
                 T == "Float32" ? Float32 :
                 T == "Int64"   ? Int64   :
                 T == "Int32"   ? Int32   :
                 T == "UInt8"   ? UInt8   :
                 error("Unsupported type $T")
        ncomp  = parse(Int, xml_attr(da_elem, "NumberOfComponents"; default="1"))
        offset = parse(Int, LightXML.attribute(da_elem, "offset"))
        pos = blob_start + offset              # position in `raw`
        nbytes = reinterpret(UInt64, raw[pos:pos+7])[1]
        nelem  = Int(nbytes ÷ sizeof(jltype))
        data   = reinterpret(jltype, raw[pos+8 : pos+8+Int(nbytes)-1])
        return copy(data), ncomp
    end

    # Points
    pts_elem = find_element(find_element(piece, "Points"), "DataArray")
    pts_flat, _ = read_array(pts_elem)
    points = reshape(pts_flat, 3, npts)

    # PointData
    pd_node = find_element(piece, "PointData")
    pointdata = Dict{String, Vector{Float64}}()
    for da in get_elements_by_tagname(pd_node, "DataArray")
        name = LightXML.attribute(da, "Name")
        arr, _ = read_array(da)
        pointdata[name] = Float64.(arr)
    end

    free(doc)
    return Dict(:points => points, :pointdata => pointdata, :npts => npts)
end

"""
    xml_attr(elem, name; default=nothing)
Read an XML attribute, returning `default` if missing.
"""
function xml_attr(e::XMLElement, name::AbstractString; default=nothing)
    val = LightXML.attribute(e, name)
    return val === nothing ? default : val
end

"""
Read a PVTU file: merge all pieces.
"""
function read_pvtu(pvtu_path::String)
    basedir = dirname(pvtu_path)
    doc  = parse_file(pvtu_path)
    root_el = LightXML.root(doc)
    pgrid = find_element(root_el, "PUnstructuredGrid")

    all_points = Vector{Matrix{Float64}}()
    all_pd     = Dict{String, Vector{Vector{Float64}}}()

    for piece_el in get_elements_by_tagname(pgrid, "Piece")
        src = LightXML.attribute(piece_el, "Source")
        vtu_path = joinpath(basedir, src)
        d = read_vtu(vtu_path)
        push!(all_points, d[:points])
        for (k, v) in d[:pointdata]
            haskey(all_pd, k) || (all_pd[k] = Vector{Vector{Float64}}())
            push!(all_pd[k], v)
        end
    end
    free(doc)

    points = hcat(all_points...)
    pointdata = Dict(k => vcat(vs...) for (k, vs) in all_pd)
    return points, pointdata
end

# ──────────────────────────────────────────────────────────────────
# 2.  Parse simulation.pvd
# ──────────────────────────────────────────────────────────────────

function parse_pvd(pvd_path::String)
    doc  = parse_file(pvd_path)
    root_el = LightXML.root(doc)
    coll = find_element(root_el, "Collection")
    timesteps = Float64[]
    files     = String[]
    for ds in get_elements_by_tagname(coll, "DataSet")
        push!(timesteps, parse(Float64, LightXML.attribute(ds, "timestep")))
        push!(files,     LightXML.attribute(ds, "file"))
    end
    free(doc)
    return timesteps, files
end

# ──────────────────────────────────────────────────────────────────
# 3.  Compute horizontally‑averaged  w'w'  profile
# ──────────────────────────────────────────────────────────────────

"""
Given points (3×N) and pointdata dict, compute w = ρw/ρ, then
horizontally average w'w' on z‑bins.

Returns (z_centers, wpwp_profile)
"""
function compute_wpwp(points, pointdata; nbins=100)
    # ρ  = pointdata["ρ"]
    ρ  = 1.0
    ρw = pointdata["ρw"]
    w  = ρw ./ ρ                           # vertical velocity

    z  = points[3, :]
    zmin, zmax = extrema(z)
    dz = (zmax - zmin) / nbins
    z_edges = range(zmin, zmax + dz/2, length=nbins+1)

    z_centers   = zeros(nbins)
    wpwp_prof   = zeros(nbins)
    count_prof  = zeros(Int, nbins)

    # first pass: bin‑averaged w̄(z)
    w_bar = zeros(nbins)
    for i in eachindex(z)
        bin = clamp(Int(floor((z[i] - zmin) / dz)) + 1, 1, nbins)
        w_bar[bin]      += w[i]
        count_prof[bin] += 1
        z_centers[bin]  += z[i]
    end
    for b in 1:nbins
        if count_prof[b] > 0
            w_bar[b]     /= count_prof[b]
            z_centers[b] /= count_prof[b]
        end
    end

    # second pass: w'w'
    fill!(wpwp_prof, 0.0)
    for i in eachindex(z)
        bin = clamp(Int(floor((z[i] - zmin) / dz)) + 1, 1, nbins)
        wp = w[i] - w_bar[bin]
        wpwp_prof[bin] += wp^2
    end
    for b in 1:nbins
        count_prof[b] > 0 && (wpwp_prof[b] /= count_prof[b])
    end

    return z_centers, wpwp_prof, count_prof
end

function compute_q(points; nbins=100, q=pointdata["ρu"])
    u  = q

    z  = points[3, :]
    zmin, zmax = extrema(z)
    dz = (zmax - zmin) / nbins

    z_centers  = zeros(nbins)
    u_bar      = zeros(nbins)
    count_prof = zeros(Int, nbins)

    for i in eachindex(z)
        bin = clamp(Int(floor((z[i] - zmin) / dz)) + 1, 1, nbins)
        u_bar[bin]      += u[i]
        count_prof[bin] += 1
        z_centers[bin]  += z[i]
    end
    for b in 1:nbins
        if count_prof[b] > 0
            u_bar[b]     /= count_prof[b]
            z_centers[b] /= count_prof[b]
        end
    end

    return z_centers, u_bar, count_prof
end


# ──────────────────────────────────────────────────────────────────
# 4.  Main
# ──────────────────────────────────────────────────────────────────

function main()
    n_last = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 5

    pvd_dir  = @__DIR__
    pvd_path = joinpath(pvd_dir, "simulation.pvd")
    isfile(pvd_path) || error("simulation.pvd not found in $pvd_dir")

    timesteps, files = parse_pvd(pvd_path)
    ntotal = length(timesteps)
    n_use  = min(n_last, ntotal)
    idx_range = (ntotal - n_use + 1):ntotal

    println("──────────────────────────────────────────")
    println("LES Statistics:  w'w'  (resolved)")
    println("  Total snapshots in PVD : $ntotal")
    println("  Averaging last $n_use snapshots")
    @printf("  Time window : %.1f s  →  %.1f s\n",
            timesteps[idx_range[1]], timesteps[idx_range[end]])
    println("──────────────────────────────────────────")

    # Accumulate time‑averaged profile
    nbins = 10000
    wpwp_avg = nothing
    u_avg    = nothing
    z_avg    = nothing
    θ_avg    = nothing
    P_avg    = nothing
    ρ_avg    = nothing
    cnt_avg  = nothing

    for (count, idx) in enumerate(idx_range)
        t    = timesteps[idx]
        file = files[idx]
        pvtu = joinpath(pvd_dir, file)
        @printf("  [%d/%d]  t = %.1f s   %s\n", count, n_use, t, file)

        points, pd = read_pvtu(pvtu)
        zc, wpwp, cnt = compute_wpwp(points, pd; nbins=nbins)
        zc, u, _      = compute_q(points; nbins=nbins, q=pd["ρu"])
        zc, θ, _      = compute_q(points; nbins=nbins, q=pd["θ_p"])
        zc, P, _      = compute_q(points; nbins=nbins, q=pd["pressure"])
        zc, ρ, _      = compute_q(points; nbins=nbins, q=pd["ρ"])

        if wpwp_avg === nothing
            wpwp_avg = copy(wpwp)
            z_avg    = copy(zc)
            u_avg    = copy(u)
            θ_avg    = copy(θ)
            P_avg    = copy(P)
            ρ_avg    = copy(ρ)
            cnt_avg  = copy(cnt)
        else
            wpwp_avg .+= wpwp
            u_avg    .+= u
            θ_avg    .+= θ
            P_avg    .+= P
            ρ_avg    .+= ρ
            cnt_avg  .+= cnt
        end
    end
    wpwp_avg ./= n_use
    u_avg    ./= n_use
    θ_avg    ./= n_use
    P_avg    ./= n_use
    ρ_avg    ./= n_use

    # ── Remove empty bins and sort by z ─────────────────────────
    mask = cnt_avg .> 0
    z_avg    = z_avg[mask]
    wpwp_avg = wpwp_avg[mask]
    u_avg    = u_avg[mask]
    θ_avg    = θ_avg[mask]
    P_avg    = P_avg[mask]
    ρ_avg    = ρ_avg[mask]

    sp = sortperm(z_avg)
    z_avg    = z_avg[sp]
    wpwp_avg = wpwp_avg[sp]
    u_avg    = u_avg[sp]
    θ_avg    = θ_avg[sp]
    P_avg    = P_avg[sp]
    ρ_avg    = ρ_avg[sp]
    zmin, zmax = extrema(z_avg)



    # ── Save profile to CSV ──────────────────────────────────────
    csv_path = joinpath(pvd_dir, "LES_profile.csv")
    PhysConst = PhysicalConst{Float64}()
    Tabs = copy(P_avg)
    open(csv_path, "w") do io
        println(io, "z,wpwp,u")
        for i in eachindex(z_avg)
            θ_d =  perfectGasLaw_ρPtoθ(PhysConst;ρ=ρ_avg[i], Press=P_avg[i])
            Tabs[i] = θ_d/(PhysConst.pref/P_avg[i])^(1/PhysConst.cpoverR)

            if z_avg[i] > 15000
                @printf(io, "%.1f %.1f %.1f 0.000 -3.00 0.00\n", z_avg[i], P_avg[i], θ_d)
            end
            # @printf(io, "%.6f %.6f %.6f %.6f %.6f\n", z_avg[i], wpwp_avg[i], u_avg[i], θ_avg[i], θ_d)
        end
    end
    println("\n  Profile saved to: $csv_path")

    # ── Plot ─────────────────────────────────────────────────────
    fig = Figure(size=(600, 800))
    ax  = Axis(fig[1, 1],
               xlabel = L"\overline{w^\prime w^\prime} \; [\mathrm{m^2/s^2}]",
               ylabel = "z  [m]",
               yticks = zmin:500:zmax,   # every 500 m
               title  = @sprintf("Resolved w'w'  (avg over last %d steps, t = %.0f–%.0f s)",
                                 n_use, timesteps[idx_range[1]], timesteps[idx_range[end]]))
    scatterlines!(ax, wpwp_avg, z_avg; linewidth=2, markersize=6, color=:blue)
    fig_path = joinpath(pvd_dir, "wpwp_profile.png")
    save(fig_path, fig; px_per_unit=2)
    println("  Figure saved to: $fig_path")
    println("──────────────────────────────────────────")

    fig = Figure(size=(600, 800))
    ax  = Axis(fig[1, 1],
               xlabel = L"\overline{u} \; [\mathrm{m/s}]",
               ylabel = "z  [m]",
               yticks = zmin:500:zmax,   # every 500 m
               title  = @sprintf("Resolved u  (avg over last %d steps, t = %.0f–%.0f s)",
                                 n_use, timesteps[idx_range[1]], timesteps[idx_range[end]]))
    scatterlines!(ax, u_avg, z_avg; linewidth=2, markersize=6, color=:black)
    fig_path = joinpath(pvd_dir, "u_profile.png")
    # fig_path = joinpath(pvd_dir, "wpwp_profile.png")
    save(fig_path, fig; px_per_unit=2)
    println("  Figure saved to: $fig_path")
    println("──────────────────────────────────────────")

    fig = Figure(size=(600, 800))
    ax  = Axis(fig[1, 1],
               xlabel = L"\overline{θ} \; [\mathrm{K}]",
               ylabel = "z  [km]",
               yticks = zmin/1000:1:zmax/1000,   # every 500 m
               title  = @sprintf("Resolved θ  (avg over last %d steps, t = %.0f–%.0f s)",
                                 n_use, timesteps[idx_range[1]], timesteps[idx_range[end]]))
    scatterlines!(ax, Tabs, z_avg./1000; linewidth=2, markersize=6, color=:black)
    fig_path = joinpath(pvd_dir, "theta_profile.png")
    # fig_path = joinpath(pvd_dir, "wpwp_profile.png")
    save(fig_path, fig; px_per_unit=2)
    println("  Figure saved to: $fig_path")
    println("──────────────────────────────────────────")

    # return z_avg, wpwp_avg
end

main()
