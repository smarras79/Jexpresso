#!/usr/bin/env julia
#
# Plot 1-D LES profile figures from les_statistics_tavg.dat and les_stress_tavg.dat
#
# Usage:
#   julia tools/plot_les_tavg.jl <output_dir>
#   julia tools/plot_les_tavg.jl output/CompEuler/LESICP6/output
#
using CairoMakie
using Printf

# ──────────────────────────────────────────────────────────────────
# Parser
# ──────────────────────────────────────────────────────────────────

"""
    read_les_tavg(filepath) -> (varnames, z, data)

Read a les_*_tavg.dat file.  Returns:
  - varnames : String vector of variable names (length nvars, excludes "z")
  - z        : Float64 vector of z-levels (length nz)
  - data     : Float64 matrix (nz × nvars)

Header format:   # time_end=...  n_samples=...  z  var1  var2  ...
Data row format: time  z  var1  var2  ...
"""
function read_les_tavg(filepath::String)
    lines = readlines(filepath)
    isempty(lines) && error("Empty file: $filepath")

    # -- parse header --
    header = lstrip(lines[1], ['#', ' '])
    tokens = split(header)
    # drop "time_end=..." and "n_samples=..." tokens
    colnames = [t for t in tokens
                if !startswith(t, "time_end=") && !startswith(t, "n_samples=")]
    # colnames[1] == "z", colnames[2:end] == variable names

    # -- parse data rows --
    rows = Vector{Vector{Float64}}()
    for line in lines[2:end]
        s = strip(line)
        (isempty(s) || s[1] == '#') && continue
        # data row: time  z  var1  var2  ...
        cols = parse.(Float64, split(s))
        push!(rows, cols[2:end])   # drop leading "time" column
    end

    isempty(rows) && error("No data rows in $filepath")
    mat = reduce(hcat, rows)'   # nz × length(colnames)

    z        = mat[:, 1]
    data     = mat[:, 2:end]
    varnames = String.(colnames[2:end])
    return varnames, z, data
end

# ──────────────────────────────────────────────────────────────────
# Plotting helpers
# ──────────────────────────────────────────────────────────────────

function nice_xlabel(name::String)
    labels = Dict(
        "u_mean"    => "ū  [m/s]",
        "v_mean"    => "v̄  [m/s]",
        "w_mean"    => "w̄  [m/s]",
        "t_mean"    => "θ̄  [K]",
        "p_mean"    => "P̄ [Pa]",   # alias if named p_mean
        "rho_mean"  => "ρ̄  [kg/m³]",
        "Press_mean"=> "P̄  [Pa]",
        "upup_res"  => "u'u'  [m²/s²]",
        "upvp_res"  => "u'v'  [m²/s²]",
        "upwp_res"  => "u'w'  [m²/s²]",
        "vpvp_res"  => "v'v'  [m²/s²]",
        "vpwp_res"  => "v'w'  [m²/s²]",
        "wpwp_res"  => "w'w'  [m²/s²]",
        "tptp_res"  => "θ'θ'  [K²]",
        "uptp_res"  => "u'θ'  [m·K/s]",
        "vptp_res"  => "v'θ'  [m·K/s]",
        "wptp_res"  => "w'θ'  [m·K/s]",
        "upup_sfs"  => "τ_uu  [m²/s²]",
        "upvp_sfs"  => "τ_uv  [m²/s²]",
        "upwp_sfs"  => "τ_uw  [m²/s²]",
        "vpvp_sfs"  => "τ_vv  [m²/s²]",
        "vpwp_sfs"  => "τ_vw  [m²/s²]",
        "wpwp_sfs"  => "τ_ww  [m²/s²]",
        "tptp_sfs"  => "τ_θθ  [K²]",
        "uptp_sfs"  => "τ_uθ  [m·K/s]",
        "vptp_sfs"  => "τ_vθ  [m·K/s]",
        "wptp_sfs"  => "τ_wθ  [m·K/s]",
    )
    get(labels, name, name)
end

"""
    plot_profiles(varnames, z, data, title, outpath; ncols=3)

Multi-panel figure: one panel per variable, variable on x-axis, z on y-axis.
"""
function plot_profiles(varnames, z, data, fig_title, outpath; ncols=3)
    nvars = length(varnames)
    nrows = cld(nvars, ncols)
    zmin, zmax = extrema(z)

    fig = Figure(size=(380 * ncols, 420 * nrows))
    Label(fig[0, :], fig_title; fontsize=16, font=:bold)

    for (k, vname) in enumerate(varnames)
        row = cld(k, ncols)
        col = mod1(k, ncols)
        ax  = Axis(fig[row, col];
                   xlabel      = nice_xlabel(vname),
                   ylabel      = "z  [m]",
                   title       = vname,
                   titlesize   = 12,
                   xlabelsize  = 11,
                   ylabelsize  = 11,
                   xticklabelsize = 9,
                   yticklabelsize = 9)
        lines!(ax, data[:, k], z; linewidth=2, color=:royalblue)
        scatter!(ax, data[:, k], z; markersize=4, color=:royalblue)
        ylims!(ax, zmin, zmax)
    end

    save(outpath, fig; px_per_unit=2)
    println("  Saved: $outpath")
end

"""
    save_individual_profiles(varnames, z, data, outdir; prefix="")

Save one PNG per variable into `outdir`. File names are `prefix_varname.png`.
"""
function save_individual_profiles(varnames, z, data, outdir; prefix="")
    mkpath(outdir)
    zmin, zmax = extrema(z)

    for (k, vname) in enumerate(varnames)
        fig = Figure(size=(480, 520))
        ax  = Axis(fig[1, 1];
                   xlabel         = nice_xlabel(vname),
                   ylabel         = "z  [m]",
                   title          = vname,
                   titlesize      = 14,
                   xlabelsize     = 12,
                   ylabelsize     = 12,
                   xticklabelsize = 10,
                   yticklabelsize = 10)
        lines!(ax, data[:, k], z; linewidth=2, color=:royalblue)
        scatter!(ax, data[:, k], z; markersize=4, color=:royalblue)
        ylims!(ax, zmin, zmax)

        fname = isempty(prefix) ? "$vname.png" : "$(prefix)_$vname.png"
        save(joinpath(outdir, fname), fig; px_per_unit=2)
        println("  Saved: $(joinpath(outdir, fname))")
    end
end

# ──────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────

function main()
    dir = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "output",
                                                  "CompEuler", "LESICP2-coarse", "output")
    dir = abspath(dir)
    isdir(dir) || error("Directory not found: $dir")

    plots_dir = joinpath(dir, "les_plots")

    # ── Mean profiles ──────────────────────────────────────────────
    stat_file = joinpath(dir, "les_statistics_tavg.dat")
    if isfile(stat_file)
        varnames, z, data = read_les_tavg(stat_file)
        println("les_statistics_tavg: $(length(varnames)) vars, $(length(z)) z-levels")
        out = joinpath(dir, "les_profiles_tavg.png")
        plot_profiles(varnames, z, data, "LES mean profiles (time-averaged)", out; ncols=3)
        save_individual_profiles(varnames, z, data, plots_dir; prefix="mean")
    else
        @warn "Not found: $stat_file"
    end

    # ── Stress components ──────────────────────────────────────────
    stress_file = joinpath(dir, "les_stress_tavg.dat")
    if isfile(stress_file)
        varnames, z, data = read_les_tavg(stress_file)
        println("les_stress_tavg: $(length(varnames)) vars, $(length(z)) z-levels")

        # Split into resolved and SFS groups
        res_idx   = findall(v -> endswith(v, "_res"), varnames)
        sfs_idx   = findall(v -> endswith(v, "_sfs"), varnames)
        other_idx = setdiff(1:length(varnames), res_idx, sfs_idx)

        if !isempty(res_idx)
            out = joinpath(dir, "les_stress_resolved_tavg.png")
            plot_profiles(varnames[res_idx], z, data[:, res_idx],
                          "Resolved Reynolds stresses (time-averaged)", out; ncols=3)
            save_individual_profiles(varnames[res_idx], z, data[:, res_idx], plots_dir; prefix="res")
        end
        if !isempty(sfs_idx)
            out = joinpath(dir, "les_stress_sfs_tavg.png")
            plot_profiles(varnames[sfs_idx], z, data[:, sfs_idx],
                          "SFS stresses (time-averaged)", out; ncols=3)
            save_individual_profiles(varnames[sfs_idx], z, data[:, sfs_idx], plots_dir; prefix="sfs")
        end
        if !isempty(other_idx)
            out = joinpath(dir, "les_stress_other_tavg.png")
            plot_profiles(varnames[other_idx], z, data[:, other_idx],
                          "Other stress variables (time-averaged)", out; ncols=3)
            save_individual_profiles(varnames[other_idx], z, data[:, other_idx], plots_dir; prefix="other")
        end
    else
        @warn "Not found: $stress_file"
    end

    println("Done.")
end

main()
