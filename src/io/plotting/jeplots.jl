using Plots #; plotlyjs()
using LaTeXStrings
using ColorSchemes
using CairoMakie
using GLMakie

include("../../kernel/mesh/mesh.jl")

function plot_curve(x, y,  title::String, fout_name::String)
    
    default(titlefont=(14, "Arial, sans-serif"),
            legendfontsize = 18,
            guidefont = (18, :darkgreen),
            tickfont = (12, :orange),
            guide = "x",
            framestyle = :zerolines, yminorgrid = true)
    
    data = scatter(x, y, title=title,
                   markersize = 5, markercolor="Blue",
                   xlabel = "x", ylabel = "q(x)",
                   legend = :none)
    
    Plots.savefig(data, fout_name)
    
end


function plot_surf()

    x = range(-3, 3, length=30)
    fig = surface(
        x, x, (x, y)->exp(-x^2 - y^2), c=:viridis, legend=:none,
        nx=50, ny=50, display_option=Plots.GR.OPTION_SHADED_MESH,  # <-- series[:extra_kwargs]
    )

    display(fig)
    
end

function plot_1d_grid(mesh::St_mesh)
    
    x = mesh.x
    npoin = length(x)
    
    plt = plot() #Clear plot
    for i=1:npoin
        display(scatter!(x, zeros(npoin), markersizes=4))
    end 
end

#
# Curves (1D) or Contours (2D) with PlotlyJS
#
function plot_results(SD::NSD_1D, x1, _, z1, title::String, fout_name::String)
    
    default(titlefont=(14, "Arial, sans-serif"),
            legendfontsize = 18,
            guidefont = (18, :darkgreen),
            tickfont = (12, :orange),
            guide = "x",
            framestyle = :zerolines, yminorgrid = true)
    
    data = scatter(x1, z1, title=title,
                   markersize = 5, markercolor="Blue",
                   xlabel = "x", ylabel = "q(x)",
                   legend = :none)
    
    Plots.savefig(data, fout_name)
    
end


"""
    This function uses the amazing package Mackie to plot arbitrarily gridded
    unstructured data to filled contour plot
"""
function plot_triangulation(SD::NSD_2D, x, y, q, title::String, fout_name::String)
    
    fig = Makie.tricontourf(x, y, q, colormap = :heat)
    save(fout_name, fig, resolution = (600, 600))
    fig
    
end
