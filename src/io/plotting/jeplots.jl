using Plots
using LaTeXStrings
using ColorSchemes
using CairoMakie
using Makie

#= CITE Mackie:
@article{DanischKrumbiegel2021,
  doi = {10.21105/joss.03349},
  url = {https://doi.org/10.21105/joss.03349},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {65},
  pages = {3349},
  author = {Simon Danisch and Julius Krumbiegel},
  title = {Makie.jl: Flexible high-performance data visualization for Julia},
  journal = {Journal of Open Source Software}
}
=#

include("../../kernel/mesh/mesh.jl")

function plot_curve(x, y,  title::String, fout_name::String)
    
    default(titlefont=(14, "Arial, sans-serif"),
            legendfontsize = 18,
            guidefont = (18, :darkgreen),
            tickfont = (12, :orange),
            guide = "x",
            framestyle = :zerolines, yminorgrid = true)
    
    data = Plots.scatter(x, y, title=title,
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
        display(Plots.scatter!(x, zeros(npoin), markersizes=4))
    end 
end

#
# Curves (1D) or Contours (2D) with PlotlyJS
#
function plot_results(SD::NSD_1D, x1, y1, z1, title::String, fout_name::String)
    
    default(titlefont=(14, "Arial, sans-serif"),
            legendfontsize = 18,
            guidefont = (18, :darkgreen),
            tickfont = (12, :orange),
            guide = "x",
            framestyle = :zerolines, yminorgrid = true)
    
    data = Plots.scatter(x1, z1, title=title,
                         markersize = 5, markercolor="Blue",
                         xlabel = "x", ylabel = "q(x)",
                         legend = :none)
    
    Plots.savefig(data, fout_name)
    
end


"""
    This function uses the amazing package Mackie to plot arbitrarily gridded
    unstructured data to filled contour plot
"""
function plot_triangulation(SD::NSD_2D, x, y, q::Array, title::String, fout_name::String)
    
    #fddor ivar=1:size(q,2)
        fig = Makie.tricontourf(x, y, q[:,1], colormap = :heat)
        save(string(fout_name), fig, resolution = (600, 600))
        fig
    #end
end
function plot_triangulation(SD::NSD_1D, x, y, q::Array, title::String, fout_name::String) nothing end
function plot_triangulation(SD::NSD_3D, x, y, q::Array, title::String, fout_name::String) nothing end
