using Plots
using Dierckx
using LaTeXStrings
using ColorSchemes
using CairoMakie
using Makie
using ImageMagick

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

function plot_curve(x, q::Array, title::String, fout_name::String)
    
    default(titlefont=(14, "Arial, sans-serif"),
            legendfontsize = 18,
            guidefont = (18, :darkgreen),
            tickfont = (12, :orange),
            guide = "x",
            framestyle = :zerolines, yminorgrid = true)
    
    data = Plots.scatter(x, q, title=title,
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
    
    plt = plot() #Clear plot
    for i=1:mesh.npoin
        display(Plots.scatter!(mesh.x[1:mesh.npoin], zeros(mesh.npoin), markersizes=4))
    end 
end

#
# Curves (1D) or Contours (2D) with PlotlyJS
#
function plot_results(SD::NSD_1D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String; iout=1, nvar=1)
    
    default(titlefont=(14, "Arial, sans-serif"),
            legendfontsize = 18,
            guidefont = (18, :darkgreen),
            tickfont = (12, :orange),
            guide = "x",
            framestyle = :zerolines, yminorgrid = true)

    npoin = mesh.npoin
    for ivar=1:nvar
        idx = (ivar - 1)*npoin
        data = Plots.scatter(mesh.x[1:npoin], q[idx+1:ivar*npoin], title=title,
                             markersize = 5, markercolor="Blue",
                             xlabel = "x", ylabel = "q(x)",
                             legend = :none)
        
        fout_name = string(OUTPUT_DIR, "/ivar", ivar, "-it", iout, ".png")
        Plots.savefig(data, fout_name)
    end
end


function plot_triangulation(SD::NSD_2D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String; iout=1, nvar=1)

"""
    This function uses the amazing package Mackie to plot arbitrarily gridded
    unstructured data to filled contour plot
"""    
    npoin = mesh.npoin
    for ivar=1:nvar
        idx = (ivar - 1)*npoin
        
        fout_name = string(OUTPUT_DIR, "/ivar", ivar, "-it", iout, ".png")
        fig, ax, sol = Makie.tricontourf(mesh.x[1:npoin], mesh.y[1:npoin], q[idx+1:ivar*npoin], colormap = :viridis)
        Colorbar(fig[1,2], colormap = :viridis)        
        save(string(fout_name), fig, resolution = (600, 600))
        fig
    end
end
function plot_triangulation(SD::NSD_1D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String; nvar=1) nothing end
function plot_triangulation(SD::NSD_3D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String; nvar=1) nothing end


function write_ascii(SD::NSD_1D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String; iout=1, nvar=1)
    
    npoin = mesh.npoin
    for ivar=1:nvar
        idx = (ivar - 1)*npoin
        
        fout_name = string(OUTPUT_DIR, "/ivar", ivar, "-it", iout, ".dat")
        open(fout_name, "w") do f
            @printf(f, " %f %f \n", mesh.x[1:npoin,1], q[idx+1:ivar*npoin])
        end
    end
end

function plot_surf3d(SD::NSD_2D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String; iout=1, nvar=1)

    xmin = minimum(mesh.x); xmax = maximum(mesh.x);
    ymin = minimum(mesh.y); ymax = maximum(mesh.y); 

    nxi = 100
    nyi = 100
    
    npoin = mesh.npoin
    for ivar=1:nvar
        idx = (ivar - 1)*npoin
        
        fout_name = string(OUTPUT_DIR, "/ivar", ivar, "-it", iout, ".png")
        spl = Spline2D(mesh.x[1:npoin], mesh.y[1:npoin], q[idx+1:ivar*npoin]; kx=3, ky=3, s=1e-4)
        xg = LinRange(xmin, xmax, nxi); yg = LinRange(ymin, ymax, nyi);
        zspl = evalgrid(spl, xg, yg);
        
        fig = Plots.surface(xg, yg, zspl'; legend=:false, xl="x", yl="y", zl=string("q", ivar)) #, title=title, titlefont=12)
        
        save(string(fout_name), fig)
        #display(fig)
    end
    
end
