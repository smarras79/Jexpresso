using Plots
using Dierckx
using LaTeXStrings
using ColorSchemes
using CairoMakie
using Makie
Makie.theme(:fonts)

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

#
# Curves (1D) or Contours (2D) with PlotlyJS
#
function plot_results(SD::NSD_1D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String; iout=1, nvar=1, PT=nothing)

    xmin = minimum(mesh.x); xmax = maximum(mesh.x);
    qmin = minimum(q);      qmax = maximum(q);
    epsi = 1.1
    npoin = floor(Int64, size(q, 1)/nvar)
   
    qout = copy(q)
    qe   = range(0,0,npoin)

    #outvar = ["ρ", "u", "e"]
    outvar = ["ρ", "u", "p"]
    if PT === CompEuler()
        #ρ
        qout[1:npoin] .= @view q[1:npoin]

        #u = ρu/ρ
        ivar = 2
        idx = (ivar - 1)*npoin
        qout[idx+1:2*npoin] .= q[idx+1:2*npoin]./q[1:npoin]
        
        ivar = 3
        idx = (ivar - 1)*npoin
        γ = 1.4

        if (outvar[3] == "e")
            #Internal energy
            qout[idx+1:3*npoin] .= (q[2*npoin+1:3*npoin] .- 0.5*q[npoin+1:2*npoin].*q[npoin+1:2*npoin]./q[1:npoin])./q[1:npoin] #internal energy: p/((γ-1)ρ)
        elseif (outvar[3] == "p")
            #Pressure
            qout[idx+1:3*npoin] .= (γ - 1.0)*(q[2*npoin+1:3*npoin] .- 0.5*q[npoin+1:2*npoin].*q[npoin+1:2*npoin]./q[1:npoin]) #Pressure
        end
    elseif PT === ShallowWater()
        Hb = zeros(npoin,nvar)
        for i=1:npoin
            x= mesh.x[i]
            Hb[i,1] = zb[i]
        end

        for ivar=1:nvar
            idx = (ivar - 1)*npoin
            qout[idx+1:ivar*npoin] .= q[idx+1:ivar*npoin] .+ Hb[:,ivar]
        end
    end
    
    for ivar=1:nvar

        idx = (ivar - 1)*npoin
        fig, ax, plt = CairoMakie.scatter(mesh.x[1:npoin], qout[idx+1:ivar*npoin];
                                          markersize = 10, markercolor="Blue",
                                          xlabel = "x", ylabel = "q(x)",
                                          fontsize = 24, fonts = (; regular = "Dejavu", weird = "Blackchancery"),  axis = (; title = string(outvar[ivar]), xlabel = "xx")
                                          )

        
        fout_name = string(OUTPUT_DIR, "/ivar", ivar, "-it", iout, ".png")        
        save(string(fout_name), fig; resolution = (600, 400))
        fig
    end
end

function plot_1d_grid(mesh::St_mesh)
    
    plt = plot() #Clear plot
    for i=1:mesh.npoin
        display(CairoMakie.scatter(mesh.x[1:mesh.npoin], zeros(mesh.npoin), markersizes=4, markercolor="Blue"))
    end 
end

function plot_triangulation(SD::NSD_2D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String; iout=1, nvar=1)

"""
    This function uses the amazing package Mackie to plot arbitrarily gridded
    unstructured data to filled contour plot
"""
    npoin = floor(Int64, size(q, 1)/nvar)
    for ivar=1:nvar
        idx = (ivar - 1)*npoin
        
        fout_name = string(OUTPUT_DIR, "/ivar", ivar, "-it", iout, ".png")
        fig, ax, sol = Makie.tricontourf(mesh.x[1:npoin], mesh.y[1:npoin], q[idx+1:ivar*npoin], colormap = :viridis)
        Colorbar(fig[1,2], colormap = :viridis)        
        save(string(fout_name), fig) #, resolution = (600, 600))
        fig
    end
end
function plot_triangulation(SD::NSD_1D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String; nvar=1) nothing end
function plot_triangulation(SD::NSD_3D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String; nvar=1) nothing end

function plot_surf3d(SD::NSD_2D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String; iout=1, nvar=1, smoothing_factor=1e-1)

    xmin = minimum(mesh.x); xmax = maximum(mesh.x);
    ymin = minimum(mesh.y); ymax = maximum(mesh.y); 

    nxi = 500
    nyi = 500
    npoin = floor(Int64, size(q, 1)/nvar)
    for ivar=1:nvar
        idx = (ivar - 1)*npoin
        
        fout_name = string(OUTPUT_DIR, "/ivar", ivar, "-it", iout, ".png")
        
        #Spline2d
        spl = Spline2D(mesh.x[1:npoin], mesh.y[1:npoin], q[idx+1:ivar*npoin]; kx=4, ky=4, s=smoothing_factor)
        xg = LinRange(xmin, xmax, nxi); yg = LinRange(ymin, ymax, nyi);
        zspl = evalgrid(spl, xg, yg);
        #End spline2d

        #figure:
        fig = Figure(resolution=(1200, 400))
        axs = [Axis3(fig[1, i]; aspect=(1, 1, 1)) for i = 1:1]
        
        hm = Makie.surface!(axs[1], xg, yg, zspl) # legend=:false, xl="x", yl="y", zl=string("q", ivar)) #, title=title, titlefont=12)
        #Colorbar(fig[1, 1], hm, height=Relative(0.5))
        
        save(string(fout_name), fig)
        #display(fig)
        fig
    end
    
end

#
# ASCII
#
function write_ascii(SD::NSD_1D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String; iout=1, nvar=1)
    
    npoin = floor(Int64, size(q, 1)/nvar)
    for ivar=1:nvar
        idx = (ivar - 1)*npoin
        
        fout_name = string(OUTPUT_DIR, "/ivar", ivar, "-it", iout, ".dat")
        open(fout_name, "w") do f
            @printf(f, " %f %f \n", mesh.x[1:npoin,1], q[idx+1:ivar*npoin])
        end
    end
end

function write_ascii(SD::NSD_2D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String; iout=1, nvar=1)
    
    npoin = floor(Int64, size(q, 1)/nvar)
    for ivar=1:nvar
        idx = (ivar - 1)*npoin
        
        fout_name = string(OUTPUT_DIR, "/ivar", ivar, "-it", iout, ".dat")
        open(fout_name, "w") do f
            @printf(f, " %f %f \n", mesh.x[1:npoin,1], mesh.y[1:npoin,1], q[idx+1:ivar*npoin])
        end
    end
end

function write_ascii(SD::NSD_3D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String; iout=1, nvar=1)
    
    npoin = floor(Int64, size(q, 1)/nvar)
    for ivar=1:nvar
        idx = (ivar - 1)*npoin
        
        fout_name = string(OUTPUT_DIR, "/ivar", ivar, "-it", iout, ".dat")
        open(fout_name, "w") do f
            @printf(f, " %f %f \n", mesh.x[1:npoin,1], mesh.y[1:npoin,1], mesh.z[1:npoin,1], q[idx+1:ivar*npoin])
        end
    end
end
