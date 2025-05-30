using Dierckx
using LaTeXStrings
using ColorSchemes
using CairoMakie
using Makie
#using Interpolations

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
journal = {Journal of Open SOurce Software}
}
=#

#
# Curves (1D) or Contours (2D) with PlotlyJS
#

function plot_initial(SD::NSD_1D, x, q, ivar, OUTPUT_DIR::String)
    
    npoin = length(q)
    fig, ax, plt = CairoMakie.scatter(x[1:npoin], q[1:npoin];
                                      markersize = 10, color="Blue",
                                      xlabel = "x", ylabel = "q(x)",
                                      fontsize = 24, fonts = (; regular = "Dejavu", weird = "Blackchancery"),
                                      axis = (; title = "u", xlabel = "x")
                                      )
    
    fout_name = string(OUTPUT_DIR, "/INIT-", ivar, ".png")
    
    save(string(fout_name), fig)
    fig
end

function plot_results(SD::NSD_1D, mesh::St_mesh, q, title::String, OUTPUT_DIR::String, outvar, inputs::Dict; iout=1, nvar=1, PT=nothing)
    
    epsi = 1.1
    npoin = mesh.npoin
        
    for ivar=1:nvar
        
        idx = (ivar - 1)*npoin
        
        CairoMakie.activate!(type = "eps")
        fig = Figure(size = (600,400),fontsize=22)
        ax = Axis(fig[1, 1], title=string(outvar[ivar]), xlabel="x")
        CairoMakie.scatter!(mesh.x[1:npoin], q[idx+1:ivar*npoin]; markersize = 10, color="Blue")
        vlines = inputs[:plot_vlines]
        hlines = inputs[:plot_hlines]
        axis = inputs[:plot_axis]
        if !(vlines == "empty")
            for i=1:size(vlines,1) 
                #vlines!(ax, [-2.5,2.5], color = :red)
                vlines!(ax,vlines[i], color = :red)
            end
        end
        if !(hlines == "empty")
            for i=1:size(hlines,1)
                #vlines!(ax, [-2.5,2.5], color = :red)
                hlines!(ax,hlines[i], color = :red)
            end
        end
        if !(axis == "empty")
            idx = (ivar-1)*2
            ylims!(ax, axis[1+idx], axis[2+idx])
        end
        fout_name = string(OUTPUT_DIR, "/ivar", ivar, "-it", iout, ".png")        
        save(string(fout_name), fig)
        fig
    end
end


function plot_results!(SD::NSD_1D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String, outvar, inputs::Dict; iout=1, nvar=1, fig=Figure(),color ="Blue",p=[],marker = :circle, PT=nothing)
    
    epsi = 1.1
    npoin = mesh.npoin
    
    for ivar=1:1
        idx = (ivar - 1)*npoin
        CairoMakie.activate!(type = "eps")
        if !(p==[]) 
            ax = Axis(fig[1, 1], title="", xlabel="")
            hidedecorations!(ax)
            push!(p,CairoMakie.scatter!(mesh.x[1:mesh.npoin_original], q[idx+1:(ivar-1)*npoin+mesh.npoin_original];marker = marker, markersize = 10, color=color))
        else
            ax = Axis(fig[1, 1], title=string(outvar[ivar]), xlabel="x")
            push!(p,CairoMakie.scatter!(mesh.x[1:mesh.npoin_original], q[idx+1:(ivar-1)*npoin+mesh.npoin_original];marker = marker, markersize = 10, color=color))
        end
        p[end].color = color
        ylims!(ax, -0.03, 0.03)
        fout_name = string(OUTPUT_DIR, "/ivar", ivar, "-it", iout, ".eps")  
        save(string(fout_name), fig)
        fig
    end
end


function plot_results(SD::NSD_1D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String, varnames; iout=1, nvar=1, PT=nothing)
    
    epsi = 1.1
    npoin = mesh.npoin
    
    qout = copy(q)
    qe   = range(0,0,npoin)

    #outvar = ["ρ", "u", "e"]
    outvar = varnames
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
                                          markersize = 10, color="Blue",
                                          xlabel = "x", ylabel = "q(x)",
                                          fontsize = 24, fonts = (; regular = "Dejavu", weird = "Blackchancery"),  axis = (; title = string(outvar[ivar]), xlabel = "xx")
                                          )

        
        fout_name = string(OUTPUT_DIR, "/ivar", ivar, "-it", iout, ".png")        
        save(string(fout_name), fig)
        fig
    end
end

function plot_1d_grid(mesh::St_mesh)
    
    plt = plot() #Clear plot
    for i=1:mesh.npoin
        display(CairoMakie.scatter(mesh.x[1:mesh.npoin], zeros(mesh.npoin), markersizes=4, color="Blue"))
    end 
end


function plot_initial(SD::NSD_2D, x::Array, q::Array, ivar, OUTPUT_DIR::String)
    nothing
end

function plot_triangulation(SD::NSD_2D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String, inputs::Dict; iout=1, nvar=1)

    """
        This function uses the amazing package Mackie to plot arbitrarily gridded
        unstructured data to filled contour plot
    """
    
    if ("Laguerre" in mesh.bdy_edge_type)
        npoin = mesh.npoin_original
    else
        npoin = floor(Int64, size(q, 1)/nvar)
    end
    npoin = mesh.npoin
    for ivar=1:nvar
        idx = (ivar - 1)*npoin
        fout_name = string(OUTPUT_DIR, "/ivar", ivar, "-it", iout, ".png")
        fig, ax, sol = Makie.tricontourf(mesh.x[1:npoin], mesh.y[1:npoin], q[idx+1:ivar*npoin], colormap = :viridis)
        
        minq = minimum(q[idx+1:ivar*npoin])
        maxq = maximum(q[idx+1:ivar*npoin])

        if (maxq > minq) 
            Lx = abs(maximum(mesh.x) - minimum(mesh.x))
            Ly = abs(maximum(mesh.y) - minimum(mesh.y))
            vlines = inputs[:plot_vlines]
            hlines = inputs[:plot_hlines]
            if !(vlines == "empty")
                for i=1:size(vlines,1)
                    vlines!(ax,vlines[i], color = :red, linestyle = :dash)
                end
            end
            if !(hlines == "empty")
                for i=1:size(hlines,1)
                    hlines!(ax,hlines[i], color = :red, linestyle = :dash)
                end
            end
            
            if (Ly > Lx)
                ax.aspect = Lx/Ly; colsize!(fig.layout, 1, Aspect(1, Lx/Ly))
            else
                #ax.aspect = Lx/Ly; #colsize!(fig.layout, 1, Aspect(1, Lx/Ly))
            end      

            Colorbar(fig[1,2], colormap = :viridis,  limits = (minq, maxq))        
            save(string(fout_name), fig) #, size = (600, 600))
            fig
        end
        
    end
end
function plot_triangulation(SD::NSD_1D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String, inputs::Dict; nvar=1) nothing end
function plot_triangulation(SD::NSD_3D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String, inputs::Dict; nvar=1) nothing end

function plot_surf3d(SD::NSD_2D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String; iout=1, nvar=1, smoothing_factor=1e-3)

    xmin = minimum(mesh.x); xmax = maximum(mesh.x);
    ymin = minimum(mesh.y); ymax = maximum(mesh.y); 

    nxi = 500
    nyi = 500
    npoin = mesh.npoin
    for ivar=1:nvar
        idx = (ivar - 1)*npoin
        
        fout_name = string(OUTPUT_DIR, "/ivar", ivar, "-it", iout, ".png")
        
        #Spline2d
        spl = Spline2D(mesh.x[1:npoin], mesh.y[1:npoin], q[idx+1:ivar*npoin]; kx=4, ky=4, s=smoothing_factor)
        xg = LinRange(xmin, xmax, nxi); yg = LinRange(ymin, ymax, nyi);
        zspl = evalgrid(spl, xg, yg);
        #End spline2d

        #figure:
        fig = Figure(size=(1200, 400))
        axs = [Axis3(fig[1, i]; aspect=(1, 1, 1), azimuth=-π/2, elevation=π/2) for i = 1:1]
        
        hm = Makie.surface!(axs[1], xg, yg, zspl) # xl="x", yl="y", zl=string("q", ivar)) #, title=title, titlefont=12)
        #Colorbar(fig[1, 1], hm, height=Relative(0.5))
        
        save(string(fout_name), fig)
        #display(fig)
        fig
    end
    
end
