using Dierckx
using Plots
using LaTeXStrings
using ColorSchemes
#using Interpolations

#
# Curves (1D) or Contours (2D) with Plots.jl
#

function plot_initial(SD::NSD_1D, x, q, ivar, OUTPUT_DIR::String)

    npoin = length(q)
    plt = Plots.scatter(x[1:npoin], q[1:npoin];
                        markersize = 5,
                        color = :blue,
                        xlabel = "x",
                        ylabel = "q(x)",
                        title = "u",
                        titlefontsize = 24,
                        guidefontsize = 18,
                        legendfontsize = 14,
                        tickfontsize = 14,
                        legend = false,
                        size = (800, 600))

    fout_name = string(OUTPUT_DIR, "/INIT-", ivar, ".png")

    Plots.savefig(plt, string(fout_name))
    plt
end

function plot_results(SD::NSD_1D, mesh::St_mesh, q, title::String, OUTPUT_DIR::String, outvar, inputs::Dict; iout=1, nvar=1, PT=nothing)

    epsi = 1.1
    npoin = mesh.npoin

    qout = reshape(q, npoin, nvar)   # 2D view (3x4) - NO allocation
    x_coords = mesh.x[1:npoin]
    sort_idx = sortperm(x_coords)
    for ivar=1:nvar

        idx = (ivar - 1)*npoin

        plt = Plots.plot(x_coords[sort_idx], qout[sort_idx, ivar];
                        line = (:blue, 2),
                        marker = (:circle, 5, :blue),
                        title = string(outvar[ivar]),
                        xlabel = "x",
                        titlefontsize = 22,
                        guidefontsize = 18,
                        legendfontsize = 14,
                        tickfontsize = 14,
                        legend = false,
                        size = (600, 400))

        vlines = inputs[:plot_vlines]
        hlines = inputs[:plot_hlines]
        axis = inputs[:plot_axis]
        if !(vlines == "empty")
            for i=1:size(vlines,1)
                Plots.vline!(plt, [vlines[i]]; color = :red, linestyle = :solid, label = "")
            end
        end
        if !(hlines == "empty")
            for i=1:size(hlines,1)
                Plots.hline!(plt, [hlines[i]]; color = :red, linestyle = :solid, label = "")
            end
        end
        if !(axis == "empty")
            idx = (ivar-1)*2
            Plots.ylims!(plt, axis[1+idx], axis[2+idx])
        end
        fout_name = string(OUTPUT_DIR, "/ivar", ivar, "-it", iout, ".png")
        Plots.savefig(plt, string(fout_name))
        plt
    end
end


function plot_results!(SD::NSD_1D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String, outvar, inputs::Dict; iout=1, nvar=1, fig=nothing, color ="blue", p=[], marker = :circle, PT=nothing)

    @print "C"

    epsi = 1.1
    npoin = mesh.npoin

    for ivar=1:1
        idx = (ivar - 1)*npoin

        if fig === nothing
            fig = Plots.plot(xlabel = "x",
                           title = string(outvar[ivar]),
                           titlefontsize = 18,
                           guidefontsize = 14,
                           legend = false)
        end

        if !(p==[])
            # Add to existing plot without decorations
            Plots.scatter!(fig, mesh.x[1:mesh.npoin_original],
                         q[idx+1:(ivar-1)*npoin+mesh.npoin_original];
                         marker = marker,
                         markersize = 5,
                         color = color,
                         label = "")
        else
            Plots.scatter!(fig, mesh.x[1:mesh.npoin_original],
                         q[idx+1:(ivar-1)*npoin+mesh.npoin_original];
                         marker = marker,
                         markersize = 5,
                         color = color,
                         label = "")
        end

        Plots.ylims!(fig, -0.03, 0.03)
        fout_name = string(OUTPUT_DIR, "/ivar", ivar, "-it", iout, ".eps")
        Plots.savefig(fig, string(fout_name))
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
        plt = Plots.scatter(mesh.x[1:npoin], qout[idx+1:ivar*npoin];
                           markersize = 5,
                           color = :blue,
                           xlabel = "x",
                           ylabel = "q(x)",
                           title = string(outvar[ivar]),
                           titlefontsize = 24,
                           guidefontsize = 18,
                           legendfontsize = 14,
                           tickfontsize = 14,
                           legend = false,
                           size = (800, 600))


        fout_name = string(OUTPUT_DIR, "/ivar", ivar, "-it", iout, ".png")
        Plots.savefig(plt, string(fout_name))
        plt
    end
end

function plot_1d_grid(mesh::St_mesh)

    plt = Plots.plot() #Clear plot
    for i=1:mesh.npoin
        display(Plots.scatter(mesh.x[1:mesh.npoin], zeros(mesh.npoin),
                             markersize = 4,
                             color = :blue,
                             legend = false))
    end
end


function plot_initial(SD::NSD_2D, x::Array, q::Array, ivar, OUTPUT_DIR::String)
    nothing
end

function plot_triangulation(SD::NSD_2D, mesh::St_mesh, q::Array, title::String, OUTPUT_DIR::String, inputs::Dict; iout=1, nvar=1)

    """
        This function uses Plots.jl to plot arbitrarily gridded
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

        # Create triangulated contour plot using Plots
        plt = Plots.scatter(mesh.x[1:npoin], mesh.y[1:npoin];
                           zcolor = q[idx+1:ivar*npoin],
                           color = :viridis,
                           marker = :circle,
                           markersize = 3,
                           markerstrokewidth = 0,
                           colorbar = true,
                           legend = false,
                           aspect_ratio = :equal,
                           size = (800, 600))

        minq = minimum(q[idx+1:ivar*npoin])
        maxq = maximum(q[idx+1:ivar*npoin])

        if (maxq > minq)
            Lx = abs(maximum(mesh.x) - minimum(mesh.x))
            Ly = abs(maximum(mesh.y) - minimum(mesh.y))
            vlines = inputs[:plot_vlines]
            hlines = inputs[:plot_hlines]
            if !(vlines == "empty")
                for i=1:size(vlines,1)
                    Plots.vline!(plt, [vlines[i]]; color = :red, linestyle = :dash, label = "")
                end
            end
            if !(hlines == "empty")
                for i=1:size(hlines,1)
                    Plots.hline!(plt, [hlines[i]]; color = :red, linestyle = :dash, label = "")
                end
            end

            Plots.clims!(plt, minq, maxq)
            Plots.savefig(plt, string(fout_name))
            plt
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
        plt = Plots.surface(xg, yg, zspl;
                           color = :viridis,
                           camera = (0, 90),  # Top-down view similar to Makie's azimuth/elevation
                           colorbar = true,
                           legend = false,
                           size = (1200, 400))

        Plots.savefig(plt, string(fout_name))
        plt
    end

end
