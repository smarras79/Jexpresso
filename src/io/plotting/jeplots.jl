#using Plots; gr()
using PlotlyJS
using LaTeXStrings
using ColorSchemes

include("../../kernel/mesh/mesh.jl")

function plot_error(x, y, title::String, legend_labels; yscale)

    xlabel = "nop"
    ylabel = "log10(ε)"
    
    plot_curve!(x, y, title::String, xlabel, ylabel, legend_labels, yscale)
    
end

function plot_curve(ξ, ψ, title)
    
    plt = plot() #Clear plot
    display(plot!(ξ, ψ, title = title, legend=false, lw = 3,
                 xtickfontsize=16, ytickfontsize=16, reuse=false,
                 xlabel="ξ",ylabel="ψ(ξ)"))
end

function scatter_curve(ξ, ψ, title)
    
    plt = plot() #Clear plot
    display(scatter!(ξ, ψ, title = title, legend=false, lw = 3,
                     xtickfontsize=16, ytickfontsize=16, reuse=false,
                     xlabel="ξ",ylabel="ψ(ξ)"))
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
# Contours with PlotlyJS
#
function jcontour(x1, y1, z1, title, fout_name)
    
    data = PlotlyJS.contour(;z=z1, x=x1, y=y1,
                            colorbar=attr(;title="",titleside="right",
                                          titlefont=attr(;size=14,
                                                         family="Arial, sans-serif")
                                          )
                            )

    layout = Layout(;title=title)
    #display(PlotlyJS.plot(data, layout)) #WARNING aspect ratio doesn't seem to work
    PlotlyJS.savefig(PlotlyJS.plot(data, layout), fout_name)
    
end
