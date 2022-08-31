using Plots; gr()
using LaTeXStrings
using ColorSchemes

include("../../mesh/mod_mesh.jl")

export plot_error
export plot_curve
export plot_1d_grid

function plot_error(x, y, title::String, legend_labels; yscale)

    xlabel = "nop"
    ylabel = "log10(ε)"
    
    plot_curve!(x, y, title::String, xlabel, ylabel, legend_labels, yscale)
    
end


function plot_curve!(x, y,
                title::String,
                xlabel::String,
                ylabel::String,
                legend_labels,
                yscale=:log10)

    # yscale is log10 by default.
    # use yscale=:none for linear
    
    fontsize = 16
    
    fig = plot(x, y, title = title, label = legend_labels, lw = 3)
    plot!(xlabel=xlabel)
    plot!(ylabel=ylabel; yscale=yscale)
    plot!(xtickfontsize=fontsize,
          ytickfontsize=fontsize, reuse=true)
    
    #fig = scatter!(x, y, linewidth=3, linestyle=:dot,markersize=3)
    #display(fig)
    
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
