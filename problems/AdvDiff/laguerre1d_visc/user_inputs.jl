#
# 1D advection-diffusion into a semi-infinite (Laguerre) layer on the right.
# A Gaussian (x0 = 27, sigma = 3) is advected at speed 1 with nu = 0.5; at
# tend = 50 it is centred at x = 77, inside the layer (x > 60). qe holds the
# exact solution at tend (see initialize.jl), so the output carries both.
# Layer nodes: x = xmax + yfac_laguerre*xi; with ngr = 31 the layer reaches
# x ~ 115, and a finer layer (smaller yfac) needs a smaller dt (explicit
# diffusion).
#
function user_inputs()
    inputs = Dict(
        :ode_solver          => SSPRK33(),
        :tend                => 50.0,
        :Δt                  => 0.002,
        :ndiagnostics_outputs => 1,
        :output_dir          => "./output/",
        :SOL_VARS_TYPE       => PERT(),
        :interpolation_nodes => "lgl",
        :nop                 => 4,
        :nop_laguerre        => 30,
        :lexact_integration  => false,
        :lsource             => false,
        :llaguerre_1d_right  => true,
        :laguerre_beta       => 1.0,
        :yfac_laguerre       => 0.5,
        :lvisc               => true,
        :μ                   => [0.5],
        :lread_gmsh          => false,
        :outformat           => "png",
        :xmin                => 0.0,
        :xmax                => 60.0,
        :nelx                => 40,
    )
    return inputs
end
