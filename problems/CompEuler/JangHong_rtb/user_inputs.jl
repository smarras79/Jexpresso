function user_inputs()
    inputs = Dict(
        :ode_solver           => CarpenterKennedy2N54(),
        :Δt                   => 0.1,
        :tinit                => 0.0,
        :tend                 => 1020.0,
        :diagnostics_at_times => (0:60:1020),
        :case                 => "rtb",
        :lsource              => true,
        :interpolation_nodes  => "lgl",
        :nop                  => 4,
        :lvisc                => true,
       # :μ                    => [0.0, 25.0, 25.0, 25.0, 0.0, 0.0],
        :μ                    => [0.0, 50.0, 50.0, 50.0, 0.0, 0.0], #default 
        :lread_gmsh           => true,
        :gmsh_filename        => "./meshes/gmsh_grids/janghong_rtb.msh",
        :outformat            => "vtk",
        :loverwrite_output    => true,
        :output_dir           => "./janghong_rtb_50",
        :loutput_pert         => true,
    )
    return inputs
end