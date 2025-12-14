# Example: Integrating Turbulence Statistics in Jexpresso
# This file demonstrates how to add turbulence statistics computation
# to your Jexpresso LES simulation

# Include the turbulence statistics module
include("../src/io/turbulence_statistics.jl")

"""
Example 1: Adding statistics to user_inputs.jl

Add these parameters to your user_inputs() Dict:
"""
function example_user_inputs_addition()
    statistics_params = Dict(
        # Enable turbulence statistics
        :l_turbulence_stats    => true,

        # Start time for statistics accumulation (after spin-up)
        :turb_stats_start_time => 100.0,

        # Statistics output frequency (in simulation time units)
        :turb_stats_output_dt  => 50.0,

        # Statistics output times
        :turb_stats_output_times => (100.0, 200.0, 500.0, 1000.0, 2000.0),

        # Reset statistics after each output (for windowed averaging)
        :turb_stats_reset_after_output => false,
    )

    return statistics_params
end

"""
Example 2: Allocating turbulence statistics in main solver setup

Add to your solver initialization (typically in the main simulation file):
"""
function example_allocation(inputs, mesh, SD, backend)

    # Get parameters
    npoin = mesh.npoin
    T = Float64  # Precision type
    l_turbulence_stats = get(inputs, :l_turbulence_stats, false)

    # Allocate turbulence statistics structure
    turb_stats = allocate_turbulence_stats(SD, npoin, T, backend;
                                          l_turbulence_stats=l_turbulence_stats)

    return turb_stats
end

"""
Example 3: Initializing statistics at start time

Add to your time integration setup:
"""
function example_initialization(inputs, turb_stats, integrator)

    if inputs[:l_turbulence_stats]
        t_start = get(inputs, :turb_stats_start_time, 0.0)

        # Initialize statistics when simulation time reaches t_start
        if integrator.t >= t_start
            initialize_turbulence_stats!(turb_stats, t_start)
        end
    end
end

"""
Example 4: Accumulating statistics during time stepping

Method A: Using DiscreteCallback (Recommended)
"""
function setup_statistics_callback(inputs, params, mesh, SD)

    if !inputs[:l_turbulence_stats]
        return nothing
    end

    # Define condition: accumulate every time step after start time
    function stats_condition(u, t, integrator)
        t_start = get(inputs, :turb_stats_start_time, 0.0)
        return t >= t_start
    end

    # Define action: accumulate statistics
    function stats_affect!(integrator)
        # Extract parameters
        turb_stats = integrator.p.turb_stats
        uaux = integrator.p.uaux
        press = integrator.p.q.press
        mueff = get(integrator.p, :mueff, nothing)  # If available
        dt = integrator.dt

        # Get dimension info
        npoin = mesh.npoin
        ndim = (SD == NSD_1D()) ? 1 : (SD == NSD_2D()) ? 2 : 3

        # Accumulate statistics
        accumulate_statistics!(turb_stats, uaux, press, mueff,
                              dt, npoin, ndim, SD)
    end

    # Create callback
    stats_cb = DiscreteCallback(stats_condition, stats_affect!,
                               save_positions=(false, false))

    return stats_cb
end

"""
Example 5: Outputting statistics at specific times

Add to your output callback:
"""
function setup_output_callback(inputs, params, mesh, SD)

    output_times = get(inputs, :turb_stats_output_times, ())

    function output_condition(u, t, integrator)
        # Check if current time is in output times
        return t in output_times
    end

    function output_affect!(integrator)
        if inputs[:l_turbulence_stats]
            # Get parameters
            turb_stats = integrator.p.turb_stats
            OUTPUT_DIR = inputs[:output_dir]
            t = integrator.t
            ndim = (SD == NSD_1D()) ? 1 : (SD == NSD_2D()) ? 2 : 3

            # Normalize accumulated statistics
            normalize_statistics!(turb_stats)

            # Print summary
            print_turbulence_stats_summary(turb_stats, ndim)

            # Write to ASCII file
            write_turbulence_stats_ascii(turb_stats, mesh, OUTPUT_DIR,
                                        "turbulence_stats", ndim, t)

            # Compute and write derived quantities
            compute_and_write_derived_stats(turb_stats, mesh, OUTPUT_DIR, ndim, t)

            # Reset for next window if requested
            if get(inputs, :turb_stats_reset_after_output, false)
                reset_turbulence_stats!(turb_stats)
            end
        end
    end

    output_cb = DiscreteCallback(output_condition, output_affect!,
                                save_positions=(false, false))

    return output_cb
end

"""
Helper function: Compute and write derived turbulence statistics
"""
function compute_and_write_derived_stats(turb_stats, mesh, OUTPUT_DIR, ndim, t)

    npoin = size(turb_stats.avvel, 1)

    # Compute RMS velocity fluctuations
    vel_rms = similar(turb_stats.avvel)
    compute_rms_fluctuations!(vel_rms, turb_stats, ndim)

    # Compute Reynolds stresses
    nstress = (ndim == 2) ? 1 : (ndim == 3) ? 3 : 1
    reynolds_stress = zeros(eltype(turb_stats.avvel), npoin, nstress)
    if ndim >= 2
        compute_reynolds_stresses!(reynolds_stress, turb_stats, ndim)
    end

    # Write derived statistics to file
    fname = @sprintf("derived_stats_t%.6f.dat", t)
    open(string(OUTPUT_DIR, "/", fname), "w") do f
        @printf(f, "# Derived Turbulence Statistics\n")
        @printf(f, "# Time: %.6f\n", t)
        @printf(f, "#\n")

        if ndim == 1
            @printf(f, "# ip | x | u'_rms\n")
            for ip = 1:npoin
                @printf(f, "%d %.6e %.6e\n",
                    ip, mesh.x[ip], vel_rms[ip,1])
            end

        elseif ndim == 2
            @printf(f, "# ip | x | y | u'_rms | v'_rms | <u'v'>\n")
            for ip = 1:npoin
                @printf(f, "%d %.6e %.6e %.6e %.6e %.6e\n",
                    ip, mesh.x[ip], mesh.y[ip],
                    vel_rms[ip,1], vel_rms[ip,2],
                    reynolds_stress[ip,1])
            end

        elseif ndim == 3
            @printf(f, "# ip | x | y | z | u'_rms | v'_rms | w'_rms | <u'v'> | <u'w'> | <v'w'>\n")
            for ip = 1:npoin
                @printf(f, "%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
                    ip, mesh.x[ip], mesh.y[ip], mesh.z[ip],
                    vel_rms[ip,1], vel_rms[ip,2], vel_rms[ip,3],
                    reynolds_stress[ip,1], reynolds_stress[ip,2], reynolds_stress[ip,3])
            end
        end
    end

    println(" # Derived statistics written to: ", OUTPUT_DIR, "/", fname)
end

"""
Example 6: Complete integration in main time loop

This shows the complete workflow:
"""
function example_complete_integration()
    println("""
    # Complete Integration Workflow:

    1. In user_inputs.jl:
       Add turbulence statistics parameters (see example_user_inputs_addition)

    2. In main solver initialization:
       # Allocate statistics
       turb_stats = allocate_turbulence_stats(SD, mesh.npoin, Float64, backend;
                                              l_turbulence_stats=inputs[:l_turbulence_stats])

       # Add to parameters struct
       params.turb_stats = turb_stats

    3. Setup callbacks:
       # Create statistics accumulation callback
       stats_cb = setup_statistics_callback(inputs, params, mesh, SD)

       # Create statistics output callback
       output_cb = setup_output_callback(inputs, params, mesh, SD)

       # Combine callbacks
       callbacks = CallbackSet(stats_cb, output_cb, other_callbacks...)

    4. Run time integration:
       sol = solve(prob, inputs[:ode_solver], callback=callbacks, ...)

    5. After simulation:
       # Final normalization and output
       if inputs[:l_turbulence_stats]
           normalize_statistics!(turb_stats)
           print_turbulence_stats_summary(turb_stats, ndim)
           write_turbulence_stats_ascii(turb_stats, mesh, OUTPUT_DIR,
                                        "final_stats", ndim, sol.t[end])
       end
    """)
end

"""
Example 7: Advanced - Spatial averaging for channel flow

For homogeneous turbulence in one direction (e.g., spanwise in channel flow):
"""
function compute_spanwise_averaged_statistics(turb_stats, mesh, ndim)

    # Assuming structured grid with spanwise direction in y
    # This is a simplified example - adapt to your mesh structure

    println(" # Computing spanwise-averaged statistics...")

    # Get unique x and z coordinates
    x_coords = unique(mesh.x)
    z_coords = (ndim == 3) ? unique(mesh.z) : [0.0]

    nx = length(x_coords)
    nz = length(z_coords)

    # Allocate arrays for averaged statistics
    avg_stats = Dict(
        :x => zeros(nx, nz),
        :z => zeros(nx, nz),
        :mean_u => zeros(nx, nz),
        :mean_v => zeros(nx, nz),
        :rms_u => zeros(nx, nz),
        :rms_v => zeros(nx, nz),
        :reynolds_stress => zeros(nx, nz)
    )

    # Average over spanwise (y) direction
    for (ix, x) in enumerate(x_coords)
        for (iz, z) in enumerate(z_coords)
            # Find all points with this x, z coordinate
            # (different y values)
            indices = findall(i -> mesh.x[i] ≈ x && mesh.z[i] ≈ z, 1:length(mesh.x))

            if !isempty(indices)
                # Compute averages
                avg_stats[:x][ix, iz] = x
                avg_stats[:z][ix, iz] = z
                avg_stats[:mean_u][ix, iz] = mean(turb_stats.avvel[indices, 1])
                avg_stats[:mean_v][ix, iz] = mean(turb_stats.avvel[indices, 2])

                # Compute RMS
                for i in indices
                    mean_u = turb_stats.avvel[i, 1]
                    mean_v = turb_stats.avvel[i, 2]
                    var_u = turb_stats.avve2[i, 1] - mean_u^2
                    var_v = turb_stats.avve2[i, 2] - mean_v^2
                    avg_stats[:rms_u][ix, iz] += sqrt(max(0, var_u))
                    avg_stats[:rms_v][ix, iz] += sqrt(max(0, var_v))

                    # Reynolds stress
                    reynolds = turb_stats.avvex[i, 1] - mean_u * mean_v
                    avg_stats[:reynolds_stress][ix, iz] += reynolds
                end

                # Normalize by number of points
                n = length(indices)
                avg_stats[:rms_u][ix, iz] /= n
                avg_stats[:rms_v][ix, iz] /= n
                avg_stats[:reynolds_stress][ix, iz] /= n
            end
        end
    end

    return avg_stats
end

"""
Example 8: Writing spatial profiles

For channel flow or boundary layer, write vertical profiles:
"""
function write_vertical_profiles(avg_stats, OUTPUT_DIR, filename, t)

    fname = string(filename, "_profiles_t", @sprintf("%.6f", t), ".dat")

    open(string(OUTPUT_DIR, "/", fname), "w") do f
        @printf(f, "# Vertical Profiles of Turbulence Statistics\n")
        @printf(f, "# Time: %.6f\n", t)
        @printf(f, "# z | <u> | <v> | u'_rms | v'_rms | <u'v'>\n")

        nz = size(avg_stats[:z], 2)
        for iz = 1:nz
            @printf(f, "%.6e %.6e %.6e %.6e %.6e %.6e\n",
                avg_stats[:z][1, iz],
                avg_stats[:mean_u][1, iz],
                avg_stats[:mean_v][1, iz],
                avg_stats[:rms_u][1, iz],
                avg_stats[:rms_v][1, iz],
                avg_stats[:reynolds_stress][1, iz])
        end
    end

    println(" # Vertical profiles written to: ", OUTPUT_DIR, "/", fname)
end

# Print usage instructions
println("="^70)
println("Turbulence Statistics Integration Examples Loaded")
println("="^70)
println()
println("Available example functions:")
println("  - example_user_inputs_addition(): Parameters to add to user_inputs()")
println("  - example_allocation(): How to allocate statistics structure")
println("  - example_initialization(): How to initialize at start time")
println("  - setup_statistics_callback(): Create accumulation callback")
println("  - setup_output_callback(): Create output callback")
println("  - example_complete_integration(): Complete workflow overview")
println("  - compute_spanwise_averaged_statistics(): Spatial averaging example")
println("  - write_vertical_profiles(): Write 1D profiles")
println()
println("Quick Start:")
println("  1. Add parameters from example_user_inputs_addition() to your user_inputs.jl")
println("  2. Follow the steps in example_complete_integration()")
println("  3. Run your simulation as usual")
println("="^70)
