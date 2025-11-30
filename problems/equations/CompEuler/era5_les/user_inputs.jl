"""
    ERA5 LES Test Case - User Inputs

    This test case demonstrates how to use ERA5 reanalysis data for
    Large Eddy Simulation weather forecasting with Jexpresso.

    Author: Jexpresso Development Team
    Date: 2025-11-30
"""

function user_inputs()
    inputs = Dict(
        #---------------------------------------------------------------------------
        # Case description
        #---------------------------------------------------------------------------
        :case                 => "era5_les",

        #---------------------------------------------------------------------------
        # Time integration
        #---------------------------------------------------------------------------
        :ode_solver           => SSPRK54(),
        :Δt                   => 1.0,        # Time step [s]
        :tinit                => 0.0,        # Initial time [s]
        :tend                 => 3600.0,     # End time [s] - 1 hour simulation
        :diagnostics_at_times => (600:600:3600),  # Output every 10 minutes

        #---------------------------------------------------------------------------
        # Integration and quadrature properties
        #---------------------------------------------------------------------------
        :interpolation_nodes  => "lgl",
        :nop                  => 4,          # Polynomial order

        #---------------------------------------------------------------------------
        # Physical parameters
        #---------------------------------------------------------------------------
        :lsource              => true,       # Enable source terms
        :lsaturation          => true,       # Enable saturation adjustment

        # Viscosity and SGS model
        :lvisc                => true,
        :visc_model           => SMAG(),     # Smagorinsky SGS model
        :ivisc_equations      => [1, 2, 3, 4, 5, 6, 7],
        # Smagorinsky constant: Cs = 0.23, input Cs^2 for momentum, Cs^2/Pr for scalars
        :μ                    => [0.1587, 0.0529, 0.0529, 0.0529, 0.1587, 0.1587, 0.1587],

        #---------------------------------------------------------------------------
        # ERA5 data configuration
        #---------------------------------------------------------------------------
        :lERA5                => true,       # Enable ERA5 data reading
        :era5_file            => "./data_files/era5_example.nc",  # ERA5 NetCDF file
        :era5_target_lat      => 40.0,       # Target latitude [degrees N]
        :era5_target_lon      => -105.0,     # Target longitude [degrees E]
        :era5_time_index      => 1,          # Time index in NetCDF file

        # Optional: Specify custom variable names if your ERA5 file uses different naming
        # :era5_var_names       => Dict(
        #     "T" => "temperature",
        #     "u" => "u_wind",
        #     "v" => "v_wind",
        #     "q" => "specific_humidity",
        #     "z" => "geopotential",
        #     "sp" => "surface_pressure"
        # ),

        # Optional: Select specific pressure levels [Pa]
        # :era5_pressure_levels => [100000.0, 95000.0, 90000.0, 85000.0, 80000.0, 75000.0],

        # Optional: Enable ERA5 nudging/relaxation forcing
        # :lERA5_forcing        => true,
        # :era5_nudging_timescale => 3600.0,  # Nudging timescale [s]

        #---------------------------------------------------------------------------
        # Mesh parameters
        #---------------------------------------------------------------------------
        :lread_gmsh           => false,      # Use simple Cartesian mesh
        :nsd                  => 3,          # 3D simulation

        # Domain: 5km × 5km × 3km LES domain
        :nelx                 => 16,
        :xmin                 => 0.0,
        :xmax                 => 5000.0,     # 5 km in x

        :nely                 => 16,
        :ymin                 => 0.0,
        :ymax                 => 5000.0,     # 5 km in y

        :nelz                 => 12,
        :zmin                 => 0.0,
        :zmax                 => 3000.0,     # 3 km in z

        #---------------------------------------------------------------------------
        # Output parameters
        #---------------------------------------------------------------------------
        :outformat            => "vtk",
        :output_dir           => "./output/",
        :loutput_pert         => false,
        :lwrite_initial       => true,

        #---------------------------------------------------------------------------
    ) #Dict

    return inputs
end
