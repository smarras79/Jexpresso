#!/usr/bin/env python3
"""
Download ERA5 Reanalysis Data for Jexpresso LES

This script downloads ERA5 pressure level data from the Copernicus Climate Data Store
for use with Jexpresso's ERA5 extension.

Requirements:
    pip install cdsapi

Setup:
    1. Register at https://cds.climate.copernicus.eu
    2. Create ~/.cdsapirc with your API key:
       url: https://cds.climate.copernicus.eu/api/v2
       key: YOUR_UID:YOUR_API_KEY

Usage:
    python download_era5.py [options]

Author: Jexpresso Development Team
Date: 2025-11-30
"""

import argparse
import cdsapi
from datetime import datetime


def download_era5_pressure_levels(
    output_file,
    year,
    month,
    day,
    time,
    area=None,
    pressure_levels=None,
    variables=None
):
    """
    Download ERA5 pressure level data.

    Parameters
    ----------
    output_file : str
        Output NetCDF filename
    year : str
        Year (e.g., '2024')
    month : str
        Month (e.g., '01')
    day : str
        Day (e.g., '15')
    time : str or list
        Time(s) in HH:MM format (e.g., '12:00' or ['00:00', '06:00', '12:00', '18:00'])
    area : list, optional
        [North, West, South, East] in degrees (e.g., [45, -110, 35, -100])
        If None, downloads global data
    pressure_levels : list, optional
        Pressure levels in hPa (e.g., ['1000', '950', '900'])
        If None, uses default levels for LES
    variables : list, optional
        List of variables to download
        If None, uses default variables for LES
    """

    # Default pressure levels for LES (surface to ~15 km)
    if pressure_levels is None:
        pressure_levels = [
            '1000', '975', '950', '925', '900', '875', '850', '825', '800',
            '775', '750', '700', '650', '600', '550', '500', '450', '400',
            '350', '300', '250', '200'
        ]

    # Default variables for Jexpresso LES
    if variables is None:
        variables = [
            'temperature',
            'u_component_of_wind',
            'v_component_of_wind',
            'vertical_velocity',
            'specific_humidity',
            'geopotential',
        ]

    # Construct request
    request = {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': variables,
        'pressure_level': pressure_levels,
        'year': year,
        'month': month,
        'day': day,
        'time': time if isinstance(time, list) else [time],
    }

    # Add spatial subsetting if specified
    if area is not None:
        request['area'] = area

    print("=" * 70)
    print("ERA5 Download Configuration")
    print("=" * 70)
    print(f"Output file:      {output_file}")
    print(f"Date:             {year}-{month}-{day}")
    print(f"Time(s):          {time}")
    print(f"Pressure levels:  {len(pressure_levels)} levels ({pressure_levels[0]} - {pressure_levels[-1]} hPa)")
    print(f"Variables:        {', '.join(variables)}")
    if area:
        print(f"Area:             N={area[0]}°, W={area[1]}°, S={area[2]}°, E={area[3]}°")
    else:
        print(f"Area:             Global")
    print("=" * 70)

    # Initialize CDS API client
    c = cdsapi.Client()

    # Submit request
    print("\nSubmitting request to CDS...")
    print("This may take several minutes to hours depending on data size.")
    print("You can close this script - download will continue on CDS servers.")
    print("Check status at: https://cds.climate.copernicus.eu/cdsapp#!/yourrequests")
    print()

    try:
        c.retrieve('reanalysis-era5-pressure-levels', request, output_file)
        print(f"\n✓ Download complete: {output_file}")
        print(f"  File size: {get_file_size(output_file)}")
    except Exception as e:
        print(f"\n✗ Download failed: {e}")
        print("\nTroubleshooting:")
        print("  - Check your API key in ~/.cdsapirc")
        print("  - Verify CDS Terms & Conditions are accepted")
        print("  - Check CDS status: https://cds.climate.copernicus.eu")
        raise


def download_surface_pressure(output_file, year, month, day, time, area=None):
    """
    Download ERA5 surface pressure (required for vertical interpolation).

    Parameters are the same as download_era5_pressure_levels.
    """

    request = {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': ['surface_pressure'],
        'year': year,
        'month': month,
        'day': day,
        'time': time if isinstance(time, list) else [time],
    }

    if area is not None:
        request['area'] = area

    print("\nDownloading surface pressure...")

    c = cdsapi.Client()
    c.retrieve('reanalysis-era5-single-levels', request, output_file)

    print(f"✓ Surface pressure downloaded: {output_file}")


def get_file_size(filepath):
    """Get human-readable file size."""
    try:
        import os
        size_bytes = os.path.getsize(filepath)
        for unit in ['B', 'KB', 'MB', 'GB']:
            if size_bytes < 1024.0:
                return f"{size_bytes:.1f} {unit}"
            size_bytes /= 1024.0
    except:
        return "unknown"


def main():
    parser = argparse.ArgumentParser(
        description='Download ERA5 data for Jexpresso LES',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download for Boulder, CO on Jan 15, 2024 at 12:00 UTC
  python download_era5.py -y 2024 -m 01 -d 15 -t 12:00 \\
      --area 45 -110 35 -100 -o era5_boulder.nc

  # Download global data for multiple times
  python download_era5.py -y 2024 -m 01 -d 15 \\
      -t 00:00 06:00 12:00 18:00 -o era5_global.nc

  # Download with custom pressure levels
  python download_era5.py -y 2024 -m 01 -d 15 -t 12:00 \\
      --levels 1000 950 900 850 800 -o era5_custom.nc
        """
    )

    parser.add_argument('-y', '--year', required=True,
                        help='Year (e.g., 2024)')
    parser.add_argument('-m', '--month', required=True,
                        help='Month (e.g., 01)')
    parser.add_argument('-d', '--day', required=True,
                        help='Day (e.g., 15)')
    parser.add_argument('-t', '--time', nargs='+', required=True,
                        help='Time(s) in HH:MM format (e.g., 12:00)')
    parser.add_argument('-o', '--output', default='era5_data.nc',
                        help='Output NetCDF file (default: era5_data.nc)')
    parser.add_argument('--area', type=float, nargs=4, metavar=('N', 'W', 'S', 'E'),
                        help='Area: North West South East (degrees)')
    parser.add_argument('--levels', nargs='+',
                        help='Pressure levels in hPa (default: standard LES levels)')
    parser.add_argument('--variables', nargs='+',
                        help='Variables to download (default: T, u, v, w, q, z)')
    parser.add_argument('--surface-pressure', action='store_true',
                        help='Also download surface pressure separately')

    args = parser.parse_args()

    # Download pressure level data
    download_era5_pressure_levels(
        output_file=args.output,
        year=args.year,
        month=args.month,
        day=args.day,
        time=args.time,
        area=args.area,
        pressure_levels=args.levels,
        variables=args.variables
    )

    # Optionally download surface pressure
    if args.surface_pressure:
        sp_file = args.output.replace('.nc', '_surface.nc')
        download_surface_pressure(
            output_file=sp_file,
            year=args.year,
            month=args.month,
            day=args.day,
            time=args.time,
            area=args.area
        )

    print("\n" + "=" * 70)
    print("Download complete!")
    print("=" * 70)
    print(f"\nNext steps:")
    print(f"  1. Verify data: ncdump -h {args.output}")
    print(f"  2. Move to Jexpresso: mv {args.output} /path/to/Jexpresso/data_files/")
    print(f"  3. Update user_inputs.jl:")
    print(f"     :era5_file => \"./data_files/{args.output}\",")
    if args.area:
        print(f"     :era5_target_lat => {(args.area[0] + args.area[2]) / 2:.1f},")
        print(f"     :era5_target_lon => {(args.area[1] + args.area[3]) / 2:.1f},")
    print()


if __name__ == '__main__':
    main()
