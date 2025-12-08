#!/usr/bin/env julia
"""
Jexpresso Dependency Installer
================================

This script ensures that all required dependencies for Jexpresso are installed
with the correct versions. It handles installation errors and provides clear
feedback about the installation process.

Usage:
    julia install_dependencies.jl

Or make it executable:
    chmod +x install_dependencies.jl
    ./install_dependencies.jl
"""

using Pkg

# ANSI color codes for terminal output
const COLOR_GREEN = "\033[32m"
const COLOR_RED = "\033[31m"
const COLOR_YELLOW = "\033[33m"
const COLOR_BLUE = "\033[34m"
const COLOR_CYAN = "\033[36m"
const COLOR_RESET = "\033[0m"

# Package versions required for Jexpresso
const REQUIRED_PACKAGES = Dict(
    "MPI" => "0.20.22",
    "MPIPreferences" => "0.1.11",
    "PackageCompiler" => "2.2.1",
    "Thermodynamics" => "0.12.7",
    "PrettyTables" => "2.4.0",
    "Crayons" => "4.1.1",
    "UnicodePlots" => "3.7.2",
    "Gridap" => "0.18.12",
    "GridapDistributed" => "0.4.7",
    "GridapGmsh" => "0.7.2",
    "GridapP4est" => "0.3.11"
)

const MIN_JULIA_VERSION = v"1.11.2"

"""
    print_header()

Print a stylized header for the installation script.
"""
function print_header()
    println("\n" * "="^70)
    println("$(COLOR_CYAN)Jexpresso Dependency Installer$(COLOR_RESET)")
    println("="^70 * "\n")
end

"""
    print_section(title::String)

Print a section header.
"""
function print_section(title::String)
    println("\n$(COLOR_BLUE)▶ $title$(COLOR_RESET)")
    println("-"^70)
end

"""
    print_success(message::String)

Print a success message in green.
"""
function print_success(message::String)
    println("$(COLOR_GREEN)✓ $message$(COLOR_RESET)")
end

"""
    print_error(message::String)

Print an error message in red.
"""
function print_error(message::String)
    println("$(COLOR_RED)✗ $message$(COLOR_RESET)")
end

"""
    print_warning(message::String)

Print a warning message in yellow.
"""
function print_warning(message::String)
    println("$(COLOR_YELLOW)⚠ $message$(COLOR_RESET)")
end

"""
    print_info(message::String)

Print an informational message.
"""
function print_info(message::String)
    println("  $message")
end

"""
    check_julia_version()

Verify that the Julia version meets the minimum requirements.
"""
function check_julia_version()
    print_section("Checking Julia Version")

    current_version = VERSION
    print_info("Current Julia version: $current_version")
    print_info("Minimum required version: $MIN_JULIA_VERSION")

    if current_version >= MIN_JULIA_VERSION
        print_success("Julia version is compatible")
        return true
    else
        print_error("Julia version is too old")
        print_warning("Please upgrade to Julia $MIN_JULIA_VERSION or higher")
        return false
    end
end

"""
    get_installed_version(package_name::String)

Get the currently installed version of a package, or nothing if not installed.
"""
function get_installed_version(package_name::String)
    try
        deps = Pkg.dependencies()
        for (uuid, dep) in deps
            if dep.name == package_name
                return dep.version
            end
        end
        return nothing
    catch e
        return nothing
    end
end

"""
    install_package(name::String, version::String; retry_count::Int=3)

Install a specific version of a package with retry logic.
"""
function install_package(name::String, version::String; retry_count::Int=3)
    for attempt in 1:retry_count
        try
            print_info("Installing $name@$version (attempt $attempt/$retry_count)...")

            # Use Pkg.add with specific version
            Pkg.add(name=name, version=version)

            # Verify installation
            installed_version = get_installed_version(name)
            if installed_version !== nothing
                installed_str = string(installed_version)
                if installed_str == version
                    print_success("Successfully installed $name@$version")
                    return true
                else
                    print_warning("Installed $name@$installed_str instead of $version")
                    # Try to pin the exact version
                    try
                        Pkg.pin(name=name, version=version)
                        installed_version = get_installed_version(name)
                        installed_str = string(installed_version)
                        if installed_str == version
                            print_success("Successfully pinned $name@$version")
                            return true
                        end
                    catch e
                        print_warning("Could not pin $name to version $version")
                    end
                end
            else
                print_error("Failed to verify installation of $name")
            end
        catch e
            print_warning("Attempt $attempt failed: $(sprint(showerror, e))")
            if attempt < retry_count
                print_info("Retrying...")
                sleep(1)  # Wait before retry
            end
        end
    end

    print_error("Failed to install $name@$version after $retry_count attempts")
    return false
end

"""
    verify_installations()

Verify that all required packages are installed with correct versions.
"""
function verify_installations()
    print_section("Verifying Installations")

    all_correct = true
    mismatched = []

    for (name, required_version) in REQUIRED_PACKAGES
        installed_version = get_installed_version(name)

        if installed_version === nothing
            print_error("$name is not installed")
            push!(mismatched, name)
            all_correct = false
        else
            installed_str = string(installed_version)
            if installed_str == required_version
                print_success("$name@$installed_str ✓")
            else
                print_warning("$name: installed $installed_str, required $required_version")
                push!(mismatched, name)
                all_correct = false
            end
        end
    end

    return all_correct, mismatched
end

"""
    activate_current_environment()

Activate the current project environment.
"""
function activate_current_environment()
    print_section("Activating Project Environment")

    project_dir = dirname(@__FILE__)
    print_info("Project directory: $project_dir")

    try
        Pkg.activate(project_dir)
        print_success("Project environment activated")
        return true
    catch e
        print_error("Failed to activate project environment: $(sprint(showerror, e))")
        return false
    end
end

"""
    install_all_packages()

Install all required packages with specific versions.
"""
function install_all_packages()
    print_section("Installing Required Packages")

    total = length(REQUIRED_PACKAGES)
    success_count = 0
    failed_packages = String[]

    for (i, (name, version)) in enumerate(sort(collect(REQUIRED_PACKAGES), by=x->x[1]))
        println("\n[$i/$total] Processing $name...")
        if install_package(name, version)
            success_count += 1
        else
            push!(failed_packages, name)
        end
    end

    println("\n" * "="^70)
    println("Installation Summary:")
    println("  Successful: $success_count/$total")
    println("  Failed: $(length(failed_packages))/$total")

    if !isempty(failed_packages)
        println("\nFailed packages:")
        for pkg in failed_packages
            println("  - $pkg")
        end
    end
    println("="^70)

    return isempty(failed_packages)
end

"""
    resolve_dependencies()

Resolve all package dependencies in the environment.
"""
function resolve_dependencies()
    print_section("Resolving Dependencies")

    try
        print_info("Running Pkg.resolve()...")
        Pkg.resolve()
        print_success("Dependencies resolved successfully")
        return true
    catch e
        print_error("Failed to resolve dependencies: $(sprint(showerror, e))")
        print_warning("You may need to manually resolve conflicts")
        return false
    end
end

"""
    precompile_packages()

Precompile all installed packages.
"""
function precompile_packages()
    print_section("Precompiling Packages")

    try
        print_info("This may take several minutes...")
        Pkg.precompile()
        print_success("All packages precompiled successfully")
        return true
    catch e
        print_error("Precompilation failed: $(sprint(showerror, e))")
        print_warning("Some packages may not work correctly")
        return false
    end
end

"""
    create_backup()

Create a backup of the current Manifest.toml if it exists.
"""
function create_backup()
    manifest_path = joinpath(dirname(@__FILE__), "Manifest.toml")

    if isfile(manifest_path)
        backup_path = joinpath(dirname(@__FILE__), "Manifest.toml.backup")
        try
            cp(manifest_path, backup_path, force=true)
            print_info("Created backup: Manifest.toml.backup")
            return true
        catch e
            print_warning("Failed to create backup: $(sprint(showerror, e))")
            return false
        end
    end
    return true
end

"""
    print_final_summary(success::Bool)

Print a final summary of the installation process.
"""
function print_final_summary(success::Bool)
    println("\n" * "="^70)
    if success
        println("$(COLOR_GREEN)Installation completed successfully!$(COLOR_RESET)")
        println("\nYou can now use Jexpresso with the correct package versions.")
        println("To get started, run:")
        println("  julia --project=. your_script.jl")
    else
        println("$(COLOR_RED)Installation completed with errors$(COLOR_RESET)")
        println("\nSome packages could not be installed with the required versions.")
        println("Please check the error messages above and try:")
        println("  1. Updating Julia to version $MIN_JULIA_VERSION or higher")
        println("  2. Running this script again")
        println("  3. Manually installing failed packages using:")
        println("     julia> using Pkg")
        println("     julia> Pkg.add(name=\"PackageName\", version=\"x.y.z\")")
    end
    println("="^70 * "\n")
end

"""
    main()

Main installation routine.
"""
function main()
    print_header()

    # Step 1: Check Julia version
    if !check_julia_version()
        print_final_summary(false)
        return 1
    end

    # Step 2: Activate project environment
    if !activate_current_environment()
        print_final_summary(false)
        return 1
    end

    # Step 3: Create backup
    create_backup()

    # Step 4: Install required packages
    install_success = install_all_packages()

    # Step 5: Resolve dependencies
    resolve_success = resolve_dependencies()

    # Step 6: Verify installations
    verify_success, mismatched = verify_installations()

    # Step 7: Precompile (optional, continue even if it fails)
    if install_success && resolve_success && verify_success
        precompile_packages()
    end

    # Final summary
    overall_success = install_success && resolve_success && verify_success
    print_final_summary(overall_success)

    return overall_success ? 0 : 1
end

# Run main function
exit(main())
