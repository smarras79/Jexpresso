#!/usr/bin/env julia

"""
CLI interface for the EquationGenerator package.

Usage:
    julia --project=. src/cli.jl equations.pdf [OPTIONS]
"""

using ArgParse

# Add current directory to load path
push!(LOAD_PATH, @__DIR__)

using EquationGenerator

function parse_commandline()
    s = ArgParseSettings(
        description = "Jexpresso Equation Generator - Generate problem directories from PDF equations",
        epilog = """
        Examples:
            julia src/cli.jl equations.pdf
            julia src/cli.jl equations.pdf -c CompEuler -o ../../problems
            julia src/cli.jl equations.pdf --dry-run
        """
    )

    @add_arg_table! s begin
        "pdf_file"
            help = "Path to PDF file containing equation specifications"
            required = true

        "--output-dir", "-o"
            help = "Output directory for generated problem (default: ../../problems)"
            default = "../../problems"

        "--category", "-c"
            help = "Problem category (e.g., CompEuler, AdvDiff). Auto-detected if not specified."
            default = nothing

        "--api-key", "-k"
            help = "Anthropic API key (or set ANTHROPIC_API_KEY env var)"
            default = nothing

        "--save-json", "-s"
            help = "Save extracted equation information to JSON file"
            action = :store_true

        "--dry-run", "-d"
            help = "Extract and analyze equations but do not generate code"
            action = :store_true
    end

    return parse_args(s)
end

function main()
    args = parse_commandline()

    try
        generate_problem(
            args["pdf_file"];
            output_dir = args["output-dir"],
            category = args["category"],
            api_key = args["api-key"],
            save_json = args["save-json"],
            dry_run = args["dry-run"]
        )
    catch e
        println(stderr, "Error: $e")
        exit(1)
    end
end

# Run main if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
