#!/usr/bin/env python3
"""
Jexpresso Equation Generator

AI-driven tool to generate problem directories from PDF equation specifications.
"""
import click
import json
from pathlib import Path
import sys

from pdf_parser import PDFParser
from equation_analyzer import EquationAnalyzer
from code_generator import CodeGenerator


@click.command()
@click.argument('pdf_file', type=click.Path(exists=True))
@click.option(
    '--output-dir',
    '-o',
    default='../../problems',
    help='Output directory for generated problem (default: ../../problems)'
)
@click.option(
    '--category',
    '-c',
    default=None,
    help='Problem category (e.g., CompEuler, AdvDiff). If not specified, inferred from equations'
)
@click.option(
    '--api-key',
    '-k',
    default=None,
    envvar='ANTHROPIC_API_KEY',
    help='Anthropic API key (or set ANTHROPIC_API_KEY env var)'
)
@click.option(
    '--save-json',
    '-s',
    is_flag=True,
    help='Save extracted equation information to JSON file'
)
@click.option(
    '--dry-run',
    '-d',
    is_flag=True,
    help='Extract and analyze equations but do not generate code'
)
def main(pdf_file, output_dir, category, api_key, save_json, dry_run):
    """
    Generate a Jexpresso problem directory from a PDF equation specification.

    PDF_FILE: Path to the PDF file containing equation definitions

    Example:
        python main.py my_equations.pdf -o problems -c CompEuler
    """
    click.echo("=" * 70)
    click.echo("Jexpresso Equation Generator")
    click.echo("=" * 70)

    # Step 1: Parse PDF
    click.echo(f"\n📄 Parsing PDF: {pdf_file}")
    try:
        parser = PDFParser(pdf_file)
        pdf_text = parser.extract_text()
        click.echo(f"✓ Extracted text from {parser.get_page_count()} pages")
    except Exception as e:
        click.echo(f"✗ Error parsing PDF: {e}", err=True)
        sys.exit(1)

    # Step 2: Analyze equations with AI
    click.echo("\n🤖 Analyzing equations with Claude AI...")
    try:
        analyzer = EquationAnalyzer(api_key=api_key)
        eq_info = analyzer.analyze_equations(pdf_text)
        analyzer.validate_equation_info(eq_info)

        click.echo(f"✓ Detected equation system: {eq_info['problem_name']}")
        click.echo(f"  - {eq_info['num_equations']} equations")
        click.echo(f"  - {eq_info.get('spatial_dimensions', 2)}D problem")
        click.echo(f"  - Variables: {', '.join(eq_info['variables'])}")

    except Exception as e:
        click.echo(f"✗ Error analyzing equations: {e}", err=True)
        sys.exit(1)

    # Save JSON if requested
    if save_json:
        json_file = Path(pdf_file).with_suffix('.json')
        with open(json_file, 'w') as f:
            json.dump(eq_info, f, indent=2)
        click.echo(f"\n💾 Saved equation info to: {json_file}")

    # Stop here if dry run
    if dry_run:
        click.echo("\n🏁 Dry run complete (no files generated)")
        _print_equation_info(eq_info)
        return

    # Step 3: Generate code
    click.echo("\n📝 Generating Julia code...")
    try:
        generator = CodeGenerator()
        problem_dir = generator.generate_problem_directory(
            eq_info=eq_info,
            output_dir=output_dir,
            problem_category=category
        )

        generator.generate_readme(eq_info, problem_dir)

        click.echo(f"\n✅ Problem directory created: {problem_dir}")
        click.echo("\n📋 Generated files:")
        for file in sorted(problem_dir.glob("*.jl")):
            click.echo(f"  - {file.name}")
        click.echo(f"  - README.md")

    except Exception as e:
        click.echo(f"\n✗ Error generating code: {e}", err=True)
        sys.exit(1)

    # Print summary
    _print_summary(eq_info, problem_dir)


def _print_equation_info(eq_info):
    """Print detailed equation information."""
    click.echo("\n" + "=" * 70)
    click.echo("Equation Information")
    click.echo("=" * 70)
    click.echo(f"\nProblem: {eq_info['problem_name']}")
    click.echo(f"Description: {eq_info.get('description', 'N/A')}")
    click.echo(f"\nEquations: {eq_info['num_equations']}")
    click.echo(f"Dimensions: {eq_info.get('spatial_dimensions', 2)}D")

    click.echo("\n Variables:")
    for i, var in enumerate(eq_info['variables'], 1):
        desc = eq_info.get('variable_descriptions', {}).get(var, '')
        click.echo(f"  {i}. {var}" + (f" - {desc}" if desc else ""))

    click.echo("\nFlux (X-direction):")
    for i, flux in enumerate(eq_info['flux_x'], 1):
        click.echo(f"  {i}. {flux}")

    if eq_info.get('flux_y'):
        click.echo("\nFlux (Y-direction):")
        for i, flux in enumerate(eq_info['flux_y'], 1):
            click.echo(f"  {i}. {flux}")

    click.echo("\nSource Terms:")
    for i, source in enumerate(eq_info['source_terms'], 1):
        click.echo(f"  {i}. {source}")


def _print_summary(eq_info, problem_dir):
    """Print summary and next steps."""
    click.echo("\n" + "=" * 70)
    click.echo("Next Steps")
    click.echo("=" * 70)
    click.echo("\n1. Review the generated files, especially:")
    click.echo(f"   - {problem_dir}/initialize.jl (initial conditions)")
    click.echo(f"   - {problem_dir}/user_bc.jl (boundary conditions)")
    click.echo(f"   - {problem_dir}/user_flux.jl (flux calculations)")
    click.echo("\n2. Update user_inputs.jl with problem-specific parameters")
    click.echo("\n3. Test with your mesh and verify results")
    click.echo("\n4. See README.md for detailed TODO list")
    click.echo("\n" + "=" * 70)


if __name__ == '__main__':
    main()
