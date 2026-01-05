"""
Code generator for creating Julia problem directories from equation information.
"""
from pathlib import Path
from typing import Dict
from jinja2 import Environment, FileSystemLoader
import os


class CodeGenerator:
    """Generates Julia code files from equation information using templates."""

    def __init__(self, template_dir: str = None):
        """
        Initialize code generator.

        Args:
            template_dir: Directory containing Jinja2 templates
        """
        if template_dir is None:
            # Default to templates directory in same location as this file
            current_dir = Path(__file__).parent
            template_dir = current_dir / "templates"

        self.template_dir = Path(template_dir)
        if not self.template_dir.exists():
            raise FileNotFoundError(f"Template directory not found: {template_dir}")

        self.env = Environment(loader=FileSystemLoader(str(self.template_dir)))

    def generate_problem_directory(
        self,
        eq_info: Dict,
        output_dir: str,
        problem_category: str = None
    ) -> Path:
        """
        Generate a complete problem directory with all Julia files.

        Args:
            eq_info: Equation information from EquationAnalyzer
            output_dir: Base output directory (e.g., "problems")
            problem_category: Category subdirectory (e.g., "CompEuler", "AdvDiff")
                             If None, uses problem_name from eq_info

        Returns:
            Path to the created problem directory
        """
        # Determine output path
        output_path = Path(output_dir)
        if problem_category is None:
            problem_category = eq_info.get("equation_type", "Custom")
            problem_category = problem_category.replace(" ", "").replace("-", "")

        problem_name = eq_info["problem_name"]
        full_path = output_path / problem_category / problem_name

        # Create directory
        full_path.mkdir(parents=True, exist_ok=True)

        # Prepare template context
        context = self._prepare_context(eq_info)

        # Generate each file
        files_to_generate = [
            "initialize.jl",
            "user_bc.jl",
            "user_flux.jl",
            "user_inputs.jl",
            "user_primitives.jl",
            "user_source.jl"
        ]

        for filename in files_to_generate:
            self._generate_file(filename, context, full_path)

        return full_path

    def _prepare_context(self, eq_info: Dict) -> Dict:
        """
        Prepare template context from equation information.

        Args:
            eq_info: Equation information dictionary

        Returns:
            Context dictionary for Jinja2 templates
        """
        # Format variables list for Julia
        qvars_list = ", ".join([f'"{var}"' for var in eq_info["variables"]])

        # Determine output variables (usually primitive form)
        if "primitive_variables" in eq_info and eq_info["primitive_variables"]:
            qoutvars = eq_info["primitive_variables"]
        else:
            # Default to same as variables
            qoutvars = eq_info["variables"]
        qoutvars_list = ", ".join([f'"{var}"' for var in qoutvars])

        # Format viscosity coefficients
        num_eqs = eq_info["num_equations"]
        viscosity = "[" + ", ".join(["0.0"] * num_eqs) + "]"

        context = {
            "problem_name": eq_info["problem_name"],
            "description": eq_info.get("description", ""),
            "num_equations": eq_info["num_equations"],
            "spatial_dim": eq_info.get("spatial_dimensions", 2),
            "variables": eq_info["variables"],
            "qvars_list": qvars_list,
            "qoutvars_list": qoutvars_list,
            "flux_x": eq_info["flux_x"],
            "flux_y": eq_info.get("flux_y", []),
            "flux_z": eq_info.get("flux_z", []),
            "source_terms": eq_info["source_terms"],
            "has_pressure": eq_info.get("has_pressure", False),
            "has_gravity": eq_info.get("has_gravity", False),
            "has_source": any(s != "0" and s != "0.0" for s in eq_info["source_terms"]),
            "viscosity_coefficients": viscosity,
            "equation_type": eq_info.get("equation_type", "custom"),
        }

        # Add enumerate function to context for templates
        context["enumerate"] = enumerate
        context["range"] = range

        return context

    def _generate_file(self, filename: str, context: Dict, output_dir: Path):
        """
        Generate a single Julia file from template.

        Args:
            filename: Name of the file to generate (e.g., "initialize.jl")
            context: Template context
            output_dir: Output directory path
        """
        template_name = f"{filename}.j2"

        try:
            template = self.env.get_template(template_name)
            content = template.render(**context)

            output_file = output_dir / filename
            with open(output_file, 'w') as f:
                f.write(content)

            print(f"✓ Generated {filename}")

        except Exception as e:
            print(f"✗ Error generating {filename}: {e}")
            raise

    def generate_readme(self, eq_info: Dict, output_dir: Path):
        """
        Generate a README file for the problem.

        Args:
            eq_info: Equation information
            output_dir: Output directory path
        """
        readme_content = f"""# {eq_info['problem_name']}

{eq_info.get('description', '')}

## Equation System

Number of equations: {eq_info['num_equations']}
Spatial dimensions: {eq_info.get('spatial_dimensions', 2)}

### Variables
{chr(10).join(f"- {var}: {eq_info.get('variable_descriptions', {}).get(var, 'N/A')}" for var in eq_info['variables'])}

### Flux Terms

#### X-direction (F)
{chr(10).join(f"{i+1}. {flux}" for i, flux in enumerate(eq_info['flux_x']))}

#### Y-direction (G)
{chr(10).join(f"{i+1}. {flux}" for i, flux in enumerate(eq_info.get('flux_y', [])))}

### Source Terms (S)
{chr(10).join(f"{i+1}. {source}" for i, source in enumerate(eq_info['source_terms']))}

## Notes

This directory was auto-generated by the Jexpresso Equation Generator tool.

### TODO
- [ ] Review and verify flux calculations
- [ ] Implement initial conditions in `initialize.jl`
- [ ] Define boundary conditions in `user_bc.jl`
- [ ] Update `user_inputs.jl` with problem-specific parameters
- [ ] Verify primitive variable transformations in `user_primitives.jl`
- [ ] Test with actual mesh and verify results
"""

        readme_file = output_dir / "README.md"
        with open(readme_file, 'w') as f:
            f.write(readme_content)

        print(f"✓ Generated README.md")
