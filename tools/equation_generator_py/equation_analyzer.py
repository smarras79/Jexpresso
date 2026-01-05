"""
Equation Analyzer using Claude API to extract structured equation information.
"""
import json
from typing import Dict, List, Optional
from anthropic import Anthropic
import os


class EquationAnalyzer:
    """Uses Claude AI to analyze and extract equation structure from text."""

    def __init__(self, api_key: Optional[str] = None):
        """
        Initialize the equation analyzer.

        Args:
            api_key: Anthropic API key (if None, reads from ANTHROPIC_API_KEY env var)
        """
        self.api_key = api_key or os.getenv("ANTHROPIC_API_KEY")
        if not self.api_key:
            raise ValueError(
                "API key required. Set ANTHROPIC_API_KEY environment variable "
                "or pass api_key parameter."
            )
        self.client = Anthropic(api_key=self.api_key)

    def analyze_equations(self, pdf_text: str) -> Dict:
        """
        Analyze equations from PDF text and extract structured information.

        Args:
            pdf_text: Text extracted from PDF

        Returns:
            Dictionary containing:
                - variables: List of variable names (e.g., ["ρ", "ρu", "ρv", "ρθ", "c1", "c2"])
                - num_equations: Number of equations
                - flux_x: List of flux terms in x-direction
                - flux_y: List of flux terms in y-direction
                - source_terms: List of source terms
                - problem_name: Suggested name for the problem
                - description: Brief description of the equation system
        """
        prompt = f"""You are analyzing a system of partial differential equations from a PDF document.
The equations are in the form:

∂_t(q) + ∂_x(F) + ∂_y(G) = S

where:
- q is the vector of conserved variables
- F is the flux in the x-direction
- G is the flux in the y-direction
- S is the source term vector

Here is the text extracted from the PDF:

{pdf_text}

Please analyze this text and extract the following information in JSON format:

{{
    "problem_name": "suggested_name_for_problem",
    "description": "brief description of the equation system",
    "num_equations": <number of equations>,
    "spatial_dimensions": <1, 2, or 3>,
    "variables": ["var1", "var2", ...],  // conserved variables in q vector
    "variable_descriptions": {{"var1": "description", ...}},
    "flux_x": ["flux_x_1", "flux_x_2", ...],  // F vector components
    "flux_y": ["flux_y_1", "flux_y_2", ...],  // G vector components (empty if 1D)
    "flux_z": ["flux_z_1", "flux_z_2", ...],  // H vector components (empty if not 3D)
    "source_terms": ["source_1", "source_2", ...],  // S vector components
    "primitive_variables": ["prim_var1", "prim_var2", ...],  // primitive variables (e.g., ρ, u, v, p)
    "has_pressure": true/false,
    "has_gravity": true/false,
    "equation_type": "euler/advection-diffusion/shallow-water/other"
}}

For example, for the compressible Euler equations with theta and tracers:
- variables: ["ρ", "ρu", "ρv", "ρθ", "c1", "c2"]
- flux_x: ["ρu", "ρu²+p", "ρvu", "ρθu", "c1u", "c2u"]
- flux_y: ["ρv", "ρuv", "ρv²+p", "ρθv", "c1v", "c2v"]
- source_terms: ["0", "0", "-ρg", "0", "0", "0"]

Extract the information and respond with ONLY the JSON object, no additional text."""

        response = self.client.messages.create(
            model="claude-sonnet-4-20250514",
            max_tokens=4096,
            messages=[{"role": "user", "content": prompt}]
        )

        # Extract JSON from response
        response_text = response.content[0].text.strip()

        # Try to parse JSON, handling potential markdown code blocks
        if response_text.startswith("```"):
            # Remove markdown code blocks
            lines = response_text.split("\n")
            response_text = "\n".join(
                line for line in lines
                if not line.strip().startswith("```")
            )

        try:
            equation_info = json.loads(response_text)
        except json.JSONDecodeError as e:
            raise ValueError(f"Failed to parse Claude response as JSON: {e}\nResponse: {response_text}")

        # Validate required fields
        required_fields = [
            "problem_name", "num_equations", "variables",
            "flux_x", "source_terms"
        ]
        for field in required_fields:
            if field not in equation_info:
                raise ValueError(f"Missing required field: {field}")

        return equation_info

    def validate_equation_info(self, eq_info: Dict) -> bool:
        """
        Validate that equation information is consistent.

        Args:
            eq_info: Equation information dictionary

        Returns:
            True if valid, raises ValueError otherwise
        """
        num_eqs = eq_info["num_equations"]
        num_vars = len(eq_info["variables"])
        num_flux_x = len(eq_info["flux_x"])
        num_source = len(eq_info["source_terms"])

        if num_vars != num_eqs:
            raise ValueError(
                f"Number of variables ({num_vars}) doesn't match "
                f"number of equations ({num_eqs})"
            )

        if num_flux_x != num_eqs:
            raise ValueError(
                f"Number of x-flux terms ({num_flux_x}) doesn't match "
                f"number of equations ({num_eqs})"
            )

        if num_source != num_eqs:
            raise ValueError(
                f"Number of source terms ({num_source}) doesn't match "
                f"number of equations ({num_eqs})"
            )

        spatial_dim = eq_info.get("spatial_dimensions", 2)
        if spatial_dim >= 2:
            num_flux_y = len(eq_info.get("flux_y", []))
            if num_flux_y != num_eqs:
                raise ValueError(
                    f"Number of y-flux terms ({num_flux_y}) doesn't match "
                    f"number of equations ({num_eqs})"
                )

        return True
