"""
    Stage 0 Test Problem: RT with Spatial AMR Detection

    This problem tests the spatial non-conformity detection added in Stage 0.
    It creates a spatial-refined mesh and verifies that the RT solver correctly
    detects and reports that spatial constraint support is not yet implemented.

    Configuration:
    - Small uniform mesh (2x1x1 elements)
    - Spatial refinement enabled (via :lRT_spatial_amr => true)
    - Angular refinement disabled for simplicity
    - Manufactured solution for testing

    Expected behavior:
    1. Mesh loads successfully (uniform)
    2. Spatial refinement is applied (creates hanging nodes)
    3. RT solver detects spatial non-conformity
    4. Clear error message displayed
    5. Test passes (detection working correctly)
"""

# Placeholder for any special initialization if needed
# Currently using defaults from user_inputs, user_bc, user_flux, user_primitives, user_source
