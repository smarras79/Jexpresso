#!/bin/bash
# Jexpresso Installation Wrapper Script
# This script provides a convenient way to run the Julia installation script

set -e  # Exit on error

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "==============================================="
echo "  Jexpresso Dependency Installer"
echo "==============================================="
echo ""

# Check if Julia is installed
if ! command -v julia &> /dev/null; then
    echo -e "${RED}Error: Julia is not installed or not in PATH${NC}"
    echo "Please install Julia 1.11.2 or higher from https://julialang.org/downloads/"
    exit 1
fi

# Check Julia version
echo "Checking Julia version..."
JULIA_VERSION=$(julia --version | grep -oP '\d+\.\d+\.\d+' | head -1)
echo "Found Julia version: $JULIA_VERSION"

# Compare versions (basic check)
REQUIRED_MAJOR=1
REQUIRED_MINOR=11
REQUIRED_PATCH=2

CURRENT_MAJOR=$(echo $JULIA_VERSION | cut -d. -f1)
CURRENT_MINOR=$(echo $JULIA_VERSION | cut -d. -f2)
CURRENT_PATCH=$(echo $JULIA_VERSION | cut -d. -f3)

if [ "$CURRENT_MAJOR" -lt "$REQUIRED_MAJOR" ] || \
   ([ "$CURRENT_MAJOR" -eq "$REQUIRED_MAJOR" ] && [ "$CURRENT_MINOR" -lt "$REQUIRED_MINOR" ]) || \
   ([ "$CURRENT_MAJOR" -eq "$REQUIRED_MAJOR" ] && [ "$CURRENT_MINOR" -eq "$REQUIRED_MINOR" ] && [ "$CURRENT_PATCH" -lt "$REQUIRED_PATCH" ]); then
    echo -e "${YELLOW}Warning: Julia version $JULIA_VERSION is older than recommended version 1.11.2${NC}"
    echo "Consider upgrading Julia for best compatibility"
    read -p "Continue anyway? (y/N) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# Run the Julia installation script
echo ""
echo "Running Julia dependency installation script..."
echo ""

julia install_dependencies.jl
EXIT_CODE=$?

echo ""
if [ $EXIT_CODE -eq 0 ]; then
    echo -e "${GREEN}Installation completed successfully!${NC}"
    echo ""
    echo "Next steps:"
    echo "  1. To use Jexpresso, run: julia --project=."
    echo "  2. See INSTALLATION.md for detailed usage instructions"
    echo "  3. Check documentation at: https://smarras79.github.io/Jexpresso/dev/"
else
    echo -e "${RED}Installation failed with errors${NC}"
    echo "Please check the error messages above and see INSTALLATION.md for troubleshooting"
    exit $EXIT_CODE
fi
