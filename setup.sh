#!/bin/bash

# STDM Package Setup Script
# This script automates the installation of the STDM package

set -e  # Exit on error

echo "============================================"
echo "STDM Package Setup"
echo "============================================"
echo ""

# Check if uv is installed
if command -v uv &> /dev/null; then
    echo "✓ uv is installed"
    USE_UV=true
else
    echo "⚠ uv not found. Will use pip instead."
    echo "  To install uv: curl -LsSf https://astral.sh/uv/install.sh | sh"
    USE_UV=false
fi

# Check Python version
echo ""
echo "Checking Python version..."
PYTHON_VERSION=$(python --version 2>&1 | awk '{print $2}')
echo "  Python version: $PYTHON_VERSION"

# Create virtual environment
echo ""
echo "Creating virtual environment..."
if [ "$USE_UV" = true ]; then
    uv venv
else
    python -m venv .venv
fi
echo "  ✓ Virtual environment created"

# Activate virtual environment
echo ""
echo "Activating virtual environment..."
if [[ "$OSTYPE" == "msys" || "$OSTYPE" == "win32" ]]; then
    source .venv/Scripts/activate
else
    source .venv/bin/activate
fi
echo "  ✓ Virtual environment activated"

# Install package
echo ""
echo "Installing STDM package..."
if [ "$USE_UV" = true ]; then
    uv pip install -e .
else
    pip install --upgrade pip
    pip install -e .
fi
echo "  ✓ Package installed"

# Optional: Install dev dependencies
read -p "Install development dependencies (pytest, jupyter, etc.)? (y/n) " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Installing development dependencies..."
    if [ "$USE_UV" = true ]; then
        uv pip install -e ".[dev]"
    else
        pip install -e ".[dev]"
    fi
    echo "  ✓ Development dependencies installed"
fi

# Verify installation
echo ""
echo "Verifying installation..."
python -c "import stdm; print(f'  ✓ STDM version: {stdm.__version__}')" || {
    echo "  ✗ Installation verification failed"
    exit 1
}

# Check if data exists
echo ""
echo "Checking for data files..."
if [ -d "data/merge-vst" ] && [ -f "data/merge-vst/vst_counts_matrix.csv" ]; then
    echo "  ✓ Merged data found"
else
    echo "  ⚠ Merged data not found in data/merge-vst/"
fi

if [ -d "data/separate-log" ]; then
    NUM_FILES=$(ls data/separate-log/*_normalized_expression.csv 2>/dev/null | wc -l)
    echo "  ✓ Found $NUM_FILES separate data files in data/separate-log/"
else
    echo "  ⚠ Separate data not found in data/separate-log/"
fi

# Create results directory
echo ""
echo "Creating results directory..."
mkdir -p results
echo "  ✓ Results directory created"

# Print next steps
echo ""
echo "============================================"
echo "Setup Complete!"
echo "============================================"
echo ""
echo "Next steps:"
echo ""
echo "1. Activate the virtual environment:"
if [[ "$OSTYPE" == "msys" || "$OSTYPE" == "win32" ]]; then
    echo "   .venv\\Scripts\\activate"
else
    echo "   source .venv/bin/activate"
fi
echo ""
echo "2. Try the quick start:"
echo "   python -c \"from stdm import DataLoader; print('STDM ready!')\""
echo ""
echo "3. Run an example:"
echo "   python examples/fit_merged_data.py"
echo ""
echo "4. Use the CLI:"
echo "   stdm-fit --help"
echo ""
echo "5. Run tests:"
echo "   pytest tests/ -v"
echo ""
echo "For more information, see README.md"
echo "============================================"

