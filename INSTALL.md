# Installation Guide for STDM

This guide provides detailed installation instructions for the STDM package.

## Prerequisites

- Python 3.9 or later
- Operating System: macOS, Linux, or Windows

## Method 1: Using uv (Recommended)

[uv](https://github.com/astral-sh/uv) is an extremely fast Python package installer and resolver, written in Rust.

### Step 1: Install uv

**macOS and Linux:**
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**Windows:**
```powershell
powershell -c "irm https://astral.sh/uv/install.ps1 | iex"
```

**Alternative (using pip):**
```bash
pip install uv
```

### Step 2: Clone the Repository

```bash
git clone https://github.com/yourusername/workflow-stdm-05.git
cd workflow-stdm-05
```

### Step 3: Create Virtual Environment

```bash
# Create a new virtual environment
uv venv

# Activate the virtual environment
# On macOS/Linux:
source .venv/bin/activate

# On Windows:
.venv\Scripts\activate
```

### Step 4: Install the Package

```bash
# Install in editable mode with all dependencies
uv pip install -e .

# Optional: Install development dependencies (for testing and development)
uv pip install -e ".[dev]"
```

### Step 5: Verify Installation

```bash
# Test the installation
python -c "import stdm; print(stdm.__version__)"

# Run the CLI
stdm-fit --help
```

## Method 2: Using pip

### Step 1: Create Virtual Environment

```bash
# Navigate to project directory
cd workflow-stdm-05

# Create virtual environment
python -m venv .venv

# Activate virtual environment
# On macOS/Linux:
source .venv/bin/activate

# On Windows:
.venv\Scripts\activate
```

### Step 2: Install Dependencies

```bash
# Upgrade pip
pip install --upgrade pip

# Install the package
pip install -e .

# Optional: Install development dependencies
pip install -e ".[dev]"
```

### Step 3: Verify Installation

```bash
python -c "import stdm; print(stdm.__version__)"
```

## Method 3: Direct Dependency Installation

If you prefer to install dependencies manually:

```bash
pip install numpy>=1.24.0
pip install pandas>=2.0.0
pip install tensorly>=0.8.0
pip install scikit-learn>=1.3.0
pip install matplotlib>=3.7.0
pip install seaborn>=0.12.0
pip install scipy>=1.10.0
```

Then install the package:
```bash
pip install -e .
```

## Development Installation

For development work (includes testing and linting tools):

```bash
# Using uv
uv pip install -e ".[dev]"

# Using pip
pip install -e ".[dev]"
```

This installs additional packages:
- pytest (testing framework)
- pytest-cov (code coverage)
- black (code formatter)
- ruff (linter)
- jupyter (notebook support)
- ipykernel (Jupyter kernel)

## Troubleshooting

### Issue: uv not found after installation

**Solution:** Add uv to your PATH:

```bash
# On macOS/Linux, add to ~/.bashrc or ~/.zshrc:
export PATH="$HOME/.cargo/bin:$PATH"

# Reload shell configuration
source ~/.bashrc  # or source ~/.zshrc
```

### Issue: TensorLy installation fails

**Solution:** Try installing with specific versions:

```bash
pip install tensorly==0.8.1
```

Or build from source:
```bash
pip install git+https://github.com/tensorly/tensorly.git
```

### Issue: NumPy compatibility errors

**Solution:** Ensure compatible versions:

```bash
pip install "numpy>=1.24.0,<2.0.0"
```

### Issue: Permission denied on macOS/Linux

**Solution:** Don't use sudo. Use virtual environments instead:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

### Issue: Missing C compiler for some dependencies

**Solution:**

**On macOS:**
```bash
xcode-select --install
```

**On Ubuntu/Debian:**
```bash
sudo apt-get install build-essential python3-dev
```

**On CentOS/RHEL:**
```bash
sudo yum install gcc python3-devel
```

### Issue: Memory errors during installation

**Solution:** Install dependencies one at a time:

```bash
pip install numpy
pip install pandas
pip install tensorly
pip install -e .
```

## Running Tests

After installation, verify everything works:

```bash
# Run all tests
pytest tests/ -v

# Run with coverage report
pytest tests/ --cov=stdm --cov-report=html

# View coverage report
open htmlcov/index.html  # macOS
xdg-open htmlcov/index.html  # Linux
```

## Running Examples

Try the example scripts:

```bash
# Make sure you have data in the data directory
ls data/merge-vst/

# Run merged data example
python examples/fit_merged_data.py

# Run separate data example
python examples/fit_separate_data.py

# Compare models
python examples/compare_models.py
```

## Uninstallation

To uninstall STDM:

```bash
pip uninstall stdm
```

To remove the virtual environment:

```bash
# Deactivate first
deactivate

# Remove the directory
rm -rf .venv
```

## Updating

To update to the latest version:

```bash
# Pull latest changes
git pull origin main

# Update dependencies
uv pip install -e . --force-reinstall

# Or with pip
pip install -e . --force-reinstall
```

## Docker Installation (Optional)

For isolated environments, you can use Docker:

```dockerfile
# Dockerfile
FROM python:3.11-slim

WORKDIR /app

COPY . .

RUN pip install -e .

CMD ["python"]
```

Build and run:
```bash
docker build -t stdm .
docker run -it -v $(pwd)/data:/app/data stdm
```

## Cloud/HPC Installation

For HPC clusters or cloud environments:

```bash
# Load Python module (if needed)
module load python/3.11

# Create virtual environment in your home directory
python -m venv ~/stdm-env
source ~/stdm-env/bin/activate

# Install with pip (uv may not be available)
pip install -e .
```

## Next Steps

After successful installation:

1. Read the [README.md](README.md) for usage examples
2. Try the [Quick Start](README.md#quick-start) guide
3. Explore the [examples/](examples/) directory
4. Check the API documentation in [README.md](README.md#api-reference)

## Getting Help

If you encounter issues:

1. Check this troubleshooting guide
2. Search existing GitHub issues
3. Open a new issue with:
   - Your operating system
   - Python version (`python --version`)
   - Error message
   - Steps to reproduce

## Support

- GitHub Issues: [Link to issues]
- Email: your.email@example.com
- Documentation: [Link to docs]

