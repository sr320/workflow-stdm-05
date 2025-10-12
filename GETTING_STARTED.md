# Getting Started with STDM

Welcome to STDM! This guide will help you get up and running with tensor decomposition analysis of your gene expression data.

## Table of Contents

1.  [Quick Setup](#quick-setup)
2.  [Your First Analysis](#your-first-analysis)
3.  [Understanding Your Results](#understanding-your-results)
4.  [Common Workflows](#common-workflows)
5.  [Troubleshooting](#troubleshooting)

## Quick Setup {#quick-setup}

### Option 1: Automated Setup (Recommended)

``` bash
cd /Users/sr320/Documents/GitHub/workflow-stdm-05
./setup.sh
```

This script will: - ✅ Detect and use `uv` (or fall back to `pip`) - ✅ Create a virtual environment - ✅ Install all dependencies - ✅ Verify the installation - ✅ Check for data files

### Option 2: Manual Setup

``` bash
# Install uv (if not already installed)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Create and activate virtual environment
uv venv
source .venv/bin/activate  # macOS/Linux
# .venv\Scripts\activate   # Windows

# Install package
uv pip install -e .
```

### Verify Installation

``` bash
# Test import
python -c "import stdm; print(f'STDM v{stdm.__version__} is ready!')"

# Test CLI
stdm-fit --help
```

## Your First Analysis {#your-first-analysis}

### Step 1: Explore Your Data

``` python
from stdm import DataLoader
from pathlib import Path

# Load your data
data_dir = Path("data/merge-vst")
loader = DataLoader(data_dir)
data = loader.load_merged_data()

print(f"Loaded: {data.shape[0]} genes × {data.shape[1]} samples")
print(f"Sample names: {data.columns[:5].tolist()}")
```

**Expected Output**:

```         
Loaded merged data: 10225 genes × 112 samples
Loaded: 10225 genes × 112 samples
Sample names: ['ACR-139-TP1', 'ACR-139-TP2', 'ACR-139-TP3', 'ACR-139-TP4', 'ACR-145-TP1']
```

### Step 2: Preprocess

``` python
# Filter low expression genes
data_filtered = loader.filter_low_expression(
    data, 
    min_expression=5.0,  # Minimum expression value
    min_samples=10       # Must be expressed in at least 10 samples
)

# Normalize (z-score across samples)
data_normalized = loader.normalize_data(
    data_filtered,
    method="zscore",
    axis=1  # Normalize each gene across samples
)

print(f"After filtering: {data_normalized.shape[0]} genes")
```

### Step 3: Build Tensor

``` python
from stdm import TensorBuilder

builder = TensorBuilder()
tensor = builder.from_merged_data(data_normalized)
metadata = builder.get_metadata()

print(f"\nTensor shape: {tensor.shape}")
print(f"Dimensions:")
print(f"  Genes: {len(metadata['genes'])}")
print(f"  Individuals: {len(metadata['individuals'])}")
print(f"  Time points: {metadata['timepoints']}")
print(f"  Species: {metadata['species']}")
```

**Expected Output**:

```         
Built tensor with shape: (9792, 28, 4, 3)
Dimensions:
  Genes: 9792
  Individuals: 28
  Time points: [1, 2, 3, 4]
  Species: ['ACR', 'POC', 'POR']
```

### Step 4: Fit Model

``` python
from stdm import CPDecomposition

# Create model
model = CPDecomposition(
    rank=10,           # Number of components
    random_state=42    # For reproducibility
)

# Fit to data
model.fit(tensor, n_iter_max=100, verbose=True)

# Get results
factors = model.get_factors()
error = model.compute_reconstruction_error(tensor)

print(f"\nReconstruction error: {error:.6f}")
```

### Step 5: Visualize

``` python
from stdm.visualization import (
    plot_components,
    plot_temporal_patterns,
    plot_species_comparison
)

# Plot all components
plot_components(
    factors,
    mode_names=['Genes', 'Individuals', 'Time', 'Species'],
    n_components=5
)

# Plot temporal patterns
plot_temporal_patterns(
    factors[2],  # Time mode factor
    metadata['timepoints']
)

# Compare species
plot_species_comparison(
    factors[3],  # Species mode factor
    metadata['species']
)
```

## Understanding Your Results {#understanding-your-results}

### What are the Factors?

After decomposition, you get 4 factor matrices:

1.  **Gene Factor** (shape: genes × rank)
    -   Shows which genes contribute to each component
    -   High values = gene is strongly associated with that component
2.  **Individual Factor** (shape: individuals × rank)
    -   Shows expression patterns for each individual
    -   Captures individual-specific variation
3.  **Time Factor** (shape: time points × rank)
    -   Shows temporal dynamics of each component
    -   Reveals patterns that change over time
4.  **Species Factor** (shape: species × rank)
    -   Shows species-specific patterns
    -   Identifies conserved vs. species-specific components

### Interpreting Components

Each component (column in the factors) represents a pattern in the data:

``` python
# Analyze component 0
print("Component 1 Analysis:")

# Temporal pattern
time_pattern = factors[2][:, 0]
print(f"\nTemporal pattern:")
for tp, val in zip(metadata['timepoints'], time_pattern):
    print(f"  TP{tp}: {val:.3f}")

# Species pattern
species_pattern = factors[3][:, 0]
print(f"\nSpecies pattern:")
for species, val in zip(metadata['species'], species_pattern):
    print(f"  {species}: {val:.3f}")

# Top genes
from stdm.utils import identify_marker_genes
markers = identify_marker_genes(
    factors[0],
    metadata['genes'],
    component_idx=0,
    top_n=10
)
print(f"\nTop 10 marker genes:")
for gene, loading in markers:
    print(f"  {gene}: {loading:.4f}")
```

## Common Workflows {#common-workflows}

### Workflow 1: Find Optimal Rank

``` python
from stdm.visualization import plot_reconstruction_error

# Test different ranks
ranks = [5, 10, 15, 20, 25, 30]
errors = {}

for rank in ranks:
    model = CPDecomposition(rank=rank, random_state=42)
    model.fit(tensor, n_iter_max=100, verbose=False)
    errors[rank] = model.compute_reconstruction_error(tensor)

# Plot results
plot_reconstruction_error(errors, save_path="results/rank_comparison.png")

# Choose rank where error plateaus (elbow method)
```

### Workflow 2: Compare CP vs Tucker

``` python
from stdm import TuckerDecomposition

# CP decomposition
cp_model = CPDecomposition(rank=10)
cp_model.fit(tensor, n_iter_max=100, verbose=False)
cp_error = cp_model.compute_reconstruction_error(tensor)

# Tucker decomposition
tucker_model = TuckerDecomposition(rank=(100, 20, 4, 3))
tucker_model.fit(tensor, n_iter_max=50, verbose=False)
tucker_error = tucker_model.compute_reconstruction_error(tensor)

print(f"CP error: {cp_error:.6f}")
print(f"Tucker error: {tucker_error:.6f}")
```

### Workflow 3: Analyze Separate Species

``` python
# Load separate species data
loader = DataLoader("data/separate-log")
species_data = loader.load_separate_data()

# Build tensor with gene alignment
builder = TensorBuilder()
tensor = builder.from_separate_data(
    species_data,
    align_genes=True  # Use common genes only
)

# Fit model
model = CPDecomposition(rank=10)
model.fit(tensor)
```

### Workflow 4: Export Results

``` python
from stdm.utils import save_results

# Save all factors and metadata
save_results(
    factors,
    metadata,
    output_dir="results/my_analysis",
    prefix="cp_rank10"
)

# This creates:
# - cp_rank10_factor_genes.csv
# - cp_rank10_factor_individuals.csv
# - cp_rank10_factor_timepoints.csv
# - cp_rank10_factor_species.csv
# - cp_rank10_metadata.json
```

## Using the Command-Line Interface

For quick analyses without writing code:

``` bash
# Basic analysis
stdm-fit \
    --data-path data/separate-log \
    --data-type merged \
    --method cp \
    --rank 15 \
    --output-dir results/cli_test

# With normalization and non-negativity
stdm-fit \
    --data-dir data/merge-vst \
    --data-type merged \
    --method cp \
    --rank 15 \
    --normalize \
    --non-negative \
    --output-dir results/cli_test_nn

# Tucker decomposition
stdm-fit \
    --data-dir data/separate-log \
    --data-type separate \
    --method tucker \
    --rank 10 \
    --normalize \
    --output-dir results/tucker_test
```

## Running Example Scripts

The `examples/` directory contains complete analysis pipelines:

``` bash
# Activate virtual environment first
source .venv/bin/activate

# Run merged data analysis
python examples/fit_merged_data.py

# Run separate species analysis  
python examples/fit_separate_data.py

# Compare different models and ranks
python examples/compare_models.py
```

## Working with Jupyter Notebooks

``` bash
# Install Jupyter (if not already installed)
uv pip install jupyter

# Start Jupyter
jupyter notebook

# Open examples/quickstart_notebook.ipynb
```

## Troubleshooting {#troubleshooting}

### Issue: "Module not found" error

``` bash
# Make sure virtual environment is activated
source .venv/bin/activate

# Reinstall package
uv pip install -e .
```

### Issue: TensorLy convergence warnings

``` python
# Increase iterations
model.fit(tensor, n_iter_max=200)

# Try different initialization
model.fit(tensor, init='svd')
```

### Issue: Out of memory

``` python
# Filter more aggressively
data = loader.filter_low_expression(data, min_expression=10.0)

# Use lower rank
model = CPDecomposition(rank=5)  # instead of rank=20
```

### Issue: Poor reconstruction

``` python
# Try non-negative decomposition
model = CPDecomposition(rank=10, non_negative=True)

# Increase rank
model = CPDecomposition(rank=20)

# Try Tucker instead of CP
model = TuckerDecomposition(rank=(100, 20, 4, 3))
```

## Next Steps

1.  **Read the full documentation**: [README.md](README.md)
2.  **Study the examples**: Explore `examples/` directory
3.  **Run tests**: `pytest tests/ -v`
4.  **Check package summary**: [PACKAGE_SUMMARY.md](PACKAGE_SUMMARY.md)
5.  **Read installation guide**: [INSTALL.md](INSTALL.md)

## Getting Help

-   **Documentation**: README.md, INSTALL.md, PACKAGE_SUMMARY.md
-   **Examples**: examples/ directory
-   **Tests**: tests/ directory for usage patterns
-   **GitHub Issues**: Report bugs or request features
-   **Email**: [your.email\@example.com](mailto:your.email@example.com){.email}

## Quick Reference

### Common Imports

``` python
# Core functionality
from stdm import DataLoader, TensorBuilder, CPDecomposition, TuckerDecomposition

# Visualization
from stdm.visualization import (
    plot_components,
    plot_temporal_patterns,
    plot_species_comparison,
    plot_gene_loadings
)

# Utilities
from stdm.utils import (
    save_results,
    compute_explained_variance,
    identify_marker_genes
)
```

### Typical Analysis Pipeline

``` python
# 1. Load data
loader = DataLoader("data/merge-vst")
data = loader.load_merged_data()

# 2. Preprocess
data = loader.filter_low_expression(data)
data = loader.normalize_data(data, method="zscore")

# 3. Build tensor
builder = TensorBuilder()
tensor = builder.from_merged_data(data)

# 4. Fit model
model = CPDecomposition(rank=10)
model.fit(tensor)

# 5. Analyze
factors = model.get_factors()
error = model.compute_reconstruction_error(tensor)

# 6. Visualize
plot_components(factors)

# 7. Save
save_results(factors, builder.get_metadata(), "results/")
```

Happy analyzing! 🧬📊
