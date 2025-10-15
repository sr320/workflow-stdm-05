# STDM: Sparse Tensor Decomposition Models for Gene Expression Data

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive Python package for fitting sparse tensor decomposition models optimized for multi-species, multi-timepoint gene expression data. STDM provides efficient implementations of CP (CANDECOMP/PARAFAC) and Tucker decomposition methods with specialized features for genomic data analysis.

## Features

- 🧬 **Optimized for Gene Expression Data**: Handles ~10k genes across multiple species and time points
- 📊 **Multiple Decomposition Methods**: CP and Tucker decomposition with non-negative and sparse variants
- 🔄 **Flexible Data Input**: Supports both merged and species-separated data formats
- 📈 **Rich Visualization**: Built-in plotting for components, temporal patterns, and species comparisons
- 🚀 **High Performance**: Leverages TensorLy for efficient tensor operations
- 📦 **Easy Installation**: Uses `uv` for fast, reproducible package management
- ⏰ **Auto-Timestamped Outputs**: Every run creates a timestamped directory with all results
- 📝 **Automatic Reporting**: Generates comprehensive README for each analysis explaining all outputs
- ✅ **Data Validation**: Validates your data format before analysis
- 🔌 **Custom Data Support**: Easy to use with your own gene expression data

## Installation

### Using uv (Recommended)

[uv](https://github.com/astral-sh/uv) is a fast Python package installer and resolver. Install it first if you haven't:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Then install STDM:

```bash
# Clone the repository
git clone https://github.com/yourusername/workflow-stdm-05.git
cd workflow-stdm-05

# Create a virtual environment and install dependencies
uv venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install the package in development mode
uv pip install -e .

# Optional: Install development dependencies
uv pip install -e ".[dev]"
```

### Using pip

Alternatively, you can use pip:

```bash
pip install -e .
```

## Quick Start

### 1. Load and Prepare Data

```python
from stdm import DataLoader, TensorBuilder

# Load merged data (all species together)
loader = DataLoader("data/merge-vst")
data = loader.load_merged_data()

# Optional: preprocess
data = loader.filter_low_expression(data, min_expression=5.0, min_samples=10)
data = loader.normalize_data(data, method="zscore")

# Build 4D tensor: Genes × Individuals × Time Points × Species
builder = TensorBuilder()
tensor = builder.from_merged_data(data)
metadata = builder.get_metadata()
```

### 2. Fit Decomposition Model

```python
from stdm import CPDecomposition

# Fit CP decomposition
model = CPDecomposition(rank=10, random_state=42)
model.fit(tensor, n_iter_max=100, verbose=True)

# Get factors and reconstruction error
factors = model.get_factors()
error = model.compute_reconstruction_error(tensor)
print(f"Reconstruction error: {error:.6f}")
```

### 3. Visualize Results

```python
from stdm.visualization import (
    plot_components,
    plot_temporal_patterns,
    plot_species_comparison,
    plot_gene_loadings
)

# Plot all components
plot_components(
    factors, 
    mode_names=['Genes', 'Individuals', 'Time', 'Species'],
    n_components=5,
    save_path="results/components.png"
)

# Plot temporal patterns
plot_temporal_patterns(
    factors[2],  # Time mode
    metadata['timepoints'],
    save_path="results/temporal.png"
)

# Plot species-specific patterns
plot_species_comparison(
    factors[3],  # Species mode
    metadata['species'],
    save_path="results/species.png"
)
```

## Data Format

### Merged Data Format

CSV file with genes as rows and samples as columns. Sample names should follow the pattern:
`SPECIES-INDIVIDUAL-TIMEPOINT` (e.g., `ACR-139-TP1`, `POR-216-TP2`)

```csv
group_id,ACR-139-TP1,ACR-139-TP2,POR-216-TP1,POC-201-TP3,...
OG_00001,5.97,6.21,4.96,8.19,...
OG_00002,7.52,8.06,8.34,10.29,...
...
```

### Separate Data Format

One CSV file per species with genes as rows and samples as columns:

```csv
group_id,ACR.100.TP1,ACR.100.TP2,ACR.100.TP3,ACR.100.TP4,...
OG_00001,2.23,7.48,10.94,12.06,...
OG_00002,2.76,0.46,0.33,1.90,...
...
```

## Usage Examples

### Example 1: Fit CP Decomposition on Merged Data

```python
from pathlib import Path
from stdm import DataLoader, TensorBuilder, CPDecomposition
from stdm.utils import save_results

# Load data
loader = DataLoader("data/merge-vst")
data = loader.load_merged_data()
data = loader.normalize_data(data, method="zscore")

# Build tensor
builder = TensorBuilder()
tensor = builder.from_merged_data(data)

# Fit model
model = CPDecomposition(rank=15, non_negative=False, random_state=42)
model.fit(tensor, n_iter_max=100, verbose=True)

# Save results
save_results(
    model.get_factors(),
    builder.get_metadata(),
    "results/cp_merged",
    prefix="cp"
)
```

### Example 2: Compare Different Ranks

```python
from stdm import CPDecomposition
from stdm.visualization import plot_reconstruction_error

ranks = [5, 10, 15, 20, 25, 30]
errors = {}

for rank in ranks:
    model = CPDecomposition(rank=rank, random_state=42)
    model.fit(tensor, n_iter_max=100, verbose=False)
    errors[rank] = model.compute_reconstruction_error(tensor)

plot_reconstruction_error(
    errors,
    title="Reconstruction Error vs Rank",
    save_path="results/rank_comparison.png"
)
```

### Example 3: Analyze Separate Species Data

```python
# Load separate data for each species
loader = DataLoader("data/separate-log")
species_data = loader.load_separate_data()

# Build tensor with gene alignment
builder = TensorBuilder()
tensor = builder.from_separate_data(species_data, align_genes=True)

# Fit Tucker decomposition
from stdm import TuckerDecomposition

tucker_rank = (100, 20, 4, 3)  # Genes, Individuals, Time, Species
model = TuckerDecomposition(rank=tucker_rank, random_state=42)
model.fit(tensor, n_iter_max=50, verbose=True)

# Analyze results
core = model.get_core_tensor()
print(f"Core tensor shape: {core.shape}")
```

### Example 4: Identify Marker Genes

```python
from stdm.utils import identify_marker_genes

# Get gene factors from decomposition
gene_factor = model.get_factors()[0]

# Identify top marker genes for component 0
markers = identify_marker_genes(
    gene_factor,
    gene_names=metadata['genes'],
    component_idx=0,
    threshold=0.5,
    top_n=50
)

print("Top marker genes for component 1:")
for gene, loading in markers[:10]:
    print(f"  {gene}: {loading:.4f}")
```

## Using Your Own Data

STDM makes it easy to analyze your own gene expression data:

### Step 1: Validate Your Data

```bash
# Check if your data format is correct
stdm-fit --data-path /path/to/your/data.csv --validate-only
```

### Step 2: Run Analysis

```bash
# Analyze with auto-timestamped output and comprehensive reporting
stdm-fit \
    --data-path /path/to/your/data.csv \
    --method cp \
    --rank 15 \
    --normalize \
    --output-dir results/my_analysis
```

### Step 3: Review Results

Each run creates a timestamped directory (e.g., `results/my_analysis/cp_20241012_143022/`) containing:
- **ANALYSIS_README.md** - Comprehensive report explaining all outputs
- **Factor matrices** - CSV files with decomposition results
- **Metadata** - JSON file with dimension labels
- **Visualizations** - PNG plots of components

**Your Data Format:**

Merged CSV:
```csv
group_id,SPECIES-IND-TP1,SPECIES-IND-TP2,...
GENE001,5.2,6.1,...
GENE002,7.3,8.2,...
```

Or separate files (one per species):
```
my_data/
├── species1_normalized_expression.csv
├── species2_normalized_expression.csv
└── species3_normalized_expression.csv
```

See `examples/custom_data_analysis.py` for a complete example.

## Command-Line Interface

STDM provides a powerful command-line interface with auto-reporting:

```bash
# Basic analysis (auto-timestamped, with report)
stdm-fit \
    --data-path data/merge-vst/vst_counts_matrix.csv \
    --method cp \
    --rank 15 \
    --normalize

# Advanced options
stdm-fit \
    --data-path data/separate-log \
    --data-type separate \
    --method tucker \
    --rank 10 \
    --non-negative \
    --normalize \
    --n-iter 200 \
    --output-dir results/tucker

# Disable auto-timestamping (overwrite results)
stdm-fit \
    --data-path your_data.csv \
    --method cp \
    --rank 10 \
    --no-timestamp

# Validate data only (no analysis)
stdm-fit --data-path your_data.csv --validate-only
```

### CLI Options

| Option | Description | Default |
|--------|-------------|---------|
| `--data-path` | Path to data file or directory | Required |
| `--data-type` | Data type: merged, separate, or auto | auto |
| `--method` | Decomposition: cp or tucker | cp |
| `--rank` | Number of components | 10 |
| `--output-dir` | Base output directory | results |
| `--no-timestamp` | Disable timestamped subdirectories | False |
| `--no-report` | Disable automatic reporting | False |
| `--normalize` | Apply z-score normalization | False |
| `--non-negative` | Use non-negative decomposition | False |
| `--n-iter` | Maximum iterations | 100 |
| `--validate-only` | Only validate data format | False |

## Understanding Your Results

### Auto-Generated Analysis Report

Every analysis generates an `ANALYSIS_README.md` that explains:

1. **Output Files** - What each CSV and JSON file contains
2. **Interpretation Guide** - How to read factor matrices
3. **Quality Metrics** - Assessment of your model fit
4. **Optimal Parameters** - Recommendations for improving results
5. **Next Steps** - Code examples for further analysis

### Example Report Structure

```
results/cp_20241012_143022/
├── ANALYSIS_README.md              ← Read this first!
├── cp_factor_genes.csv             ← Gene loadings (10225 × 10)
├── cp_factor_individuals.csv       ← Individual patterns (28 × 10)
├── cp_factor_timepoints.csv        ← Temporal dynamics (4 × 10)
├── cp_factor_species.csv           ← Species patterns (3 × 10)
├── cp_metadata.json                ← Dimension labels
└── cp_components.png               ← Visualization
```

### Reading Factor Matrices

```python
import pandas as pd

# Load gene factor
genes = pd.read_csv("results/*/cp_factor_genes.csv", index_col=0)

# Get top genes for component 0
top_genes = genes["component_0"].abs().nlargest(50)
print(top_genes)

# Export for GO enrichment
top_genes.index.to_series().to_csv("component_0_genes.txt", index=False)
```

## Example Scripts

The `examples/` directory contains complete analysis scripts:

- **`fit_merged_data.py`**: Comprehensive analysis of merged data with CP and Tucker decomposition
- **`fit_separate_data.py`**: Analysis of species-separated data
- **`compare_models.py`**: Compare different ranks and methods
- **`custom_data_analysis.py`**: 🆕 Complete guide for using your own data with auto-reporting
- **`quickstart_notebook.ipynb`**: Interactive Jupyter tutorial

Run examples:

```bash
# Fit models on merged data
python examples/fit_merged_data.py

# Analyze separate data
python examples/fit_separate_data.py

# Compare models and ranks
python examples/compare_models.py

# 🆕 Use your own data (demonstrates all new features)
python examples/custom_data_analysis.py
```

## API Reference

### Core Classes

#### `DataLoader`

Load and preprocess gene expression data.

**Methods:**
- `load_merged_data(filename)`: Load merged CSV file
- `load_separate_data(species_files)`: Load separate species files
- `normalize_data(data, method)`: Normalize expression values
- `filter_low_expression(data, min_expression, min_samples)`: Filter low-expression genes
- `parse_sample_names(sample_names, delimiter)`: Parse sample metadata

#### `TensorBuilder`

Construct tensors from gene expression data.

**Methods:**
- `from_merged_data(data)`: Build tensor from merged data
- `from_separate_data(species_data)`: Build tensor from separate data
- `get_metadata()`: Get tensor dimension labels
- `unfold(mode)`: Unfold tensor along specified mode

#### `CPDecomposition`

CP (CANDECOMP/PARAFAC) decomposition.

**Parameters:**
- `rank`: Number of components
- `non_negative`: Enforce non-negativity constraint
- `l2_reg`: L2 regularization for sparsity
- `random_state`: Random seed

**Methods:**
- `fit(tensor, n_iter_max, tol, verbose)`: Fit the model
- `reconstruct()`: Reconstruct tensor from factors
- `get_factors()`: Get factor matrices
- `get_gene_components(top_n)`: Get top genes per component

#### `TuckerDecomposition`

Tucker decomposition.

**Parameters:**
- `rank`: Tuple of ranks for each mode
- `non_negative`: Enforce non-negativity constraint
- `random_state`: Random seed

**Methods:**
- `fit(tensor, n_iter_max, tol, verbose)`: Fit the model
- `reconstruct()`: Reconstruct tensor from factors
- `get_factors()`: Get factor matrices
- `get_core_tensor()`: Get core tensor

### Visualization Functions

- `plot_components(factors, mode_names, n_components)`: Plot factor loadings
- `plot_factor_heatmap(factor, labels, title)`: Heatmap of factor matrix
- `plot_reconstruction_error(errors)`: Plot error vs rank
- `plot_temporal_patterns(time_factor, timepoints)`: Plot temporal dynamics
- `plot_species_comparison(species_factor, species_names)`: Compare species patterns
- `plot_gene_loadings(gene_factor, gene_names, component_idx)`: Plot top gene loadings

### New Features (v0.1.0) 🆕

#### `DataValidator`

Validate user-provided gene expression data.

**Methods**:
```python
validator = DataValidator(strict=False)

# Validate merged CSV file
is_valid, messages = validator.validate_merged_data("data.csv")

# Validate separate species files
is_valid, messages = validator.validate_separate_data("data_dir/")

# Print validation report
validator.print_report()
```

#### `AnalysisReport`

Generate comprehensive analysis reports automatically.

**Usage**:
```python
from stdm import AnalysisReport

reporter = AnalysisReport(output_dir)
report_path = reporter.generate_run_report(
    method="CP",
    tensor_shape=tensor.shape,
    rank=10,
    reconstruction_error=0.25,
    explained_variance=0.82,
    runtime_seconds=120,
    metadata=metadata,
    factors=factors
)
# Creates detailed ANALYSIS_README.md in output_dir
```

#### `create_timestamped_output_dir`

Create timestamped directories for each analysis run.

**Usage**:
```python
from stdm import create_timestamped_output_dir

# Creates directory like: results/cp_20241012_143022/
output_dir = create_timestamped_output_dir("results", prefix="cp")
```

### Utility Functions

- `save_results(factors, metadata, output_dir)`: Save decomposition results
- `load_results(input_dir)`: Load saved results
- `compute_explained_variance(tensor, reconstructed)`: Calculate explained variance
- `normalize_factors(factors)`: Normalize factor matrices
- `identify_marker_genes(gene_factor, gene_names, component_idx)`: Find marker genes
- `cross_validate_rank(tensor, ranks, decomposition_class)`: Cross-validation for rank selection

## Testing

Run tests with pytest:

```bash
# Install test dependencies
uv pip install -e ".[dev]"

# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ --cov=stdm --cov-report=html
```

## Project Structure

```
workflow-stdm-05/
├── src/stdm/              # Main package
│   ├── __init__.py
│   ├── data_loader.py     # Data loading utilities
│   ├── tensor_builder.py  # Tensor construction
│   ├── decomposition.py   # CP and Tucker models
│   ├── visualization.py   # Plotting functions
│   ├── utils.py          # Utility functions
│   └── cli.py            # Command-line interface
├── examples/              # Example scripts
│   ├── fit_merged_data.py
│   ├── fit_separate_data.py
│   └── compare_models.py
├── tests/                # Unit tests
│   └── test_basic.py
├── data/                 # Gene expression data
│   ├── merge-vst/
│   ├── separate-log/
│   └── separate-syn/
├── pyproject.toml        # Package configuration
├── README.md            # This file
└── .python-version      # Python version specification
```

## Performance Tips

1. **Filtering**: Remove low-expression genes to reduce tensor size
2. **Normalization**: Z-score normalization often works best for decomposition
3. **Initialization**: Try different random seeds for reproducibility
4. **Rank Selection**: Use cross-validation or elbow method to choose optimal rank
5. **Iterations**: Start with fewer iterations for exploration, increase for final analysis
6. **Memory**: For large tensors, consider Tucker decomposition or reduce rank

## Troubleshooting

### Installation Issues

```bash
# If uv installation fails, try updating
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clear uv cache if needed
uv cache clean
```

### Memory Errors

```python
# Reduce tensor size by filtering genes
data = loader.filter_low_expression(data, min_expression=10.0)

# Or use lower rank
model = CPDecomposition(rank=5)  # Instead of rank=20
```

### Convergence Issues

```python
# Increase iterations
model.fit(tensor, n_iter_max=200)

# Try different initialization
model.fit(tensor, init='svd')  # Instead of 'random'

# Use non-negative decomposition for count data
model = CPDecomposition(rank=10, non_negative=True)
```

## Citation

If you use STDM in your research, please cite:

```bibtex
@software{stdm2024,
  title = {STDM: Sparse Tensor Decomposition Models for Gene Expression Data},
  author = {Your Name},
  year = {2024},
  url = {https://github.com/yourusername/workflow-stdm-05}
}
```

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- Built with [TensorLy](https://tensorly.org/) for tensor operations
- Inspired by research in multi-species genomic analysis
- Developed for coral reef genomics studies

## Support

For questions, issues, or feature requests:
- Open an issue on GitHub
- Email: your.email@example.com
- Documentation: [Link to docs]

## Changelog

### Version 0.1.0 (2024-10-12)

- Initial release
- CP and Tucker decomposition implementations
- Support for merged and separate data formats
- Comprehensive visualization suite
- CLI interface
- Example scripts and documentation
