# STDM Package Summary

## Package Overview

**STDM (Sparse Tensor Decomposition Models)** is a comprehensive Python package designed for analyzing multi-species, multi-timepoint gene expression data using tensor decomposition methods.

### Key Features

✅ **Data Loading & Preprocessing**
- Support for merged and species-separated CSV data
- Automatic sample name parsing (SPECIES-ID-TPX format)
- Data normalization (z-score, min-max, log2)
- Low-expression gene filtering

✅ **Tensor Construction**
- 4D tensors: Genes × Individuals × Time Points × Species
- Automatic gene alignment across species
- Metadata tracking for all dimensions

✅ **Decomposition Methods**
- **CP Decomposition**: CANDECOMP/PARAFAC with optional non-negativity
- **Tucker Decomposition**: Higher-order SVD with core tensor
- Sparse regularization options
- Reconstruction error metrics

✅ **Visualization Suite**
- Component loadings across all modes
- Temporal pattern plots
- Species comparison charts
- Gene loading heatmaps
- Reconstruction error curves

✅ **Analysis Tools**
- Marker gene identification
- Cross-validation for rank selection
- Model comparison utilities
- Result export (CSV + JSON)

## Package Structure

```
workflow-stdm-05/
├── src/stdm/                    # Main package source
│   ├── __init__.py             # Package initialization
│   ├── data_loader.py          # Data loading utilities
│   ├── tensor_builder.py       # Tensor construction
│   ├── decomposition.py        # CP & Tucker models
│   ├── visualization.py        # Plotting functions
│   ├── utils.py                # Helper utilities
│   └── cli.py                  # Command-line interface
│
├── examples/                    # Example scripts
│   ├── fit_merged_data.py      # Merged data analysis
│   ├── fit_separate_data.py    # Separate species analysis
│   ├── compare_models.py       # Model comparison
│   └── quickstart_notebook.ipynb # Interactive tutorial
│
├── tests/                       # Unit tests
│   └── test_basic.py           # Basic functionality tests
│
├── data/                        # Gene expression data
│   ├── merge-vst/              # Merged species data
│   │   └── vst_counts_matrix.csv
│   ├── separate-log/           # Log-transformed separate data
│   │   ├── apul_normalized_expression.csv
│   │   ├── peve_normalized_expression.csv
│   │   └── ptua_normalized_expression.csv
│   └── separate-syn/           # Syn-normalized separate data
│       ├── apul_normalized_expression.csv
│       ├── peve_normalized_expression.csv
│       └── ptua_normalized_expression.csv
│
├── pyproject.toml              # Package configuration
├── README.md                   # Main documentation
├── INSTALL.md                  # Installation guide
├── LICENSE                     # MIT license
├── .python-version             # Python version spec
└── .gitignore                  # Git ignore rules
```

## Core Components

### 1. DataLoader Class

**Purpose**: Load and preprocess gene expression data

**Key Methods**:
```python
loader = DataLoader("data/merge-vst")

# Load data
data = loader.load_merged_data()
species_data = loader.load_separate_data()

# Preprocess
data = loader.normalize_data(data, method="zscore")
data = loader.filter_low_expression(data, min_expression=5.0)

# Parse metadata
species, individuals, timepoints = loader.parse_sample_names(samples)
```

### 2. TensorBuilder Class

**Purpose**: Construct multi-dimensional tensors from data

**Key Methods**:
```python
builder = TensorBuilder()

# Build from merged data (4D tensor)
tensor = builder.from_merged_data(data)

# Build from separate data
tensor = builder.from_separate_data(species_data, align_genes=True)

# Access metadata
metadata = builder.get_metadata()
# Returns: {'genes': [...], 'individuals': [...], 'timepoints': [...], 'species': [...]}

# Unfold tensor
matrix = builder.unfold(mode=0)  # Unfold along genes
```

### 3. CPDecomposition Class

**Purpose**: Fit CP decomposition models

**Key Parameters**:
- `rank`: Number of components (int)
- `non_negative`: Enforce non-negativity (bool)
- `l2_reg`: Sparsity regularization (float)
- `random_state`: Random seed (int)

**Key Methods**:
```python
model = CPDecomposition(rank=10, random_state=42)
model.fit(tensor, n_iter_max=100, verbose=True)

# Get results
factors = model.get_factors()  # List of factor matrices
error = model.compute_reconstruction_error(tensor)
reconstructed = model.reconstruct()

# Analyze genes
top_genes = model.get_gene_components(top_n=50)
```

### 4. TuckerDecomposition Class

**Purpose**: Fit Tucker decomposition models

**Key Parameters**:
- `rank`: Tuple of ranks for each mode
- `non_negative`: Enforce non-negativity (bool)
- `random_state`: Random seed (int)

**Key Methods**:
```python
model = TuckerDecomposition(rank=(100, 20, 4, 3), random_state=42)
model.fit(tensor, n_iter_max=50, verbose=True)

# Get results
factors = model.get_factors()
core = model.get_core_tensor()
error = model.compute_reconstruction_error(tensor)
```

### 5. Visualization Functions

**Available Plots**:
```python
from stdm.visualization import (
    plot_components,
    plot_factor_heatmap,
    plot_reconstruction_error,
    plot_temporal_patterns,
    plot_species_comparison,
    plot_gene_loadings
)

# Plot all components
plot_components(factors, mode_names=['Genes', 'Individuals', 'Time', 'Species'])

# Plot temporal dynamics
plot_temporal_patterns(time_factor, timepoints)

# Compare species
plot_species_comparison(species_factor, species_names)

# Top genes for a component
plot_gene_loadings(gene_factor, gene_names, component_idx=0)
```

### 6. Utility Functions

**Key Utilities**:
```python
from stdm.utils import (
    save_results,
    load_results,
    compute_explained_variance,
    normalize_factors,
    identify_marker_genes,
    compare_components,
    cross_validate_rank
)

# Save/load results
save_results(factors, metadata, "results/", prefix="cp")
factors, metadata = load_results("results/", prefix="cp")

# Analysis
explained_var = compute_explained_variance(tensor, reconstructed)
markers = identify_marker_genes(gene_factor, gene_names, component_idx=0)
similarity = compare_components(factor1, factor2, comp1_idx, comp2_idx)

# Model selection
errors = cross_validate_rank(tensor, ranks=[5,10,15,20], CPDecomposition)
```

## Data Format Requirements

### Merged Data

**File**: Single CSV with all species combined

**Format**:
- Rows: Genes (e.g., OG_00001, OG_00002, ...)
- Columns: Samples with naming pattern `SPECIES-INDIVIDUAL-TPX`
  - Example: `ACR-139-TP1`, `POR-216-TP2`, `POC-201-TP3`
- Values: Normalized expression values (float)

**Example**:
```csv
group_id,ACR-139-TP1,ACR-139-TP2,POR-216-TP1,POC-201-TP3
OG_00001,5.97,6.21,4.96,8.19
OG_00002,7.52,8.06,8.34,10.29
```

### Separate Data

**Files**: One CSV per species

**Format**:
- Rows: Genes (e.g., OG_00001, OG_00002, ...)
- Columns: Samples with naming pattern `SPECIES.INDIVIDUAL.TPX`
  - Example: `ACR.100.TP1`, `ACR.100.TP2`, `ACR.100.TP3`
- Values: Normalized expression values (float)

**Example** (apul_normalized_expression.csv):
```csv
group_id,ACR.100.TP1,ACR.100.TP2,ACR.100.TP3,ACR.100.TP4
OG_00001,2.23,7.48,10.94,12.06
OG_00002,2.76,0.46,0.33,1.90
```

## Installation

### Quick Install with uv (Recommended)

```bash
# Install uv
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clone and install
git clone https://github.com/yourusername/workflow-stdm-05.git
cd workflow-stdm-05
uv venv
source .venv/bin/activate
uv pip install -e .
```

### Install with pip

```bash
cd workflow-stdm-05
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

## Quick Start

### Python API

```python
from stdm import DataLoader, TensorBuilder, CPDecomposition

# Load and preprocess
loader = DataLoader("data/merge-vst")
data = loader.load_merged_data()
data = loader.normalize_data(data, method="zscore")

# Build tensor
builder = TensorBuilder()
tensor = builder.from_merged_data(data)

# Fit model
model = CPDecomposition(rank=10)
model.fit(tensor)

# Analyze
factors = model.get_factors()
error = model.compute_reconstruction_error(tensor)
```

### Command-Line Interface

```bash
stdm-fit \
    --data-dir data/merge-vst \
    --data-type merged \
    --method cp \
    --rank 15 \
    --output-dir results/
```

## Use Cases

### 1. Multi-Species Gene Expression Analysis

**Goal**: Identify conserved and species-specific expression patterns

**Approach**:
- Load merged data from all species
- Fit CP decomposition with appropriate rank
- Analyze species factor to identify species-specific components
- Extract marker genes for each pattern

### 2. Temporal Dynamics

**Goal**: Discover temporal gene expression patterns

**Approach**:
- Build tensor with time as one dimension
- Fit decomposition model
- Visualize temporal factor loadings
- Identify genes with similar temporal profiles

### 3. Cross-Species Comparison

**Goal**: Compare gene expression across species at different time points

**Approach**:
- Load separate species data
- Align genes across species
- Fit Tucker decomposition for flexibility
- Compare factors across species

### 4. Marker Gene Discovery

**Goal**: Find genes characteristic of specific conditions

**Approach**:
- Fit decomposition model
- Identify components of interest
- Extract top genes from gene factor
- Export for functional enrichment analysis

## Testing

```bash
# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ --cov=stdm --cov-report=html

# Run specific test
pytest tests/test_basic.py::TestCPDecomposition -v
```

## Performance Considerations

| Data Size | Genes | Samples | Recommended Rank | Expected Time |
|-----------|-------|---------|------------------|---------------|
| Small     | 1k    | 50      | 5-10            | < 1 min       |
| Medium    | 5k    | 100     | 10-20           | 5-10 min      |
| Large     | 10k   | 200     | 15-30           | 15-30 min     |

**Tips**:
- Filter low-expression genes before decomposition
- Start with lower ranks for exploration
- Use Tucker for very large tensors
- Consider non-negative decomposition for count data

## Dependencies

**Core**:
- numpy >= 1.24.0
- pandas >= 2.0.0
- tensorly >= 0.8.0
- scipy >= 1.10.0

**Visualization**:
- matplotlib >= 3.7.0
- seaborn >= 0.12.0

**ML**:
- scikit-learn >= 1.3.0

**Development** (optional):
- pytest >= 7.0.0
- black >= 23.0.0
- ruff >= 0.1.0
- jupyter >= 1.0.0

## Citation

```bibtex
@software{stdm2024,
  title = {STDM: Sparse Tensor Decomposition Models for Gene Expression Data},
  author = {Your Name},
  year = {2024},
  url = {https://github.com/yourusername/workflow-stdm-05}
}
```

## Support & Resources

- **Documentation**: README.md
- **Installation Guide**: INSTALL.md
- **Examples**: examples/ directory
- **Tests**: tests/ directory
- **Issues**: GitHub Issues
- **License**: MIT

## Version History

- **v0.1.0** (2024-10-12): Initial release with CP and Tucker decomposition, comprehensive visualization, and example scripts

