# STDM: Sparse Tensor Decomposition Models for Gene Expression Data

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Python package for analyzing multi-species, multi-timepoint gene expression data using tensor decomposition methods. Optimized for ~10k genes across multiple species and time points.

## Quick Install

### Option 1: Automated Setup (Recommended)
```bash
cd workflow-stdm-05
./setup.sh
```

### Option 2: Manual Install
```bash
# Install uv (fast Python package manager)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Setup environment
git clone https://github.com/yourusername/workflow-stdm-05.git
cd workflow-stdm-05
uv venv
source .venv/bin/activate
uv pip install -e .
```

## Quick Start: Run on Your Data

STDM works with 3 types of datasets in the `data/` directory:

### Dataset 1: Merged Data (All Species Together)
```bash
# Basic analysis
stdm-fit --data-path data/merge-vst/vst_counts_matrix.csv --method cp --rank 15

# With normalization and auto-reporting
stdm-fit --data-path data/merge-vst/vst_counts_matrix.csv --method cp --rank 15 --normalize
```

### Dataset 2: Separate Species (Log Normalized)
```bash
# Auto-detects separate species format
stdm-fit --data-path data/separate-log --data-type separate --method cp --rank 10

# Tucker decomposition
stdm-fit --data-path data/separate-log --data-type separate --method tucker --rank 10
```

### Dataset 3: Separate Species (Syn Normalized)
```bash
# Auto-detects separate species format
stdm-fit --data-path data/separate-syn --data-type separate --method cp --rank 10
```

## Python API Examples

### Load and Analyze Merged Data
```python
from stdm import DataLoader, TensorBuilder, CPDecomposition

# Load merged data (all species together)
loader = DataLoader("data/merge-vst")
data = loader.load_merged_data()
data = loader.normalize_data(data, method="zscore")

# Build tensor and fit model
builder = TensorBuilder()
tensor = builder.from_merged_data(data)
model = CPDecomposition(rank=15, random_state=42)
model.fit(tensor, n_iter_max=100, verbose=True)

# Get results
factors = model.get_factors()
error = model.compute_reconstruction_error(tensor)
print(f"Reconstruction error: {error:.6f}")
```

### Load and Analyze Separate Species Data
```python
# Load separate species data
loader = DataLoader("data/separate-log")
species_data = loader.load_separate_data()

# Build tensor with gene alignment
builder = TensorBuilder()
tensor = builder.from_separate_data(species_data, align_genes=True)

# Fit Tucker model (better for separate data)
from stdm import TuckerDecomposition
model = TuckerDecomposition(rank=(100, 20, 4, 3), random_state=42)
model.fit(tensor, n_iter_max=50, verbose=True)
```

### Visualize Results
```python
from stdm.visualization import plot_components, plot_temporal_patterns, plot_species_comparison

metadata = builder.get_metadata()
factors = model.get_factors()

# Plot all components
plot_components(factors, mode_names=['Genes', 'Individuals', 'Time', 'Species'])

# Plot temporal patterns
plot_temporal_patterns(factors[2], metadata['timepoints'])

# Compare species patterns
plot_species_comparison(factors[3], metadata['species'])
```

## Understanding Results

Each run creates a timestamped directory (e.g., `results/cp_20241012_143022/`) with:
- **ANALYSIS_README.md** - Detailed report explaining all outputs
- **Factor matrices** - CSV files with decomposition results
- **Visualizations** - PNG plots of components

## Data Formats

### Merged Format (data/merge-vst/)
Single CSV with all species:
```csv
group_id,ACR-139-TP1,ACR-139-TP2,POR-216-TP1,POC-201-TP3,...
OG_00001,5.97,6.21,4.96,8.19,...
```

### Separate Format (data/separate-log/, data/separate-syn/)
One CSV per species:
```csv
group_id,ACR.100.TP1,ACR.100.TP2,ACR.100.TP3,...
OG_00001,2.23,7.48,10.94,12.06,...
```

## Example Scripts

```bash
# Run complete analysis on merged data
python examples/fit_merged_data.py

# Analyze separate species data
python examples/fit_separate_data.py

# Compare different models and ranks
python examples/compare_models.py
```

## Advanced Features

- **Auto-timestamped outputs** - Never overwrite results
- **Automatic analysis reports** - Comprehensive documentation for each run
- **Custom data validation** - Check your data format before analysis
- **Multiple decomposition methods** - CP and Tucker with various constraints
- **Rich visualizations** - Components, temporal patterns, species comparisons

## Need Help?

- **Full documentation**: See the backed up files (README-old.md, etc.)
- **Interactive tutorial**: `examples/quickstart_notebook.ipynb`
- **API reference**: Check the source code docstrings
- **Troubleshooting**: See the detailed guides in backed up documentation

## Citation

```bibtex
@software{stdm2024,
  title = {STDM: Sparse Tensor Decomposition Models for Gene Expression Data},
  author = {Your Name},
  year = {2024},
  url = {https://github.com/yourusername/workflow-stdm-05}
}
```
