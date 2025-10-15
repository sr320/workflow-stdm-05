# STDM: Sparse Tensor Decomposition for Gene Expression Data

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Analyze multi-species, multi-timepoint gene expression data and get:**
1. 📈 **Line plots** of component expression over time
2. 📝 **Gene lists** for each component

Perfect for temporal RNA-seq experiments across multiple species or conditions.

---

## Quick Install & Run

```bash
# 1. Clone and install
git clone https://github.com/yourusername/workflow-stdm-05.git
cd workflow-stdm-05
./setup.sh

# 2. Run analysis (uses example data)
python examples/quickstart.py

# 3. Check results/
#    📈 component_expression_over_time.png
#    📝 component_genes_summary.txt
```

**That's it!** Your results are in the timestamped `results/quickstart_*` directory.

---

## What You Get

### 1. Component Expression Over Time 📈

A clear line plot showing how each component changes across your time points:

![Component Expression](https://via.placeholder.com/800x400/4A90E2/FFFFFF?text=Component+Expression+Over+Time)

Each line is a distinct expression pattern discovered in your data. Use this to:
- Identify temporal trends
- Find co-regulated gene modules
- Compare expression dynamics

### 2. Top Genes Per Component 📝

Three formats for easy analysis:

**`component_genes_summary.txt`** - Quick overview:
```
================================================================================
COMPONENT 1
================================================================================
Top 10 genes (total: 50):

   1. OG_00686              (loading:   5.0069)
   2. OG_00683              (loading:   6.7286)
   3. OG_00685              (loading:   7.3699)
   ...
```

**`component_genes.csv`** - Import into Excel/R:
```csv
Component,Rank,Gene,Loading,AbsLoading
1,1,OG_00686,5.006909,5.006909
1,2,OG_00683,6.728631,6.728631
...
```

**`component_X_genes.txt`** - Full list per component

---

## Using Your Own Data

### Data Format

Your data should be a CSV with genes as rows and samples as columns:

```csv
group_id,ACR-139-TP1,ACR-139-TP2,POR-216-TP1,POC-201-TP3,...
OG_00001,5.97,6.21,4.96,8.19,...
OG_00002,3.45,7.82,5.11,6.23,...
```

**Column naming convention:** `{SPECIES}-{INDIVIDUAL}-TP{TIMEPOINT}`
- Examples: `ACR-139-TP1`, `POR-216-TP2`, `POC-45-TP4`

### Run on Your Data

**Option 1: Command Line (Recommended)**
```bash
stdm-fit --data-path your_data.csv --method cp --rank 10
```

**Option 2: Python Script**

Edit `examples/quickstart.py`:
```python
DATA_PATH = Path("path/to/your_data.csv")
N_COMPONENTS = 10  # Adjust as needed
```

Then run:
```bash
python examples/quickstart.py
```

**Option 3: Python API**
```python
from stdm import DataLoader, TensorBuilder, CPDecomposition
from stdm import plot_temporal_expression, extract_component_genes

# Load your data
loader = DataLoader("path/to/data/directory")
data = loader.load_merged_data("your_file.csv")

# Build tensor
builder = TensorBuilder()
tensor = builder.from_merged_data(data)
metadata = builder.get_metadata()

# Fit model
model = CPDecomposition(rank=10, random_state=42)
model.fit(tensor, n_iter_max=100)

# Get results
factors = model.get_factors()

# Generate outputs
plot_temporal_expression(
    factors[2],  # Time factor
    metadata['timepoints'],
    save_path="temporal_plot.png"
)

extract_component_genes(
    factors[0],  # Gene factor
    metadata['genes'],
    output_dir="results",
    top_n=50
)
```

---

## Understanding Your Results

### Components = Expression Patterns

Each component represents a group of co-expressed genes with a specific temporal pattern:

| Component | Temporal Pattern | Interpretation Example |
|-----------|-----------------|------------------------|
| Component 1 | Increases steadily | Developmental program |
| Component 2 | Peaks at TP2 | Stress response |
| Component 3 | Decreases over time | Early development genes |
| Component 4 | Flat across time | Housekeeping genes |

### Gene Loadings

The **loading** value tells you how strongly a gene contributes to that component:

- **High absolute loading (>5)**: Gene is a key marker for this pattern
- **Medium loading (2-5)**: Gene participates in this pattern  
- **Low loading (<2)**: Gene is not strongly associated

Use the top genes for:
- GO enrichment analysis
- Pathway analysis
- Experimental validation
- Literature mining

---

## Choosing the Number of Components

Start with **rank=10** and adjust based on your data:

```bash
# Try different ranks
stdm-fit --data-path data.csv --rank 5
stdm-fit --data-path data.csv --rank 10
stdm-fit --data-path data.csv --rank 15
```

**Guidelines:**
- **Small datasets (<1k genes)**: 5-10 components
- **Medium datasets (1-5k genes)**: 10-15 components
- **Large datasets (>5k genes)**: 15-25 components

Look at the **reconstruction error** in the output:
- Lower is better
- < 0.3 = excellent fit
- 0.3-0.5 = good fit
- \> 0.5 = consider increasing rank

---

## Command-Line Options

```bash
stdm-fit [OPTIONS]

Required:
  --data-path PATH         Path to your CSV file or directory

Optional:
  --method {cp,tucker}     Decomposition method (default: cp)
  --rank INT               Number of components (default: 10)
  --normalize              Apply z-score normalization
  --non-negative           Enforce non-negative factorization
  --n-iter INT             Max iterations (default: 100)
  --output-dir PATH        Output directory (default: results)
  --no-report              Skip detailed analysis report
  --validate-only          Check data format without running analysis

Examples:
  # Basic analysis
  stdm-fit --data-path data/my_data.csv --rank 10
  
  # With normalization
  stdm-fit --data-path data/my_data.csv --rank 15 --normalize
  
  # Tucker decomposition (for complex data)
  stdm-fit --data-path data/my_data.csv --method tucker --rank 20
```

---

## Example Data

Three datasets are included in `data/`:

1. **`data/merge-vst/`** - All species merged (easiest to start with)
   ```bash
   stdm-fit --data-path data/merge-vst/vst_counts_matrix.csv --rank 10
   ```

2. **`data/separate-log/`** - Species separate, log-normalized
   ```bash
   stdm-fit --data-path data/separate-log --data-type separate --rank 10
   ```

3. **`data/separate-syn/`** - Species separate, syn-normalized
   ```bash
   stdm-fit --data-path data/separate-syn --data-type separate --rank 10
   ```

---

## Output Files Explained

Every analysis creates a timestamped directory (e.g., `results/cp_20251015_143022/`) with:

### 🎯 Key Outputs (Start Here!)
- **`component_expression_over_time.png`** - Temporal line plot
- **`component_genes_summary.txt`** - Quick gene overview
- **`component_genes.csv`** - Full gene table

### 📊 Detailed Results
- **`ANALYSIS_README.md`** - Comprehensive analysis report
- **`*_factor_genes.csv`** - Full gene factor matrix
- **`*_factor_timepoints.csv`** - Time factor matrix
- **`*_factor_individuals.csv`** - Individual factor matrix
- **`*_factor_species.csv`** - Species factor matrix
- **`*_metadata.json`** - Dimension labels and parameters

### 🔬 Advanced Visualizations
- **`*_components.png`** - All factors visualized
- Individual component files for detailed analysis

---

## Troubleshooting

### "No module named 'stdm'"
```bash
# Make sure you're in the right directory and installed
cd workflow-stdm-05
./setup.sh
source .venv/bin/activate
```

### "ValueError: Could not parse sample name"
Your column names must follow the pattern: `SPECIES-INDIVIDUAL-TPTIMEPOINT`
- ✅ Good: `ACR-139-TP1`, `POR-45-TP2`
- ❌ Bad: `sample1`, `ACR_139_T1`

### "Reconstruction error is high (>0.5)"
Try:
1. Increase number of components: `--rank 15` or `--rank 20`
2. Increase iterations: `--n-iter 200`
3. Use Tucker decomposition: `--method tucker`
4. Normalize your data: `--normalize`

### Need more help?
```bash
# Validate your data format
stdm-fit --data-path your_data.csv --validate-only

# Check examples
python examples/quickstart.py
python examples/fit_merged_data.py
```

---

## Advanced Usage

### Multiple Datasets

Analyze separate species files:
```bash
# Put files in a directory: data/my_species/
#   - species1_expression.csv
#   - species2_expression.csv
#   - species3_expression.csv

stdm-fit --data-path data/my_species --data-type separate --rank 10
```

### Batch Processing

```python
# Process multiple ranks
for rank in [5, 10, 15, 20]:
    os.system(f"stdm-fit --data-path data.csv --rank {rank}")
```

### Custom Analysis

See `examples/` for more:
- `fit_merged_data.py` - Full pipeline with custom options
- `fit_separate_data.py` - Multi-species analysis
- `compare_models.py` - Compare CP vs Tucker
- `custom_data_analysis.py` - Advanced customization

---

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

---

## References

- **CP Decomposition**: Harshman, R. A. (1970). Foundations of the PARAFAC procedure
- **Tucker Decomposition**: Tucker, L. R. (1966). Some mathematical notes on three-mode factor analysis
- **TensorLy**: Kossaifi et al. (2019). TensorLy: Tensor Learning in Python

---

## License

MIT License - see [LICENSE](LICENSE) file for details.

---

## Questions or Issues?

- 📖 Check the [detailed report](results/*/ANALYSIS_README.md) generated with each run
- 💡 See [examples/](examples/) for working code
- 🐛 Open an issue on GitHub with your error message and data format

**Happy analyzing! 🧬📊**
