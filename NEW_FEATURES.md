# New Features in STDM v0.1.0

## Summary of Enhancements

Three major features have been added to make STDM more user-friendly and production-ready:

1. **🎯 Auto-Timestamped Outputs** - Never overwrite results again
2. **📝 Automatic Analysis Reports** - Comprehensive README generated for every run
3. **✅ Custom Data Support** - Easy validation and analysis of your own data

---

## 1. Auto-Timestamped Outputs ⏰

### What It Does

Every analysis run now creates a unique timestamped directory, ensuring you never accidentally overwrite previous results.

### Example

```bash
# Run analysis multiple times
stdm-fit --data-path data/merge-vst/vst_counts_matrix.csv --method cp --rank 10
stdm-fit --data-path data/merge-vst/vst_counts_matrix.csv --method cp --rank 15
stdm-fit --data-path data/merge-vst/vst_counts_matrix.csv --method tucker --rank 10
```

**Creates:**
```
results/
├── cp_20241012_143022/      # First run (rank=10)
├── cp_20241012_144530/      # Second run (rank=15)
└── tucker_20241012_150145/  # Third run (Tucker)
```

### Benefits

- ✅ Compare results across different parameters
- ✅ Track analysis history
- ✅ No risk of overwriting important results
- ✅ Each run has unique identifier

### Disable If Needed

```bash
# Use --no-timestamp to overwrite results
stdm-fit --data-path data.csv --no-timestamp
```

---

## 2. Automatic Analysis Reports 📝

### What It Does

Every analysis generates a comprehensive `ANALYSIS_README.md` that explains:

1. **What each output file contains**
2. **How to interpret your results**
3. **Quality assessment of your model**
4. **Optimal parameter recommendations**
5. **Next steps and code examples**

### Generated Report Structure

**ANALYSIS_README.md** contains:

```markdown
# Tensor Decomposition Analysis Report

Generated: 2024-10-12 14:30:22
Analysis ID: 20241012_143022

## Quick Summary
- Method, rank, performance metrics

## Output Files Explained
- Detailed description of each CSV and JSON file
- What each column/row represents
- How to use the data

## How to Interpret Your Results
- Understanding factor values
- Reading components
- Key findings from your data

## Optimal Parameters & Recommendations
- Quality assessment
- Suggested improvements
- Parameter tuning guidelines

## Next Steps
- Code examples for further analysis
- How to extract marker genes
- Validation approaches
```

### Example Content

The report explains things like:

**For Gene Factor:**
```
cp_factor_genes.csv
- Size: Genes × Components
- Content: Gene loadings for each component
- Use: Identify which genes contribute to each pattern
- Interpretation:
  • High positive values = gene strongly associated
  • Values near zero = gene not involved
  • Can be used for marker gene identification
```

**Quality Metrics:**
```
Reconstruction Error: 0.256123
- Excellent fit (< 0.3)
- Captures most major patterns

Explained Variance: 0.8234 (82.34%)
- Good - captures most major patterns
```

**Recommendations:**
```
Your Data: Current rank (10) seems appropriate

Suggested Experiments:
1. Test multiple ranks: [5, 10, 15, 20]
2. Try Tucker decomposition for more flexibility
3. Apply more stringent gene filtering
```

### Benefits

- ✅ Self-documenting analyses
- ✅ No need to remember what each file contains
- ✅ Optimal parameters identified automatically
- ✅ Reproducible research

### Disable If Needed

```bash
stdm-fit --data-path data.csv --no-report
```

---

## 3. Custom Data Support ✅

### What It Does

New data validation and easy integration for your own gene expression data.

### Step 1: Validate Your Data

```bash
# Check if your data format is correct
stdm-fit --data-path /path/to/your/data.csv --validate-only
```

**Output:**
```
============================================================
DATA VALIDATION REPORT
============================================================
✅ File readable as CSV
📊 Data shape: 8523 genes × 96 samples
✅ No missing values
✅ All values are numeric
✅ Sample names follow expected format: SPECIES-ID-TPX
📈 Value range: [0.15, 18.42]
ℹ️  Large values detected. Consider log-transformation...

📋 Validation Summary:
   Genes: 8523
   Samples: 96
   Data type: Numeric
   Missing values: 0
============================================================
```

### Step 2: Run Analysis

```bash
stdm-fit \
    --data-path /path/to/your/data.csv \
    --method cp \
    --rank 15 \
    --normalize
```

### Your Data Format

**Merged CSV:**
```csv
group_id,SPP1-01-TP1,SPP1-01-TP2,SPP2-05-TP1,...
GENE001,5.2,6.1,4.8,...
GENE002,7.3,8.2,6.9,...
```

**Key Requirements:**
- First column: Gene IDs
- Other columns: Samples
- Sample format: `SPECIES-INDIVIDUAL-TIMEPOINT`
- All numeric values
- No missing data (or will be filled with 0)

**Separate Files:**
```
my_data/
├── species1_normalized_expression.csv
├── species2_normalized_expression.csv
└── species3_normalized_expression.csv
```

### Data Validator Features

```python
from stdm import DataValidator

validator = DataValidator(strict=False)

# Validate your data
is_valid, messages = validator.validate_merged_data("your_data.csv")

# Print detailed report
validator.print_report()

# Check specific issues
for msg in messages:
    print(msg)
```

**Checks:**
- ✅ File format and readability
- ✅ Data shape and size
- ✅ Missing values
- ✅ Data types (numeric)
- ✅ Sample name formats
- ✅ Value ranges
- ✅ Zero variance genes
- ✅ Gene alignment (for separate files)

### Benefits

- ✅ Catch data issues before analysis
- ✅ Clear error messages
- ✅ Format validation
- ✅ Works with any gene expression data

---

## Complete Example Workflow

### Using CLI

```bash
# 1. Validate your data
stdm-fit --data-path your_data.csv --validate-only

# 2. Run analysis (auto-timestamped, with report)
stdm-fit \
    --data-path your_data.csv \
    --method cp \
    --rank 15 \
    --normalize \
    --output-dir results/my_project

# 3. Check results
ls results/my_project/cp_20241012_*/
cat results/my_project/cp_20241012_*/ANALYSIS_README.md
```

### Using Python API

```python
from stdm import (
    DataValidator,
    DataLoader,
    TensorBuilder,
    CPDecomposition,
    create_timestamped_output_dir,
    AnalysisReport
)
from stdm.utils import save_results, compute_explained_variance
import time

# 1. Validate data
validator = DataValidator()
is_valid, msgs = validator.validate_merged_data("your_data.csv")
validator.print_report()

if not is_valid:
    print("Fix data issues and try again")
    exit()

# 2. Create timestamped output directory
output_dir = create_timestamped_output_dir("results", prefix="my_analysis")
print(f"Results will be saved to: {output_dir}")

# 3. Load and process
start_time = time.time()
loader = DataLoader(".")
data = loader.load_merged_data("your_data.csv")
data = loader.normalize_data(data, method="zscore")

# 4. Build tensor
builder = TensorBuilder()
tensor = builder.from_merged_data(data)
metadata = builder.get_metadata()

# 5. Fit model
model = CPDecomposition(rank=15, random_state=42)
model.fit(tensor, n_iter_max=100, verbose=True)

# 6. Compute metrics
runtime = time.time() - start_time
factors = model.get_factors()
error = model.compute_reconstruction_error(tensor)
reconstructed = model.reconstruct()
explained_var = compute_explained_variance(tensor, reconstructed)

# 7. Save results
save_results(factors, metadata, output_dir, prefix="cp")

# 8. Generate automatic report
reporter = AnalysisReport(output_dir)
report_path = reporter.generate_run_report(
    method="CP",
    tensor_shape=tensor.shape,
    rank=15,
    reconstruction_error=error,
    explained_variance=explained_var,
    runtime_seconds=runtime,
    metadata=metadata,
    factors=factors
)

print(f"\n✅ Analysis complete!")
print(f"📁 Results: {output_dir}")
print(f"📄 Report: {report_path}")
```

---

## New Example Script

**`examples/custom_data_analysis.py`**

Complete demonstration of all new features. Run it to see:
- Data validation in action
- Timestamped outputs
- Automatic report generation
- How to use your own data

```bash
python examples/custom_data_analysis.py
```

---

## Updated CLI Options

```bash
stdm-fit --help
```

**New Options:**
- `--data-path` - Path to data file or directory (replaces `--data-dir`)
- `--data-type auto` - Auto-detect merged vs separate (new default)
- `--no-timestamp` - Disable timestamped directories
- `--no-report` - Disable automatic report generation
- `--validate-only` - Only validate data format
- `--n-iter` - Control max iterations

---

## Benefits Summary

| Feature | Before | After |
|---------|--------|-------|
| **Output Management** | Manual directory creation, risk of overwrite | Auto-timestamped, safe comparison |
| **Result Documentation** | Manual README writing | Automatic comprehensive reports |
| **Data Validation** | Trial and error | Validate before running |
| **Custom Data** | Complex setup | Simple 3-step process |
| **Parameter Guidance** | External research needed | Recommendations in report |

---

## Migration Guide

**Old workflow:**
```bash
stdm-fit --data-dir data/merge-vst --data-type merged ...
```

**New workflow (still works):**
```bash
stdm-fit --data-path data/merge-vst/vst_counts_matrix.csv ...
```

Or use automatic detection:
```bash
stdm-fit --data-path data/merge-vst/vst_counts_matrix.csv  # Auto-detects merged
stdm-fit --data-path data/separate-log/                     # Auto-detects separate
```

**All old commands still work!** New features are additions, not breaking changes.

---

## Documentation

- **README.md** - Updated with new features
- **GETTING_STARTED.md** - Includes new workflows
- **examples/custom_data_analysis.py** - Complete example
- **Auto-generated ANALYSIS_README.md** - Per-run documentation

---

## Questions?

**Q: Will old results still work?**
A: Yes! No breaking changes, just new features.

**Q: Can I disable timestamping?**
A: Yes, use `--no-timestamp` flag.

**Q: Can I disable reports?**
A: Yes, use `--no-report` flag.

**Q: What if my data format is different?**
A: Run with `--validate-only` to see specific issues.

**Q: How do I use my own data?**
A: See `examples/custom_data_analysis.py` or use `stdm-fit --data-path your_data.csv --validate-only`

---

**Version**: 0.1.0
**Date**: October 12, 2024
**Status**: ✅ Production Ready

All new features are fully tested and documented!

