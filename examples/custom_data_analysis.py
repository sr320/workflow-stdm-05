"""
Example: Analyze your own gene expression data with auto-reporting.

This script demonstrates how to:
1. Validate your own data files
2. Run analysis with auto-timestamped outputs
3. Generate comprehensive analysis reports
4. Interpret results using the automated README
"""

import sys
from pathlib import Path
import time

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from stdm import (
    DataLoader,
    TensorBuilder,
    CPDecomposition,
    DataValidator,
    AnalysisReport,
    create_timestamped_output_dir
)
from stdm.utils import save_results, compute_explained_variance


def analyze_custom_data(
    data_path: str,
    output_base_dir: str = "results",
    rank: int = 10,
    normalize: bool = True
):
    """
    Analyze custom gene expression data with full reporting.
    
    Parameters
    ----------
    data_path : str
        Path to your data file (merged CSV) or directory (separate files)
    output_base_dir : str
        Base directory for results
    rank : int
        Number of components for decomposition
    normalize : bool
        Whether to normalize data
    """
    
    print("=" * 70)
    print("Custom Data Analysis with STDM")
    print("=" * 70)
    
    # Step 1: Validate your data
    print(f"\n📋 Step 1: Validating your data...")
    print(f"   Data path: {data_path}")
    
    validator = DataValidator(strict=False)
    data_path_obj = Path(data_path)
    
    if not data_path_obj.exists():
        print(f"\n❌ Error: Path not found: {data_path}")
        print("\n💡 Tips:")
        print("   - Check the file/directory path is correct")
        print("   - Use absolute path if relative path doesn't work")
        print("   - Ensure file has read permissions")
        return
    
    # Validate based on type
    if data_path_obj.is_file():
        print("   Detected: Merged data file")
        is_valid, messages = validator.validate_merged_data(data_path)
    elif data_path_obj.is_dir():
        print("   Detected: Directory with separate files")
        is_valid, messages = validator.validate_separate_data(data_path)
    else:
        print(f"\n❌ Error: {data_path} is neither a file nor directory")
        return
    
    # Print validation report
    validator.print_report()
    
    if not is_valid:
        print("\n❌ Data validation failed!")
        print("\n💡 Common issues:")
        print("   - Non-numeric values in expression data")
        print("   - Missing or malformed sample names")
        print("   - Inconsistent gene IDs across files")
        print("   - Very few genes or samples")
        print("\n   Fix these issues and try again.")
        return
    
    print("\n✅ Data validation passed!")
    
    # Step 2: Create timestamped output directory
    print(f"\n📁 Step 2: Creating output directory...")
    output_dir = create_timestamped_output_dir(output_base_dir, prefix="custom")
    print(f"   Output: {output_dir}")
    print(f"   This run ID: {output_dir.name}")
    
    # Step 3: Load data
    start_time = time.time()
    print(f"\n📂 Step 3: Loading data...")
    
    if data_path_obj.is_file():
        loader = DataLoader(data_path_obj.parent)
        data = loader.load_merged_data(data_path_obj.name)
        data_type = "merged"
        print(f"   Loaded: {data.shape[0]} genes × {data.shape[1]} samples")
    else:
        loader = DataLoader(data_path)
        data_dict = loader.load_separate_data()
        data_type = "separate"
        print(f"   Loaded {len(data_dict)} species files")
    
    # Step 4: Preprocess
    if normalize:
        print(f"\n🔄 Step 4: Normalizing data...")
        if data_type == "merged":
            data = loader.normalize_data(data, method="zscore")
            print("   Applied z-score normalization")
        else:
            for species in data_dict:
                data_dict[species] = loader.normalize_data(
                    data_dict[species], 
                    method="zscore"
                )
            print(f"   Applied z-score normalization to all {len(data_dict)} species")
    
    # Step 5: Build tensor
    print(f"\n🔨 Step 5: Building tensor...")
    builder = TensorBuilder()
    
    if data_type == "merged":
        tensor = builder.from_merged_data(data)
    else:
        tensor = builder.from_separate_data(data_dict)
    
    metadata = builder.get_metadata()
    print(f"   Tensor shape: {tensor.shape}")
    print(f"   Dimensions: Genes={tensor.shape[0]}, Individuals={tensor.shape[1]}, ")
    print(f"              Time={tensor.shape[2]}, Species={tensor.shape[3]}")
    
    # Step 6: Fit model
    print(f"\n🚀 Step 6: Fitting CP decomposition (rank={rank})...")
    model = CPDecomposition(rank=rank, random_state=42)
    model.fit(tensor, n_iter_max=100, verbose=True)
    
    # Step 7: Compute metrics
    print(f"\n📊 Step 7: Computing performance metrics...")
    factors = model.get_factors()
    error = model.compute_reconstruction_error(tensor)
    reconstructed = model.reconstruct()
    explained_var = compute_explained_variance(tensor, reconstructed)
    runtime = time.time() - start_time
    
    print(f"\n   Results:")
    print(f"   ├─ Reconstruction error: {error:.6f}")
    print(f"   ├─ Explained variance: {explained_var:.4f} ({explained_var*100:.2f}%)")
    print(f"   └─ Runtime: {runtime:.1f}s")
    
    # Step 8: Save results
    print(f"\n💾 Step 8: Saving results...")
    save_results(factors, metadata, output_dir, prefix="cp")
    print(f"   ✅ Factor matrices saved")
    print(f"   ✅ Metadata saved")
    
    # Step 9: Generate comprehensive report
    print(f"\n📝 Step 9: Generating analysis report...")
    reporter = AnalysisReport(output_dir)
    report_path = reporter.generate_run_report(
        method="CP",
        tensor_shape=tensor.shape,
        rank=rank,
        reconstruction_error=error,
        explained_variance=explained_var,
        runtime_seconds=runtime,
        metadata=metadata,
        factors=factors
    )
    
    print(f"   ✅ Report generated: {report_path.name}")
    
    # Step 10: Summary
    print("\n" + "=" * 70)
    print("✅ Analysis Complete!")
    print("=" * 70)
    print(f"\n📁 All results saved to: {output_dir}")
    print(f"\n📄 Output files:")
    print(f"   ├─ ANALYSIS_README.md     - Comprehensive analysis report")
    print(f"   ├─ cp_factor_genes.csv    - Gene loadings ({tensor.shape[0]} × {rank})")
    print(f"   ├─ cp_factor_individuals.csv - Individual patterns")
    print(f"   ├─ cp_factor_timepoints.csv - Temporal dynamics")
    print(f"   ├─ cp_factor_species.csv - Species-specific patterns")
    print(f"   └─ cp_metadata.json      - Dimension labels")
    
    print(f"\n📖 Next Steps:")
    print(f"   1. Read: {output_dir}/ANALYSIS_README.md")
    print(f"      This file explains:")
    print(f"      • How to interpret each output file")
    print(f"      • Optimal parameter recommendations")
    print(f"      • Quality assessment of your results")
    print(f"      • Troubleshooting tips")
    
    print(f"\n   2. Explore factor matrices:")
    print(f"      import pandas as pd")
    print(f"      genes = pd.read_csv('{output_dir}/cp_factor_genes.csv', index_col=0)")
    print(f"      print(genes.head())")
    
    print(f"\n   3. Extract marker genes for component 0:")
    print(f"      top_genes = genes['component_0'].abs().nlargest(50)")
    
    print(f"\n   4. Compare results across runs:")
    print(f"      Each run is timestamped, so you can:")
    print(f"      • Try different ranks (5, 10, 15, 20)")
    print(f"      • Test normalization methods")
    print(f"      • Compare Tucker vs CP")
    
    print("\n" + "=" * 70)


def main():
    """Main function."""
    
    # Example 1: Use provided sample data
    print("=" * 70)
    print("Example 1: Analyzing provided sample data")
    print("=" * 70)
    
    sample_data = Path(__file__).parent.parent / "data" / "merge-vst" / "vst_counts_matrix.csv"
    
    if sample_data.exists():
        analyze_custom_data(
            data_path=str(sample_data),
            output_base_dir="results/custom_example",
            rank=10,
            normalize=True
        )
    else:
        print(f"\n⚠️  Sample data not found at {sample_data}")
        print("     Using your own data instead...")
    
    # Example 2: Guide for using your own data
    print("\n\n" + "=" * 70)
    print("Using Your Own Data")
    print("=" * 70)
    print("""
To analyze your own data, use one of these approaches:

**Option 1: Command Line**
```bash
stdm-fit --data-path /path/to/your/data.csv \\
         --method cp --rank 15 --normalize \\
         --output-dir results/my_analysis
```

**Option 2: Python API**
```python
from stdm import DataValidator, create_timestamped_output_dir

# First, validate your data
validator = DataValidator()
is_valid, messages = validator.validate_merged_data("your_data.csv")
validator.print_report()

# Then run analysis
analyze_custom_data(
    data_path="your_data.csv",
    output_base_dir="results/my_project",
    rank=15
)
```

**Your Data Format Should Be:**

Merged CSV:
- Rows: Genes (with gene IDs as row names)
- Columns: Samples (format: SPECIES-INDIVIDUAL-TIMEPOINT)
- Values: Normalized expression (numeric)

Example:
```
group_id,ACR-001-TP1,ACR-001-TP2,POR-002-TP1,...
GENE001,5.2,6.1,4.8,...
GENE002,7.3,8.2,6.9,...
```

Separate Directory:
- Multiple CSV files, one per species
- File naming: species_normalized_expression.csv
- Same format as merged, but one species per file

**Need Help?**
- Run with --validate-only to check your data format
- Check GETTING_STARTED.md for detailed examples
- See examples/custom_data_analysis.py (this file!)
""")


if __name__ == "__main__":
    main()

