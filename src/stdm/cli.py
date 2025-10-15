"""
Command-line interface for STDM package with auto-reporting and timestamping.
"""

import argparse
import sys
import time
from pathlib import Path
from stdm import (
    DataLoader, TensorBuilder, CPDecomposition, TuckerDecomposition,
    AnalysisReport, create_timestamped_output_dir, DataValidator
)
from stdm.utils import save_results, compute_explained_variance, extract_component_genes
from stdm.visualization import plot_components, plot_temporal_expression


def main():
    """Main CLI entry point with auto-reporting."""
    parser = argparse.ArgumentParser(
        description="Sparse Tensor Decomposition Models for Gene Expression Data",
        epilog="Generates timestamped outputs and comprehensive analysis reports."
    )
    
    parser.add_argument(
        "--data-path",
        type=str,
        required=True,
        help="Path to data file (merged) or directory (separate)"
    )
    
    parser.add_argument(
        "--data-type",
        type=str,
        choices=["merged", "separate", "auto"],
        default="auto",
        help="Type of data: merged, separate, or auto-detect (default: auto)"
    )
    
    parser.add_argument(
        "--method",
        type=str,
        choices=["cp", "tucker"],
        default="cp",
        help="Decomposition method: cp or tucker (default: cp)"
    )
    
    parser.add_argument(
        "--rank",
        type=int,
        default=10,
        help="Rank for CP or first mode rank for Tucker (default: 10)"
    )
    
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results",
        help="Base output directory (default: results)"
    )
    
    parser.add_argument(
        "--no-timestamp",
        action="store_true",
        help="Disable timestamped subdirectories"
    )
    
    parser.add_argument(
        "--no-report",
        action="store_true",
        help="Disable automatic report generation"
    )
    
    parser.add_argument(
        "--non-negative",
        action="store_true",
        help="Use non-negative decomposition"
    )
    
    parser.add_argument(
        "--normalize",
        action="store_true",
        help="Normalize data before decomposition (z-score)"
    )
    
    parser.add_argument(
        "--n-iter",
        type=int,
        default=100,
        help="Maximum iterations for decomposition (default: 100)"
    )
    
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Only validate data without running analysis"
    )
    
    args = parser.parse_args()
    
    print("=" * 70)
    print("STDM: Sparse Tensor Decomposition Models")
    print("=" * 70)
    
    # Validate data
    print(f"\n📋 Validating data from {args.data_path}...")
    validator = DataValidator(strict=False)
    
    data_path = Path(args.data_path)
    if data_path.is_file():
        is_valid, messages = validator.validate_merged_data(data_path)
    elif data_path.is_dir():
        is_valid, messages = validator.validate_separate_data(data_path)
    else:
        print(f"❌ Path not found: {args.data_path}")
        sys.exit(1)
    
    validator.print_report()
    
    if not is_valid:
        print("❌ Data validation failed. Please fix errors and try again.")
        sys.exit(1)
    
    if args.validate_only:
        print("✅ Validation complete. Use without --validate-only to run analysis.")
        sys.exit(0)
    
    # Create output directory
    if args.no_timestamp:
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        output_dir = create_timestamped_output_dir(args.output_dir, prefix=args.method)
    
    print(f"\n📁 Output directory: {output_dir}")
    
    # Start timing
    start_time = time.time()
    
    # Load data
    print(f"\n📂 Loading data...")
    loader = DataLoader(args.data_path if data_path.is_dir() else data_path.parent)
    
    if data_path.is_file():
        data = loader.load_merged_data(data_path.name)
        data_type = "merged"
    else:
        data_dict = loader.load_separate_data()
        data_type = "separate"
    
    # Normalize if requested
    if args.normalize:
        print("\n🔄 Normalizing data (z-score)...")
        if data_type == "merged":
            data = loader.normalize_data(data, method="zscore")
        else:
            data_dict = {
                species: loader.normalize_data(df, method="zscore")
                for species, df in data_dict.items()
            }
    
    # Build tensor
    print("\n🔨 Building tensor...")
    builder = TensorBuilder()
    
    if data_type == "merged":
        tensor = builder.from_merged_data(data)
    else:
        tensor = builder.from_separate_data(data_dict)
    
    metadata = builder.get_metadata()
    
    # Fit decomposition
    print(f"\n🚀 Fitting {args.method.upper()} decomposition (rank={args.rank})...")
    
    if args.method == "cp":
        model = CPDecomposition(
            rank=args.rank,
            non_negative=args.non_negative,
            random_state=42
        )
    else:
        # Use proportional ranks for Tucker
        shape = tensor.shape
        tucker_rank = (
            min(args.rank * 2, shape[0]),
            min(args.rank, shape[1]),
            min(4, shape[2]),
            min(3, shape[3])
        )
        model = TuckerDecomposition(
            rank=tucker_rank,
            non_negative=args.non_negative,
            random_state=42
        )
        print(f"   Tucker ranks: {tucker_rank}")
    
    model.fit(tensor, n_iter_max=args.n_iter, verbose=True)
    
    # Compute metrics
    runtime = time.time() - start_time
    factors = model.get_factors()
    error = model.compute_reconstruction_error(tensor)
    reconstructed = model.reconstruct()
    explained_var = compute_explained_variance(tensor, reconstructed)
    
    print(f"\n📊 Results:")
    print(f"   Reconstruction error: {error:.6f}")
    print(f"   Explained variance: {explained_var:.4f} ({explained_var*100:.2f}%)")
    print(f"   Runtime: {runtime:.1f}s")
    
    # Save results
    print(f"\n💾 Saving results...")
    save_results(factors, metadata, output_dir, prefix=args.method)
    
    # Generate visualizations
    print("\n📈 Generating visualizations...")
    try:
        mode_names = ['Genes', 'Individuals', 'Time Points', 'Species']
        plot_components(
            factors,
            mode_names=mode_names,
            n_components=min(5, args.rank),
            save_path=output_dir / f"{args.method}_components.png"
        )
        
        # Generate temporal expression plot (key output)
        if len(factors) >= 3 and len(metadata.get('timepoints', [])) > 0:
            plot_temporal_expression(
                factors[2],  # Time factor
                metadata['timepoints'],
                n_components=args.rank,
                save_path=output_dir / "component_expression_over_time.png",
                title="Component Expression Across Time Points"
            )
    except Exception as e:
        print(f"   ⚠️  Could not generate plots: {e}")
    
    # Extract gene lists (key output)
    print("\n📝 Extracting top genes per component...")
    try:
        if len(factors) >= 1 and len(metadata.get('genes', [])) > 0:
            extract_component_genes(
                factors[0],  # Gene factor
                metadata['genes'],
                output_dir,
                top_n=50
            )
    except Exception as e:
        print(f"   ⚠️  Could not extract genes: {e}")
    
    # Generate report
    if not args.no_report:
        print("\n📝 Generating analysis report...")
        reporter = AnalysisReport(output_dir)
        rank = args.rank if args.method == "cp" else tucker_rank
        report_path = reporter.generate_run_report(
            method=args.method.upper(),
            tensor_shape=tensor.shape,
            rank=rank,
            reconstruction_error=error,
            explained_variance=explained_var,
            runtime_seconds=runtime,
            metadata=metadata,
            factors=factors
        )
        print(f"   📄 Report saved: {report_path}")
    
    print("\n" + "=" * 70)
    print("✅ Analysis completed successfully!")
    print(f"📁 All results saved to: {output_dir}")
    print("=" * 70)
    print("\n🎯 KEY OUTPUTS:")
    print(f"  📈 Temporal plot: component_expression_over_time.png")
    print(f"  📝 Gene lists: component_genes_summary.txt & component_genes.csv")
    print("\n📚 Additional outputs:")
    print(f"  1. Analysis report: ANALYSIS_README.md")
    print(f"  2. Factor matrices: *_factor_*.csv")
    print(f"  3. Visualizations: *.png")
    print("=" * 70)


if __name__ == "__main__":
    main()

