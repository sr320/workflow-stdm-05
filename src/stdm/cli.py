"""
Command-line interface for STDM package.
"""

import argparse
import sys
from pathlib import Path
from stdm import DataLoader, TensorBuilder, CPDecomposition, TuckerDecomposition
from stdm.utils import save_results
from stdm.visualization import plot_components, plot_reconstruction_error


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Sparse Tensor Decomposition Models for Gene Expression Data"
    )
    
    parser.add_argument(
        "--data-dir",
        type=str,
        required=True,
        help="Path to data directory"
    )
    
    parser.add_argument(
        "--data-type",
        type=str,
        choices=["merged", "separate"],
        required=True,
        help="Type of data: merged or separate"
    )
    
    parser.add_argument(
        "--method",
        type=str,
        choices=["cp", "tucker"],
        default="cp",
        help="Decomposition method: cp or tucker"
    )
    
    parser.add_argument(
        "--rank",
        type=int,
        default=10,
        help="Rank for CP decomposition or first mode rank for Tucker"
    )
    
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results",
        help="Output directory for results"
    )
    
    parser.add_argument(
        "--non-negative",
        action="store_true",
        help="Use non-negative decomposition"
    )
    
    parser.add_argument(
        "--normalize",
        action="store_true",
        help="Normalize data before decomposition"
    )
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("STDM: Sparse Tensor Decomposition Models")
    print("=" * 60)
    
    # Load data
    print(f"\nLoading data from {args.data_dir}...")
    loader = DataLoader(args.data_dir)
    
    if args.data_type == "merged":
        data = loader.load_merged_data()
    else:
        data_dict = loader.load_separate_data()
    
    # Normalize if requested
    if args.normalize:
        print("\nNormalizing data...")
        if args.data_type == "merged":
            data = loader.normalize_data(data, method="zscore")
        else:
            data_dict = {
                species: loader.normalize_data(df, method="zscore")
                for species, df in data_dict.items()
            }
    
    # Build tensor
    print("\nBuilding tensor...")
    builder = TensorBuilder()
    
    if args.data_type == "merged":
        tensor = builder.from_merged_data(data)
    else:
        tensor = builder.from_separate_data(data_dict)
    
    metadata = builder.get_metadata()
    
    # Fit decomposition
    print(f"\nFitting {args.method.upper()} decomposition...")
    
    if args.method == "cp":
        model = CPDecomposition(
            rank=args.rank,
            non_negative=args.non_negative
        )
    else:
        # Use proportional ranks for Tucker
        shape = tensor.shape
        tucker_rank = (
            min(args.rank, shape[0]),
            min(args.rank // 2, shape[1]),
            min(4, shape[2]),
            min(3, shape[3])
        )
        model = TuckerDecomposition(
            rank=tucker_rank,
            non_negative=args.non_negative
        )
    
    model.fit(tensor)
    
    # Save results
    print(f"\nSaving results to {args.output_dir}...")
    factors = model.get_factors()
    save_results(factors, metadata, args.output_dir, prefix=args.method)
    
    print("\n" + "=" * 60)
    print("Decomposition completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()

