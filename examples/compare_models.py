"""
Example: Compare different tensor decomposition models and ranks.

This script demonstrates how to:
1. Test multiple decomposition ranks
2. Compare CP vs Tucker decomposition
3. Evaluate model performance using cross-validation
4. Select optimal model parameters
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import numpy as np
import matplotlib.pyplot as plt
from stdm import DataLoader, TensorBuilder, CPDecomposition, TuckerDecomposition
from stdm.visualization import plot_reconstruction_error
from stdm.utils import compute_explained_variance


def compare_ranks(tensor, ranks, metadata, output_dir):
    """
    Compare different ranks for CP decomposition.
    
    Parameters
    ----------
    tensor : np.ndarray
        Input tensor
    ranks : list of int
        List of ranks to test
    metadata : dict
        Tensor metadata
    output_dir : Path
        Output directory
    """
    print("\n" + "=" * 70)
    print("Comparing CP Decomposition Ranks")
    print("=" * 70)
    
    errors = {}
    explained_vars = {}
    
    for rank in ranks:
        print(f"\nTesting rank {rank}...")
        model = CPDecomposition(rank=rank, random_state=42)
        model.fit(tensor, n_iter_max=100, verbose=False)
        
        error = model.compute_reconstruction_error(tensor)
        reconstructed = model.reconstruct()
        explained_var = compute_explained_variance(tensor, reconstructed)
        
        errors[rank] = error
        explained_vars[rank] = explained_var
        
        print(f"  Reconstruction error: {error:.6f}")
        print(f"  Explained variance: {explained_var:.4f}")
    
    # Plot results
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot reconstruction error
    ranks_sorted = sorted(errors.keys())
    ax1.plot(ranks_sorted, [errors[r] for r in ranks_sorted], 'o-', linewidth=2, markersize=8)
    ax1.set_xlabel('Rank', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Reconstruction Error', fontsize=12, fontweight='bold')
    ax1.set_title('CP Decomposition: Reconstruction Error', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # Plot explained variance
    ax2.plot(ranks_sorted, [explained_vars[r] for r in ranks_sorted], 'o-', 
             linewidth=2, markersize=8, color='green')
    ax2.set_xlabel('Rank', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Explained Variance', fontsize=12, fontweight='bold')
    ax2.set_title('CP Decomposition: Explained Variance', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    save_path = output_dir / "rank_comparison.png"
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"\nSaved rank comparison plot to {save_path}")
    plt.close()
    
    # Find optimal rank (elbow method)
    # Simple heuristic: find where improvement decreases significantly
    errors_list = [errors[r] for r in ranks_sorted]
    improvements = [errors_list[i] - errors_list[i+1] for i in range(len(errors_list)-1)]
    
    if improvements:
        # Find where improvement drops below 50% of max improvement
        max_improvement = max(improvements)
        threshold = 0.5 * max_improvement
        
        for i, imp in enumerate(improvements):
            if imp < threshold:
                optimal_rank = ranks_sorted[i]
                break
        else:
            optimal_rank = ranks_sorted[-1]
        
        print(f"\nSuggested optimal rank: {optimal_rank}")
    
    return errors, explained_vars


def compare_methods(tensor, rank, metadata, output_dir):
    """
    Compare CP and Tucker decomposition methods.
    
    Parameters
    ----------
    tensor : np.ndarray
        Input tensor
    rank : int
        Rank for CP, base rank for Tucker
    metadata : dict
        Tensor metadata
    output_dir : Path
        Output directory
    """
    print("\n" + "=" * 70)
    print("Comparing CP vs Tucker Decomposition")
    print("=" * 70)
    
    results = {}
    
    # CP decomposition
    print(f"\n1. Fitting CP decomposition (rank={rank})...")
    cp_model = CPDecomposition(rank=rank, random_state=42)
    cp_model.fit(tensor, n_iter_max=100, verbose=False)
    
    cp_error = cp_model.compute_reconstruction_error(tensor)
    cp_reconstructed = cp_model.reconstruct()
    cp_explained_var = compute_explained_variance(tensor, cp_reconstructed)
    
    results['CP'] = {
        'error': cp_error,
        'explained_var': cp_explained_var,
        'n_parameters': sum(f.size for f in cp_model.get_factors())
    }
    
    print(f"  Reconstruction error: {cp_error:.6f}")
    print(f"  Explained variance: {cp_explained_var:.4f}")
    print(f"  Number of parameters: {results['CP']['n_parameters']}")
    
    # Tucker decomposition
    shape = tensor.shape
    tucker_rank = (
        min(rank * 2, shape[0]),
        min(rank, shape[1]),
        min(4, shape[2]),
        min(3, shape[3])
    )
    
    print(f"\n2. Fitting Tucker decomposition (rank={tucker_rank})...")
    tucker_model = TuckerDecomposition(rank=tucker_rank, random_state=42)
    tucker_model.fit(tensor, n_iter_max=50, verbose=False)
    
    tucker_error = tucker_model.compute_reconstruction_error(tensor)
    tucker_reconstructed = tucker_model.reconstruct()
    tucker_explained_var = compute_explained_variance(tensor, tucker_reconstructed)
    
    tucker_factors = tucker_model.get_factors()
    tucker_core = tucker_model.get_core_tensor()
    results['Tucker'] = {
        'error': tucker_error,
        'explained_var': tucker_explained_var,
        'n_parameters': sum(f.size for f in tucker_factors) + tucker_core.size
    }
    
    print(f"  Reconstruction error: {tucker_error:.6f}")
    print(f"  Explained variance: {tucker_explained_var:.4f}")
    print(f"  Number of parameters: {results['Tucker']['n_parameters']}")
    
    # Summary comparison
    print("\n" + "=" * 70)
    print("Method Comparison Summary")
    print("=" * 70)
    print(f"\n{'Method':<15} {'Error':<15} {'Explained Var':<15} {'Parameters':<15}")
    print("-" * 60)
    for method, res in results.items():
        print(f"{method:<15} {res['error']:<15.6f} {res['explained_var']:<15.4f} {res['n_parameters']:<15}")
    
    return results


def main():
    """Main function for model comparison."""
    
    # Configuration
    DATA_DIR = Path(__file__).parent.parent / "data" / "merge-vst"
    OUTPUT_DIR = Path(__file__).parent.parent / "results" / "comparison"
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("Tensor Decomposition Model Comparison")
    print("=" * 70)
    
    # Load and prepare data
    print("\n1. Loading and preparing data...")
    loader = DataLoader(DATA_DIR)
    data = loader.load_merged_data()
    
    # Filter and normalize
    data = loader.filter_low_expression(data, min_expression=5.0, min_samples=10)
    data = loader.normalize_data(data, method="zscore", axis=1)
    
    # Build tensor
    builder = TensorBuilder()
    tensor = builder.from_merged_data(data)
    metadata = builder.get_metadata()
    
    # Compare different ranks
    ranks_to_test = [5, 10, 15, 20, 25, 30]
    print(f"\n2. Testing ranks: {ranks_to_test}")
    errors, explained_vars = compare_ranks(tensor, ranks_to_test, metadata, OUTPUT_DIR)
    
    # Compare methods at optimal rank
    optimal_rank = 15  # You can adjust based on rank comparison results
    print(f"\n3. Comparing methods at rank {optimal_rank}...")
    method_results = compare_methods(tensor, optimal_rank, metadata, OUTPUT_DIR)
    
    print("\n" + "=" * 70)
    print("Model comparison complete! Results saved to:", OUTPUT_DIR)
    print("=" * 70)


if __name__ == "__main__":
    main()

