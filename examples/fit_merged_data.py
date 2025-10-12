"""
Example: Fit tensor decomposition models on merged gene expression data.

This script demonstrates how to:
1. Load merged gene expression data
2. Build a 4D tensor (Genes × Individuals × Time × Species)
3. Fit CP and Tucker decompositions
4. Visualize and save results
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import numpy as np
from stdm import (
    DataLoader,
    TensorBuilder,
    CPDecomposition,
    TuckerDecomposition,
    plot_components,
    plot_factor_heatmap,
    plot_reconstruction_error,
    plot_temporal_patterns,
    plot_species_comparison,
    plot_gene_loadings,
)
from stdm.utils import save_results, compute_explained_variance


def main():
    """Main function to fit models on merged data."""
    
    # Configuration
    DATA_DIR = Path(__file__).parent.parent / "data" / "merge-vst"
    OUTPUT_DIR = Path(__file__).parent.parent / "results" / "merged"
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("Tensor Decomposition on Merged Gene Expression Data")
    print("=" * 70)
    
    # Step 1: Load data
    print("\n1. Loading merged data...")
    loader = DataLoader(DATA_DIR)
    data = loader.load_merged_data()
    
    # Optional: Filter low expression genes
    print("\n2. Filtering low expression genes...")
    data = loader.filter_low_expression(data, min_expression=5.0, min_samples=10)
    
    # Optional: Normalize data
    print("\n3. Normalizing data (z-score)...")
    data = loader.normalize_data(data, method="zscore", axis=1)
    
    # Step 2: Build tensor
    print("\n4. Building 4D tensor...")
    builder = TensorBuilder()
    tensor = builder.from_merged_data(data)
    metadata = builder.get_metadata()
    
    print(f"\nTensor statistics:")
    print(f"  Shape: {tensor.shape}")
    print(f"  Min: {tensor.min():.3f}")
    print(f"  Max: {tensor.max():.3f}")
    print(f"  Mean: {tensor.mean():.3f}")
    print(f"  Std: {tensor.std():.3f}")
    
    # Step 3: Fit CP decomposition
    print("\n5. Fitting CP decomposition...")
    cp_rank = 10
    cp_model = CPDecomposition(
        rank=cp_rank,
        non_negative=False,  # Set to True if your data requires non-negativity
        random_state=42
    )
    cp_model.fit(tensor, n_iter_max=100, verbose=True)
    
    cp_error = cp_model.compute_reconstruction_error(tensor)
    cp_reconstructed = cp_model.reconstruct()
    cp_explained_var = compute_explained_variance(tensor, cp_reconstructed)
    
    print(f"\nCP Decomposition Results:")
    print(f"  Rank: {cp_rank}")
    print(f"  Reconstruction error: {cp_error:.6f}")
    print(f"  Explained variance: {cp_explained_var:.4f}")
    
    # Step 4: Fit Tucker decomposition
    print("\n6. Fitting Tucker decomposition...")
    tucker_rank = (100, 20, 4, 3)  # Genes, Individuals, Time, Species
    tucker_model = TuckerDecomposition(
        rank=tucker_rank,
        non_negative=False,
        random_state=42
    )
    tucker_model.fit(tensor, n_iter_max=50, verbose=True)
    
    tucker_error = tucker_model.compute_reconstruction_error(tensor)
    tucker_reconstructed = tucker_model.reconstruct()
    tucker_explained_var = compute_explained_variance(tensor, tucker_reconstructed)
    
    print(f"\nTucker Decomposition Results:")
    print(f"  Rank: {tucker_rank}")
    print(f"  Reconstruction error: {tucker_error:.6f}")
    print(f"  Explained variance: {tucker_explained_var:.4f}")
    
    # Step 5: Save results
    print("\n7. Saving results...")
    cp_factors = cp_model.get_factors()
    tucker_factors = tucker_model.get_factors()
    
    save_results(cp_factors, metadata, OUTPUT_DIR / "cp", prefix="cp")
    save_results(tucker_factors, metadata, OUTPUT_DIR / "tucker", prefix="tucker")
    
    # Step 6: Visualize results
    print("\n8. Creating visualizations...")
    
    # Plot CP components
    mode_names = ['Genes', 'Individuals', 'Time Points', 'Species']
    plot_components(
        cp_factors,
        mode_names=mode_names,
        n_components=5,
        save_path=OUTPUT_DIR / "cp_components.png"
    )
    
    # Plot temporal patterns
    if len(metadata['timepoints']) > 0:
        plot_temporal_patterns(
            cp_factors[2],  # Time mode
            metadata['timepoints'],
            n_components=5,
            save_path=OUTPUT_DIR / "temporal_patterns.png"
        )
    
    # Plot species comparison
    if len(metadata['species']) > 0:
        plot_species_comparison(
            cp_factors[3],  # Species mode
            metadata['species'],
            n_components=5,
            save_path=OUTPUT_DIR / "species_comparison.png"
        )
    
    # Plot top gene loadings for first component
    if len(metadata['genes']) > 0:
        plot_gene_loadings(
            cp_factors[0],  # Gene mode
            gene_names=metadata['genes'],
            component_idx=0,
            top_n=30,
            save_path=OUTPUT_DIR / "top_genes_comp1.png"
        )
    
    # Compare reconstruction errors
    errors = {"CP": cp_error, "Tucker": tucker_error}
    print("\n9. Comparing models...")
    print(f"\nModel Comparison:")
    print(f"  CP reconstruction error: {cp_error:.6f} (explained var: {cp_explained_var:.4f})")
    print(f"  Tucker reconstruction error: {tucker_error:.6f} (explained var: {tucker_explained_var:.4f})")
    
    print("\n" + "=" * 70)
    print("Analysis complete! Results saved to:", OUTPUT_DIR)
    print("=" * 70)


if __name__ == "__main__":
    main()

