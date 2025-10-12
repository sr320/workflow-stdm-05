"""
Example: Fit tensor decomposition models on separate species gene expression data.

This script demonstrates how to:
1. Load separate gene expression data for each species
2. Build a 4D tensor aligning genes across species
3. Fit CP and Tucker decompositions
4. Compare patterns across species
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
    plot_components,
    plot_temporal_patterns,
    plot_species_comparison,
)
from stdm.utils import save_results, compute_explained_variance


def analyze_subdirectory(subdir_name: str, output_prefix: str):
    """
    Analyze data from a specific subdirectory.
    
    Parameters
    ----------
    subdir_name : str
        Name of the subdirectory (e.g., 'separate-log', 'separate-syn')
    output_prefix : str
        Prefix for output files
    """
    # Configuration
    DATA_DIR = Path(__file__).parent.parent / "data" / subdir_name
    OUTPUT_DIR = Path(__file__).parent.parent / "results" / output_prefix
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print("\n" + "=" * 70)
    print(f"Analyzing: {subdir_name}")
    print("=" * 70)
    
    # Step 1: Load separate species data
    print("\n1. Loading separate species data...")
    loader = DataLoader(DATA_DIR)
    species_data = loader.load_separate_data()
    
    if len(species_data) == 0:
        print(f"No data found in {DATA_DIR}")
        return
    
    print(f"\nLoaded {len(species_data)} species:")
    for species, df in species_data.items():
        print(f"  {species}: {df.shape[0]} genes × {df.shape[1]} samples")
    
    # Step 2: Optional preprocessing
    print("\n2. Preprocessing data...")
    # Filter low expression genes for each species
    for species in species_data:
        species_data[species] = loader.filter_low_expression(
            species_data[species], 
            min_expression=3.0, 
            min_samples=3
        )
    
    # Normalize each species dataset
    for species in species_data:
        species_data[species] = loader.normalize_data(
            species_data[species], 
            method="zscore", 
            axis=1
        )
    
    # Step 3: Build tensor
    print("\n3. Building 4D tensor from separate data...")
    builder = TensorBuilder()
    tensor = builder.from_separate_data(species_data, align_genes=True)
    metadata = builder.get_metadata()
    
    print(f"\nTensor statistics:")
    print(f"  Shape: {tensor.shape}")
    print(f"  Non-zero entries: {np.count_nonzero(tensor)}")
    print(f"  Sparsity: {1 - np.count_nonzero(tensor) / tensor.size:.2%}")
    
    # Step 4: Fit CP decomposition
    print("\n4. Fitting CP decomposition...")
    cp_rank = 10
    cp_model = CPDecomposition(
        rank=cp_rank,
        non_negative=False,
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
    
    # Step 5: Save results
    print("\n5. Saving results...")
    cp_factors = cp_model.get_factors()
    save_results(cp_factors, metadata, OUTPUT_DIR, prefix="cp")
    
    # Step 6: Visualize results
    print("\n6. Creating visualizations...")
    
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
            cp_factors[2],
            metadata['timepoints'],
            n_components=5,
            save_path=OUTPUT_DIR / "temporal_patterns.png"
        )
    
    # Plot species comparison
    if len(metadata['species']) > 0:
        plot_species_comparison(
            cp_factors[3],
            metadata['species'],
            n_components=5,
            save_path=OUTPUT_DIR / "species_comparison.png"
        )
    
    print(f"\nResults saved to: {OUTPUT_DIR}")


def main():
    """Main function to analyze all separate data directories."""
    
    print("=" * 70)
    print("Tensor Decomposition on Separate Species Data")
    print("=" * 70)
    
    # Analyze log-transformed data
    analyze_subdirectory("separate-log", "separate_log")
    
    # Analyze other normalized data
    analyze_subdirectory("separate-syn", "separate_syn")
    
    print("\n" + "=" * 70)
    print("All analyses complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()

