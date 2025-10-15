"""
Quickstart: Get Component Expression Over Time and Gene Lists

This is the simplest way to run STDM analysis and get the two key outputs:
1. Line plot showing component expression across time points
2. List of top genes in each component

Just run this script and check the results/ directory!
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from stdm import (
    DataLoader,
    TensorBuilder,
    CPDecomposition,
    plot_temporal_expression,
    extract_component_genes,
    create_timestamped_output_dir
)


def main():
    """Run quickstart analysis."""
    
    print("\n" + "=" * 70)
    print("STDM QUICKSTART: Component Expression & Gene Lists")
    print("=" * 70)
    
    # === CONFIGURATION ===
    # Change these paths to match your data
    DATA_PATH = Path(__file__).parent.parent / "data" / "merge-vst" / "vst_counts_matrix.csv"
    N_COMPONENTS = 10  # Number of components to extract
    TOP_GENES = 50     # Top genes per component
    
    # Create output directory with timestamp
    output_dir = create_timestamped_output_dir("results", prefix="quickstart")
    print(f"\n📁 Output directory: {output_dir}\n")
    
    # === STEP 1: LOAD DATA ===
    print("📂 Step 1/4: Loading data...")
    loader = DataLoader(DATA_PATH.parent)
    data = loader.load_merged_data(DATA_PATH.name)
    print(f"   Loaded {len(data)} genes")
    
    # === STEP 2: BUILD TENSOR ===
    print("\n🔨 Step 2/4: Building tensor...")
    builder = TensorBuilder()
    tensor = builder.from_merged_data(data)
    metadata = builder.get_metadata()
    print(f"   Tensor shape: {tensor.shape}")
    print(f"   Dimensions: {len(metadata['genes'])} genes × "
          f"{len(metadata['individuals'])} individuals × "
          f"{len(metadata['timepoints'])} timepoints × "
          f"{len(metadata['species'])} species")
    
    # === STEP 3: FIT MODEL ===
    print(f"\n🚀 Step 3/4: Fitting CP decomposition with {N_COMPONENTS} components...")
    model = CPDecomposition(rank=N_COMPONENTS, random_state=42)
    model.fit(tensor, n_iter_max=100, verbose=False)
    
    error = model.compute_reconstruction_error(tensor)
    print(f"   ✓ Model fitted successfully")
    print(f"   Reconstruction error: {error:.4f}")
    
    # Get factors
    factors = model.get_factors()
    gene_factor = factors[0]      # Genes × Components
    time_factor = factors[2]      # Time × Components
    
    # === STEP 4: GENERATE OUTPUTS ===
    print("\n📊 Step 4/4: Generating outputs...")
    
    # OUTPUT 1: Temporal expression line plot
    plot_path = output_dir / "component_expression_over_time.png"
    plot_temporal_expression(
        time_factor,
        metadata['timepoints'],
        n_components=N_COMPONENTS,
        save_path=plot_path,
        title="Component Expression Across Time Points"
    )
    print(f"   ✓ Temporal plot: {plot_path.name}")
    
    # OUTPUT 2: Gene lists
    extract_component_genes(
        gene_factor,
        metadata['genes'],
        output_dir,
        top_n=TOP_GENES
    )
    
    # === SUMMARY ===
    print("\n" + "=" * 70)
    print("✅ ANALYSIS COMPLETE!")
    print("=" * 70)
    print(f"\n📁 Results saved to: {output_dir}\n")
    print("Key outputs:")
    print(f"  1. 📈 Temporal plot: component_expression_over_time.png")
    print(f"  2. 📝 Gene lists:")
    print(f"     - component_genes_summary.txt (overview)")
    print(f"     - component_genes.csv (full data)")
    print(f"     - component_*_genes.txt (per-component files)")
    print("\n" + "=" * 70)
    
    # Print preview of results
    print("\n📊 PREVIEW: Component Expression Over Time")
    print("-" * 70)
    print(f"{'Timepoint':<12}", end="")
    for i in range(min(5, N_COMPONENTS)):
        print(f"Comp {i+1:>2}  ", end="")
    print()
    print("-" * 70)
    
    for t_idx, timepoint in enumerate(metadata['timepoints']):
        print(f"{str(timepoint):<12}", end="")
        for i in range(min(5, N_COMPONENTS)):
            print(f"{time_factor[t_idx, i]:>8.3f}", end="")
        print()
    
    print("\n📝 PREVIEW: Top 5 Genes in Component 1")
    print("-" * 70)
    component_0 = gene_factor[:, 0]
    top_5_indices = np.argsort(np.abs(component_0))[-5:][::-1]
    
    import numpy as np  # Import here for the preview
    
    for rank, idx in enumerate(top_5_indices, 1):
        gene_name = metadata['genes'][idx]
        loading = component_0[idx]
        print(f"  {rank}. {gene_name:<20} (loading: {loading:>8.4f})")
    
    print(f"\n  See {output_dir}/component_genes_summary.txt for all components")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    import numpy as np  # Needed for main execution
    main()

