"""
Utility functions for tensor operations and analysis.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple
from pathlib import Path
import json


def save_results(
    factors: List[np.ndarray],
    metadata: Dict,
    output_dir: str,
    prefix: str = "decomposition"
):
    """
    Save decomposition results to disk.
    
    Parameters
    ----------
    factors : list of np.ndarray
        Factor matrices from decomposition.
    metadata : dict
        Metadata dictionary with dimension labels.
    output_dir : str
        Output directory path.
    prefix : str, default="decomposition"
        Prefix for output files.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Save each factor as CSV
    mode_names = ['genes', 'individuals', 'timepoints', 'species']
    for i, factor in enumerate(factors):
        mode_name = mode_names[i] if i < len(mode_names) else f"mode_{i}"
        filename = output_path / f"{prefix}_factor_{mode_name}.csv"
        
        # Create DataFrame with appropriate labels
        if mode_name in metadata and metadata[mode_name] is not None:
            index_labels = metadata[mode_name]
        else:
            index_labels = [f"{mode_name}_{j}" for j in range(factor.shape[0])]
        
        column_labels = [f"component_{j}" for j in range(factor.shape[1])]
        
        df = pd.DataFrame(factor, index=index_labels, columns=column_labels)
        df.to_csv(filename)
        print(f"Saved {mode_name} factor to {filename}")
    
    # Save metadata as JSON
    metadata_file = output_path / f"{prefix}_metadata.json"
    # Convert numpy types to native Python types for JSON serialization
    metadata_serializable = {}
    for key, value in metadata.items():
        if isinstance(value, np.ndarray):
            metadata_serializable[key] = value.tolist()
        elif isinstance(value, (list, tuple)) and len(value) > 0 and isinstance(value[0], np.integer):
            metadata_serializable[key] = [int(x) for x in value]
        else:
            metadata_serializable[key] = value
    
    with open(metadata_file, 'w') as f:
        json.dump(metadata_serializable, f, indent=2)
    print(f"Saved metadata to {metadata_file}")


def load_results(input_dir: str, prefix: str = "decomposition") -> Tuple[List[np.ndarray], Dict]:
    """
    Load decomposition results from disk.
    
    Parameters
    ----------
    input_dir : str
        Input directory path.
    prefix : str, default="decomposition"
        Prefix of result files.
    
    Returns
    -------
    tuple of (factors, metadata)
        Loaded factor matrices and metadata.
    """
    input_path = Path(input_dir)
    
    # Load metadata
    metadata_file = input_path / f"{prefix}_metadata.json"
    with open(metadata_file, 'r') as f:
        metadata = json.load(f)
    
    # Load factors
    factors = []
    mode_names = ['genes', 'individuals', 'timepoints', 'species']
    
    for mode_name in mode_names:
        filename = input_path / f"{prefix}_factor_{mode_name}.csv"
        if filename.exists():
            df = pd.read_csv(filename, index_col=0)
            factors.append(df.values)
            print(f"Loaded {mode_name} factor from {filename}")
        else:
            break
    
    return factors, metadata


def compute_explained_variance(
    tensor: np.ndarray,
    reconstructed: np.ndarray
) -> float:
    """
    Compute the explained variance ratio.
    
    Parameters
    ----------
    tensor : np.ndarray
        Original tensor.
    reconstructed : np.ndarray
        Reconstructed tensor.
    
    Returns
    -------
    float
        Explained variance ratio (0 to 1).
    """
    ss_total = np.sum((tensor - np.mean(tensor)) ** 2)
    ss_residual = np.sum((tensor - reconstructed) ** 2)
    explained_var = 1 - (ss_residual / ss_total)
    return explained_var


def normalize_factors(factors: List[np.ndarray]) -> List[np.ndarray]:
    """
    Normalize factor matrices to unit norm.
    
    Parameters
    ----------
    factors : list of np.ndarray
        Factor matrices to normalize.
    
    Returns
    -------
    list of np.ndarray
        Normalized factor matrices.
    """
    normalized = []
    for factor in factors:
        norms = np.linalg.norm(factor, axis=0, keepdims=True)
        norms[norms == 0] = 1  # Avoid division by zero
        normalized.append(factor / norms)
    return normalized


def identify_marker_genes(
    gene_factor: np.ndarray,
    gene_names: List[str],
    component_idx: int,
    threshold: float = 0.5,
    top_n: Optional[int] = None
) -> List[Tuple[str, float]]:
    """
    Identify marker genes for a specific component.
    
    Parameters
    ----------
    gene_factor : np.ndarray
        Gene mode factor matrix.
    gene_names : list of str
        List of gene names.
    component_idx : int
        Component index.
    threshold : float, default=0.5
        Minimum absolute loading threshold.
    top_n : int, optional
        Maximum number of genes to return.
    
    Returns
    -------
    list of (gene_name, loading)
        List of marker genes and their loadings.
    """
    component = gene_factor[:, component_idx]
    
    # Filter by threshold
    mask = np.abs(component) >= threshold
    indices = np.where(mask)[0]
    
    # Sort by absolute loading
    sorted_indices = indices[np.argsort(np.abs(component[indices]))[::-1]]
    
    # Limit to top_n if specified
    if top_n is not None:
        sorted_indices = sorted_indices[:top_n]
    
    markers = [(gene_names[i], component[i]) for i in sorted_indices]
    return markers


def compare_components(
    factor1: np.ndarray,
    factor2: np.ndarray,
    component_idx1: int,
    component_idx2: int
) -> float:
    """
    Compute similarity between two components using cosine similarity.
    
    Parameters
    ----------
    factor1 : np.ndarray
        First factor matrix.
    factor2 : np.ndarray
        Second factor matrix.
    component_idx1 : int
        Component index in first factor.
    component_idx2 : int
        Component index in second factor.
    
    Returns
    -------
    float
        Cosine similarity between the two components.
    """
    comp1 = factor1[:, component_idx1]
    comp2 = factor2[:, component_idx2]
    
    # Ensure same length
    min_len = min(len(comp1), len(comp2))
    comp1 = comp1[:min_len]
    comp2 = comp2[:min_len]
    
    # Compute cosine similarity
    dot_product = np.dot(comp1, comp2)
    norm1 = np.linalg.norm(comp1)
    norm2 = np.linalg.norm(comp2)
    
    if norm1 == 0 or norm2 == 0:
        return 0.0
    
    similarity = dot_product / (norm1 * norm2)
    return similarity


def cross_validate_rank(
    tensor: np.ndarray,
    ranks: List[int],
    decomposition_class,
    n_folds: int = 5,
    random_state: int = 42,
    **decomposition_kwargs
) -> Dict[int, float]:
    """
    Perform cross-validation to select optimal rank.
    
    Parameters
    ----------
    tensor : np.ndarray
        Input tensor.
    ranks : list of int
        List of ranks to test.
    decomposition_class : class
        Decomposition class (CPDecomposition or TuckerDecomposition).
    n_folds : int, default=5
        Number of cross-validation folds.
    random_state : int, default=42
        Random seed.
    **decomposition_kwargs
        Additional arguments for decomposition.
    
    Returns
    -------
    dict
        Dictionary mapping rank to average reconstruction error.
    """
    np.random.seed(random_state)
    
    # Create random mask for cross-validation
    mask = np.ones(tensor.shape, dtype=bool)
    flat_indices = np.arange(tensor.size)
    np.random.shuffle(flat_indices)
    
    fold_size = len(flat_indices) // n_folds
    
    errors = {rank: [] for rank in ranks}
    
    for fold in range(n_folds):
        # Create test mask
        test_indices = flat_indices[fold * fold_size:(fold + 1) * fold_size]
        test_mask = np.zeros(tensor.size, dtype=bool)
        test_mask[test_indices] = True
        test_mask = test_mask.reshape(tensor.shape)
        
        # Create training tensor (set test entries to zero)
        train_tensor = tensor.copy()
        train_tensor[test_mask] = 0
        
        for rank in ranks:
            # Fit decomposition
            model = decomposition_class(rank, random_state=random_state, **decomposition_kwargs)
            model.fit(train_tensor, verbose=False)
            
            # Compute error on test set
            reconstructed = model.reconstruct()
            test_error = np.linalg.norm(
                (tensor - reconstructed)[test_mask]
            ) / np.linalg.norm(tensor[test_mask])
            
            errors[rank].append(test_error)
    
    # Average errors across folds
    avg_errors = {rank: np.mean(errors[rank]) for rank in ranks}
    
    return avg_errors

