"""
Utilities for constructing tensors from gene expression data.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Union


class TensorBuilder:
    """
    Build multi-dimensional tensors from gene expression data.
    
    Supports constructing tensors with dimensions:
    - Genes × Samples × Time Points
    - Genes × Individuals × Time Points × Species
    """
    
    def __init__(self):
        """Initialize the tensor builder."""
        self.tensor = None
        self.metadata = {}
    
    def from_merged_data(
        self,
        data: pd.DataFrame,
        parse_sample_names: bool = True,
        delimiter: str = "-"
    ) -> np.ndarray:
        """
        Build a tensor from merged gene expression data.
        
        Creates a 4D tensor: Genes × Individuals × Time Points × Species
        
        Parameters
        ----------
        data : pd.DataFrame
            Gene expression matrix with samples as columns.
        parse_sample_names : bool, default=True
            Whether to parse sample names to extract metadata.
        delimiter : str, default="-"
            Delimiter used in sample names.
        
        Returns
        -------
        np.ndarray
            4D tensor of gene expression data.
        """
        from stdm.data_loader import DataLoader
        
        # Parse sample names
        sample_names = data.columns.tolist()
        species_list, individual_list, timepoint_list = DataLoader.parse_sample_names(
            sample_names, delimiter
        )
        
        # Get unique values
        genes = data.index.tolist()
        unique_species = sorted(set(species_list))
        unique_individuals = sorted(set(individual_list))
        unique_timepoints = sorted(set(timepoint_list))
        
        # Create tensor
        n_genes = len(genes)
        n_individuals = len(unique_individuals)
        n_timepoints = len(unique_timepoints)
        n_species = len(unique_species)
        
        tensor = np.zeros((n_genes, n_individuals, n_timepoints, n_species))
        
        # Fill tensor
        for i, sample in enumerate(sample_names):
            gene_expr = data[sample].values
            species_idx = unique_species.index(species_list[i])
            individual_idx = unique_individuals.index(individual_list[i])
            timepoint_idx = unique_timepoints.index(timepoint_list[i])
            
            tensor[:, individual_idx, timepoint_idx, species_idx] = gene_expr
        
        # Store metadata
        self.metadata = {
            "genes": genes,
            "individuals": unique_individuals,
            "timepoints": unique_timepoints,
            "species": unique_species,
            "shape": tensor.shape,
            "sample_names": sample_names,
            "species_list": species_list,
            "individual_list": individual_list,
            "timepoint_list": timepoint_list,
        }
        
        self.tensor = tensor
        print(f"Built tensor with shape: {tensor.shape}")
        print(f"  Genes: {n_genes}")
        print(f"  Individuals: {n_individuals}")
        print(f"  Time points: {n_timepoints}")
        print(f"  Species: {n_species}")
        
        return tensor
    
    def from_separate_data(
        self,
        species_data: Dict[str, pd.DataFrame],
        align_genes: bool = True,
        fill_missing: float = 0.0
    ) -> np.ndarray:
        """
        Build a tensor from separate species data.
        
        Creates a 4D tensor: Genes × Individuals × Time Points × Species
        
        Parameters
        ----------
        species_data : dict
            Dictionary mapping species names to expression DataFrames.
        align_genes : bool, default=True
            Whether to align genes across species (use intersection).
        fill_missing : float, default=0.0
            Value to use for missing data.
        
        Returns
        -------
        np.ndarray
            4D tensor of gene expression data.
        """
        from stdm.data_loader import DataLoader
        
        # Find common genes
        if align_genes:
            gene_sets = [set(df.index) for df in species_data.values()]
            common_genes = sorted(list(set.intersection(*gene_sets)))
            print(f"Found {len(common_genes)} common genes across all species")
        else:
            all_genes = set()
            for df in species_data.values():
                all_genes.update(df.index)
            common_genes = sorted(list(all_genes))
            print(f"Using {len(common_genes)} total genes (union across species)")
        
        # Parse metadata for each species
        species_metadata = {}
        all_individuals = set()
        all_timepoints = set()
        
        for species_name, df in species_data.items():
            sample_names = df.columns.tolist()
            _, individuals, timepoints = DataLoader.parse_sample_names(sample_names, ".")
            species_metadata[species_name] = {
                "sample_names": sample_names,
                "individuals": individuals,
                "timepoints": timepoints,
            }
            all_individuals.update(individuals)
            all_timepoints.update(timepoints)
        
        # Sort unique values
        unique_species = sorted(species_data.keys())
        unique_individuals = sorted(all_individuals)
        unique_timepoints = sorted(all_timepoints)
        
        # Create tensor
        n_genes = len(common_genes)
        n_individuals = len(unique_individuals)
        n_timepoints = len(unique_timepoints)
        n_species = len(unique_species)
        
        tensor = np.full(
            (n_genes, n_individuals, n_timepoints, n_species), 
            fill_missing, 
            dtype=np.float32
        )
        
        # Fill tensor
        for species_idx, species_name in enumerate(unique_species):
            df = species_data[species_name]
            metadata = species_metadata[species_name]
            
            for i, sample in enumerate(metadata["sample_names"]):
                if sample not in df.columns:
                    continue
                
                individual = metadata["individuals"][i]
                timepoint = metadata["timepoints"][i]
                
                try:
                    individual_idx = unique_individuals.index(individual)
                    timepoint_idx = unique_timepoints.index(timepoint)
                    
                    # Get gene expression for common genes
                    for gene_idx, gene in enumerate(common_genes):
                        if gene in df.index:
                            tensor[gene_idx, individual_idx, timepoint_idx, species_idx] = \
                                df.loc[gene, sample]
                except ValueError:
                    continue
        
        # Store metadata
        self.metadata = {
            "genes": common_genes,
            "individuals": unique_individuals,
            "timepoints": unique_timepoints,
            "species": unique_species,
            "shape": tensor.shape,
        }
        
        self.tensor = tensor
        print(f"Built tensor with shape: {tensor.shape}")
        print(f"  Genes: {n_genes}")
        print(f"  Individuals: {n_individuals}")
        print(f"  Time points: {n_timepoints}")
        print(f"  Species: {n_species}")
        
        return tensor
    
    def get_metadata(self) -> Dict:
        """
        Get tensor metadata.
        
        Returns
        -------
        dict
            Metadata dictionary containing dimension labels and other info.
        """
        return self.metadata
    
    def unfold(self, mode: int) -> np.ndarray:
        """
        Unfold (matricize) the tensor along a specified mode.
        
        Parameters
        ----------
        mode : int
            Mode along which to unfold (0=genes, 1=individuals, 2=time, 3=species).
        
        Returns
        -------
        np.ndarray
            Unfolded matrix.
        """
        if self.tensor is None:
            raise ValueError("No tensor has been built yet")
        
        return self._unfold_tensor(self.tensor, mode)
    
    @staticmethod
    def _unfold_tensor(tensor: np.ndarray, mode: int) -> np.ndarray:
        """
        Helper function to unfold a tensor.
        
        Parameters
        ----------
        tensor : np.ndarray
            Input tensor.
        mode : int
            Mode along which to unfold.
        
        Returns
        -------
        np.ndarray
            Unfolded matrix.
        """
        shape = tensor.shape
        n_dims = len(shape)
        
        if mode < 0 or mode >= n_dims:
            raise ValueError(f"Invalid mode {mode} for tensor with {n_dims} dimensions")
        
        # Move the unfolding mode to the front
        modes = [mode] + [i for i in range(n_dims) if i != mode]
        tensor_permuted = np.transpose(tensor, modes)
        
        # Reshape to matrix
        matrix = tensor_permuted.reshape(shape[mode], -1)
        return matrix

