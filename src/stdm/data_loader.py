"""
Data loading and preprocessing utilities for gene expression data.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union


class DataLoader:
    """
    Load and preprocess gene expression data from CSV files.
    
    Handles both merged data (all species together) and separate data 
    (one file per species).
    """
    
    def __init__(self, data_dir: Union[str, Path]):
        """
        Initialize the data loader.
        
        Parameters
        ----------
        data_dir : str or Path
            Path to the directory containing gene expression data.
        """
        self.data_dir = Path(data_dir)
        if not self.data_dir.exists():
            raise ValueError(f"Data directory does not exist: {self.data_dir}")
    
    def load_merged_data(self, filename: str = "vst_counts_matrix.csv") -> pd.DataFrame:
        """
        Load merged gene expression data containing all species.
        
        Parameters
        ----------
        filename : str, default="vst_counts_matrix.csv"
            Name of the merged data file.
        
        Returns
        -------
        pd.DataFrame
            Gene expression matrix with genes as rows and samples as columns.
        """
        filepath = self.data_dir / filename
        if not filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")
        
        df = pd.read_csv(filepath, index_col=0)
        print(f"Loaded merged data: {df.shape[0]} genes × {df.shape[1]} samples")
        return df
    
    def load_separate_data(
        self, 
        species_files: Optional[Dict[str, str]] = None
    ) -> Dict[str, pd.DataFrame]:
        """
        Load separate gene expression data for each species.
        
        Parameters
        ----------
        species_files : dict, optional
            Mapping from species name to filename. If None, will attempt to 
            auto-discover files with pattern "*_normalized_expression.csv".
        
        Returns
        -------
        dict
            Dictionary mapping species names to their expression DataFrames.
        """
        if species_files is None:
            # Auto-discover species files
            species_files = {}
            for file in self.data_dir.glob("*_normalized_expression.csv"):
                species_name = file.stem.replace("_normalized_expression", "")
                species_files[species_name] = file.name
        
        data = {}
        for species, filename in species_files.items():
            filepath = self.data_dir / filename
            if not filepath.exists():
                print(f"Warning: File not found for species {species}: {filepath}")
                continue
            
            df = pd.read_csv(filepath, index_col=0)
            data[species] = df
            print(f"Loaded {species}: {df.shape[0]} genes × {df.shape[1]} samples")
        
        return data
    
    @staticmethod
    def parse_sample_names(
        sample_names: List[str], 
        delimiter: str = "-"
    ) -> Tuple[List[str], List[str], List[int]]:
        """
        Parse sample names to extract species, individual IDs, and time points.
        
        Assumes sample format: SPECIES-ID-TPX (e.g., ACR-139-TP1, POR-216-TP2)
        
        Parameters
        ----------
        sample_names : list of str
            List of sample names to parse.
        delimiter : str, default="-"
            Delimiter used in sample names.
        
        Returns
        -------
        tuple of (species_list, individual_list, timepoint_list)
            Parsed components of sample names.
        """
        species = []
        individuals = []
        timepoints = []
        
        for name in sample_names:
            parts = name.split(delimiter)
            if len(parts) >= 3:
                species.append(parts[0])
                individuals.append(f"{parts[0]}-{parts[1]}")
                # Extract time point number (TP1 -> 1, TP2 -> 2, etc.)
                tp_str = parts[2].replace("TP", "")
                timepoints.append(int(tp_str))
            else:
                print(f"Warning: Could not parse sample name: {name}")
                species.append("Unknown")
                individuals.append("Unknown")
                timepoints.append(0)
        
        return species, individuals, timepoints
    
    @staticmethod
    def normalize_data(
        data: pd.DataFrame, 
        method: str = "zscore", 
        axis: int = 1
    ) -> pd.DataFrame:
        """
        Normalize gene expression data.
        
        Parameters
        ----------
        data : pd.DataFrame
            Gene expression matrix.
        method : str, default="zscore"
            Normalization method: "zscore", "minmax", or "log2".
        axis : int, default=1
            Axis along which to normalize (0=genes, 1=samples).
        
        Returns
        -------
        pd.DataFrame
            Normalized gene expression matrix.
        """
        if method == "zscore":
            mean = data.mean(axis=axis, keepdims=True)
            std = data.std(axis=axis, keepdims=True)
            normalized = (data - mean) / (std + 1e-8)
        elif method == "minmax":
            min_val = data.min(axis=axis, keepdims=True)
            max_val = data.max(axis=axis, keepdims=True)
            normalized = (data - min_val) / (max_val - min_val + 1e-8)
        elif method == "log2":
            normalized = np.log2(data + 1)
        else:
            raise ValueError(f"Unknown normalization method: {method}")
        
        return pd.DataFrame(normalized, index=data.index, columns=data.columns)
    
    @staticmethod
    def filter_low_expression(
        data: pd.DataFrame, 
        min_expression: float = 1.0, 
        min_samples: int = 5
    ) -> pd.DataFrame:
        """
        Filter out genes with low expression across samples.
        
        Parameters
        ----------
        data : pd.DataFrame
            Gene expression matrix.
        min_expression : float, default=1.0
            Minimum expression threshold.
        min_samples : int, default=5
            Minimum number of samples that must exceed threshold.
        
        Returns
        -------
        pd.DataFrame
            Filtered gene expression matrix.
        """
        mask = (data > min_expression).sum(axis=1) >= min_samples
        filtered = data[mask]
        print(f"Filtered {data.shape[0] - filtered.shape[0]} low-expression genes")
        print(f"Remaining genes: {filtered.shape[0]}")
        return filtered

