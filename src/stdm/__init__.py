"""
STDM: Sparse Tensor Decomposition Models for Gene Expression Data

A Python package for fitting sparse tensor decomposition models optimized for 
gene expression data from multiple species and time points.
"""

__version__ = "0.1.0"
__author__ = "Your Name"

from stdm.data_loader import DataLoader
from stdm.tensor_builder import TensorBuilder
from stdm.decomposition import CPDecomposition, TuckerDecomposition
from stdm.visualization import plot_components, plot_reconstruction_error, plot_factor_heatmap
from stdm.auto_report import AnalysisReport, create_timestamped_output_dir
from stdm.data_validator import DataValidator, validate_and_load

__all__ = [
    "DataLoader",
    "TensorBuilder",
    "CPDecomposition",
    "TuckerDecomposition",
    "plot_components",
    "plot_reconstruction_error",
    "plot_factor_heatmap",
    "AnalysisReport",
    "create_timestamped_output_dir",
    "DataValidator",
    "validate_and_load",
]

