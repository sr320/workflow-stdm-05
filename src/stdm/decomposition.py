"""
Sparse tensor decomposition models for gene expression data.
"""

import numpy as np
import tensorly as tl
from tensorly.decomposition import parafac, tucker, non_negative_parafac, non_negative_tucker
from typing import Dict, List, Optional, Tuple, Union
import warnings


class BaseDecomposition:
    """
    Base class for tensor decomposition models.
    """
    
    def __init__(self, rank: Union[int, Tuple[int, ...]], random_state: int = 42):
        """
        Initialize the decomposition model.
        
        Parameters
        ----------
        rank : int or tuple of int
            Rank for the decomposition. For CP, this is a single integer.
            For Tucker, this is a tuple specifying rank for each mode.
        random_state : int, default=42
            Random seed for reproducibility.
        """
        self.rank = rank
        self.random_state = random_state
        self.factors = None
        self.core = None
        self.reconstruction_error = None
        self.fitted = False
    
    def fit(self, tensor: np.ndarray, **kwargs):
        """
        Fit the decomposition model to the tensor.
        
        Parameters
        ----------
        tensor : np.ndarray
            Input tensor to decompose.
        **kwargs
            Additional arguments passed to the decomposition method.
        """
        raise NotImplementedError("Subclasses must implement fit()")
    
    def reconstruct(self) -> np.ndarray:
        """
        Reconstruct the tensor from the decomposition.
        
        Returns
        -------
        np.ndarray
            Reconstructed tensor.
        """
        raise NotImplementedError("Subclasses must implement reconstruct()")
    
    def compute_reconstruction_error(self, tensor: np.ndarray) -> float:
        """
        Compute the reconstruction error.
        
        Parameters
        ----------
        tensor : np.ndarray
            Original tensor.
        
        Returns
        -------
        float
            Relative reconstruction error.
        """
        if not self.fitted:
            raise ValueError("Model must be fitted before computing reconstruction error")
        
        reconstructed = self.reconstruct()
        error = np.linalg.norm(tensor - reconstructed) / np.linalg.norm(tensor)
        self.reconstruction_error = error
        return error
    
    def get_factors(self) -> Union[List[np.ndarray], Tuple]:
        """
        Get the factor matrices.
        
        Returns
        -------
        list or tuple
            Factor matrices from the decomposition.
        """
        if not self.fitted:
            raise ValueError("Model must be fitted before accessing factors")
        return self.factors


class CPDecomposition(BaseDecomposition):
    """
    CP (CANDECOMP/PARAFAC) tensor decomposition.
    
    Decomposes a tensor into a sum of rank-1 tensors:
    X ≈ Σ λᵣ (aᵣ ⊗ bᵣ ⊗ cᵣ ⊗ ...)
    
    Optimized for gene expression data with sparsity and non-negativity constraints.
    """
    
    def __init__(
        self, 
        rank: int, 
        random_state: int = 42,
        non_negative: bool = False,
        l2_reg: float = 0.0,
        sparse_component: Optional[List[int]] = None
    ):
        """
        Initialize CP decomposition.
        
        Parameters
        ----------
        rank : int
            Number of components.
        random_state : int, default=42
            Random seed for reproducibility.
        non_negative : bool, default=False
            Whether to enforce non-negativity constraints.
        l2_reg : float, default=0.0
            L2 regularization parameter for sparsity.
        sparse_component : list of int, optional
            Indices of modes to apply sparsity regularization.
        """
        super().__init__(rank, random_state)
        self.non_negative = non_negative
        self.l2_reg = l2_reg
        self.sparse_component = sparse_component
        self.weights = None
    
    def fit(
        self, 
        tensor: np.ndarray, 
        init: str = "random",
        n_iter_max: int = 100,
        tol: float = 1e-6,
        verbose: bool = True,
        **kwargs
    ):
        """
        Fit CP decomposition to the tensor.
        
        Parameters
        ----------
        tensor : np.ndarray
            Input tensor to decompose.
        init : str, default="random"
            Initialization method: "random" or "svd".
        n_iter_max : int, default=100
            Maximum number of iterations.
        tol : float, default=1e-6
            Convergence tolerance.
        verbose : bool, default=True
            Whether to print progress.
        **kwargs
            Additional arguments passed to the decomposition.
        """
        tl.set_backend('numpy')
        np.random.seed(self.random_state)
        
        if verbose:
            print(f"Fitting CP decomposition with rank={self.rank}")
            print(f"  Tensor shape: {tensor.shape}")
            print(f"  Non-negative: {self.non_negative}")
            print(f"  L2 regularization: {self.l2_reg}")
        
        # Handle missing values (zeros or NaNs)
        mask = ~np.isnan(tensor)
        if not np.all(mask):
            if verbose:
                print(f"  Detected missing values: {np.sum(~mask)} entries")
            tensor = np.nan_to_num(tensor, nan=0.0)
        
        try:
            if self.non_negative:
                # Non-negative CP decomposition
                result = non_negative_parafac(
                    tensor,
                    rank=self.rank,
                    init=init,
                    n_iter_max=n_iter_max,
                    tol=tol,
                    verbose=verbose,
                    random_state=self.random_state,
                    **kwargs
                )
            else:
                # Standard CP decomposition
                result = parafac(
                    tensor,
                    rank=self.rank,
                    init=init,
                    n_iter_max=n_iter_max,
                    tol=tol,
                    verbose=verbose,
                    random_state=self.random_state,
                    **kwargs
                )
            
            # Extract weights and factors
            if isinstance(result, tuple):
                self.weights, self.factors = result
            else:
                self.weights = result[0]
                self.factors = result[1]
            
            # Apply sparsity regularization if specified
            if self.l2_reg > 0 and self.sparse_component is not None:
                self._apply_sparsity()
            
            self.fitted = True
            
            if verbose:
                error = self.compute_reconstruction_error(tensor)
                print(f"  Reconstruction error: {error:.6f}")
        
        except Exception as e:
            print(f"Error during CP decomposition: {e}")
            raise
    
    def reconstruct(self) -> np.ndarray:
        """
        Reconstruct the tensor from CP factors.
        
        Returns
        -------
        np.ndarray
            Reconstructed tensor.
        """
        if not self.fitted:
            raise ValueError("Model must be fitted before reconstruction")
        
        from tensorly.cp_tensor import cp_to_tensor
        return cp_to_tensor((self.weights, self.factors))
    
    def _apply_sparsity(self):
        """Apply L2 regularization to promote sparsity in selected components."""
        for mode_idx in self.sparse_component:
            if mode_idx < len(self.factors):
                factor = self.factors[mode_idx]
                # Soft thresholding
                threshold = self.l2_reg * np.std(factor)
                self.factors[mode_idx] = np.sign(factor) * np.maximum(
                    np.abs(factor) - threshold, 0
                )
    
    def get_gene_components(self, top_n: int = 50) -> Dict[int, List[Tuple[int, float]]]:
        """
        Get top genes for each component.
        
        Parameters
        ----------
        top_n : int, default=50
            Number of top genes to return per component.
        
        Returns
        -------
        dict
            Dictionary mapping component index to list of (gene_idx, weight) tuples.
        """
        if not self.fitted:
            raise ValueError("Model must be fitted first")
        
        gene_factor = self.factors[0]  # Genes are typically the first mode
        components = {}
        
        for r in range(self.rank):
            component_weights = gene_factor[:, r]
            top_indices = np.argsort(np.abs(component_weights))[-top_n:][::-1]
            components[r] = [(idx, component_weights[idx]) for idx in top_indices]
        
        return components


class TuckerDecomposition(BaseDecomposition):
    """
    Tucker tensor decomposition.
    
    Decomposes a tensor into a core tensor and factor matrices:
    X ≈ G ×₁ A ×₂ B ×₃ C ×₄ D
    
    Provides more flexibility than CP decomposition with a dense core tensor.
    """
    
    def __init__(
        self, 
        rank: Tuple[int, ...], 
        random_state: int = 42,
        non_negative: bool = False
    ):
        """
        Initialize Tucker decomposition.
        
        Parameters
        ----------
        rank : tuple of int
            Rank for each mode (e.g., (100, 20, 4, 3) for genes, individuals, time, species).
        random_state : int, default=42
            Random seed for reproducibility.
        non_negative : bool, default=False
            Whether to enforce non-negativity constraints.
        """
        super().__init__(rank, random_state)
        self.non_negative = non_negative
    
    def fit(
        self, 
        tensor: np.ndarray, 
        init: str = "random",
        n_iter_max: int = 100,
        tol: float = 1e-5,
        verbose: bool = True,
        **kwargs
    ):
        """
        Fit Tucker decomposition to the tensor.
        
        Parameters
        ----------
        tensor : np.ndarray
            Input tensor to decompose.
        init : str, default="random"
            Initialization method.
        n_iter_max : int, default=100
            Maximum number of iterations.
        tol : float, default=1e-5
            Convergence tolerance.
        verbose : bool, default=True
            Whether to print progress.
        **kwargs
            Additional arguments.
        """
        tl.set_backend('numpy')
        np.random.seed(self.random_state)
        
        if verbose:
            print(f"Fitting Tucker decomposition with rank={self.rank}")
            print(f"  Tensor shape: {tensor.shape}")
            print(f"  Non-negative: {self.non_negative}")
        
        # Handle missing values
        tensor = np.nan_to_num(tensor, nan=0.0)
        
        try:
            if self.non_negative:
                # Non-negative Tucker decomposition
                result = non_negative_tucker(
                    tensor,
                    rank=self.rank,
                    init=init,
                    n_iter_max=n_iter_max,
                    tol=tol,
                    verbose=verbose,
                    random_state=self.random_state,
                    **kwargs
                )
            else:
                # Standard Tucker decomposition
                result = tucker(
                    tensor,
                    rank=self.rank,
                    init=init,
                    n_iter_max=n_iter_max,
                    tol=tol,
                    verbose=verbose,
                    random_state=self.random_state,
                    **kwargs
                )
            
            # Extract core and factors
            if isinstance(result, tuple):
                self.core, self.factors = result
            else:
                self.core = result[0]
                self.factors = result[1]
            
            self.fitted = True
            
            if verbose:
                error = self.compute_reconstruction_error(tensor)
                print(f"  Reconstruction error: {error:.6f}")
                print(f"  Core tensor shape: {self.core.shape}")
        
        except Exception as e:
            print(f"Error during Tucker decomposition: {e}")
            raise
    
    def reconstruct(self) -> np.ndarray:
        """
        Reconstruct the tensor from Tucker factors.
        
        Returns
        -------
        np.ndarray
            Reconstructed tensor.
        """
        if not self.fitted:
            raise ValueError("Model must be fitted before reconstruction")
        
        from tensorly.tucker_tensor import tucker_to_tensor
        return tucker_to_tensor((self.core, self.factors))
    
    def get_core_tensor(self) -> np.ndarray:
        """
        Get the core tensor.
        
        Returns
        -------
        np.ndarray
            Core tensor from Tucker decomposition.
        """
        if not self.fitted:
            raise ValueError("Model must be fitted first")
        return self.core

