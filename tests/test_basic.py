"""
Basic tests for STDM package functionality.
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import numpy as np
import pandas as pd
import pytest
from stdm import DataLoader, TensorBuilder, CPDecomposition, TuckerDecomposition
from stdm.utils import normalize_factors, compute_explained_variance


class TestDataLoader:
    """Test DataLoader functionality."""
    
    def test_parse_sample_names(self):
        """Test parsing of sample names."""
        sample_names = ["ACR-139-TP1", "POR-216-TP2", "POC-201-TP3"]
        species, individuals, timepoints = DataLoader.parse_sample_names(sample_names)
        
        assert len(species) == 3
        assert species[0] == "ACR"
        assert species[1] == "POR"
        assert species[2] == "POC"
        
        assert len(timepoints) == 3
        assert timepoints[0] == 1
        assert timepoints[1] == 2
        assert timepoints[2] == 3
    
    def test_normalize_data_zscore(self):
        """Test z-score normalization."""
        data = pd.DataFrame(np.random.randn(100, 50))
        normalized = DataLoader.normalize_data(data, method="zscore", axis=1)
        
        # Check that each row has mean ~0 and std ~1
        assert np.allclose(normalized.mean(axis=1), 0, atol=1e-10)
        assert np.allclose(normalized.std(axis=1), 1, atol=1e-10)
    
    def test_normalize_data_minmax(self):
        """Test min-max normalization."""
        data = pd.DataFrame(np.random.randn(100, 50))
        normalized = DataLoader.normalize_data(data, method="minmax", axis=1)
        
        # Check that each row is in [0, 1]
        assert normalized.min().min() >= 0
        assert normalized.max().max() <= 1
    
    def test_filter_low_expression(self):
        """Test filtering of low expression genes."""
        data = pd.DataFrame(np.random.randn(100, 50))
        # Set some rows to low values
        data.iloc[:20, :] = 0.1
        
        filtered = DataLoader.filter_low_expression(data, min_expression=1.0, min_samples=5)
        
        # Should filter out low expression genes
        assert filtered.shape[0] < data.shape[0]


class TestTensorBuilder:
    """Test TensorBuilder functionality."""
    
    def test_from_synthetic_merged_data(self):
        """Test building tensor from synthetic merged data."""
        # Create synthetic data
        genes = [f"Gene_{i}" for i in range(100)]
        samples = []
        for species in ["ACR", "POR", "POC"]:
            for ind in range(5):
                for tp in range(1, 5):
                    samples.append(f"{species}-{ind:03d}-TP{tp}")
        
        data = pd.DataFrame(
            np.random.randn(len(genes), len(samples)),
            index=genes,
            columns=samples
        )
        
        builder = TensorBuilder()
        tensor = builder.from_merged_data(data)
        
        # Check tensor shape
        assert tensor.ndim == 4
        assert tensor.shape[0] == 100  # genes
        assert tensor.shape[2] == 4    # time points
        assert tensor.shape[3] == 3    # species
    
    def test_from_separate_data(self):
        """Test building tensor from separate species data."""
        genes = [f"Gene_{i}" for i in range(50)]
        
        species_data = {}
        for species in ["species_a", "species_b"]:
            samples = [f"{species.upper()}.{i}.TP{tp}" for i in range(3) for tp in range(1, 5)]
            species_data[species] = pd.DataFrame(
                np.random.randn(len(genes), len(samples)),
                index=genes,
                columns=samples
            )
        
        builder = TensorBuilder()
        tensor = builder.from_separate_data(species_data)
        
        # Check tensor shape
        assert tensor.ndim == 4
        assert tensor.shape[0] == 50  # genes
        assert tensor.shape[3] == 2   # species
    
    def test_unfold(self):
        """Test tensor unfolding."""
        tensor = np.random.randn(10, 5, 4, 3)
        builder = TensorBuilder()
        builder.tensor = tensor
        
        # Unfold along different modes
        for mode in range(4):
            unfolded = builder.unfold(mode)
            assert unfolded.shape[0] == tensor.shape[mode]
            assert unfolded.size == tensor.size


class TestCPDecomposition:
    """Test CP decomposition."""
    
    def test_fit_small_tensor(self):
        """Test CP decomposition on small tensor."""
        tensor = np.random.randn(20, 10, 4, 3)
        rank = 5
        
        model = CPDecomposition(rank=rank, random_state=42)
        model.fit(tensor, n_iter_max=50, verbose=False)
        
        assert model.fitted
        assert len(model.factors) == 4
        assert model.factors[0].shape == (20, rank)
        assert model.factors[1].shape == (10, rank)
        assert model.factors[2].shape == (4, rank)
        assert model.factors[3].shape == (3, rank)
    
    def test_reconstruct(self):
        """Test tensor reconstruction from CP decomposition."""
        tensor = np.random.randn(15, 8, 4, 3)
        rank = 5
        
        model = CPDecomposition(rank=rank, random_state=42)
        model.fit(tensor, n_iter_max=50, verbose=False)
        
        reconstructed = model.reconstruct()
        
        assert reconstructed.shape == tensor.shape
        
        # Check reconstruction error is reasonable
        error = np.linalg.norm(tensor - reconstructed) / np.linalg.norm(tensor)
        assert error < 1.0  # Should be better than random
    
    def test_non_negative_cp(self):
        """Test non-negative CP decomposition."""
        tensor = np.abs(np.random.randn(15, 8, 4, 3))  # Ensure non-negative input
        rank = 5
        
        model = CPDecomposition(rank=rank, non_negative=True, random_state=42)
        model.fit(tensor, n_iter_max=50, verbose=False)
        
        # Check all factors are non-negative
        for factor in model.factors:
            assert np.all(factor >= 0)


class TestTuckerDecomposition:
    """Test Tucker decomposition."""
    
    def test_fit_small_tensor(self):
        """Test Tucker decomposition on small tensor."""
        tensor = np.random.randn(20, 10, 4, 3)
        rank = (10, 5, 4, 3)
        
        model = TuckerDecomposition(rank=rank, random_state=42)
        model.fit(tensor, n_iter_max=30, verbose=False)
        
        assert model.fitted
        assert len(model.factors) == 4
        assert model.core is not None
        assert model.core.shape == rank
    
    def test_reconstruct(self):
        """Test tensor reconstruction from Tucker decomposition."""
        tensor = np.random.randn(15, 8, 4, 3)
        rank = (10, 5, 4, 3)
        
        model = TuckerDecomposition(rank=rank, random_state=42)
        model.fit(tensor, n_iter_max=30, verbose=False)
        
        reconstructed = model.reconstruct()
        
        assert reconstructed.shape == tensor.shape
        
        # Check reconstruction error
        error = np.linalg.norm(tensor - reconstructed) / np.linalg.norm(tensor)
        assert error < 1.0


class TestUtils:
    """Test utility functions."""
    
    def test_normalize_factors(self):
        """Test factor normalization."""
        factors = [np.random.randn(20, 5), np.random.randn(10, 5)]
        normalized = normalize_factors(factors)
        
        # Check that each column has unit norm
        for factor in normalized:
            norms = np.linalg.norm(factor, axis=0)
            assert np.allclose(norms, 1.0)
    
    def test_compute_explained_variance(self):
        """Test explained variance computation."""
        tensor = np.random.randn(10, 8, 4)
        # Perfect reconstruction
        explained_var = compute_explained_variance(tensor, tensor)
        assert np.isclose(explained_var, 1.0)
        
        # Poor reconstruction
        random_tensor = np.random.randn(*tensor.shape)
        explained_var = compute_explained_variance(tensor, random_tensor)
        assert explained_var < 1.0


def run_tests():
    """Run all tests."""
    pytest.main([__file__, "-v"])


if __name__ == "__main__":
    run_tests()

