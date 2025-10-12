"""
Data validation for user-provided gene expression data.

Ensures data meets requirements before analysis.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import warnings


class DataValidator:
    """
    Validate user-provided gene expression data.
    
    Checks format, structure, and quality of input data.
    """
    
    def __init__(self, strict: bool = False):
        """
        Initialize validator.
        
        Parameters
        ----------
        strict : bool, default=False
            If True, fail on warnings. If False, issue warnings but continue.
        """
        self.strict = strict
        self.validation_report = []
    
    def validate_merged_data(
        self,
        data_path: Union[str, Path]
    ) -> Tuple[bool, List[str]]:
        """
        Validate merged gene expression CSV file.
        
        Parameters
        ----------
        data_path : str or Path
            Path to CSV file.
        
        Returns
        -------
        tuple of (is_valid, messages)
            Whether data is valid and list of validation messages.
        """
        self.validation_report = []
        data_path = Path(data_path)
        
        # Check file exists
        if not data_path.exists():
            self.validation_report.append(f"❌ File not found: {data_path}")
            return False, self.validation_report
        
        # Check file extension
        if data_path.suffix not in ['.csv', '.tsv', '.txt']:
            self.validation_report.append(
                f"⚠️  Unexpected file extension: {data_path.suffix}. "
                f"Expected .csv, .tsv, or .txt"
            )
        
        try:
            # Try to load data
            data = pd.read_csv(data_path, index_col=0, nrows=5)
            self.validation_report.append("✅ File readable as CSV")
            
        except Exception as e:
            self.validation_report.append(f"❌ Cannot read file: {e}")
            return False, self.validation_report
        
        # Load full data for validation
        try:
            data = pd.read_csv(data_path, index_col=0)
        except Exception as e:
            self.validation_report.append(f"❌ Cannot load full file: {e}")
            return False, self.validation_report
        
        # Check data shape
        n_genes, n_samples = data.shape
        self.validation_report.append(f"📊 Data shape: {n_genes} genes × {n_samples} samples")
        
        if n_genes < 100:
            msg = f"⚠️  Very few genes ({n_genes}). Expected at least 100."
            self.validation_report.append(msg)
            if self.strict:
                return False, self.validation_report
        
        if n_samples < 10:
            msg = f"⚠️  Very few samples ({n_samples}). Expected at least 10."
            self.validation_report.append(msg)
            if self.strict:
                return False, self.validation_report
        
        # Check for missing values
        n_missing = data.isna().sum().sum()
        if n_missing > 0:
            pct_missing = 100 * n_missing / data.size
            msg = f"⚠️  Missing values: {n_missing} ({pct_missing:.2f}%). Will be filled with 0."
            self.validation_report.append(msg)
        else:
            self.validation_report.append("✅ No missing values")
        
        # Check data types
        if not all(data.dtypes.apply(lambda x: np.issubdtype(x, np.number))):
            self.validation_report.append("❌ Data contains non-numeric values")
            return False, self.validation_report
        else:
            self.validation_report.append("✅ All values are numeric")
        
        # Check sample names format
        sample_names = data.columns.tolist()
        valid_formats = self._check_sample_name_format(sample_names)
        
        if valid_formats['merged']:
            self.validation_report.append(
                "✅ Sample names follow expected format: SPECIES-ID-TPX"
            )
        elif valid_formats['separate']:
            self.validation_report.append(
                "✅ Sample names follow expected format: SPECIES.ID.TPX"
            )
        else:
            self.validation_report.append(
                "⚠️  Sample names don't follow standard format. "
                "Expected: SPECIES-ID-TPX or SPECIES.ID.TPX"
            )
            self.validation_report.append(
                f"   Examples: {sample_names[:3]}"
            )
        
        # Check value ranges
        min_val, max_val = data.min().min(), data.max().max()
        self.validation_report.append(f"📈 Value range: [{min_val:.2f}, {max_val:.2f}]")
        
        if min_val < 0:
            self.validation_report.append(
                "⚠️  Negative values detected. "
                "Ensure data is properly normalized if using non-negative decomposition."
            )
        
        if max_val > 100:
            self.validation_report.append(
                "ℹ️  Large values detected. Consider log-transformation if data is count-based."
            )
        
        # Check for zero variance features
        zero_var = (data.std(axis=1) == 0).sum()
        if zero_var > 0:
            self.validation_report.append(
                f"⚠️  {zero_var} genes have zero variance (constant across samples). "
                f"Consider filtering."
            )
        
        # Summary
        self.validation_report.append("\n📋 Validation Summary:")
        self.validation_report.append(f"   Genes: {n_genes}")
        self.validation_report.append(f"   Samples: {n_samples}")
        self.validation_report.append(f"   Data type: Numeric")
        self.validation_report.append(f"   Missing values: {n_missing}")
        
        return True, self.validation_report
    
    def validate_separate_data(
        self,
        data_dir: Union[str, Path],
        file_pattern: str = "*_normalized_expression.csv"
    ) -> Tuple[bool, List[str]]:
        """
        Validate separate species data files.
        
        Parameters
        ----------
        data_dir : str or Path
            Directory containing species files.
        file_pattern : str, default="*_normalized_expression.csv"
            Glob pattern for data files.
        
        Returns
        -------
        tuple of (is_valid, messages)
            Whether data is valid and list of validation messages.
        """
        self.validation_report = []
        data_dir = Path(data_dir)
        
        # Check directory exists
        if not data_dir.exists():
            self.validation_report.append(f"❌ Directory not found: {data_dir}")
            return False, self.validation_report
        
        if not data_dir.is_dir():
            self.validation_report.append(f"❌ Not a directory: {data_dir}")
            return False, self.validation_report
        
        # Find data files
        data_files = list(data_dir.glob(file_pattern))
        
        if len(data_files) == 0:
            self.validation_report.append(
                f"❌ No files matching pattern '{file_pattern}' found in {data_dir}"
            )
            return False, self.validation_report
        
        self.validation_report.append(
            f"✅ Found {len(data_files)} data files"
        )
        
        # Validate each file
        all_gene_sets = []
        for file_path in data_files:
            species = file_path.stem.replace("_normalized_expression", "")
            self.validation_report.append(f"\n📁 Validating {species}:")
            
            try:
                data = pd.read_csv(file_path, index_col=0)
                self.validation_report.append(
                    f"   ✅ {data.shape[0]} genes × {data.shape[1]} samples"
                )
                all_gene_sets.append(set(data.index))
                
            except Exception as e:
                self.validation_report.append(f"   ❌ Error loading: {e}")
                return False, self.validation_report
        
        # Check gene overlap
        if len(all_gene_sets) > 1:
            common_genes = set.intersection(*all_gene_sets)
            union_genes = set.union(*all_gene_sets)
            
            self.validation_report.append(f"\n🧬 Gene Alignment:")
            self.validation_report.append(f"   Common genes: {len(common_genes)}")
            self.validation_report.append(f"   Total unique genes: {len(union_genes)}")
            self.validation_report.append(
                f"   Overlap: {100*len(common_genes)/len(union_genes):.1f}%"
            )
            
            if len(common_genes) < 100:
                self.validation_report.append(
                    f"   ⚠️  Very few common genes. Check gene ID formats."
                )
                if self.strict:
                    return False, self.validation_report
        
        return True, self.validation_report
    
    def _check_sample_name_format(self, sample_names: List[str]) -> Dict[str, bool]:
        """
        Check if sample names follow expected formats.
        
        Parameters
        ----------
        sample_names : list of str
            List of sample names.
        
        Returns
        -------
        dict
            Dictionary with 'merged' and 'separate' format validity.
        """
        # Check for merged format: SPECIES-ID-TPX
        merged_format = 0
        for name in sample_names[:min(10, len(sample_names))]:
            parts = name.split('-')
            if len(parts) >= 3 and 'TP' in parts[-1].upper():
                merged_format += 1
        
        # Check for separate format: SPECIES.ID.TPX
        separate_format = 0
        for name in sample_names[:min(10, len(sample_names))]:
            parts = name.split('.')
            if len(parts) >= 3 and 'TP' in parts[-1].upper():
                separate_format += 1
        
        return {
            'merged': merged_format >= 5,
            'separate': separate_format >= 5
        }
    
    def print_report(self):
        """Print validation report."""
        print("\n" + "="*60)
        print("DATA VALIDATION REPORT")
        print("="*60)
        for message in self.validation_report:
            print(message)
        print("="*60 + "\n")


def validate_and_load(
    data_path: Union[str, Path],
    data_type: str = "auto",
    strict: bool = False
) -> Tuple[bool, Optional[pd.DataFrame], List[str]]:
    """
    Validate and load gene expression data.
    
    Parameters
    ----------
    data_path : str or Path
        Path to data file or directory.
    data_type : str, default="auto"
        Type of data: "merged", "separate", or "auto" to detect.
    strict : bool, default=False
        Whether to fail on warnings.
    
    Returns
    -------
    tuple of (is_valid, data, messages)
        Validation result, loaded data (if valid), and messages.
    """
    validator = DataValidator(strict=strict)
    data_path = Path(data_path)
    
    # Auto-detect data type
    if data_type == "auto":
        if data_path.is_file():
            data_type = "merged"
        elif data_path.is_dir():
            data_type = "separate"
        else:
            return False, None, [f"❌ Cannot determine data type for: {data_path}"]
    
    # Validate
    if data_type == "merged":
        is_valid, messages = validator.validate_merged_data(data_path)
        if is_valid:
            data = pd.read_csv(data_path, index_col=0)
            return True, data, messages
        else:
            return False, None, messages
    
    elif data_type == "separate":
        is_valid, messages = validator.validate_separate_data(data_path)
        return is_valid, None, messages
    
    else:
        return False, None, [f"❌ Unknown data type: {data_type}"]

