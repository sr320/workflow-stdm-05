"""
Visualization utilities for tensor decomposition results.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Optional, Tuple, Union
import warnings

warnings.filterwarnings('ignore')
sns.set_style('whitegrid')


def plot_components(
    factors: List[np.ndarray],
    mode_names: Optional[List[str]] = None,
    n_components: int = 5,
    figsize: Tuple[int, int] = (15, 10),
    save_path: Optional[str] = None
):
    """
    Plot factor loadings for each mode.
    
    Parameters
    ----------
    factors : list of np.ndarray
        List of factor matrices from decomposition.
    mode_names : list of str, optional
        Names for each mode (e.g., ['Genes', 'Individuals', 'Time', 'Species']).
    n_components : int, default=5
        Number of components to visualize.
    figsize : tuple, default=(15, 10)
        Figure size.
    save_path : str, optional
        Path to save the figure.
    """
    n_modes = len(factors)
    if mode_names is None:
        mode_names = [f"Mode {i}" for i in range(n_modes)]
    
    fig, axes = plt.subplots(n_modes, n_components, figsize=figsize)
    
    if n_modes == 1:
        axes = axes.reshape(1, -1)
    if n_components == 1:
        axes = axes.reshape(-1, 1)
    
    for mode_idx, (factor, mode_name) in enumerate(zip(factors, mode_names)):
        rank = min(factor.shape[1], n_components)
        
        for comp_idx in range(rank):
            ax = axes[mode_idx, comp_idx]
            component = factor[:, comp_idx]
            
            # Plot component weights
            ax.plot(component, linewidth=1.5, alpha=0.7)
            ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5, alpha=0.3)
            
            if mode_idx == 0:
                ax.set_title(f'Component {comp_idx + 1}', fontsize=10, fontweight='bold')
            if comp_idx == 0:
                ax.set_ylabel(mode_name, fontsize=10, fontweight='bold')
            
            ax.tick_params(labelsize=8)
            ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved component plot to {save_path}")
    
    plt.show()


def plot_factor_heatmap(
    factor: np.ndarray,
    labels: Optional[List[str]] = None,
    title: str = "Factor Heatmap",
    figsize: Tuple[int, int] = (12, 8),
    cmap: str = "RdBu_r",
    save_path: Optional[str] = None
):
    """
    Plot a heatmap of a factor matrix.
    
    Parameters
    ----------
    factor : np.ndarray
        Factor matrix to visualize.
    labels : list of str, optional
        Row labels for the heatmap.
    title : str, default="Factor Heatmap"
        Plot title.
    figsize : tuple, default=(12, 8)
        Figure size.
    cmap : str, default="RdBu_r"
        Colormap to use.
    save_path : str, optional
        Path to save the figure.
    """
    plt.figure(figsize=figsize)
    
    # Center colormap at zero
    vmax = np.max(np.abs(factor))
    
    sns.heatmap(
        factor,
        cmap=cmap,
        center=0,
        vmin=-vmax,
        vmax=vmax,
        yticklabels=labels if labels is not None else False,
        xticklabels=[f"C{i+1}" for i in range(factor.shape[1])],
        cbar_kws={'label': 'Loading'},
        linewidths=0
    )
    
    plt.title(title, fontsize=14, fontweight='bold', pad=20)
    plt.xlabel('Component', fontsize=12, fontweight='bold')
    plt.ylabel('Features', fontsize=12, fontweight='bold')
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved heatmap to {save_path}")
    
    plt.show()


def plot_reconstruction_error(
    errors: Union[List[float], Dict[int, float]],
    xlabel: str = "Rank",
    ylabel: str = "Reconstruction Error",
    title: str = "Reconstruction Error vs Rank",
    figsize: Tuple[int, int] = (10, 6),
    save_path: Optional[str] = None
):
    """
    Plot reconstruction error as a function of rank.
    
    Parameters
    ----------
    errors : list or dict
        Reconstruction errors. If dict, keys are ranks and values are errors.
    xlabel : str, default="Rank"
        X-axis label.
    ylabel : str, default="Reconstruction Error"
        Y-axis label.
    title : str
        Plot title.
    figsize : tuple, default=(10, 6)
        Figure size.
    save_path : str, optional
        Path to save the figure.
    """
    plt.figure(figsize=figsize)
    
    if isinstance(errors, dict):
        ranks = sorted(errors.keys())
        error_values = [errors[r] for r in ranks]
    else:
        ranks = list(range(1, len(errors) + 1))
        error_values = errors
    
    plt.plot(ranks, error_values, 'o-', linewidth=2, markersize=8, alpha=0.7)
    plt.xlabel(xlabel, fontsize=12, fontweight='bold')
    plt.ylabel(ylabel, fontsize=12, fontweight='bold')
    plt.title(title, fontsize=14, fontweight='bold', pad=20)
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved error plot to {save_path}")
    
    plt.show()


def plot_temporal_patterns(
    time_factor: np.ndarray,
    timepoints: List[int],
    n_components: int = 5,
    figsize: Tuple[int, int] = (12, 6),
    save_path: Optional[str] = None
):
    """
    Plot temporal patterns from time mode factor.
    
    Parameters
    ----------
    time_factor : np.ndarray
        Time mode factor matrix (timepoints × components).
    timepoints : list of int
        List of time point labels.
    n_components : int, default=5
        Number of components to plot.
    figsize : tuple, default=(12, 6)
        Figure size.
    save_path : str, optional
        Path to save the figure.
    """
    plt.figure(figsize=figsize)
    
    n_comp = min(time_factor.shape[1], n_components)
    
    for i in range(n_comp):
        plt.plot(
            timepoints, 
            time_factor[:, i], 
            'o-', 
            linewidth=2, 
            markersize=8, 
            label=f'Component {i+1}',
            alpha=0.7
        )
    
    plt.xlabel('Time Point', fontsize=12, fontweight='bold')
    plt.ylabel('Loading', fontsize=12, fontweight='bold')
    plt.title('Temporal Patterns', fontsize=14, fontweight='bold', pad=20)
    plt.legend(loc='best', fontsize=10)
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved temporal pattern plot to {save_path}")
    
    plt.show()


def plot_species_comparison(
    species_factor: np.ndarray,
    species_names: List[str],
    n_components: int = 5,
    figsize: Tuple[int, int] = (10, 6),
    save_path: Optional[str] = None
):
    """
    Plot species-specific patterns from species mode factor.
    
    Parameters
    ----------
    species_factor : np.ndarray
        Species mode factor matrix (species × components).
    species_names : list of str
        List of species names.
    n_components : int, default=5
        Number of components to plot.
    figsize : tuple, default=(10, 6)
        Figure size.
    save_path : str, optional
        Path to save the figure.
    """
    n_comp = min(species_factor.shape[1], n_components)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    x = np.arange(len(species_names))
    width = 0.8 / n_comp
    
    for i in range(n_comp):
        offset = width * (i - n_comp / 2 + 0.5)
        ax.bar(
            x + offset, 
            species_factor[:, i], 
            width, 
            label=f'Component {i+1}',
            alpha=0.7
        )
    
    ax.set_xlabel('Species', fontsize=12, fontweight='bold')
    ax.set_ylabel('Loading', fontsize=12, fontweight='bold')
    ax.set_title('Species-Specific Patterns', fontsize=14, fontweight='bold', pad=20)
    ax.set_xticks(x)
    ax.set_xticklabels(species_names, fontsize=10)
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')
    ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved species comparison plot to {save_path}")
    
    plt.show()


def plot_gene_loadings(
    gene_factor: np.ndarray,
    gene_names: Optional[List[str]] = None,
    component_idx: int = 0,
    top_n: int = 20,
    figsize: Tuple[int, int] = (10, 8),
    save_path: Optional[str] = None
):
    """
    Plot top gene loadings for a specific component.
    
    Parameters
    ----------
    gene_factor : np.ndarray
        Gene mode factor matrix (genes × components).
    gene_names : list of str, optional
        List of gene names/IDs.
    component_idx : int, default=0
        Component index to visualize.
    top_n : int, default=20
        Number of top genes to show.
    figsize : tuple, default=(10, 8)
        Figure size.
    save_path : str, optional
        Path to save the figure.
    """
    component = gene_factor[:, component_idx]
    
    # Get top genes by absolute loading
    top_indices = np.argsort(np.abs(component))[-top_n:]
    top_loadings = component[top_indices]
    
    if gene_names is not None:
        top_genes = [gene_names[i] for i in top_indices]
    else:
        top_genes = [f"Gene {i}" for i in top_indices]
    
    # Sort by loading value
    sort_idx = np.argsort(top_loadings)
    top_loadings = top_loadings[sort_idx]
    top_genes = [top_genes[i] for i in sort_idx]
    
    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    colors = ['red' if x < 0 else 'blue' for x in top_loadings]
    
    ax.barh(range(len(top_genes)), top_loadings, color=colors, alpha=0.7)
    ax.set_yticks(range(len(top_genes)))
    ax.set_yticklabels(top_genes, fontsize=9)
    ax.set_xlabel('Loading', fontsize=12, fontweight='bold')
    ax.set_title(f'Top {top_n} Genes - Component {component_idx + 1}', 
                 fontsize=14, fontweight='bold', pad=20)
    ax.axvline(x=0, color='k', linestyle='-', linewidth=1)
    ax.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved gene loadings plot to {save_path}")
    
    plt.show()


def plot_temporal_expression(
    time_factor: np.ndarray,
    timepoints: List[Union[int, str]],
    n_components: Optional[int] = None,
    figsize: Tuple[int, int] = (14, 8),
    save_path: Optional[str] = None,
    title: str = "Component Expression Over Time"
):
    """
    Create a clear line plot showing component expression across time points.
    
    This is the primary visualization for temporal patterns in the data.
    Each line represents one component's activity level over time.
    
    Parameters
    ----------
    time_factor : np.ndarray
        Time mode factor matrix (timepoints × components).
    timepoints : list
        List of time point labels.
    n_components : int, optional
        Number of components to plot. If None, plots all components.
    figsize : tuple, default=(14, 8)
        Figure size.
    save_path : str, optional
        Path to save the figure.
    title : str
        Plot title.
    """
    plt.figure(figsize=figsize)
    
    n_comp = time_factor.shape[1] if n_components is None else min(time_factor.shape[1], n_components)
    
    # Use a colormap for better distinction
    colors = plt.cm.tab10(np.linspace(0, 1, n_comp))
    
    for i in range(n_comp):
        plt.plot(
            timepoints, 
            time_factor[:, i], 
            'o-', 
            linewidth=2.5, 
            markersize=10, 
            label=f'Component {i+1}',
            alpha=0.8,
            color=colors[i]
        )
    
    plt.xlabel('Time Point', fontsize=14, fontweight='bold')
    plt.ylabel('Component Expression', fontsize=14, fontweight='bold')
    plt.title(title, fontsize=16, fontweight='bold', pad=20)
    plt.legend(loc='best', fontsize=11, framealpha=0.9)
    plt.grid(True, alpha=0.3, linestyle='--')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    
    # Add zero line for reference
    plt.axhline(y=0, color='black', linestyle='-', linewidth=0.8, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"✅ Saved temporal expression plot to {save_path}")
    else:
        plt.show()
    
    plt.close()

