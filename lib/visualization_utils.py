#!/usr/bin/env python3
"""
Visualization Utilities

Functions for creating color palettes, Venn diagrams, and heatmaps.

"""

import polars as pl
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle
from pathlib import Path
from typing import List, Optional, Tuple, Union
from itertools import product


class ColorUtils:
    """Color palette utilities."""
    
    @staticmethod
    def create_colormap(low: str, mid: Optional[str], high: str, n: int = 256) -> LinearSegmentedColormap:
        """
        Create a color gradient.
        
        Args:
            low: Low color (e.g., 'green', '#00FF00')
            mid: Optional middle color
            high: High color (e.g., 'red', '#FF0000')
            n: Number of colors
            
        Returns:
            LinearSegmentedColormap
        """
        if mid is None:
            colors = [low, high]
            cmap = LinearSegmentedColormap.from_list('custom', colors, N=n)
        else:
            colors = [low, mid, high]
            cmap = LinearSegmentedColormap.from_list('custom', colors, N=n)
        
        return cmap
    
    @staticmethod
    def greenred(n: int = 256) -> LinearSegmentedColormap:
        """
        Create green-black-red colormap.
        
        Args:
            n: Number of colors
            
        Returns:
            LinearSegmentedColormap
        """
        return ColorUtils.create_colormap('green', 'black', 'red', n)


class VennDiagram:
    """Venn diagram utilities."""
    
    @staticmethod
    def venn_counts(data: pl.DataFrame, include: str = 'both') -> pl.DataFrame:
        """
        Count items in Venn diagram regions.
        
        Args:
            data: DataFrame with boolean or numeric columns
            include: 'both', 'up', or 'down'
            
        Returns:
            DataFrame with counts for each region
        """
        if include == 'up':
            data = (data > 0).cast(pl.Int32)
        elif include == 'down':
            data = (data < 0).cast(pl.Int32)
        else:  # both
            data = (data.with_columns(pl.all().abs()) > 0).cast(pl.Int32)
        
        n_sets = data.shape[1]
        if n_sets > 3:
            raise ValueError("Can't create Venn diagram for more than 3 sets")
        
        # Generate all possible combinations
        combinations = list(product([0, 1], repeat=n_sets))
        
        # Count items in each region
        counts = []
        for combo in combinations:
            mask = pl.lit(True)
            for i, val in enumerate(combo):
                col = data.columns[i]
                if val == 1:
                    mask = mask & (data[col] == 1)
                else:
                    mask = mask & (data[col] == 0)
            
            count = data.filter(mask).shape[0]
            counts.append(list(combo) + [count])
        
        # Create result DataFrame
        result_columns = data.columns + ['Counts']
        return pl.DataFrame(counts, schema=result_columns, orient="row")
    
    @staticmethod
    def plot_venn2(counts: pl.DataFrame, names: List[str], 
                   output_file: Optional[Path] = None, figsize: Tuple[int, int] = (8, 8)):
        """
        Plot 2-set Venn diagram.
        
        Args:
            counts: Venn counts DataFrame
            names: Set names
            output_file: Optional output file path
            figsize: Figure size
        """
        fig, ax = plt.subplots(figsize=figsize)
        ax.set_xlim(-4, 4)
        ax.set_ylim(-4, 4)
        ax.set_aspect('equal')
        ax.axis('off')
        
        # Draw circles
        circle1 = Circle((-1, 0), 1.5, fill=False, edgecolor='blue', linewidth=2)
        circle2 = Circle((1, 0), 1.5, fill=False, edgecolor='red', linewidth=2)
        ax.add_patch(circle1)
        ax.add_patch(circle2)
        
        # Add labels
        ax.text(-1.5, 2, names[0], fontsize=14, ha='center', weight='bold')
        ax.text(1.5, 2, names[1], fontsize=14, ha='center', weight='bold')
        
        # Add counts
        counts_list = counts['Counts'].to_list()
        ax.text(-1.5, 0, str(counts_list[2]), fontsize=12, ha='center')  # Left only
        ax.text(1.5, 0, str(counts_list[1]), fontsize=12, ha='center')   # Right only
        ax.text(0, 0, str(counts_list[3]), fontsize=12, ha='center')     # Overlap
        ax.text(0, -2.5, str(counts_list[0]), fontsize=10, ha='center')  # Neither
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Venn diagram saved to {output_file}")
        else:
            plt.show()
        
        plt.close()
    
    @staticmethod
    def plot_venn3(counts: pl.DataFrame, names: List[str],
                   output_file: Optional[Path] = None, figsize: Tuple[int, int] = (10, 10)):
        """
        Plot 3-set Venn diagram.
        
        Args:
            counts: Venn counts DataFrame
            names: Set names
            output_file: Optional output file path
            figsize: Figure size
        """
        fig, ax = plt.subplots(figsize=figsize)
        ax.set_xlim(-4, 4)
        ax.set_ylim(-4, 4)
        ax.set_aspect('equal')
        ax.axis('off')
        
        # Draw circles
        r = 1.5
        sqrt3 = np.sqrt(3)
        
        circle1 = Circle((-1, 1/sqrt3), r, fill=False, edgecolor='blue', linewidth=2)
        circle2 = Circle((1, 1/sqrt3), r, fill=False, edgecolor='red', linewidth=2)
        circle3 = Circle((0, -2/sqrt3), r, fill=False, edgecolor='green', linewidth=2)
        
        ax.add_patch(circle1)
        ax.add_patch(circle2)
        ax.add_patch(circle3)
        
        # Add labels
        ax.text(-1.5, 2.5, names[0], fontsize=14, ha='center', weight='bold')
        ax.text(1.5, 2.5, names[1], fontsize=14, ha='center', weight='bold')
        ax.text(0, -3.2, names[2], fontsize=14, ha='center', weight='bold')
        
        # Add counts (simplified positioning)
        counts_list = counts['Counts'].to_list()
        positions = [
            (0, -3, counts_list[0]),      # None
            (-1.5, 1, counts_list[4]),    # A only
            (1.5, 1, counts_list[2]),     # B only  
            (0, -1.8, counts_list[1]),    # C only
            (-0.7, -0.3, counts_list[5]), # A & C
            (0.7, -0.3, counts_list[3]),  # B & C
            (0, 1, counts_list[6]),       # A & B
            (0, 0, counts_list[7])        # A & B & C
        ]
        
        for x, y, count in positions:
            ax.text(x, y, str(count), fontsize=11, ha='center')
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Venn diagram saved to {output_file}")
        else:
            plt.show()
        
        plt.close()


class HeatmapPlotter:
    """Heatmap plotting utilities."""
    
    @staticmethod
    def plot_heatmap(data: Union[pl.DataFrame, np.ndarray],
                     row_labels: Optional[List[str]] = None,
                     col_labels: Optional[List[str]] = None,
                     cmap: str = 'RdYlGn_r',
                     vmin: Optional[float] = None,
                     vmax: Optional[float] = None,
                     center: Optional[float] = 0,
                     figsize: Tuple[int, int] = (10, 8),
                     output_file: Optional[Path] = None,
                     title: Optional[str] = None,
                     **kwargs):
        """
        Create a heatmap.
        
        Args:
            data: Data matrix (DataFrame or numpy array)
            row_labels: Row labels
            col_labels: Column labels
            cmap: Colormap name
            vmin: Minimum value for colormap
            vmax: Maximum value for colormap
            center: Center value for diverging colormap
            figsize: Figure size
            output_file: Optional output file path
            title: Plot title
            **kwargs: Additional arguments for seaborn.heatmap
        """
        # Convert to numpy if DataFrame
        if isinstance(data, pl.DataFrame):
            if row_labels is None:
                row_labels = data.columns
            data_array = data.to_numpy()
        else:
            data_array = data
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Plot heatmap
        sns.heatmap(
            data_array,
            ax=ax,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            center=center,
            xticklabels=col_labels if col_labels else False,
            yticklabels=row_labels if row_labels else False,
            cbar_kws={'label': 'Value'},
            **kwargs
        )
        
        if title:
            ax.set_title(title, fontsize=14, weight='bold')
        
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Heatmap saved to {output_file}")
        else:
            plt.show()
        
        plt.close()
    
    @staticmethod
    def plot_clustered_heatmap(data: Union[pl.DataFrame, np.ndarray],
                               row_labels: Optional[List[str]] = None,
                               col_labels: Optional[List[str]] = None,
                               cmap: str = 'RdYlGn_r',
                               vmin: Optional[float] = None,
                               vmax: Optional[float] = None,
                               figsize: Tuple[int, int] = (12, 10),
                               output_file: Optional[Path] = None,
                               title: Optional[str] = None,
                               **kwargs):
        """
        Create a clustered heatmap with dendrograms.
        
        Args:
            data: Data matrix
            row_labels: Row labels
            col_labels: Column labels
            cmap: Colormap name
            vmin: Minimum value
            vmax: Maximum value
            figsize: Figure size
            output_file: Output file path
            title: Plot title
            **kwargs: Additional seaborn.clustermap arguments
        """
        # Convert to pandas for seaborn.clustermap
        if isinstance(data, pl.DataFrame):
            import pandas as pd
            data_pd = data.to_pandas()
            if row_labels:
                data_pd.index = row_labels
            if col_labels:
                data_pd.columns = col_labels
        else:
            import pandas as pd
            data_pd = pd.DataFrame(
                data,
                index=row_labels,
                columns=col_labels
            )
        data_pd.set_index("gene_id", drop=True, inplace=True)

        # Create clustered heatmap
        g = sns.clustermap(
            data_pd,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            figsize=figsize,
            cbar_kws={'label': 'Value'},
            **kwargs
        )
        
        if title:
            g.fig.suptitle(title, fontsize=14, weight='bold', y=0.98)
        
        if output_file:
            g.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Clustered heatmap saved to {output_file}")
        else:
            plt.show()
        
        plt.close()
