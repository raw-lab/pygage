#!/usr/bin/env python3
"""
Visualization Utilities

Functions for creating color palettes, Venn diagrams, and heatmaps.

"""

import polars as pl
import argparse
from pathlib import Path

from pygage.visualization_utils import VennDiagram, HeatmapPlotter

def main():
    """Command-line interface for visualization utilities."""
    parser = argparse.ArgumentParser(
        description='Visualization utilities for gene set analysis'
    )
    
    subparsers = parser.add_subparsers(dest='plot_type', help='Plot type')
    
    # Venn diagram
    venn_parser = subparsers.add_parser('venn', help='Create Venn diagram')
    venn_parser.add_argument('--input', type=Path, required=True, help='Input CSV/TSV file')
    venn_parser.add_argument('--names', nargs='+', required=True, help='Set names')
    venn_parser.add_argument('--include', choices=['both', 'up', 'down'], default='both')
    venn_parser.add_argument('--output', type=Path, required=True, help='Output image file')
    
    # Heatmap
    heatmap_parser = subparsers.add_parser('heatmap', help='Create heatmap')
    heatmap_parser.add_argument('--input', type=Path, required=True, help='Input CSV/TSV file')
    heatmap_parser.add_argument('--output', type=Path, required=True, help='Output image file')
    heatmap_parser.add_argument('--cluster', action='store_true', help='Add clustering')
    heatmap_parser.add_argument('--cmap', default='RdYlGn_r', help='Colormap')
    heatmap_parser.add_argument('--title', help='Plot title')
    
    args = parser.parse_args()
    
    if args.plot_type == 'venn':
        # Read data
        if args.input.suffix == '.csv':
            data = pl.read_csv(args.input)
        else:
            data = pl.read_csv(args.input, separator='\t')
        
        # Create Venn counts
        venn = VennDiagram()
        counts = venn.venn_counts(data, include=args.include)
        
        # Plot
        n_sets = data.shape[1]
        if n_sets == 2:
            venn.plot_venn2(counts, args.names, args.output)
        elif n_sets == 3:
            venn.plot_venn3(counts, args.names, args.output)
        else:
            print(f"Error: Can only plot 2 or 3 sets, got {n_sets}")
    
    elif args.plot_type == 'heatmap':
        # Read data
        if args.input.suffix == '.csv':
            data = pl.read_csv(args.input)
        else:
            data = pl.read_csv(args.input, separator='\t')
        
        # Plot
        plotter = HeatmapPlotter()
        if args.cluster:
            plotter.plot_clustered_heatmap(
                data,
                cmap=args.cmap,
                output_file=args.output,
                title=args.title
            )
        else:
            plotter.plot_heatmap(
                data,
                cmap=args.cmap,
                output_file=args.output,
                title=args.title
            )
    
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
