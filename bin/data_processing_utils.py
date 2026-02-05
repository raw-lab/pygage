#!/usr/bin/env python3
"""
Data Processing and Expression Analysis Utilities

Functions for data transformation, normalization, and gene extraction.

"""

import polars as pl
import numpy as np
import argparse
from pathlib import Path
from typing import Optional, List, Tuple, Union
from visualization_utils import HeatmapPlotter
import matplotlib.pyplot as plt


class DataTransformer:
    """Data transformation utilities."""
    
    @staticmethod
    def row_normalize(data: pl.DataFrame) -> pl.DataFrame:
        """
        Normalize matrix rows to z-scores.
        
        Args:
            data: DataFrame to normalize
            
        Returns:
            Row-normalized DataFrame
        """
        # Get numeric columns
        numeric_cols = data.select(pl.col(pl.NUMERIC_DTYPES)).columns
        
        # Calculate row-wise mean and std
        result = data.clone()
        
        for col in numeric_cols:
            values = data[col].to_numpy()
            mean = np.nanmean(values)
            std = np.nanstd(values)
            
            if std > 0:
                normalized = (values - mean) / std
            else:
                normalized = values - mean
            
            result = result.with_columns(pl.Series(col, normalized))
        
        return result
    
    @staticmethod
    def prepare_paired_data(data: pl.DataFrame,
                           ref_indices: List[int],
                           samp_indices: List[int],
                           comparison: str = 'paired',
                           use_fold: bool = True) -> pl.DataFrame:
        """
        Prepare paired or unpaired expression data.
        
        Args:
            data: Expression DataFrame
            ref_indices: Reference sample column indices
            samp_indices: Sample column indices
            comparison: 'paired' or 'unpaired'
            use_fold: Use fold change (log ratio)
            
        Returns:
            Processed DataFrame
        """
        ref_cols = [data.columns[i] for i in ref_indices]
        samp_cols = [data.columns[i] for i in samp_indices]
        
        if comparison == 'paired':
            if len(ref_cols) != len(samp_cols):
                raise ValueError("Paired comparison requires equal number of ref and samp columns")
            
            result_data = []
            for ref_col, samp_col in zip(ref_cols, samp_cols):
                if use_fold:
                    # Log fold change
                    diff = (
                        data[samp_col].log() - data[ref_col].log()
                    )
                else:
                    diff = data[samp_col] - data[ref_col]
                
                result_data.append(diff.alias(f"{samp_col}_vs_{ref_col}"))
            
            result = pl.concat(result_data, how='horizontal')
            
        elif comparison == 'unpaired':
            result_data = []
            
            for samp_col in samp_cols:
                for ref_col in ref_cols:
                    if use_fold:
                        diff = (
                            data[samp_col].log() - data[ref_col].log()
                        )
                    else:
                        diff = data[samp_col] - data[ref_col]
                    
                    result_data.append(diff.alias(f"{samp_col}_vs_{ref_col}"))
            
            result = pl.concat(result_data, how='horizontal')
        
        else:
            raise ValueError(f"Unknown comparison type: {comparison}")
        
        return result


class GeneExtractor:
    """Gene extraction utilities."""
    
    @staticmethod
    def extract_essential_genes(gene_set: List[str],
                               expression_data: pl.DataFrame,
                               gene_col: str = 'gene_id',
                               threshold: float = 1.0,
                               rank_by_abs: bool = False) -> pl.DataFrame:
        """
        Extract essential genes from a gene set.
        
        Args:
            gene_set: List of gene IDs
            expression_data: Expression DataFrame
            gene_col: Gene ID column name
            threshold: Z-score threshold for essential genes
            rank_by_abs: Rank by absolute value
            
        Returns:
            Filtered DataFrame with essential genes
        """
        # Filter to gene set
        filtered = expression_data.filter(
            pl.col(gene_col).is_in(gene_set)
        )
        
        if filtered.shape[0] == 0:
            return pl.DataFrame()
        
        # Calculate mean expression across samples
        numeric_cols = filtered.select(pl.col(pl.NUMERIC_DTYPES)).columns
        
        if len(numeric_cols) == 0:
            return filtered
        
        # Calculate row means
        means = filtered.select(numeric_cols).mean_horizontal()
        filtered = filtered.with_columns(means.alias('mean_expression'))
        
        # Rank by mean (or absolute mean)
        if rank_by_abs:
            filtered = filtered.with_columns(
                pl.col('mean_expression').abs().alias('abs_mean')
            ).sort('abs_mean', descending=True)
        else:
            filtered = filtered.sort('mean_expression', descending=True)
        
        # Calculate z-scores relative to all genes
        global_mean = expression_data.select(numeric_cols).mean_horizontal().mean()
        global_std = expression_data.select(numeric_cols).mean_horizontal().std()
        
        filtered = filtered.with_columns(
            ((pl.col('mean_expression') - global_mean) / global_std).alias('z_score')
        )
        
        # Filter by threshold
        essential = filtered.filter(
            pl.col('z_score').abs() > threshold
        )
        
        return essential.drop(['mean_expression', 'abs_mean', 'z_score'], strict=False)


class GeneDataExporter:
    """Gene data export and visualization."""
    
    @staticmethod
    def export_gene_data(genes: List[str],
                        expression_data: pl.DataFrame,
                        gene_col: str = 'gene_id',
                        output_file: Optional[Path] = None,
                        create_heatmap: bool = False,
                        heatmap_output: Optional[Path] = None,
                        normalize: bool = True):
        """
        Export and visualize gene expression data.
        
        Args:
            genes: List of gene IDs to export
            expression_data: Expression DataFrame
            gene_col: Gene ID column name
            output_file: Output file for data (CSV/TSV)
            create_heatmap: Generate heatmap
            heatmap_output: Heatmap output file
            normalize: Normalize rows for heatmap
        """
        # Filter to genes
        gene_data = expression_data.filter(
            pl.col(gene_col).is_in(genes)
        )
        
        if gene_data.shape[0] == 0:
            print(f"Warning: No genes found in expression data")
            return
        
        # Export data
        if output_file:
            if output_file.suffix == '.csv':
                gene_data.write_csv(output_file)
            else:
                gene_data.write_csv(output_file, separator='\t')
            print(f"Gene data exported to {output_file}")
        
        # Create heatmap
        if create_heatmap:
            # Get numeric columns
            numeric_cols = gene_data.select(pl.col(pl.NUMERIC_DTYPES)).columns
            
            if len(numeric_cols) > 0:
                # Prepare data for heatmap
                heatmap_data = gene_data.select([gene_col] + numeric_cols)
                
                # Normalize if requested
                if normalize:
                    transformer = DataTransformer()
                    heatmap_data = transformer.row_normalize(heatmap_data)
                
                # Get labels
                row_labels = heatmap_data[gene_col].to_list()
                col_labels = numeric_cols
                
                # Create heatmap
                plotter = HeatmapPlotter()
                plotter.plot_heatmap(
                    heatmap_data.select(numeric_cols),
                    row_labels=row_labels,
                    col_labels=col_labels,
                    output_file=heatmap_output,
                    title='Gene Expression Heatmap',
                    center=0 if normalize else None
                )
    
    @staticmethod
    def create_scatterplot(expression_data: pl.DataFrame,
                          ref_col: str,
                          samp_col: str,
                          gene_col: str = 'gene_id',
                          genes: Optional[List[str]] = None,
                          output_file: Optional[Path] = None,
                          title: Optional[str] = None):
        """
        Create scatterplot comparing reference vs sample.
        
        Args:
            expression_data: Expression DataFrame
            ref_col: Reference column name
            samp_col: Sample column name
            gene_col: Gene ID column
            genes: Optional list of genes to highlight
            output_file: Output file path
            title: Plot title
        """
        # Filter to genes if specified
        if genes:
            plot_data = expression_data.filter(
                pl.col(gene_col).is_in(genes)
            )
        else:
            plot_data = expression_data
        
        # Get data
        x = plot_data[ref_col].to_numpy()
        y = plot_data[samp_col].to_numpy()
        
        # Create plot
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(x, y, alpha=0.5, s=20)
        
        # Add diagonal line
        min_val = min(x.min(), y.min())
        max_val = max(x.max(), y.max())
        ax.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5, linewidth=2)
        
        ax.set_xlabel(f'{ref_col} (Control)', fontsize=12)
        ax.set_ylabel(f'{samp_col} (Experiment)', fontsize=12)
        
        if title:
            ax.set_title(title, fontsize=14, weight='bold')
        else:
            ax.set_title('Expression Comparison', fontsize=14, weight='bold')
        
        ax.set_aspect('equal')
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Scatterplot saved to {output_file}")
        else:
            plt.show()
        
        plt.close()


def main():
    """Command-line interface for data processing."""
    parser = argparse.ArgumentParser(
        description='Data processing utilities for gene expression analysis'
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Command')
    
    # Normalize
    norm_parser = subparsers.add_parser('normalize', help='Normalize expression data')
    norm_parser.add_argument('--input', type=Path, required=True, help='Input file')
    norm_parser.add_argument('--output', type=Path, required=True, help='Output file')
    
    # Extract essential genes
    extract_parser = subparsers.add_parser('extract', help='Extract essential genes')
    extract_parser.add_argument('--input', type=Path, required=True, help='Expression data')
    extract_parser.add_argument('--genes', type=Path, required=True, help='Gene list file')
    extract_parser.add_argument('--output', type=Path, required=True, help='Output file')
    extract_parser.add_argument('--gene-col', default='gene_id', help='Gene ID column')
    extract_parser.add_argument('--threshold', type=float, default=1.0, help='Z-score threshold')
    
    # Export and visualize
    export_parser = subparsers.add_parser('export', help='Export gene data')
    export_parser.add_argument('--input', type=Path, required=True, help='Expression data')
    export_parser.add_argument('--genes', type=Path, required=True, help='Gene list file')
    export_parser.add_argument('--output', type=Path, required=True, help='Output data file')
    export_parser.add_argument('--heatmap', type=Path, help='Heatmap output file')
    export_parser.add_argument('--gene-col', default='gene_id', help='Gene ID column')
    export_parser.add_argument('--normalize', action='store_true', help='Normalize for heatmap')
    
    args = parser.parse_args()
    
    if args.command == 'normalize':
        # Read data
        if args.input.suffix == '.csv':
            data = pl.read_csv(args.input)
        else:
            data = pl.read_csv(args.input, separator='\t')
        
        # Normalize
        transformer = DataTransformer()
        normalized = transformer.row_normalize(data)
        
        # Write output
        if args.output.suffix == '.csv':
            normalized.write_csv(args.output)
        else:
            normalized.write_csv(args.output, separator='\t')
        
        print(f"Normalized data written to {args.output}")
    
    elif args.command == 'extract':
        # Read data
        if args.input.suffix == '.csv':
            data = pl.read_csv(args.input)
        else:
            data = pl.read_csv(args.input, separator='\t')
        
        # Read gene list
        with open(args.genes) as f:
            gene_list = [line.strip() for line in f if line.strip()]
        
        # Extract essential genes
        extractor = GeneExtractor()
        essential = extractor.extract_essential_genes(
            gene_list,
            data,
            gene_col=args.gene_col,
            threshold=args.threshold
        )
        
        # Write output
        if args.output.suffix == '.csv':
            essential.write_csv(args.output)
        else:
            essential.write_csv(args.output, separator='\t')
        
        print(f"Extracted {essential.shape[0]} essential genes to {args.output}")
    
    elif args.command == 'export':
        # Read data
        if args.input.suffix == '.csv':
            data = pl.read_csv(args.input)
        else:
            data = pl.read_csv(args.input, separator='\t')
        
        # Read gene list
        with open(args.genes) as f:
            gene_list = [line.strip() for line in f if line.strip()]
        
        # Export
        exporter = GeneDataExporter()
        exporter.export_gene_data(
            gene_list,
            data,
            gene_col=args.gene_col,
            output_file=args.output,
            create_heatmap=args.heatmap is not None,
            heatmap_output=args.heatmap,
            normalize=args.normalize
        )
    
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
