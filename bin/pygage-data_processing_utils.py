#!/usr/bin/env python3
"""
Data Processing and Expression Analysis Utilities

Functions for data transformation, normalization, and gene extraction.

"""

import polars as pl
import argparse
from pathlib import Path

from pygage.data_processing_utils import DataTransformer, GeneExtractor, GeneDataExporter

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
