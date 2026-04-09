#!/usr/bin/env python3
"""
GAGE Core Analysis Functions

Core functions for Generally Applicable Gene Set Enrichment (GAGE) analysis.

"""

import polars as pl
import argparse
from pathlib import Path
import json

from pygage.core import GAGEPreparation, GAGEAnalysis


def main():
    """Command-line interface for GAGE analysis."""
    parser = argparse.ArgumentParser(
        description='GAGE (Generally Applicable Gene Set Enrichment) analysis'
    )
    
    parser.add_argument(
        '--expression',
        type=Path,
        required=True,
        help='Expression data file (CSV/TSV)'
    )
    
    parser.add_argument(
        '--gene-sets',
        type=Path,
        required=True,
        help='Gene sets JSON file'
    )
    
    parser.add_argument(
        '--gene-col',
        default='gene_id',
        help='Gene ID column name'
    )
    
    parser.add_argument(
        '--ref-indices',
        nargs='+',
        type=int,
        help='Reference sample column indices'
    )
    
    parser.add_argument(
        '--samp-indices',
        nargs='+',
        type=int,
        help='Sample column indices'
    )
    
    parser.add_argument(
        '--comparison',
        choices=['paired', 'unpaired', '1ongroup', 'as.group'],
        default='paired',
        help='Comparison type'
    )
    
    parser.add_argument(
        '--test-method',
        choices=['t-test', 'z-test', 'ks-test'],
        default='t-test',
        help='Statistical test method'
    )
    
    parser.add_argument(
        '--cutoff',
        type=float,
        default=0.1,
        help='Q-value cutoff for significance'
    )
    
    parser.add_argument(
        '--output',
        type=Path,
        required=True,
        help='Output directory for results'
    )
    
    args = parser.parse_args()
    
    # Read expression data
    if args.expression.suffix == '.csv':
        expr_data = pl.read_csv(args.expression)
    else:
        expr_data = pl.read_csv(args.expression, separator='\t')
    
    print(f"Loaded expression data: {expr_data.shape}")
    
    # Read gene sets
    with open(args.gene_sets) as f:
        gene_sets_data = json.load(f)
    
    # Extract gene sets dictionary
    if 'gene_sets' in gene_sets_data:
        gene_sets = gene_sets_data['gene_sets']
    else:
        gene_sets = gene_sets_data
    
    print(f"Loaded {len(gene_sets)} gene sets")
    
    # Prepare expression data
    prep = GAGEPreparation()
    prepared = prep.prepare_expression(
        expr_data,
        ref_indices=args.ref_indices,
        samp_indices=args.samp_indices,
        comparison=args.comparison
    )

    # Append the gene-col back to the dataframe
    prepared = prepared.with_columns(expr_data.get_column(args.gene_col)) 
    
    # Run GAGE
    
    gage = GAGEAnalysis()
    results = gage.run_gage(
        prepared,
        gene_sets,
        gene_col=args.gene_col,
        test_method=args.test_method
    )
    
    # Filter significant
    significant = gage.filter_significant(cutoff=args.cutoff)
    
    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)
    
    # Write results
    for key, df in significant.items():
        output_file = args.output / f"{key}.tsv"
        df.write_csv(output_file, separator='\t')
        print(f"Wrote {df.shape[0]} rows to {output_file}")
    
    # Summary
    print("\nSummary:")
    print(f"  Up-regulated gene sets: {significant['greater'].shape[0]}")
    if 'less' in significant:
        print(f"  Down-regulated gene sets: {significant['less'].shape[0]}")


if __name__ == '__main__':
    main()
