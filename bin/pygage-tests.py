#!/usr/bin/env python3
"""
Statistical Tests for GAGE Analysis

Functions for performing various statistical tests on gene sets.

"""

import polars as pl
import argparse
from pathlib import Path
import json

from pygage.tests import GeneSetTests

def main():
    """Command-line interface for gene set statistical tests."""
    parser = argparse.ArgumentParser(
        description='Statistical tests for gene set analysis'
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
        '--method',
        choices=['ks-test', 't-test', 'z-test'],
        default='t-test',
        help='Statistical test method'
    )
    
    parser.add_argument(
        '--min-size',
        type=int,
        default=10,
        help='Minimum gene set size'
    )
    
    parser.add_argument(
        '--max-size',
        type=int,
        default=500,
        help='Maximum gene set size'
    )
    
    parser.add_argument(
        '--output',
        type=Path,
        required=True,
        help='Output file for results'
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
    
    if 'gene_sets' in gene_sets_data:
        gene_sets = gene_sets_data['gene_sets']
    else:
        gene_sets = gene_sets_data
    
    print(f"Loaded {len(gene_sets)} gene sets")
    
    # Run test
    tester = GeneSetTests()
    set_size_range = (args.min_size, args.max_size)
    
    if args.method == 'ks-test':
        results = tester.kolmogorov_smirnov_test(
            expr_data,
            gene_sets,
            gene_col=args.gene_col,
            set_size_range=set_size_range
        )
    elif args.method == 't-test':
        results = tester.t_test(
            expr_data,
            gene_sets,
            gene_col=args.gene_col,
            set_size_range=set_size_range
        )
    else:  # z-test
        results = tester.z_test(
            expr_data,
            gene_sets,
            gene_col=args.gene_col,
            set_size_range=set_size_range
        )
    
    # Write results
    results_df = results['results']
    
    if args.output.suffix == '.csv':
        results_df.write_csv(args.output)
    else:
        results_df.write_csv(args.output, separator='\t')
    
    print(f"Test method: {results['method']}")
    print(f"Tested {results_df.shape[0]} gene sets")
    print(f"Results written to {args.output}")
    
    # Summary statistics
    if 'p_greater' in results_df.columns:
        sig_up = results_df.filter(pl.col('p_greater') < 0.05).shape[0]
        print(f"Significant up-regulated (p < 0.05): {sig_up}")
    
    if 'p_less' in results_df.columns:
        sig_down = results_df.filter(pl.col('p_less') < 0.05).shape[0]
        print(f"Significant down-regulated (p < 0.05): {sig_down}")


if __name__ == '__main__':
    main()
