#!/usr/bin/env python3
"""
GAGE Results Analysis and Comparison

Functions for comparing, filtering, and grouping GAGE results.

"""

import polars as pl
import argparse
from pathlib import Path
import json

from pygage.results_analysis import ResultsComparator, SignificanceFilter, GeneSetGrouper

def main():
    """Command-line interface for results analysis."""
    parser = argparse.ArgumentParser(
        description='Analyze and compare GAGE results'
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Command')
    
    # Compare command
    compare_parser = subparsers.add_parser('compare', help='Compare multiple results')
    compare_parser.add_argument(
        '--inputs',
        nargs='+',
        type=Path,
        required=True,
        help='Input result files'
    )
    compare_parser.add_argument(
        '--names',
        nargs='+',
        required=True,
        help='Sample names'
    )
    compare_parser.add_argument(
        '--cutoff',
        type=float,
        default=0.1,
        help='Q-value cutoff'
    )
    compare_parser.add_argument(
        '--output',
        type=Path,
        required=True,
        help='Output file'
    )
    compare_parser.add_argument(
        '--venn',
        type=Path,
        help='Create Venn diagram (output image file)'
    )
    
    # Filter command
    filter_parser = subparsers.add_parser('filter', help='Filter significant results')
    filter_parser.add_argument(
        '--greater',
        type=Path,
        required=True,
        help='Greater results file'
    )
    filter_parser.add_argument(
        '--less',
        type=Path,
        help='Less results file (optional)'
    )
    filter_parser.add_argument(
        '--cutoff',
        type=float,
        default=0.1,
        help='Q-value cutoff'
    )
    filter_parser.add_argument(
        '--output-dir',
        type=Path,
        required=True,
        help='Output directory'
    )
    
    # Group command
    group_parser = subparsers.add_parser('group', help='Group overlapping gene sets')
    group_parser.add_argument(
        '--results',
        type=Path,
        required=True,
        help='GAGE results file'
    )
    group_parser.add_argument(
        '--gene-sets',
        type=Path,
        required=True,
        help='Gene sets JSON file'
    )
    group_parser.add_argument(
        '--expression',
        type=Path,
        required=True,
        help='Expression data file'
    )
    group_parser.add_argument(
        '--output',
        type=Path,
        required=True,
        help='Output JSON file'
    )
    
    args = parser.parse_args()
    
    if args.command == 'compare':
        comparator = ResultsComparator()
        
        combined = comparator.compare_results(
            args.inputs,
            args.names,
            q_cutoff=args.cutoff,
            output_file=args.output
        )
        
        print(f"\nComparison summary:")
        print(f"  Total gene sets: {combined.shape[0]}")
        print(f"  Significant in all: {combined.filter(pl.col('hits') == len(args.names)).shape[0]}")
        
        if args.venn and len(args.inputs) in [2, 3]:
            comparator.create_venn_comparison(
                args.inputs,
                args.names,
                q_cutoff=args.cutoff,
                output_file=args.venn
            )
    
    elif args.command == 'filter':
        # Read results
        if args.greater.suffix == '.csv':
            greater = pl.read_csv(args.greater)
        else:
            greater = pl.read_csv(args.greater, separator='\t')
        
        results = {'greater': greater}
        
        if args.less:
            if args.less.suffix == '.csv':
                less = pl.read_csv(args.less)
            else:
                less = pl.read_csv(args.less, separator='\t')
            results['less'] = less
        
        # Filter
        filterer = SignificanceFilter()
        filtered = filterer.filter_significant(results, cutoff=args.cutoff)
        
        # Write output
        args.output_dir.mkdir(parents=True, exist_ok=True)
        
        for key, df in filtered.items():
            if df.shape[0] > 0:
                output_file = args.output_dir / f"{key}_significant.tsv"
                df.write_csv(output_file, separator='\t')
                print(f"Wrote {df.shape[0]} {key} gene sets to {output_file}")
    
    elif args.command == 'group':
        # Read results
        if args.results.suffix == '.csv':
            results = pl.read_csv(args.results)
        else:
            results = pl.read_csv(args.results, separator='\t')
        
        # Read gene sets
        with open(args.gene_sets) as f:
            gene_sets_data = json.load(f)
        
        if 'gene_sets' in gene_sets_data:
            gene_sets = gene_sets_data['gene_sets']
        else:
            gene_sets = gene_sets_data
        
        # Read expression data
        if args.expression.suffix == '.csv':
            expr_data = pl.read_csv(args.expression)
        else:
            expr_data = pl.read_csv(args.expression, separator='\t')
        
        # Group
        grouper = GeneSetGrouper()
        groups = grouper.group_gene_sets(
            results,
            gene_sets,
            expr_data,
            output_file=args.output
        )
        
        print(f"\nFound {len(groups)} gene set groups")
        for group_name, members in groups.items():
            print(f"  {group_name}: {len(members)} gene sets")
    
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
