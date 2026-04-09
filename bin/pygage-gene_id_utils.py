#!/usr/bin/env python3
"""
Gene ID Conversion Utilities

Functions for converting between Entrez Gene IDs and official gene symbols
for human genes using the egSymb mapping data.

"""

import polars as pl
import argparse
from pathlib import Path

from pygage.gene_id_utils import GeneIDConverter

def main():
    """Command-line interface for gene ID conversion."""
    parser = argparse.ArgumentParser(
        description='Convert between Entrez Gene IDs and gene symbols'
    )
    
    parser.add_argument(
        'input_file',
        type=Path,
        help='Input file with gene IDs or symbols (one per line or CSV/TSV)'
    )
    
    parser.add_argument(
        '--mapping',
        type=Path,
        required=True,
        help='Mapping file (CSV/TSV) with columns: entrez_id, symbol'
    )
    
    parser.add_argument(
        '--direction',
        choices=['eg2sym', 'sym2eg'],
        required=True,
        help='Conversion direction: eg2sym (ID to Symbol) or sym2eg (Symbol to ID)'
    )
    
    parser.add_argument(
        '--output',
        type=Path,
        help='Output file (default: stdout)'
    )
    
    parser.add_argument(
        '--column',
        type=str,
        help='Column name if input is CSV/TSV (otherwise assumes one ID per line)'
    )
    
    args = parser.parse_args()
    
    # Initialize converter
    converter = GeneIDConverter(args.mapping)
    
    # Read input
    if args.input_file.suffix in ['.csv', '.tsv', '.txt']:
        if args.column:
            if args.input_file.suffix == '.csv':
                input_df = pl.read_csv(args.input_file)
            else:
                input_df = pl.read_csv(args.input_file, separator='\t')
            input_values = input_df[args.column].to_list()
        else:
            # Assume single column
            with open(args.input_file) as f:
                input_values = [line.strip() for line in f if line.strip()]
    else:
        with open(args.input_file) as f:
            input_values = [line.strip() for line in f if line.strip()]
    
    # Convert
    if args.direction == 'eg2sym':
        output_df = converter.eg2sym(input_values)
    else:
        output_df = converter.sym2eg(input_values)
    
    # Write output
    if args.output:
        if args.output.suffix == '.csv':
            output_df.write_csv(args.output)
        else:
            output_df.write_csv(args.output, separator='\t')
        print(f"Results written to {args.output}")
    else:
        print(output_df)


if __name__ == '__main__':
    main()
