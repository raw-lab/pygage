#!/usr/bin/env python3
"""
Pathway and Gene Set Database Utilities

Functions for retrieving gene sets from KEGG pathways and Gene Ontology databases.

"""

import argparse
from pathlib import Path
import json

from pygage.pathway_database_utils import KEGGPathwayRetriever, GOGeneSetRetriever

def main():
    """Command-line interface for pathway retrieval."""
    parser = argparse.ArgumentParser(
        description='Retrieve pathway gene sets from KEGG or GO databases'
    )
    
    subparsers = parser.add_subparsers(dest='database', help='Database to query')
    
    # KEGG subcommand
    kegg_parser = subparsers.add_parser('kegg', help='Retrieve KEGG pathways')
    kegg_parser.add_argument(
        '--species',
        default='hsa',
        help='Species code (default: hsa for human)'
    )
    kegg_parser.add_argument(
        '--id-type',
        choices=['kegg', 'entrez'],
        default='kegg',
        help='Gene ID type (default: kegg)'
    )
    kegg_parser.add_argument(
        '--output',
        type=Path,
        required=True,
        help='Output JSON file for gene sets'
    )
    
    # GO subcommand
    go_parser = subparsers.add_parser('go', help='Retrieve GO gene sets')
    go_parser.add_argument(
        '--annotation-file',
        type=Path,
        required=True,
        help='GO annotation file (GAF format)'
    )
    go_parser.add_argument(
        '--species',
        default='human',
        help='Species name (default: human)'
    )
    go_parser.add_argument(
        '--output',
        type=Path,
        required=True,
        help='Output JSON file for gene sets'
    )
    
    args = parser.parse_args()

    if args.database == 'kegg':
        retriever = KEGGPathwayRetriever()
        print("Retrieving pathway genes...")
        results = retriever.get_pathway_genes(
            species=args.species,
            id_type=args.id_type
        )
        
        print("Saving pathway genes...")
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"Retrieved {len(results['gene_sets'])} KEGG pathways")
        print(f"Results written to {args.output}")
    
    elif args.database == 'go':
        retriever = GOGeneSetRetriever()
        results = retriever.get_go_gene_sets(
            species=args.species,
            annotation_file=args.annotation_file
        )
        
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"Retrieved {len(results['gene_sets'])} GO terms")
        print(f"Results written to {args.output}")
    
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
