#!/usr/bin/env python3
"""
Pathway and Gene Set Database Utilities

Functions for retrieving gene sets from KEGG pathways and Gene Ontology databases.

"""

import polars as pl
import argparse
import requests
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import json
from collections import defaultdict


class KEGGPathwayRetriever:
    """Retrieve KEGG pathway gene sets."""
    
    KEGG_REST_BASE = "https://rest.kegg.jp"
    
    def __init__(self):
        self.species_info = None
        self.pathways = None
    
    def get_species_code(self, species: str = "hsa") -> Dict[str, str]:
        """
        Get KEGG species code information.
        
        Args:
            species: Species code (e.g., 'hsa' for human)
            
        Returns:
            Dictionary with species information
        """
        url = f"{self.KEGG_REST_BASE}/list/organism"
        response = requests.get(url)
        
        if response.status_code != 200:
            raise RuntimeError(f"Failed to fetch KEGG species list: {response.status_code}")
        
        # Parse organism list
        for line in response.text.strip().split('\n'):
            parts = line.split('\t')
            if len(parts) >= 2 and parts[1] == species:
                self.species_info = {
                    'code': parts[1],
                    'name': parts[2] if len(parts) > 2 else '',
                    'full_name': parts[3] if len(parts) > 3 else ''
                }
                return self.species_info
        
        raise ValueError(f"Unknown species: {species}")
    
    def get_pathway_genes(self, species: str = "hsa", id_type: str = "kegg") -> Dict[str, Dict]:
        """
        Retrieve KEGG gene sets for a species.
        
        Args:
            species: Species code (default: "hsa" for human)
            id_type: ID type - "kegg" or "entrez"
            
        Returns:
            Dictionary containing:
                - gene_sets: Dict[pathway_id, List[gene_ids]]
                - pathway_names: Dict[pathway_id, pathway_name]
                - categories: Dict with signal, metabolic, disease indices
        """
        # Get species info
        species_info = self.get_species_code(species)
        
        # Get pathway list for species
        url = f"{self.KEGG_REST_BASE}/list/pathway/{species}"
        response = requests.get(url)
        
        if response.status_code != 200:
            raise RuntimeError(f"Failed to fetch pathways: {response.status_code}")
        
        pathway_names = {}
        for line in response.text.strip().split('\n'):
            parts = line.split('\t')
            if len(parts) >= 2:
                pathway_id = parts[0].replace(f'path:{species}', '')
                pathway_name = parts[1]
                pathway_names[pathway_id] = pathway_name
        
        # Get genes for each pathway
        gene_sets = {}
        for pathway_id in pathway_names.keys():
            url = f"{self.KEGG_REST_BASE}/link/{species}/pathway:{species}{pathway_id}"
            response = requests.get(url)
            
            if response.status_code == 200:
                genes = []
                for line in response.text.strip().split('\n'):
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        gene = parts[1].replace(f'{species}:', '')
                        genes.append(gene)
                gene_sets[pathway_id] = genes
        
        # Convert to Entrez IDs if requested
        if id_type == "entrez":
            gene_sets = self._convert_to_entrez(species, gene_sets)
        
        return {
            'gene_sets': gene_sets,
            'pathway_names': pathway_names,
            'categories': self._categorize_pathways(pathway_names)
        }
    
    def _convert_to_entrez(self, species: str, gene_sets: Dict[str, List[str]]) -> Dict[str, List[str]]:
        """Convert KEGG gene IDs to Entrez IDs."""
        # Get conversion mapping
        url = f"{self.KEGG_REST_BASE}/conv/{species}/ncbi-geneid"
        response = requests.get(url)
        
        if response.status_code != 200:
            print(f"Warning: Could not convert to Entrez IDs")
            return gene_sets
        
        # Build conversion map
        conversion_map = {}
        for line in response.text.strip().split('\n'):
            parts = line.split('\t')
            if len(parts) >= 2:
                kegg_id = parts[0].replace(f'{species}:', '')
                entrez_id = parts[1].replace('ncbi-geneid:', '')
                conversion_map[kegg_id] = entrez_id
        
        # Convert gene sets
        converted_sets = {}
        for pathway_id, genes in gene_sets.items():
            converted_genes = [conversion_map.get(g, g) for g in genes]
            converted_sets[pathway_id] = converted_genes
        
        return converted_sets
    
    def _categorize_pathways(self, pathway_names: Dict[str, str]) -> Dict[str, List[str]]:
        """Categorize pathways into signal, metabolic, disease."""
        categories = {
            'signal': [],
            'metabolic': [],
            'disease': []
        }
        
        for pathway_id, name in pathway_names.items():
            name_lower = name.lower()
            if any(term in name_lower for term in ['signal', 'pathway', 'receptor']):
                categories['signal'].append(pathway_id)
            elif any(term in name_lower for term in ['metabol', 'biosynthesis', 'degradation']):
                categories['metabolic'].append(pathway_id)
            elif any(term in name_lower for term in ['disease', 'cancer', 'infection']):
                categories['disease'].append(pathway_id)
        
        return categories


class GOGeneSetRetriever:
    """Retrieve Gene Ontology gene sets."""
    
    def __init__(self):
        self.gene_sets = None
    
    def get_go_gene_sets(self, species: str = "human", annotation_file: Optional[Path] = None) -> Dict:
        """
        Retrieve Gene Ontology gene sets.
        
        Args:
            species: Species name (e.g., "human", "mouse")
            annotation_file: Optional GO annotation file (GAF format)
            
        Returns:
            Dictionary containing:
                - gene_sets: Dict[go_term, List[gene_ids]]
                - go_categories: Dict[category, List[go_terms]]
                - go_names: Dict[go_term, description]
        """
        if annotation_file is None:
            raise ValueError(
                "GO annotation file required. Download from "
                "http://geneontology.org/docs/download-go-annotations/"
            )
        
        # Parse GO annotation file (GAF format)
        gene_sets = defaultdict(list)
        go_names = {}
        
        with open(annotation_file) as f:
            for line in f:
                if line.startswith('!'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 13:
                    continue
                
                gene_id = parts[1]  # Gene ID
                go_id = parts[4]     # GO ID
                go_name = parts[9]   # Gene name/description
                
                gene_sets[go_id].append(gene_id)
                if go_id not in go_names:
                    go_names[go_id] = go_name
        
        # Categorize by GO domain (BP, MF, CC)
        go_categories = self._categorize_go_terms(list(gene_sets.keys()))
        
        return {
            'gene_sets': dict(gene_sets),
            'go_categories': go_categories,
            'go_names': go_names
        }
    
    def _categorize_go_terms(self, go_terms: List[str]) -> Dict[str, List[str]]:
        """Categorize GO terms by domain (BP, MF, CC)."""
        # This is simplified - in practice would query GO database
        categories = {
            'BP': [],  # Biological Process
            'MF': [],  # Molecular Function
            'CC': []   # Cellular Component
        }
        
        # Note: Would need GO OBO file to properly categorize
        # This is a placeholder implementation
        
        return categories


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
        results = retriever.get_pathway_genes(
            species=args.species,
            id_type=args.id_type
        )
        
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
