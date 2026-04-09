#!/usr/bin/env python3
"""
Gene ID Conversion Utilities

Functions for converting between Entrez Gene IDs and official gene symbols
for human genes using the egSymb mapping data.

"""

import polars as pl
import argparse
from pathlib import Path
from typing import Optional, Union, List


class GeneIDConverter:
    """Class for converting between Entrez Gene IDs and gene symbols."""
    
    def __init__(self, mapping_file: Optional[Path] = None):
        """
        Initialize the converter with a mapping file.
        
        Args:
            mapping_file: Path to mapping file with columns [entrez_id, symbol]
                         If None, will look for default 'egSymb.csv' or 'egSymb.tsv'
        """
        self.mapping_df = None
        if mapping_file is None:
            mapping_file = Path(__file__).parent / "data" / "egSymb.tsv"
        self.load_mapping(mapping_file)
    
    def load_mapping(self, mapping_file: Path):
        """
        Load gene ID to symbol mapping.
        
        Args:
            mapping_file: Path to mapping file (CSV or TSV)
        """

        if mapping_file.suffix == '.csv':
            self.mapping_df = pl.read_csv(mapping_file, infer_schema=False)
        elif mapping_file.suffix in ['.tsv', '.txt']:
            self.mapping_df = pl.read_csv(mapping_file, infer_schema=False, separator='\t')
        else:
            raise ValueError(f"Unsupported file format: {mapping_file.suffix}")
        
        # Ensure columns are named correctly
        if self.mapping_df.shape[1] >= 2:
            self.mapping_df = self.mapping_df.select([
                pl.col(self.mapping_df.columns[0]).alias('entrez_id'),
                pl.col(self.mapping_df.columns[1]).alias('symbol')
            ])
    
    def eg2sym(self, entrez_ids: Union[List[str], List[int], pl.Series]) -> List[Optional[str]]:
        """
        Convert Entrez Gene IDs to official gene symbols.
        
        Args:
            entrez_ids: List or Series of Entrez Gene IDs
            
        Returns:
            List of gene symbols (None for missing IDs)
        """

        if self.mapping_df is None:
            raise ValueError("Mapping data not loaded. Call load_mapping() first.")
        
        # Convert to list if needed
        if isinstance(entrez_ids, pl.Series):
            entrez_ids = entrez_ids.to_list()
        
        # Convert to strings for consistent matching
        entrez_ids = [str(x) for x in entrez_ids]
        
        # Create lookup DataFrame
        lookup_df = pl.DataFrame({'entrez_id': entrez_ids})
        
        # Join with mapping
        result = lookup_df.join(
            self.mapping_df.select(['entrez_id', 'symbol']),
            on='entrez_id',
            how='left'
        )
        result = result.select([
                pl.col(result.columns[0]).alias('input'),
                pl.col(result.columns[1]).alias('output')
            ])
        
        return result
    
    def sym2eg(self, symbols: Union[List[str], pl.Series]) -> List[Optional[str]]:
        """
        Convert official gene symbols to Entrez Gene IDs.
        
        Args:
            symbols: List or Series of gene symbols
            
        Returns:
            List of Entrez Gene IDs (None for missing symbols)
        """
        if self.mapping_df is None:
            raise ValueError("Mapping data not loaded. Call load_mapping() first.")
        
        # Convert to list if needed
        if isinstance(symbols, pl.Series):
            symbols = symbols.to_list()
        
        # Create lookup DataFrame
        lookup_df = pl.DataFrame({'symbol': symbols})
        
        # Join with mapping
        result = lookup_df.join(
            self.mapping_df.select(['symbol', 'entrez_id']),
            on='symbol',
            how='left'
        )
        result = result.select([
                pl.col(result.columns[0]).alias('input'),
                pl.col(result.columns[1]).alias('output')
            ])
        
        return result
