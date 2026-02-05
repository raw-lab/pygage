#!/usr/bin/env python3
"""
Statistical Tests for GAGE Analysis

Functions for performing various statistical tests on gene sets.

"""

import polars as pl
import numpy as np
from scipy import stats
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import json


class GeneSetTests:
    """Statistical testing for gene sets."""
    
    @staticmethod
    def kolmogorov_smirnov_test(expression_data: pl.DataFrame,
                                gene_sets: Dict[str, List[str]],
                                gene_col: str = 'gene_id',
                                set_size_range: Tuple[int, int] = (10, 500),
                                same_dir: bool = True) -> Dict[str, any]:
        """
        Kolmogorov-Smirnov test for gene sets.
        
        Args:
            expression_data: Ranked expression DataFrame
            gene_sets: Dictionary of gene sets
            gene_col: Gene ID column name
            set_size_range: (min, max) gene set size
            same_dir: Test for same direction (up and down separately)
            
        Returns:
            Dictionary with results, p_results, ps_results, and set sizes
        """
        numeric_cols = expression_data.select(pl.col(pl.NUMERIC_DTYPES)).columns
        n_samples = len(numeric_cols)
        n_genes = expression_data.shape[0]
        
        # Rank transform each column
        ranked_data = expression_data.clone()
        for col in numeric_cols:
            ranks = stats.rankdata(expression_data[col].to_numpy())
            ranked_data = ranked_data.with_columns(
                pl.Series(col, ranks)
            )
        
        # All gene ranks
        all_ranks = ranked_data.select(numeric_cols).mean_horizontal().to_numpy()
        
        results = []
        
        for set_name, gene_list in gene_sets.items():
            # Filter to genes in set
            set_data = ranked_data.filter(
                pl.col(gene_col).is_in(gene_list)
            )
            
            set_size = set_data.shape[0]
            
            # Check size constraints
            if set_size < set_size_range[0] or set_size > set_size_range[1]:
                continue
            
            # Get ranks for genes in set
            set_ranks = set_data.select(numeric_cols).mean_horizontal().to_numpy()
            
            # Complement ranks (genes not in set)
            complement_ranks = np.setdiff1d(all_ranks, set_ranks)
            
            # KS test for down-regulation (less)
            ks_stat_less, p_less = stats.ks_2samp(
                set_ranks,
                complement_ranks,
                alternative='less'
            )
            
            # KS test for up-regulation (greater)
            ks_stat_greater, p_greater = stats.ks_2samp(
                set_ranks,
                complement_ranks,
                alternative='greater'
            )
            
            max_stat = max(ks_stat_less, ks_stat_greater)
            
            result = {
                'gene_set': set_name,
                'set_size': set_size,
                'statistic': max_stat,
                'p_greater': p_greater,
                'p_less': p_less if same_dir else None
            }
            
            results.append(result)
        
        return {
            'results': pl.DataFrame(results),
            'method': 'ks-test'
        }
    
    @staticmethod
    def t_test(expression_data: pl.DataFrame,
              gene_sets: Dict[str, List[str]],
              gene_col: str = 'gene_id',
              set_size_range: Tuple[int, int] = (10, 500),
              same_dir: bool = True) -> Dict[str, any]:
        """
        t-test for gene sets.
        
        Args:
            expression_data: Expression DataFrame
            gene_sets: Dictionary of gene sets
            gene_col: Gene ID column name
            set_size_range: (min, max) gene set size
            same_dir: Test for same direction
            
        Returns:
            Dictionary with results, p_results, ps_results, and set sizes
        """
        numeric_cols = expression_data.select(pl.col(pl.NUMERIC_DTYPES)).columns
        
        # Calculate global statistics
        global_means = []
        global_vars = []
        
        for col in numeric_cols:
            values = expression_data[col].to_numpy()
            global_means.append(np.nanmean(values))
            global_vars.append(np.nanvar(values, ddof=1))
        
        global_mean = np.mean(global_means)
        global_var = np.mean(global_vars)
        
        results = []
        
        for set_name, gene_list in gene_sets.items():
            # Filter to genes in set
            set_data = expression_data.filter(
                pl.col(gene_col).is_in(gene_list)
            )
            
            set_size = set_data.shape[0]
            
            # Check size constraints
            if set_size < set_size_range[0] or set_size > set_size_range[1]:
                continue
            
            # Calculate set statistics
            set_means = []
            set_vars = []
            
            for col in numeric_cols:
                values = set_data[col].to_numpy()
                set_means.append(np.nanmean(values))
                set_vars.append(np.nanvar(values, ddof=1))
            
            # Welch's t-test components
            set_mean = np.mean(set_means)
            set_var = np.mean(set_vars)
            
            # Standard error
            se = np.sqrt(set_var / set_size + global_var / set_size)
            
            # t-statistic
            if se > 0:
                t_stat = (set_mean - global_mean) / se
                
                # Degrees of freedom (Welch-Satterthwaite)
                df = ((set_var / set_size + global_var / set_size) ** 2) / \
                     ((set_var / set_size) ** 2 / (set_size - 1) + 
                      (global_var / set_size) ** 2 / (set_size - 1))
                
                # P-values
                p_greater = 1 - stats.t.cdf(t_stat, df)
                p_less = stats.t.cdf(t_stat, df)
            else:
                t_stat = 0
                p_greater = 0.5
                p_less = 0.5
            
            result = {
                'gene_set': set_name,
                'set_size': set_size,
                'statistic': t_stat,
                'p_greater': p_greater,
                'p_less': p_less if same_dir else None
            }
            
            results.append(result)
        
        return {
            'results': pl.DataFrame(results),
            'method': 't-test'
        }
    
    @staticmethod
    def z_test(expression_data: pl.DataFrame,
              gene_sets: Dict[str, List[str]],
              gene_col: str = 'gene_id',
              set_size_range: Tuple[int, int] = (10, 500),
              same_dir: bool = True) -> Dict[str, any]:
        """
        z-test for gene sets.
        
        Args:
            expression_data: Expression DataFrame
            gene_sets: Dictionary of gene sets
            gene_col: Gene ID column name
            set_size_range: (min, max) gene set size
            same_dir: Test for same direction
            
        Returns:
            Dictionary with results, p_results, ps_results, and set sizes
        """
        numeric_cols = expression_data.select(pl.col(pl.NUMERIC_DTYPES)).columns
        
        # Calculate global statistics
        global_means = []
        global_stds = []
        
        for col in numeric_cols:
            values = expression_data[col].to_numpy()
            global_means.append(np.nanmean(values))
            global_stds.append(np.nanstd(values, ddof=1))
        
        global_mean = np.mean(global_means)
        global_std = np.mean(global_stds)
        
        results = []
        
        for set_name, gene_list in gene_sets.items():
            # Filter to genes in set
            set_data = expression_data.filter(
                pl.col(gene_col).is_in(gene_list)
            )
            
            set_size = set_data.shape[0]
            
            # Check size constraints
            if set_size < set_size_range[0] or set_size > set_size_range[1]:
                continue
            
            # Calculate set mean
            set_means = []
            for col in numeric_cols:
                values = set_data[col].to_numpy()
                set_means.append(np.nanmean(values))
            
            set_mean = np.mean(set_means)
            
            # z-statistic
            if global_std > 0:
                z_stat = (set_mean - global_mean) / (global_std / np.sqrt(set_size))
                
                # P-values
                p_greater = 1 - stats.norm.cdf(z_stat)
                p_less = stats.norm.cdf(z_stat)
            else:
                z_stat = 0
                p_greater = 0.5
                p_less = 0.5
            
            result = {
                'gene_set': set_name,
                'set_size': set_size,
                'statistic': z_stat,
                'p_greater': p_greater,
                'p_less': p_less if same_dir else None
            }
            
            results.append(result)
        
        return {
            'results': pl.DataFrame(results),
            'method': 'z-test'
        }


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
