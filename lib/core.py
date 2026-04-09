#!/usr/bin/env python3
"""
GAGE Core Analysis Functions

Core functions for Generally Applicable Gene Set Enrichment (GAGE) analysis.

"""

import polars as pl
import numpy as np
from scipy import stats
from typing import Dict, List, Optional, Tuple


class GAGEPreparation:
    """GAGE data preparation utilities."""
    
    @staticmethod
    def prepare_expression(data: pl.DataFrame,
                          ref_indices: Optional[List[int]] = None,
                          samp_indices: Optional[List[int]] = None,
                          comparison: str = 'paired',
                          same_dir: bool = True,
                          use_fold: bool = True) -> pl.DataFrame:
        """
        Prepare expression data for GAGE analysis.
        
        Args:
            data: Expression DataFrame
            ref_indices: Reference sample column indices
            samp_indices: Sample column indices
            comparison: 'paired', 'unpaired', '1ongroup', or 'as.group'
            same_dir: Test for same direction
            use_fold: Use fold change (log ratio)
            
        Returns:
            Prepared expression DataFrame
        """
        if ref_indices is None and samp_indices is not None:
            all_indices = set(range(len(data.columns)))
            ref_indices = list(all_indices - set(samp_indices))
        elif samp_indices is None and ref_indices is not None:
            all_indices = set(range(len(data.columns)))
            samp_indices = list(all_indices - set(ref_indices))
        
        if ref_indices is not None and samp_indices is not None:
            ref_cols = [data.columns[i] for i in ref_indices]
            samp_cols = [data.columns[i] for i in samp_indices]
            
            if comparison == 'paired':
                if len(ref_cols) != len(samp_cols):
                    raise ValueError("Paired comparison requires equal number of samples")
                
                result_cols = []
                for ref_col, samp_col in zip(ref_cols, samp_cols):
                    if use_fold:
                        diff = data[samp_col].log() - data[ref_col].log()
                    else:
                        diff = data[samp_col] - data[ref_col]
                    result_cols.append(diff.alias(f"diff_{samp_col}"))
                
                result = pl.DataFrame(result_cols)
            
            elif comparison == 'unpaired':
                result_cols = []
                for samp_col in samp_cols:
                    for ref_col in ref_cols:
                        if use_fold:
                            diff = data[samp_col].log() - data[ref_col].log()
                        else:
                            diff = data[samp_col] - data[ref_col]
                        result_cols.append(diff.alias(f"{samp_col}_vs_{ref_col}"))
                
                result = pl.DataFrame(result_cols)
            
            elif comparison == '1ongroup':
                ref_mean = data.select(ref_cols).mean_horizontal()
                result_cols = []
                for samp_col in samp_cols:
                    if use_fold:
                        diff = data[samp_col].log() - ref_mean.log()
                    else:
                        diff = data[samp_col] - ref_mean
                    result_cols.append(diff.alias(f"{samp_col}_vs_mean"))
                
                result = pl.DataFrame(result_cols)
            
            elif comparison == 'as.group':
                samp_mean = data.select(samp_cols).mean_horizontal()
                ref_mean = data.select(ref_cols).mean_horizontal()
                
                if use_fold:
                    diff = samp_mean.log() - ref_mean.log()
                else:
                    diff = samp_mean - ref_mean
                
                result = pl.DataFrame({'group_diff': diff})
            
            else:
                raise ValueError(f"Unknown comparison type: {comparison}")
        else:
            result = data
        
        # Apply same direction transformation
        if not same_dir:
            numeric_cols = result.select(pl.col(pl.NUMERIC_DTYPES)).columns
            for col in numeric_cols:
                result = result.with_columns(
                    pl.col(col).abs().alias(col)
                )
        
        return result


class GAGEAnalysis:
    """Main GAGE analysis class."""
    
    def __init__(self):
        self.results = None
    
    def run_gage(self,
                 expression_data: pl.DataFrame,
                 gene_sets: Dict[str, List[str]],
                 gene_col: str = 'gene_id',
                 set_size_range: Tuple[int, int] = (10, 500),
                 same_dir: bool = True,
                 test_method: str = 't-test',
                 fdr_method: str = 'BH') -> Dict[str, pl.DataFrame]:
        """
        Run GAGE analysis.
        
        Args:
            expression_data: Prepared expression DataFrame
            gene_sets: Dictionary of gene set name -> list of gene IDs
            gene_col: Gene ID column name
            set_size_range: (min, max) gene set size
            same_dir: Test for same direction
            test_method: Statistical test ('t-test', 'z-test', 'ks-test')
            fdr_method: FDR adjustment method
            
        Returns:
            Dictionary with 'greater', 'less' (if same_dir), and 'stats' DataFrames
        """
        results = []
        
        # Get numeric columns (expression values)
        numeric_cols = expression_data.select(pl.col(pl.NUMERIC_DTYPES)).columns
        
        if len(numeric_cols) == 0:
            raise ValueError("No numeric columns found in expression data")
        
        # Calculate global statistics
        global_mean = expression_data.select(numeric_cols).mean_horizontal().mean()
        global_std = expression_data.select(numeric_cols).mean_horizontal().std()
        
        # Test each gene set
        for set_name, gene_list in gene_sets.items():
            # Filter to genes in set
            set_data = expression_data.filter(
                pl.col(gene_col).is_in(gene_list)
            )
            
            set_size = set_data.shape[0]
            
            # Check size constraints
            if set_size < set_size_range[0] or set_size > set_size_range[1]:
                continue
            
            # Calculate statistics
            set_values = set_data.select(numeric_cols).mean_horizontal().to_numpy()
            
            if test_method == 't-test':
                t_stat, p_greater = stats.ttest_1samp(set_values, global_mean, alternative='greater')
                _, p_less = stats.ttest_1samp(set_values, global_mean, alternative='less')
            elif test_method == 'z-test':
                z_stat = (np.mean(set_values) - global_mean) / (global_std / np.sqrt(set_size))
                p_greater = 1 - stats.norm.cdf(z_stat)
                p_less = stats.norm.cdf(z_stat)
            else:
                # KS test
                _, p_value = stats.ks_2samp(set_values, 
                                           expression_data.select(numeric_cols).mean_horizontal().to_numpy())
                p_greater = p_value / 2
                p_less = p_value / 2
            
            stat_mean = np.mean(set_values)
            
            results.append({
                'gene_set': set_name,
                'set_size': set_size,
                'stat_mean': stat_mean,
                'p_greater': p_greater,
                'p_less': p_less if same_dir else None
            })
        
        # Create results DataFrame
        results_df = pl.DataFrame(results)
        
        # Apply FDR correction
        if fdr_method == 'BH':
            # Benjamini-Hochberg
            n = len(results_df)
            
            # Greater
            greater_df = results_df.sort('p_greater')
            greater_ranks = np.arange(1, n + 1)
            q_greater = np.minimum.accumulate(
                (n / greater_ranks[::-1]) * greater_df['p_greater'].to_numpy()[::-1]
            )[::-1]
            greater_df = greater_df.with_columns(
                pl.Series('q_greater', np.minimum(q_greater, 1.0))
            )
            
            if same_dir:
                # Less
                less_df = results_df.sort('p_less')
                less_ranks = np.arange(1, n + 1)
                q_less = np.minimum.accumulate(
                    (n / less_ranks[::-1]) * less_df['p_less'].to_numpy()[::-1]
                )[::-1]
                less_df = less_df.with_columns(
                    pl.Series('q_less', np.minimum(q_less, 1.0))
                )
            
        # Prepare output
        output = {
            'greater': greater_df.sort('p_greater'),
            'stats': results_df.select(['gene_set', 'stat_mean', 'set_size'])
        }
        
        if same_dir:
            output['less'] = less_df.sort('p_less')
        
        self.results = output
        return output
    
    def filter_significant(self,
                          cutoff: float = 0.1,
                          use_q: bool = True) -> Dict[str, pl.DataFrame]:
        """
        Filter significant gene sets.
        
        Args:
            cutoff: P-value or Q-value cutoff
            use_q: Use Q-value (True) or P-value (False)
            
        Returns:
            Filtered results dictionary
        """
        if self.results is None:
            raise ValueError("No results available. Run run_gage() first.")
        
        col_name = 'q_greater' if use_q else 'p_greater'
        
        filtered = {}
        
        # Filter greater
        filtered['greater'] = self.results['greater'].filter(
            pl.col(col_name) < cutoff
        )
        
        # Filter less if available
        if 'less' in self.results:
            col_name_less = 'q_less' if use_q else 'p_less'
            filtered['less'] = self.results['less'].filter(
                pl.col(col_name_less) < cutoff
            )
        
        filtered['stats'] = self.results['stats']
        
        return filtered
