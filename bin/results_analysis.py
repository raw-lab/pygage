#!/usr/bin/env python3
"""
GAGE Results Analysis and Comparison

Functions for comparing, filtering, and grouping GAGE results.

"""

import polars as pl
import numpy as np
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from visualization_utils import VennDiagram, HeatmapPlotter
import json


class ResultsComparator:
    """Compare differential expression results across multiple datasets."""
    
    @staticmethod
    def compare_results(result_files: List[Path],
                       sample_names: List[str],
                       q_cutoff: float = 0.1,
                       output_file: Optional[Path] = None) -> pl.DataFrame:
        """
        Compare results across multiple datasets.
        
        Args:
            result_files: List of result file paths
            sample_names: Names for each dataset
            q_cutoff: Q-value cutoff for significance
            output_file: Optional output file
            
        Returns:
            Combined results DataFrame
        """
        if len(result_files) != len(sample_names):
            raise ValueError("Number of files must match number of sample names")
        
        # Read all results
        all_results = []
        for file_path, name in zip(result_files, sample_names):
            if file_path.suffix == '.csv':
                df = pl.read_csv(file_path)
            else:
                df = pl.read_csv(file_path, separator='\t')
            
            # Rename columns
            df = df.rename({
                'stat_mean': f'{name}_stat',
                'q_val': f'{name}_q' if 'q_val' in df.columns else f'{name}_p',
                'q_greater': f'{name}_q' if 'q_greater' in df.columns else f'{name}_p',
            })
            
            all_results.append(df)
        
        # Merge all results
        combined = all_results[0]
        gene_set_col = 'gene_set' if 'gene_set' in combined.columns else combined.columns[0]
        
        for df in all_results[1:]:
            combined = combined.join(
                df,
                on=gene_set_col,
                how='outer'
            )
        
        # Calculate hits (number of significant results)
        q_cols = [col for col in combined.columns if col.endswith('_q') or col.endswith('_p')]
        
        # Count non-significant results
        nsig_expr = pl.lit(0)
        for col in q_cols:
            nsig_expr = nsig_expr + (pl.col(col) > q_cutoff).cast(pl.Int32)
        
        combined = combined.with_columns(
            (len(q_cols) - nsig_expr).alias('hits')
        )
        
        # Sort by hits and p-values
        combined = combined.sort('hits', descending=True)
        
        # Write output
        if output_file:
            if output_file.suffix == '.csv':
                combined.write_csv(output_file)
            else:
                combined.write_csv(output_file, separator='\t')
            print(f"Combined results written to {output_file}")
        
        return combined
    
    @staticmethod
    def create_venn_comparison(result_files: List[Path],
                              sample_names: List[str],
                              q_cutoff: float = 0.1,
                              output_file: Optional[Path] = None):
        """
        Create Venn diagram comparing significant results.
        
        Args:
            result_files: List of result file paths (2-3 files)
            sample_names: Names for each dataset
            q_cutoff: Q-value cutoff
            output_file: Output image file
        """
        if len(result_files) > 3:
            raise ValueError("Can only create Venn diagrams for 2-3 comparisons")
        
        # Read results and create binary matrix
        sig_matrix_data = {}
        
        for file_path, name in zip(result_files, sample_names):
            if file_path.suffix == '.csv':
                df = pl.read_csv(file_path)
            else:
                df = pl.read_csv(file_path, separator='\t')
            
            # Get q-value column
            q_col = 'q_val' if 'q_val' in df.columns else 'q_greater'
            
            if q_col not in df.columns:
                raise ValueError(f"No q-value column found in {file_path}")
            
            # Create binary significance indicator
            gene_set_col = 'gene_set' if 'gene_set' in df.columns else df.columns[0]
            sig_matrix_data[name] = df.select([
                gene_set_col,
                (pl.col(q_col) < q_cutoff).cast(pl.Int32).alias(name)
            ])
        
        # Merge all
        combined = sig_matrix_data[sample_names[0]]
        gene_set_col = combined.columns[0]
        
        for name in sample_names[1:]:
            combined = combined.join(
                sig_matrix_data[name],
                on=gene_set_col,
                how='outer'
            ).fill_null(0)
        
        # Create Venn diagram
        venn_data = combined.select(sample_names)
        venn = VennDiagram()
        counts = venn.venn_counts(venn_data, include='both')
        
        if len(sample_names) == 2:
            venn.plot_venn2(counts, sample_names, output_file)
        else:
            venn.plot_venn3(counts, sample_names, output_file)


class GeneSetGrouper:
    """Group overlapping essential gene sets."""
    
    @staticmethod
    def group_gene_sets(results: pl.DataFrame,
                       gene_sets: Dict[str, List[str]],
                       expression_data: pl.DataFrame,
                       gene_col: str = 'gene_id',
                       p_cutoff: float = 0.01,
                       overlap_cutoff: float = 1e-10,
                       output_file: Optional[Path] = None) -> Dict[str, List[str]]:
        """
        Group overlapping gene sets based on shared genes.
        
        Args:
            results: GAGE results DataFrame
            gene_sets: Dictionary of gene sets
            expression_data: Expression DataFrame
            gene_col: Gene ID column name
            p_cutoff: P-value cutoff for significance
            overlap_cutoff: P-value cutoff for overlap significance
            output_file: Optional output file
            
        Returns:
            Dictionary of gene set groups
        """
        from scipy.stats import hypergeom
        
        # Filter to significant gene sets
        p_col = 'p_val' if 'p_val' in results.columns else 'p_greater'
        sig_sets = results.filter(pl.col(p_col) < p_cutoff)
        
        if sig_sets.shape[0] < 2:
            print("Less than 2 significant gene sets, no grouping needed")
            return {}
        
        set_names = sig_sets['gene_set'].to_list()
        
        # Calculate pairwise overlaps
        n_sets = len(set_names)
        overlap_matrix = np.ones((n_sets, n_sets))
        
        # Total gene universe
        total_genes = expression_data.shape[0]
        
        for i in range(n_sets):
            for j in range(i + 1, n_sets):
                set_i = set(gene_sets[set_names[i]])
                set_j = set(gene_sets[set_names[j]])
                
                overlap = len(set_i & set_j)
                
                if overlap > 0:
                    # Hypergeometric test
                    p_overlap = hypergeom.sf(
                        overlap - 1,
                        total_genes,
                        len(set_i),
                        len(set_j)
                    )
                    
                    overlap_matrix[i, j] = p_overlap
                    overlap_matrix[j, i] = p_overlap
        
        # Create adjacency matrix (significant overlaps)
        adjacency = (overlap_matrix < overlap_cutoff).astype(int)
        np.fill_diagonal(adjacency, 0)
        
        # Find connected components (simple implementation)
        visited = set()
        groups = []
        
        def dfs(node, component):
            visited.add(node)
            component.append(node)
            for neighbor in range(n_sets):
                if adjacency[node, neighbor] and neighbor not in visited:
                    dfs(neighbor, component)
        
        for i in range(n_sets):
            if i not in visited:
                component = []
                dfs(i, component)
                groups.append([set_names[idx] for idx in component])
        
        # Create output
        group_dict = {}
        for i, group in enumerate(groups):
            group_dict[f"Group_{i+1}"] = group
        
        if output_file:
            with open(output_file, 'w') as f:
                json.dump(group_dict, f, indent=2)
            print(f"Gene set groups written to {output_file}")
        
        return group_dict


class SignificanceFilter:
    """Filter significant gene sets from GAGE results."""
    
    @staticmethod
    def filter_significant(results: Dict[str, pl.DataFrame],
                          cutoff: float = 0.1,
                          use_q: bool = True,
                          dual_sig: int = 2) -> Dict[str, pl.DataFrame]:
        """
        Filter significant gene sets.
        
        Args:
            results: Dictionary with 'greater' and optionally 'less' DataFrames
            cutoff: Significance cutoff
            use_q: Use q-value (True) or p-value (False)
            dual_sig: Dual significance mode:
                      0 - Exclusive (significant in one direction only)
                      1 - Prefer better p-value
                      2 - Both directions allowed
            
        Returns:
            Filtered results dictionary
        """
        col_name = 'q_val' if use_q else 'p_val'
        if col_name not in results['greater'].columns:
            col_name = 'q_greater' if use_q else 'p_greater'
        
        filtered = {}
        
        if 'less' in results and results['less'] is not None:
            col_name_less = col_name.replace('greater', 'less')
            
            greater_sig = results['greater'].filter(pl.col(col_name) < cutoff)
            less_sig = results['less'].filter(pl.col(col_name_less) < cutoff)
            
            if dual_sig == 0:
                # Exclusive - remove genes significant in both
                greater_only = greater_sig.filter(
                    ~pl.col('gene_set').is_in(less_sig['gene_set'])
                )
                less_only = less_sig.filter(
                    ~pl.col('gene_set').is_in(greater_sig['gene_set'])
                )
                filtered['greater'] = greater_only
                filtered['less'] = less_only
            
            elif dual_sig == 1:
                # Prefer better p-value
                all_sig = pl.concat([
                    greater_sig.with_columns(pl.lit('greater').alias('direction')),
                    less_sig.with_columns(pl.lit('less').alias('direction'))
                ])
                
                # Keep best p-value for each gene set
                best = all_sig.sort('gene_set', col_name).unique(subset=['gene_set'], keep='first')
                
                filtered['greater'] = best.filter(pl.col('direction') == 'greater')
                filtered['less'] = best.filter(pl.col('direction') == 'less')
            
            else:  # dual_sig == 2
                # Both allowed
                filtered['greater'] = greater_sig
                filtered['less'] = less_sig
        
        else:
            # Only greater direction
            filtered['greater'] = results['greater'].filter(pl.col(col_name) < cutoff)
        
        # Copy stats
        if 'stats' in results:
            filtered['stats'] = results['stats']
        
        return filtered


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
