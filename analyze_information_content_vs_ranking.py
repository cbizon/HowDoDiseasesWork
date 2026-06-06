#!/usr/bin/env python3

import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from collections import defaultdict
import argparse

class InformationContentAnalyzer:
    def __init__(self):
        self.information_content = {}

    def extract_information_content(self, nodes_file: str = '/Users/bizon/Projects/ROBOKOP/graph/nodes.jsonl'):
        """Extract information content values from ROBOKOP nodes file"""
        cache_file = 'robokop_information_content_cache.json'

        try:
            print("Loading information content from cache...")
            with open(cache_file, 'r') as f:
                self.information_content = json.load(f)
            print(f"Loaded {len(self.information_content):,} terms with information content from cache")
        except FileNotFoundError:
            print("Cache not found. Extracting information content from nodes file...")
            print(f"Processing {nodes_file}...")

            processed = 0
            found_ic = 0

            with open(nodes_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    if line_num % 1000000 == 0:
                        print(f"  Processed {line_num:,} nodes, found {found_ic:,} with information_content")

                    try:
                        node = json.loads(line.strip())
                        curie = node.get('id', '')

                        if curie and 'information_content' in node:
                            ic_value = node['information_content']
                            if isinstance(ic_value, (int, float)) and not np.isnan(ic_value):
                                self.information_content[curie] = ic_value
                                found_ic += 1

                        processed += 1

                    except json.JSONDecodeError:
                        continue

            print(f"Extracted information content for {found_ic:,} terms from {processed:,} total nodes")

            # Cache the results
            with open(cache_file, 'w') as f:
                json.dump(self.information_content, f)
            print(f"Cached information content to {cache_file}")

    def load_ranking_analysis(self, ranking_file: str = 'biological_process_ranking_analysis.tsv'):
        """Load the biological process ranking analysis"""
        print(f"Loading ranking analysis from {ranking_file}...")
        df = pd.read_csv(ranking_file, sep='\t')
        print(f"Loaded {len(df):,} disease-term pairs")
        return df

    def analyze_found_vs_not_found(self, df):
        """Compare information content distributions for found vs not found terms"""
        print("\nAnalyzing information content: Found vs Not Found terms...")

        # Add information content to dataframe
        df['Information_Content'] = df['ROBOKOP_Process_CURIE'].map(self.information_content)

        # Separate found vs not found
        found_mask = df['Enrichment_Rank'] != 'Not_Found'
        found_df = df[found_mask & df['Information_Content'].notna()]
        not_found_df = df[~found_mask & df['Information_Content'].notna()]

        print(f"Terms with information content:")
        print(f"  Found in enrichment: {len(found_df):,}")
        print(f"  Not found in enrichment: {len(not_found_df):,}")
        print(f"  Missing information content: {df['Information_Content'].isna().sum():,}")

        if len(found_df) == 0 or len(not_found_df) == 0:
            print("Cannot perform comparison - one group is empty")
            return df

        # Statistical test
        found_ic = found_df['Information_Content'].values
        not_found_ic = not_found_df['Information_Content'].values

        statistic, p_value = stats.mannwhitneyu(found_ic, not_found_ic, alternative='two-sided')

        print(f"\nInformation Content Statistics:")
        print(f"  Found - Mean: {np.mean(found_ic):.3f}, Median: {np.median(found_ic):.3f}, Std: {np.std(found_ic):.3f}")
        print(f"  Not Found - Mean: {np.mean(not_found_ic):.3f}, Median: {np.median(not_found_ic):.3f}, Std: {np.std(not_found_ic):.3f}")
        print(f"  Mann-Whitney U test p-value: {p_value:.2e}")

        # Create histograms
        plt.figure(figsize=(12, 6))

        plt.subplot(1, 2, 1)
        plt.hist(found_ic, bins=50, alpha=0.7, label=f'Found (n={len(found_ic):,})', color='green')
        plt.hist(not_found_ic, bins=50, alpha=0.7, label=f'Not Found (n={len(not_found_ic):,})', color='red')
        plt.xlabel('Information Content')
        plt.ylabel('Frequency')
        plt.title('Information Content Distribution')
        plt.legend()
        plt.grid(True, alpha=0.3)

        plt.subplot(1, 2, 2)
        plt.boxplot([found_ic, not_found_ic], tick_labels=['Found', 'Not Found'])
        plt.ylabel('Information Content')
        plt.title(f'Information Content Box Plot\np-value: {p_value:.2e}')
        plt.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig('information_content_found_vs_not_found.png', dpi=300, bbox_inches='tight')
        plt.show()

        return df

    def analyze_ic_vs_rank(self, df):
        """Analyze information content vs enrichment rank for found terms"""
        print("\nAnalyzing Information Content vs Enrichment Rank...")

        # Filter to found terms with numeric ranks and information content
        found_mask = df['Enrichment_Rank'] != 'Not_Found'
        df_found = df[found_mask & df['Information_Content'].notna()].copy()

        # Convert rank to numeric
        df_found['Enrichment_Rank_Numeric'] = pd.to_numeric(df_found['Enrichment_Rank'], errors='coerce')
        df_found = df_found[df_found['Enrichment_Rank_Numeric'].notna()]

        print(f"Terms for rank analysis: {len(df_found):,}")

        if len(df_found) == 0:
            print("No terms available for rank analysis")
            return

        ic_values = df_found['Information_Content'].values
        ranks = df_found['Enrichment_Rank_Numeric'].values

        # Calculate correlation
        correlation, p_value = stats.spearmanr(ic_values, ranks)

        print(f"Spearman correlation between Information Content and Rank: {correlation:.4f} (p={p_value:.2e})")

        # Create scatter plot
        plt.figure(figsize=(12, 8))

        plt.subplot(2, 2, 1)
        plt.scatter(ic_values, ranks, alpha=0.6, s=20)
        plt.xlabel('Information Content')
        plt.ylabel('Enrichment Rank')
        plt.title(f'Information Content vs Enrichment Rank\nSpearman r = {correlation:.4f}, p = {p_value:.2e}')
        plt.grid(True, alpha=0.3)

        # Log scale version
        plt.subplot(2, 2, 2)
        plt.scatter(ic_values, ranks, alpha=0.6, s=20)
        plt.xlabel('Information Content')
        plt.ylabel('Enrichment Rank (log scale)')
        plt.yscale('log')
        plt.title('Information Content vs Enrichment Rank (Log Scale)')
        plt.grid(True, alpha=0.3)

        # Binned analysis
        plt.subplot(2, 2, 3)
        # Create IC bins
        n_bins = 10
        ic_bins = np.percentile(ic_values, np.linspace(0, 100, n_bins + 1))
        # Remove duplicate bin edges
        ic_bins = np.unique(ic_bins)
        if len(ic_bins) < 2:
            print("Not enough unique IC values for binning")
            return
        df_found['IC_Bin'] = pd.cut(df_found['Information_Content'], ic_bins, include_lowest=True)

        bin_stats = df_found.groupby('IC_Bin')['Enrichment_Rank_Numeric'].agg(['mean', 'median', 'count'])
        bin_centers = [interval.mid for interval in bin_stats.index]

        plt.plot(bin_centers, bin_stats['mean'], 'o-', label='Mean Rank')
        plt.plot(bin_centers, bin_stats['median'], 's-', label='Median Rank')
        plt.xlabel('Information Content (Bin Centers)')
        plt.ylabel('Average Enrichment Rank')
        plt.title('Binned Analysis: IC vs Average Rank')
        plt.legend()
        plt.grid(True, alpha=0.3)

        # Distribution of ranks by IC quartiles
        plt.subplot(2, 2, 4)
        ic_quartiles = pd.qcut(df_found['Information_Content'], 4, labels=['Q1', 'Q2', 'Q3', 'Q4'])
        quartile_data = [df_found[ic_quartiles == q]['Enrichment_Rank_Numeric'].values for q in ['Q1', 'Q2', 'Q3', 'Q4']]

        plt.boxplot(quartile_data, tick_labels=['Q1\n(Low IC)', 'Q2', 'Q3', 'Q4\n(High IC)'])
        plt.ylabel('Enrichment Rank')
        plt.title('Rank Distribution by IC Quartiles')
        plt.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig('information_content_vs_rank_analysis.png', dpi=300, bbox_inches='tight')
        plt.show()

        # Print quartile statistics
        print(f"\nRank statistics by Information Content quartiles:")
        for q in ['Q1', 'Q2', 'Q3', 'Q4']:
            q_data = df_found[ic_quartiles == q]
            q_ranks = q_data['Enrichment_Rank_Numeric']
            print(f"  {q}: Mean rank = {q_ranks.mean():.1f}, Median rank = {q_ranks.median():.1f}, n = {len(q_ranks):,}")

def main():
    parser = argparse.ArgumentParser(description='Analyze information content vs enrichment ranking')
    parser.add_argument('--nodes-file', default='/Users/bizon/Projects/ROBOKOP/graph/nodes.jsonl',
                        help='Path to ROBOKOP nodes file')
    parser.add_argument('--ranking-file', default='biological_process_ranking_analysis.tsv',
                        help='Path to ranking analysis file')

    args = parser.parse_args()

    analyzer = InformationContentAnalyzer()

    # Extract information content from nodes
    analyzer.extract_information_content(args.nodes_file)

    # Load ranking analysis
    df = analyzer.load_ranking_analysis(args.ranking_file)

    # Analyze found vs not found
    df = analyzer.analyze_found_vs_not_found(df)

    # Analyze IC vs rank for found terms
    analyzer.analyze_ic_vs_rank(df)

    print("\nAnalysis complete!")

if __name__ == "__main__":
    main()