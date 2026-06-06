#!/usr/bin/env python3

import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LogNorm
from collections import defaultdict
import argparse

def load_enrichment_results(jsonl_file='fast_enrichment_results.jsonl'):
    """Load enrichment results from JSONL file"""
    print(f"Loading enrichment results from {jsonl_file}...")
    results = []

    with open(jsonl_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line_num % 1000 == 0:
                print(f"  Processed {line_num:,} diseases")

            try:
                result = json.loads(line.strip())
                results.append(result)
            except json.JSONDecodeError:
                continue

    print(f"Loaded enrichment results for {len(results):,} diseases")
    return results

def extract_counts_by_pvalue(enrichment_results, p_value_cutoff=1e-5):
    """Extract gene and term counts for each disease and category at given p-value cutoff"""
    print(f"Extracting counts with p-value cutoff: {p_value_cutoff:.0e}")

    data = []
    categories = ['biolink:BiologicalProcess', 'biolink:MolecularActivity', 'biolink:Pathway']

    for result in enrichment_results:
        disease_curie = result['disease']['curie']
        disease_name = result['disease']['name']
        gene_count = result['gene_count']

        for category in categories:
            if category in result['enrichment_results']:
                # Count terms below p-value cutoff
                terms_below_cutoff = [
                    term for term in result['enrichment_results'][category]
                    if term['p_value'] <= p_value_cutoff
                ]
                term_count = len(terms_below_cutoff)

                # Only include diseases with at least one significant term
                if term_count > 0:
                    data.append({
                        'disease_curie': disease_curie,
                        'disease_name': disease_name,
                        'gene_count': gene_count,
                        'term_count': term_count,
                        'category': category.replace('biolink:', ''),
                        'p_value_cutoff': p_value_cutoff
                    })

    df = pd.DataFrame(data)
    print(f"Found {len(df):,} disease-category pairs with terms below p={p_value_cutoff:.0e}")

    # Print summary by category
    for category in df['category'].unique():
        cat_df = df[df['category'] == category]
        print(f"  {category}: {len(cat_df):,} diseases")

    return df

def create_2d_histogram(df, p_value_cutoff=1e-5, output_file=None):
    """Create 2D histograms similar to the ROBOKOP version but for enrichment results"""

    if len(df) == 0:
        print("No data to plot!")
        return

    categories = df['category'].unique()
    n_cats = len(categories)

    # Calculate correlations for each category
    correlations = {}
    for category in categories:
        cat_df = df[df['category'] == category]
        if len(cat_df) > 1:
            correlation = np.corrcoef(cat_df['gene_count'], cat_df['term_count'])[0, 1]
            correlations[category] = correlation
        else:
            correlations[category] = np.nan

    # Set up the plot
    if n_cats == 4:  # Include "All Combined"
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        axes = axes.flatten()
    else:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        axes = axes.flatten()

    # Color maps for each category
    cmaps = ['Blues', 'Reds', 'Greens', 'Purples']

    plot_idx = 0

    # Plot individual categories
    for i, category in enumerate(sorted(categories)):
        cat_df = df[df['category'] == category]

        if len(cat_df) == 0:
            continue

        ax = axes[plot_idx]

        # Create 2D histogram
        gene_counts = cat_df['gene_count'].values
        term_counts = cat_df['term_count'].values

        # Determine bins
        gene_max = max(gene_counts) + 5
        term_max = max(term_counts) + 5

        gene_bins = np.arange(0, gene_max, max(1, gene_max // 50))
        term_bins = np.arange(0, term_max, max(1, term_max // 50))

        # Create histogram
        hist, gene_edges, term_edges = np.histogram2d(
            gene_counts, term_counts, bins=[gene_bins, term_bins]
        )

        # Plot with logarithmic color scale
        # Add small value to avoid log(0)
        hist_log = hist.T + 1e-10
        im = ax.imshow(
            hist_log,
            origin='lower',
            extent=[gene_edges[0], gene_edges[-1], term_edges[0], term_edges[-1]],
            cmap=cmaps[i % len(cmaps)],
            aspect='auto',
            norm=LogNorm(vmin=hist_log[hist_log > 0].min(), vmax=hist_log.max())
        )

        corr = correlations.get(category, np.nan)
        ax.set_title(f'{category}, r = {corr:.3f}')
        ax.set_xlabel('Gene Count')
        ax.set_ylabel('Term Count')
        ax.grid(True, alpha=0.3)

        # Add colorbar
        plt.colorbar(im, ax=ax, label='Count')

        plot_idx += 1

    # Create "All Combined" plot if we have multiple categories
    if len(categories) > 1 and plot_idx < 4:
        ax = axes[plot_idx]

        # Combine all data
        all_gene_counts = df['gene_count'].values
        all_term_counts = df['term_count'].values

        # Create 2D histogram
        gene_max = max(all_gene_counts) + 5
        term_max = max(all_term_counts) + 5

        gene_bins = np.arange(0, gene_max, max(1, gene_max // 50))
        term_bins = np.arange(0, term_max, max(1, term_max // 50))

        hist, gene_edges, term_edges = np.histogram2d(
            all_gene_counts, all_term_counts, bins=[gene_bins, term_bins]
        )

        # Plot with logarithmic color scale
        # Add small value to avoid log(0)
        hist_log = hist.T + 1e-10
        im = ax.imshow(
            hist_log,
            origin='lower',
            extent=[gene_edges[0], gene_edges[-1], term_edges[0], term_edges[-1]],
            cmap='viridis',
            aspect='auto',
            norm=LogNorm(vmin=hist_log[hist_log > 0].min(), vmax=hist_log.max())
        )

        all_corr = np.corrcoef(all_gene_counts, all_term_counts)[0, 1]
        ax.set_title(f'All Categories Combined, r = {all_corr:.3f}')
        ax.set_xlabel('Gene Count')
        ax.set_ylabel('Term Count')
        ax.grid(True, alpha=0.3)

        # Add colorbar
        plt.colorbar(im, ax=ax, label='Count')

        plot_idx += 1

    # Hide unused subplots
    for i in range(plot_idx, 4):
        axes[i].set_visible(False)

    plt.suptitle(f'Gene Count vs Term Count: Enrichment Results (p ≤ {p_value_cutoff:.0e})',
                 fontsize=16, y=0.98)
    plt.tight_layout()

    # Save plot
    if output_file is None:
        output_file = f'enrichment_gene_vs_term_counts_p{p_value_cutoff:.0e}.png'

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.show()

    print(f"Plot saved as {output_file}")

    return df

def print_summary_stats(df, p_value_cutoff):
    """Print summary statistics"""
    print(f"\nSummary Statistics (p ≤ {p_value_cutoff:.0e}):")
    print("=" * 50)

    overall_stats = df.groupby('category').agg({
        'gene_count': ['count', 'mean', 'std', 'min', 'max'],
        'term_count': ['mean', 'std', 'min', 'max']
    }).round(1)

    print(overall_stats)

    # Correlation summary
    print(f"\nCorrelations by Category:")
    print("-" * 30)
    for category in df['category'].unique():
        cat_df = df[df['category'] == category]
        if len(cat_df) > 1:
            correlation = np.corrcoef(cat_df['gene_count'], cat_df['term_count'])[0, 1]
            print(f"{category}: {correlation:.3f}")

def main():
    parser = argparse.ArgumentParser(description='Create gene vs term count 2D histograms from enrichment results')
    parser.add_argument('--input', default='fast_enrichment_results.jsonl',
                        help='Input JSONL file with enrichment results')
    parser.add_argument('--p-cutoff', type=float, default=1e-5,
                        help='P-value cutoff for including terms (default: 1e-5)')
    parser.add_argument('--output',
                        help='Output plot filename (auto-generated if not provided)')

    args = parser.parse_args()

    # Load enrichment results
    enrichment_results = load_enrichment_results(args.input)

    # Extract counts with p-value cutoff
    df = extract_counts_by_pvalue(enrichment_results, args.p_cutoff)

    if len(df) == 0:
        print(f"No significant terms found at p ≤ {args.p_cutoff:.0e}")
        return

    # Create 2D histogram plot
    create_2d_histogram(df, args.p_cutoff, args.output)

    # Print summary statistics
    print_summary_stats(df, args.p_cutoff)

if __name__ == "__main__":
    main()