#!/usr/bin/env python3
"""
Create Gene Count Histogram for Stratification Analysis

This script creates a histogram of gene counts per disease to help determine
appropriate breakpoints for stratified precision@k and recall@k analysis.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def load_gene_counts():
    """Load disease gene counts from TSV file"""
    print("Loading disease gene counts...")
    df = pd.read_csv('disease_gene_counts.tsv', sep='\t')
    print(f"Loaded {len(df)} diseases")

    # Basic statistics
    total_diseases = len(df)
    diseases_with_genes = len(df[df['gene_count'] > 0])
    diseases_with_3plus_genes = len(df[df['gene_count'] >= 3])

    print(f"Total diseases: {total_diseases:,}")
    print(f"Diseases with genes: {diseases_with_genes:,} ({diseases_with_genes/total_diseases*100:.1f}%)")
    print(f"Diseases with ≥3 genes: {diseases_with_3plus_genes:,} ({diseases_with_3plus_genes/total_diseases*100:.1f}%)")

    return df

def create_gene_count_histogram(df):
    """Create histogram of gene counts"""

    # Filter to diseases with genes
    df_with_genes = df[df['gene_count'] > 0]
    gene_counts = df_with_genes['gene_count']

    # Create figure with multiple subplots for different views
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))

    # Plot 1: Full distribution (log scale)
    ax = axes[0, 0]
    ax.hist(gene_counts, bins=100, alpha=0.7, color='blue', edgecolor='black')
    ax.set_xlabel('Gene Count')
    ax.set_ylabel('Number of Diseases')
    ax.set_title('Distribution of Gene Counts per Disease (All)')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)

    # Add statistics text
    stats_text = f'Mean: {gene_counts.mean():.1f}\nMedian: {gene_counts.median():.0f}\nMax: {gene_counts.max()}'
    ax.text(0.7, 0.8, stats_text, transform=ax.transAxes,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    # Plot 2: Zoomed in (0-100 genes)
    ax = axes[0, 1]
    ax.hist(gene_counts[gene_counts <= 100], bins=50, alpha=0.7, color='green', edgecolor='black')
    ax.set_xlabel('Gene Count')
    ax.set_ylabel('Number of Diseases')
    ax.set_title('Distribution of Gene Counts (0-100 genes)')
    ax.grid(True, alpha=0.3)

    # Plot 3: Diseases with ≥3 genes only
    ax = axes[1, 0]
    eligible_counts = gene_counts[gene_counts >= 3]
    ax.hist(eligible_counts, bins=50, alpha=0.7, color='red', edgecolor='black')
    ax.set_xlabel('Gene Count')
    ax.set_ylabel('Number of Diseases')
    ax.set_title('Distribution of Gene Counts (≥3 genes only)')
    ax.grid(True, alpha=0.3)

    # Add percentiles for stratification
    percentiles = [25, 50, 75, 90, 95]
    perc_values = np.percentile(eligible_counts, percentiles)
    for p, val in zip(percentiles, perc_values):
        ax.axvline(val, color='orange', linestyle='--', alpha=0.8,
                  label=f'{p}th percentile: {val:.0f}')
    ax.legend()

    # Plot 4: Cumulative distribution
    ax = axes[1, 1]
    sorted_counts = np.sort(eligible_counts)
    cumulative_pct = np.arange(1, len(sorted_counts) + 1) / len(sorted_counts) * 100
    ax.plot(sorted_counts, cumulative_pct, color='purple', linewidth=2)
    ax.set_xlabel('Gene Count')
    ax.set_ylabel('Cumulative Percentage')
    ax.set_title('Cumulative Distribution (≥3 genes)')
    ax.grid(True, alpha=0.3)

    # Add suggested breakpoints
    suggested_breakpoints = [3, 10, 25, 50, 100]
    for bp in suggested_breakpoints:
        if bp <= sorted_counts.max():
            pct = np.searchsorted(sorted_counts, bp) / len(sorted_counts) * 100
            ax.axvline(bp, color='red', linestyle='--', alpha=0.7)
            ax.text(bp, pct + 5, f'{bp} genes\n({pct:.1f}%)',
                   ha='center', fontsize=9,
                   bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))

    plt.tight_layout()
    plt.savefig('gene_count_histogram.png', dpi=300, bbox_inches='tight')
    print("Gene count histogram saved to gene_count_histogram.png")
    plt.show()

    return eligible_counts

def suggest_stratification_breakpoints(gene_counts):
    """Suggest breakpoints for stratification based on distribution"""
    print("\n=== Stratification Analysis ===")

    # Calculate key statistics
    percentiles = [10, 25, 50, 75, 90, 95]
    perc_values = np.percentile(gene_counts, percentiles)

    print("Gene count percentiles (diseases with ≥3 genes):")
    for p, val in zip(percentiles, perc_values):
        print(f"  {p:2d}th percentile: {val:5.0f} genes")

    # Suggest stratification groups
    print("\nSuggested stratification groups:")

    # Option 1: Quartile-based
    breakpoints_1 = [3, int(perc_values[1]), int(perc_values[2]), int(perc_values[3]), float('inf')]
    print(f"Option 1 (Quartile-based): {breakpoints_1[:-1]} genes (last group: >{breakpoints_1[-2]})")

    # Option 2: Log-scale inspired
    breakpoints_2 = [3, 10, 25, 100, float('inf')]
    print(f"Option 2 (Log-scale): {breakpoints_2[:-1]} genes (last group: >{breakpoints_2[-2]})")

    # Option 3: Even distribution
    n_groups = 5
    group_size = len(gene_counts) // n_groups
    breakpoints_3 = []
    for i in range(n_groups):
        idx = min(i * group_size, len(gene_counts) - 1)
        breakpoints_3.append(int(np.sort(gene_counts)[idx]))
    breakpoints_3.append(float('inf'))
    print(f"Option 3 (Even distribution): {breakpoints_3[:-1]} genes (last group: >{breakpoints_3[-2]})")

    # Show group sizes for each option
    print("\nGroup sizes for each option:")
    for opt_num, breakpoints in enumerate([breakpoints_1, breakpoints_2, breakpoints_3], 1):
        print(f"Option {opt_num}:")
        for i in range(len(breakpoints) - 1):
            min_genes = breakpoints[i]
            max_genes = breakpoints[i + 1]
            if max_genes == float('inf'):
                count = len(gene_counts[gene_counts >= min_genes])
                print(f"  Group {i+1}: ≥{min_genes} genes: {count} diseases ({count/len(gene_counts)*100:.1f}%)")
            else:
                count = len(gene_counts[(gene_counts >= min_genes) & (gene_counts < max_genes)])
                print(f"  Group {i+1}: {min_genes}-{max_genes-1} genes: {count} diseases ({count/len(gene_counts)*100:.1f}%)")

    return breakpoints_2  # Return log-scale option as default

def main():
    """Main function"""
    df = load_gene_counts()
    eligible_counts = create_gene_count_histogram(df)
    suggested_breakpoints = suggest_stratification_breakpoints(eligible_counts)

    print(f"\nRecommended breakpoints: {suggested_breakpoints[:-1]} (last group: >{suggested_breakpoints[-2]})")
    print("These breakpoints provide reasonable group sizes for stratified analysis.")

if __name__ == "__main__":
    main()