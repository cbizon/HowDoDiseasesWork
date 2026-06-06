#!/usr/bin/env python3
"""
Analyze Gene Count vs Term Count Relationship for Diseases

This script creates 2D histograms showing the relationship between:
- Number of genes connected to a disease (from our enrichment analysis)
- Number of terms (BiologicalProcess/MolecularActivity/Pathway) connected to the disease (from ROBOKOP)

This helps us understand if diseases with more genes tend to have more known mechanisms.
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from typing import Dict, List, Tuple
import argparse

class GeneTermAnalyzer:
    """Class to analyze gene count vs term count relationships"""

    def __init__(self):
        """Initialize the analyzer"""
        self.enrichment_data = []
        self.robokop_term_counts = defaultdict(lambda: defaultdict(int))  # disease_curie -> category -> count
        self.categories = ['biolink:BiologicalProcess', 'biolink:MolecularActivity', 'biolink:Pathway']

    def load_enrichment_results(self, jsonl_file: str) -> None:
        """
        Load enrichment analysis results to get gene counts per disease

        Args:
            jsonl_file: Path to enrichment results JSONL file
        """
        print(f"Loading enrichment results from {jsonl_file}...")

        with open(jsonl_file, 'r') as f:
            for line in f:
                try:
                    result = json.loads(line.strip())
                    self.enrichment_data.append(result)
                except json.JSONDecodeError as e:
                    continue

        print(f"Loaded {len(self.enrichment_data)} disease enrichment results")

    def load_robokop_term_counts(self, jsonl_file: str) -> None:
        """
        Load ROBOKOP disease-term edges and count terms per disease by category

        Args:
            jsonl_file: Path to ROBOKOP disease-term edges JSONL file
        """
        print(f"Loading ROBOKOP term counts from {jsonl_file}...")

        with open(jsonl_file, 'r') as f:
            for line in f:
                try:
                    edge = json.loads(line.strip())
                    disease_curie = edge['source_curie']
                    term_category = edge['target_type']

                    # Count terms by category for each disease
                    if term_category in self.categories:
                        self.robokop_term_counts[disease_curie][term_category] += 1
                except (json.JSONDecodeError, KeyError) as e:
                    continue

        # Calculate totals
        total_diseases = len(self.robokop_term_counts)
        category_stats = {}

        for category in self.categories:
            diseases_with_terms = sum(1 for disease in self.robokop_term_counts
                                    if self.robokop_term_counts[disease][category] > 0)
            total_terms = sum(self.robokop_term_counts[disease][category]
                            for disease in self.robokop_term_counts)
            category_stats[category] = (diseases_with_terms, total_terms)

        print(f"Loaded term counts for {total_diseases} diseases from ROBOKOP")
        for category in self.categories:
            diseases, terms = category_stats[category]
            print(f"  {category.replace('biolink:', '')}: {diseases} diseases, {terms} terms")

    def prepare_data_for_plotting(self) -> Tuple[Dict[str, List[int]], Dict[str, List[int]]]:
        """
        Prepare gene counts and term counts for plotting

        Returns:
            Tuple of (gene_counts_by_category, term_counts_by_category)
        """
        print("Preparing data for plotting...")

        gene_counts = {cat: [] for cat in self.categories}
        term_counts = {cat: [] for cat in self.categories}

        # Find diseases that appear in both datasets
        matched_diseases = 0

        for result in self.enrichment_data:
            disease_curie = result['disease']['curie']
            gene_count = result.get('gene_count', 0)

            # Check if this disease has ROBOKOP term data
            if disease_curie in self.robokop_term_counts:
                matched_diseases += 1

                # For each category, add gene count and corresponding term count
                for category in self.categories:
                    term_count = self.robokop_term_counts[disease_curie][category]

                    # Only include if there are terms in this category or genes
                    if term_count > 0 or gene_count > 0:
                        gene_counts[category].append(gene_count)
                        term_counts[category].append(term_count)

        print(f"Found {matched_diseases} diseases in both datasets")

        # Print statistics
        for category in self.categories:
            count = len(gene_counts[category])
            if count > 0:
                avg_genes = np.mean(gene_counts[category])
                avg_terms = np.mean(term_counts[category])
                print(f"  {category.replace('biolink:', '')}: {count} data points, "
                      f"avg genes: {avg_genes:.1f}, avg terms: {avg_terms:.1f}")

        return gene_counts, term_counts

    def create_2d_histograms(self, gene_counts: Dict[str, List[int]],
                           term_counts: Dict[str, List[int]],
                           output_file: str = 'gene_vs_term_counts_2d_hist.png') -> None:
        """
        Create 2D histograms showing gene count vs term count relationship

        Args:
            gene_counts: Dictionary mapping category to list of gene counts
            term_counts: Dictionary mapping category to list of term counts
            output_file: Output filename for plot
        """
        print(f"Creating 2D histograms...")

        # Set up the subplot layout
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Gene Count vs Term Count: 2D Histograms by Category', fontsize=16)

        # Color maps for each category
        cmaps = ['Blues', 'Reds', 'Greens']

        # Individual category plots
        for i, category in enumerate(self.categories):
            ax = axes[i//2, i%2]

            if len(gene_counts[category]) > 0:
                # Create 2D histogram
                h = ax.hist2d(gene_counts[category], term_counts[category],
                            bins=30, cmap=cmaps[i], alpha=0.8)

                # Add colorbar
                plt.colorbar(h[3], ax=ax, label='Count')

                # Calculate and display correlation
                correlation = np.corrcoef(gene_counts[category], term_counts[category])[0, 1]

                ax.set_xlabel('Gene Count', fontsize=12)
                ax.set_ylabel('Term Count', fontsize=12)
                ax.set_title(f'{category.replace("biolink:", "")}\\nr = {correlation:.3f}', fontsize=12)
                ax.grid(True, alpha=0.3)
            else:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(category.replace('biolink:', ''), fontsize=12)

        # Combined plot in bottom right
        ax = axes[1, 1]

        # Combine all data points
        all_genes = []
        all_terms = []
        for category in self.categories:
            all_genes.extend(gene_counts[category])
            all_terms.extend(term_counts[category])

        if all_genes:
            h = ax.hist2d(all_genes, all_terms, bins=40, cmap='viridis', alpha=0.8)
            plt.colorbar(h[3], ax=ax, label='Count')

            # Calculate overall correlation
            overall_correlation = np.corrcoef(all_genes, all_terms)[0, 1]

            ax.set_xlabel('Gene Count', fontsize=12)
            ax.set_ylabel('Term Count', fontsize=12)
            ax.set_title(f'All Categories Combined\\nr = {overall_correlation:.3f}', fontsize=12)
            ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"2D histograms saved to {output_file}")
        plt.show()

    def create_summary_statistics(self, gene_counts: Dict[str, List[int]],
                                term_counts: Dict[str, List[int]],
                                output_file: str = 'gene_term_correlation_stats.tsv') -> None:
        """
        Create summary statistics table

        Args:
            gene_counts: Dictionary mapping category to list of gene counts
            term_counts: Dictionary mapping category to list of term counts
            output_file: Output TSV filename
        """
        print(f"Creating summary statistics...")

        with open(output_file, 'w') as f:
            f.write("Category\\tData_Points\\tGenes_Mean\\tGenes_Std\\tTerms_Mean\\tTerms_Std\\tCorrelation\\n")

            for category in self.categories:
                if len(gene_counts[category]) > 0:
                    genes = np.array(gene_counts[category])
                    terms = np.array(term_counts[category])

                    correlation = np.corrcoef(genes, terms)[0, 1]

                    f.write(f"{category.replace('biolink:', '')}\\t"
                           f"{len(genes)}\\t"
                           f"{genes.mean():.1f}\\t"
                           f"{genes.std():.1f}\\t"
                           f"{terms.mean():.1f}\\t"
                           f"{terms.std():.1f}\\t"
                           f"{correlation:.3f}\\n")
                else:
                    f.write(f"{category.replace('biolink:', '')}\\t0\\t0\\t0\\t0\\t0\\tN/A\\n")

            # Overall statistics
            all_genes = []
            all_terms = []
            for category in self.categories:
                all_genes.extend(gene_counts[category])
                all_terms.extend(term_counts[category])

            if all_genes:
                genes = np.array(all_genes)
                terms = np.array(all_terms)
                correlation = np.corrcoef(genes, terms)[0, 1]

                f.write(f"All_Combined\\t"
                       f"{len(genes)}\\t"
                       f"{genes.mean():.1f}\\t"
                       f"{genes.std():.1f}\\t"
                       f"{terms.mean():.1f}\\t"
                       f"{terms.std():.1f}\\t"
                       f"{correlation:.3f}\\n")

        print(f"Summary statistics saved to {output_file}")

    def run_analysis(self) -> None:
        """Run complete gene vs term count analysis"""
        # Prepare data
        gene_counts, term_counts = self.prepare_data_for_plotting()

        # Create visualizations
        self.create_2d_histograms(gene_counts, term_counts)

        # Create summary statistics
        self.create_summary_statistics(gene_counts, term_counts)

def main():
    """Main function to run gene vs term count analysis"""

    parser = argparse.ArgumentParser(description='Analyze relationship between gene counts and term counts for diseases')
    parser.add_argument('--enrichment', '-e', default='fast_enrichment_results.jsonl',
                       help='Enrichment results JSONL file (default: fast_enrichment_results.jsonl)')
    parser.add_argument('--robokop', '-r', default='robokop_disease_term_edges.jsonl',
                       help='ROBOKOP disease-term edges JSONL file (default: robokop_disease_term_edges.jsonl)')
    parser.add_argument('--output-plot', '-o', default='gene_vs_term_counts_2d_hist.png',
                       help='Output plot filename (default: gene_vs_term_counts_2d_hist.png)')
    parser.add_argument('--output-stats', '-s', default='gene_term_correlation_stats.tsv',
                       help='Output statistics filename (default: gene_term_correlation_stats.tsv)')

    args = parser.parse_args()

    # Initialize analyzer
    analyzer = GeneTermAnalyzer()

    # Load data
    analyzer.load_enrichment_results(args.enrichment)
    analyzer.load_robokop_term_counts(args.robokop)

    # Run analysis
    analyzer.run_analysis()

    print(f"\\nAnalysis complete!")
    print(f"2D histogram saved to: {args.output_plot}")
    print(f"Statistics saved to: {args.output_stats}")

if __name__ == "__main__":
    main()