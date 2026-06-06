#!/usr/bin/env python3
"""
Stratified Precision@K and Recall@K Analysis

This script performs stratified analysis of enrichment results by gene count,
using tertile breakpoints (33rd and 66th percentiles) to create three equal-sized groups.
"""

import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from collections import defaultdict
from typing import Dict, List, Tuple
import argparse

class StratifiedAnalyzer:
    """Stratified precision/recall analyzer"""

    def __init__(self):
        """Initialize analyzer"""
        self.enrichment_data = []
        self.ground_truth = defaultdict(lambda: defaultdict(set))
        self.disease_gene_counts = {}
        self.categories = ['biolink:BiologicalProcess', 'biolink:MolecularActivity', 'biolink:Pathway']

        # Tertile breakpoints (calculated from gene count distribution)
        self.stratification_groups = {
            'Low (3-5 genes)': (3, 5),
            'Medium (6-19 genes)': (6, 19),
            'High (≥20 genes)': (20, float('inf'))
        }

    def load_disease_gene_counts(self, tsv_file: str):
        """Load disease gene counts for stratification"""
        print(f"Loading disease gene counts from {tsv_file}...")

        df = pd.read_csv(tsv_file, sep='\t')

        # Create mapping from disease curie to gene count
        for _, row in df.iterrows():
            self.disease_gene_counts[row['disease_curie']] = row['gene_count']

        print(f"Loaded gene counts for {len(self.disease_gene_counts)} diseases")

        # Show stratification group sizes
        print("\nStratification groups:")
        for group_name, (min_genes, max_genes) in self.stratification_groups.items():
            if max_genes == float('inf'):
                count = sum(1 for count in self.disease_gene_counts.values() if count >= min_genes)
            else:
                count = sum(1 for count in self.disease_gene_counts.values()
                           if min_genes <= count <= max_genes)
            print(f"  {group_name}: {count} diseases")

    def load_enrichment_results(self, jsonl_file: str):
        """Load enrichment results"""
        print(f"Loading enrichment results from {jsonl_file}...")

        with open(jsonl_file, 'r') as f:
            for line in f:
                try:
                    result = json.loads(line.strip())
                    self.enrichment_data.append(result)
                except json.JSONDecodeError:
                    continue

        print(f"Loaded {len(self.enrichment_data)} enrichment results")

    def load_ground_truth(self, jsonl_file: str):
        """Load ground truth data"""
        print(f"Loading ground truth from {jsonl_file}...")

        with open(jsonl_file, 'r') as f:
            for line in f:
                try:
                    edge = json.loads(line.strip())
                    disease_curie = edge['source_curie']
                    term_curie = edge['target_curie']
                    term_category = edge['target_type']

                    if term_category in self.categories:
                        self.ground_truth[disease_curie][term_category].add(term_curie)
                except (json.JSONDecodeError, KeyError):
                    continue

        print(f"Loaded ground truth for {len(self.ground_truth)} diseases")

    def get_disease_group(self, disease_curie: str) -> str:
        """Get stratification group for a disease"""
        gene_count = self.disease_gene_counts.get(disease_curie, 0)

        for group_name, (min_genes, max_genes) in self.stratification_groups.items():
            if max_genes == float('inf'):
                if gene_count >= min_genes:
                    return group_name
            else:
                if min_genes <= gene_count <= max_genes:
                    return group_name

        return None  # Disease not in any group (shouldn't happen for diseases with ≥3 genes)

    def calculate_stratified_precision_recall(self, category: str, max_k: int = 50) -> Dict[str, Tuple[List[int], List[float], List[float]]]:
        """Calculate precision@k and recall@k for each stratification group"""
        print(f"Calculating stratified Precision@K and Recall@K for {category}...")

        # Group results by stratification group
        grouped_results = defaultdict(list)

        for result in self.enrichment_data:
            disease_curie = result['disease']['curie']
            group = self.get_disease_group(disease_curie)

            if group and disease_curie in self.ground_truth and category in self.ground_truth[disease_curie]:
                grouped_results[group].append(result)

        # Calculate metrics for each group
        group_metrics = {}
        k_values = list(range(1, max_k + 1))

        for group_name, group_results in grouped_results.items():
            print(f"  Processing {group_name}: {len(group_results)} diseases")

            precision_scores = []
            recall_scores = []

            for k in k_values:
                total_precision = 0
                total_recall = 0
                num_diseases = 0

                for result in group_results:
                    disease_curie = result['disease']['curie']

                    # Skip if no ground truth
                    if (disease_curie not in self.ground_truth or
                        category not in self.ground_truth[disease_curie] or
                        not self.ground_truth[disease_curie][category]):
                        continue

                    # Get enriched terms sorted by p-value
                    enrichments = result.get('enrichment_results', {}).get(category, [])
                    if not enrichments:
                        continue

                    enrichments.sort(key=lambda x: x['p_value'])
                    predicted_terms = [t['curie'] for t in enrichments[:k]]

                    if not predicted_terms:
                        continue

                    true_terms = self.ground_truth[disease_curie][category]
                    correct_predictions = set(predicted_terms).intersection(true_terms)

                    # Calculate precision@k and recall@k
                    precision = len(correct_predictions) / len(predicted_terms)
                    recall = len(correct_predictions) / len(true_terms) if true_terms else 0

                    total_precision += precision
                    total_recall += recall
                    num_diseases += 1

                avg_precision = total_precision / num_diseases if num_diseases > 0 else 0
                avg_recall = total_recall / num_diseases if num_diseases > 0 else 0

                precision_scores.append(avg_precision)
                recall_scores.append(avg_recall)

            group_metrics[group_name] = (k_values, precision_scores, recall_scores)

        return group_metrics

    def plot_stratified_results(self, all_results: Dict[str, Dict[str, Tuple]], output_file: str = 'stratified_precision_recall.png'):
        """Create stratified precision/recall plots"""

        # Create comprehensive plot
        fig, axes = plt.subplots(3, 2, figsize=(16, 18))

        colors = {'Low (3-5 genes)': 'blue', 'Medium (6-19 genes)': 'red', 'High (≥20 genes)': 'green'}

        # Plot for each category
        for cat_idx, category in enumerate(self.categories):
            category_short = category.replace('biolink:', '')

            if category in all_results:
                # Precision plot
                ax_prec = axes[cat_idx, 0]
                for group_name, (k_values, precision_scores, recall_scores) in all_results[category].items():
                    ax_prec.plot(k_values, precision_scores,
                               color=colors[group_name],
                               label=group_name,
                               linewidth=2, marker='o', markersize=2)

                ax_prec.set_xlabel('K')
                ax_prec.set_ylabel('Precision@K')
                ax_prec.set_title(f'{category_short}: Precision@K by Gene Count')
                ax_prec.legend()
                ax_prec.grid(True, alpha=0.3)

                # Recall plot
                ax_rec = axes[cat_idx, 1]
                for group_name, (k_values, precision_scores, recall_scores) in all_results[category].items():
                    ax_rec.plot(k_values, recall_scores,
                               color=colors[group_name],
                               label=group_name,
                               linewidth=2, marker='o', markersize=2)

                ax_rec.set_xlabel('K')
                ax_rec.set_ylabel('Recall@K')
                ax_rec.set_title(f'{category_short}: Recall@K by Gene Count')
                ax_rec.legend()
                ax_rec.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Stratified plots saved to {output_file}")
        plt.show()

    def print_stratified_summary(self, all_results: Dict[str, Dict[str, Tuple]]):
        """Print summary statistics for stratified analysis"""
        print("\n=== Stratified Analysis Summary ===")

        for category in self.categories:
            category_short = category.replace('biolink:', '')
            print(f"\n{category_short}:")

            if category in all_results:
                for group_name, (k_values, precision_scores, recall_scores) in all_results[category].items():
                    if precision_scores and recall_scores:
                        p1 = precision_scores[0] if len(precision_scores) > 0 else 0
                        p10 = precision_scores[9] if len(precision_scores) > 9 else 0
                        p20 = precision_scores[19] if len(precision_scores) > 19 else 0

                        r1 = recall_scores[0] if len(recall_scores) > 0 else 0
                        r10 = recall_scores[9] if len(recall_scores) > 9 else 0
                        r20 = recall_scores[19] if len(recall_scores) > 19 else 0

                        print(f"  {group_name}:")
                        print(f"    P@1={p1:.4f}, P@10={p10:.4f}, P@20={p20:.4f}")
                        print(f"    R@1={r1:.4f}, R@10={r10:.4f}, R@20={r20:.4f}")

    def run_stratified_analysis(self, max_k: int = 50) -> Dict[str, Dict[str, Tuple]]:
        """Run complete stratified analysis"""
        all_results = {}

        for category in self.categories:
            category_results = self.calculate_stratified_precision_recall(category, max_k)
            all_results[category] = category_results

        return all_results

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Stratified precision/recall analysis by gene count')
    parser.add_argument('--enrichment', '-e', default='fast_enrichment_results.jsonl',
                       help='Enrichment results JSONL file')
    parser.add_argument('--ground-truth', '-g', default='robokop_disease_term_edges_with_subclass.jsonl',
                       help='Ground truth JSONL file')
    parser.add_argument('--gene-counts', '-c', default='disease_gene_counts.tsv',
                       help='Disease gene counts TSV file')
    parser.add_argument('--max-k', '-k', type=int, default=50,
                       help='Maximum K value')

    args = parser.parse_args()

    # Initialize analyzer
    analyzer = StratifiedAnalyzer()
    analyzer.load_disease_gene_counts(args.gene_counts)
    analyzer.load_enrichment_results(args.enrichment)
    analyzer.load_ground_truth(args.ground_truth)

    # Run stratified analysis
    results = analyzer.run_stratified_analysis(args.max_k)

    # Create plots and summary
    analyzer.plot_stratified_results(results)
    analyzer.print_stratified_summary(results)

    print("\nStratified analysis complete!")

if __name__ == "__main__":
    main()