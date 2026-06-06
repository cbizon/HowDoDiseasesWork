#!/usr/bin/env python3
"""
Precision@K and Recall@K Analysis for Enrichment Results

This script creates detailed precision and recall curves to complement the hits@k analysis,
showing how well our enrichment predictions perform across different values of K.
"""

import json
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from typing import Dict, List, Tuple
import argparse

class PrecisionRecallAnalyzer:
    """Class to analyze precision and recall at different K values"""

    def __init__(self):
        """Initialize the analyzer"""
        self.enrichment_data = []
        self.ground_truth = defaultdict(lambda: defaultdict(set))
        self.categories = ['biolink:BiologicalProcess', 'biolink:MolecularActivity', 'biolink:Pathway']

    def load_enrichment_results(self, jsonl_file: str) -> None:
        """Load enrichment analysis results"""
        print(f"Loading enrichment results from {jsonl_file}...")

        with open(jsonl_file, 'r') as f:
            for line in f:
                try:
                    result = json.loads(line.strip())
                    self.enrichment_data.append(result)
                except json.JSONDecodeError:
                    continue

        print(f"Loaded {len(self.enrichment_data)} disease enrichment results")

    def load_ground_truth(self, jsonl_file: str) -> None:
        """Load ground truth from ROBOKOP disease-term edges"""
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

    def calculate_precision_recall_at_k(self, category: str, max_k: int = 50) -> Tuple[List[int], List[float], List[float]]:
        """
        Calculate precision@k and recall@k for a specific category

        Returns:
            Tuple of (k_values, precision_scores, recall_scores)
        """
        print(f"Calculating Precision@K and Recall@K for {category}...")

        k_values = list(range(1, max_k + 1))
        precision_scores = []
        recall_scores = []

        for k in k_values:
            total_precision = 0
            total_recall = 0
            num_diseases = 0

            for result in self.enrichment_data:
                disease_curie = result['disease']['curie']

                # Skip if no ground truth for this disease and category
                if (disease_curie not in self.ground_truth or
                    category not in self.ground_truth[disease_curie] or
                    not self.ground_truth[disease_curie][category]):
                    continue

                # Get enriched terms for this category, sorted by p-value
                enrichments = result.get('enrichment_results', {}).get(category, [])
                if not enrichments:
                    continue

                # Sort by p-value and get top k
                enrichments.sort(key=lambda x: x['p_value'])
                predicted_terms = [t['curie'] for t in enrichments[:k]]

                if not predicted_terms:
                    continue

                true_terms = self.ground_truth[disease_curie][category]

                # Calculate precision@k and recall@k
                correct_predictions = set(predicted_terms).intersection(true_terms)

                precision = len(correct_predictions) / len(predicted_terms)
                recall = len(correct_predictions) / len(true_terms) if true_terms else 0

                total_precision += precision
                total_recall += recall
                num_diseases += 1

            avg_precision = total_precision / num_diseases if num_diseases > 0 else 0
            avg_recall = total_recall / num_diseases if num_diseases > 0 else 0

            precision_scores.append(avg_precision)
            recall_scores.append(avg_recall)

        return k_values, precision_scores, recall_scores

    def plot_precision_recall_curves(self, results: Dict[str, Tuple[List[int], List[float], List[float]]],
                                   output_file: str = 'precision_recall_at_k.png') -> None:
        """Create precision and recall curves"""

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        colors = ['blue', 'red', 'green']

        # Precision@K plot
        for i, category in enumerate(self.categories):
            if category in results:
                k_values, precision_scores, recall_scores = results[category]
                ax1.plot(k_values, precision_scores,
                        color=colors[i],
                        label=category.replace('biolink:', ''),
                        marker='o', markersize=2, linewidth=2)

        ax1.set_xlabel('K (Number of Top Predictions)', fontsize=12)
        ax1.set_ylabel('Precision@K', fontsize=12)
        ax1.set_title('Precision@K vs K', fontsize=14)
        ax1.legend(fontsize=10)
        ax1.grid(True, alpha=0.3)

        # Set appropriate y-limits based on data
        all_precision_values = []
        for cat in results:
            if results[cat][1]:  # precision scores exist
                all_precision_values.extend(results[cat][1])
        max_precision = max(all_precision_values) if all_precision_values else 0.1
        ax1.set_ylim(0, max(0.05, max_precision * 1.1))

        # Recall@K plot
        for i, category in enumerate(self.categories):
            if category in results:
                k_values, precision_scores, recall_scores = results[category]
                ax2.plot(k_values, recall_scores,
                        color=colors[i],
                        label=category.replace('biolink:', ''),
                        marker='o', markersize=2, linewidth=2)

        ax2.set_xlabel('K (Number of Top Predictions)', fontsize=12)
        ax2.set_ylabel('Recall@K', fontsize=12)
        ax2.set_title('Recall@K vs K', fontsize=14)
        ax2.legend(fontsize=10)
        ax2.grid(True, alpha=0.3)

        # Set appropriate y-limits based on data
        all_recall_values = []
        for cat in results:
            if results[cat][2]:  # recall scores exist
                all_recall_values.extend(results[cat][2])
        max_recall = max(all_recall_values) if all_recall_values else 0.1
        ax2.set_ylim(0, max(0.1, max_recall * 1.1))

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Precision/Recall plots saved to {output_file}")
        plt.show()

    def create_combined_metrics_plot(self, pr_results: Dict, hits_results: Dict,
                                   output_file: str = 'combined_metrics_at_k.png') -> None:
        """Create a comprehensive plot showing Hits@K, Precision@K, and Recall@K"""

        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        colors = ['blue', 'red', 'green']

        # Combined plot: all metrics for BiologicalProcess
        ax = axes[0, 0]
        category = 'biolink:BiologicalProcess'
        if category in pr_results and category in hits_results:
            k_values, precision_scores, recall_scores = pr_results[category]
            _, hits_scores = hits_results[category]

            ax.plot(k_values, hits_scores, color='blue', label='Hits@K', linewidth=2)
            ax.plot(k_values, precision_scores, color='red', label='Precision@K', linewidth=2)
            ax.plot(k_values, recall_scores, color='green', label='Recall@K', linewidth=2)

        ax.set_xlabel('K', fontsize=12)
        ax.set_ylabel('Score', fontsize=12)
        ax.set_title('BiologicalProcess: All Metrics', fontsize=12)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, 1)

        # MolecularActivity metrics
        ax = axes[0, 1]
        category = 'biolink:MolecularActivity'
        if category in pr_results and category in hits_results:
            k_values, precision_scores, recall_scores = pr_results[category]
            _, hits_scores = hits_results[category]

            ax.plot(k_values, hits_scores, color='blue', label='Hits@K', linewidth=2)
            ax.plot(k_values, precision_scores, color='red', label='Precision@K', linewidth=2)
            ax.plot(k_values, recall_scores, color='green', label='Recall@K', linewidth=2)

        ax.set_xlabel('K', fontsize=12)
        ax.set_ylabel('Score', fontsize=12)
        ax.set_title('MolecularActivity: All Metrics', fontsize=12)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, 0.1)  # Lower scale for this category

        # Pathway metrics
        ax = axes[1, 0]
        category = 'biolink:Pathway'
        if category in pr_results and category in hits_results:
            k_values, precision_scores, recall_scores = pr_results[category]
            _, hits_scores = hits_results[category]

            ax.plot(k_values, hits_scores, color='blue', label='Hits@K', linewidth=2)
            ax.plot(k_values, precision_scores, color='red', label='Precision@K', linewidth=2)
            ax.plot(k_values, recall_scores, color='green', label='Recall@K', linewidth=2)

        ax.set_xlabel('K', fontsize=12)
        ax.set_ylabel('Score', fontsize=12)
        ax.set_title('Pathway: All Metrics', fontsize=12)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, 0.05)  # Lower scale for this category

        # Summary comparison at K=10
        ax = axes[1, 1]
        categories_short = ['BiologicalProcess', 'MolecularActivity', 'Pathway']
        k10_hits = [hits_results[f'biolink:{cat}'][1][9] if f'biolink:{cat}' in hits_results else 0 for cat in categories_short]
        k10_precision = [pr_results[f'biolink:{cat}'][1][9] if f'biolink:{cat}' in pr_results else 0 for cat in categories_short]
        k10_recall = [pr_results[f'biolink:{cat}'][2][9] if f'biolink:{cat}' in pr_results else 0 for cat in categories_short]

        x = np.arange(len(categories_short))
        width = 0.25

        ax.bar(x - width, k10_hits, width, label='Hits@10', color='blue', alpha=0.7)
        ax.bar(x, k10_precision, width, label='Precision@10', color='red', alpha=0.7)
        ax.bar(x + width, k10_recall, width, label='Recall@10', color='green', alpha=0.7)

        ax.set_xlabel('Category', fontsize=12)
        ax.set_ylabel('Score', fontsize=12)
        ax.set_title('Metrics Comparison at K=10', fontsize=12)
        ax.set_xticks(x)
        ax.set_xticklabels(categories_short, rotation=45)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3, axis='y')

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Combined metrics plot saved to {output_file}")
        plt.show()

    def run_analysis(self, max_k: int = 50) -> Dict:
        """Run complete precision/recall analysis"""
        results = {}

        for category in self.categories:
            k_values, precision_scores, recall_scores = self.calculate_precision_recall_at_k(category, max_k)
            results[category] = (k_values, precision_scores, recall_scores)

        return results

def load_hits_results(tsv_file: str) -> Dict:
    """Load hits@k results from detailed TSV file"""
    hits_results = {}

    try:
        with open(tsv_file, 'r') as f:
            lines = f.readlines()

        # Parse header
        header = lines[1].strip().split('\t')  # Skip first empty line
        categories = ['biolink:BiologicalProcess', 'biolink:MolecularActivity', 'biolink:Pathway']

        # Parse data
        k_values = []
        bp_scores = []
        ma_scores = []
        pw_scores = []

        for line in lines[2:]:  # Skip header lines
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    k_values.append(int(parts[0]))
                    bp_scores.append(float(parts[1]))
                    ma_scores.append(float(parts[2]))
                    pw_scores.append(float(parts[3]))

        hits_results['biolink:BiologicalProcess'] = (k_values, bp_scores)
        hits_results['biolink:MolecularActivity'] = (k_values, ma_scores)
        hits_results['biolink:Pathway'] = (k_values, pw_scores)

    except FileNotFoundError:
        print(f"Could not find {tsv_file}")

    return hits_results

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Analyze precision and recall at K values')
    parser.add_argument('--enrichment', '-e', default='fast_enrichment_results.jsonl',
                       help='Enrichment results JSONL file')
    parser.add_argument('--ground-truth', '-g', default='robokop_disease_term_edges_with_subclass.jsonl',
                       help='Ground truth JSONL file')
    parser.add_argument('--max-k', '-k', type=int, default=100,
                       help='Maximum K value (default: 100)')
    parser.add_argument('--hits-file', default='enrichment_hits_at_k_detailed.tsv',
                       help='Hits@K results TSV file')

    args = parser.parse_args()

    # Initialize analyzer
    analyzer = PrecisionRecallAnalyzer()
    analyzer.load_enrichment_results(args.enrichment)
    analyzer.load_ground_truth(args.ground_truth)

    # Run precision/recall analysis
    pr_results = analyzer.run_analysis(args.max_k)

    # Load hits@k results
    hits_results = load_hits_results(args.hits_file)

    # Create plots
    analyzer.plot_precision_recall_curves(pr_results)

    if hits_results:
        analyzer.create_combined_metrics_plot(pr_results, hits_results)

    print("\\nPrecision/Recall analysis complete!")

if __name__ == "__main__":
    main()