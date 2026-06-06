#!/usr/bin/env python3
"""
Compare Enrichment Results vs ROBOKOP Ground Truth
This script performs hits@k analysis comparing our enrichment analysis results
with ROBOKOP knowledge graph as ground truth.

Since we don't have the complete ROBOKOP extraction yet, this script demonstrates
the analysis approach using mock ground truth data. Once ROBOKOP data is available,
we can replace the mock data with real ground truth.

The analysis calculates:
- Hits@K for each category (BiologicalProcess, MolecularActivity, Pathway)
- Plots showing hits@K vs K for each category
- Overall precision/recall metrics
"""

import json
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, List, Set, Tuple
import argparse
from collections import defaultdict

class EnrichmentEvaluator:
    """Class to evaluate enrichment results against ground truth"""

    def __init__(self):
        """Initialize the evaluator"""
        self.enrichment_data = []
        self.ground_truth = defaultdict(lambda: defaultdict(set))  # disease_curie -> category -> set of term_curies
        self.categories = ['biolink:BiologicalProcess', 'biolink:MolecularActivity', 'biolink:Pathway']

    def load_enrichment_results(self, jsonl_file: str) -> None:
        """
        Load enrichment analysis results from JSONL file

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
                    print(f"Error parsing line: {e}")
                    continue

        print(f"Loaded {len(self.enrichment_data)} disease enrichment results")

    def load_ground_truth_robokop(self, jsonl_file: str) -> None:
        """
        Load ground truth from ROBOKOP disease-term edges

        Args:
            jsonl_file: Path to ROBOKOP disease-term edges JSONL file
        """
        print(f"Loading ROBOKOP ground truth from {jsonl_file}...")

        try:
            with open(jsonl_file, 'r') as f:
                for line in f:
                    try:
                        edge = json.loads(line.strip())
                        disease_curie = edge['source_curie']
                        term_curie = edge['target_curie']
                        term_category = edge['target_type']  # Get the category of the term

                        # Only add if it's one of our target categories
                        if term_category in self.categories:
                            self.ground_truth[disease_curie][term_category].add(term_curie)
                    except (json.JSONDecodeError, KeyError) as e:
                        continue

            # Calculate totals
            total_diseases = len(self.ground_truth)
            category_stats = {}
            total_edges = 0

            for category in self.categories:
                category_edges = sum(len(self.ground_truth[disease][category])
                                   for disease in self.ground_truth)
                category_stats[category] = category_edges
                total_edges += category_edges

            print(f"Loaded ground truth: {total_diseases} diseases, {total_edges} disease-term edges")
            for category in self.categories:
                print(f"  {category.replace('biolink:', '')}: {category_stats[category]} edges")

        except FileNotFoundError:
            print(f"Ground truth file {jsonl_file} not found. Creating mock ground truth for demonstration.")
            self._create_mock_ground_truth()

    def _create_mock_ground_truth(self) -> None:
        """
        Create mock ground truth data for demonstration purposes
        This simulates what real ROBOKOP ground truth would look like
        """
        print("Creating mock ground truth data...")

        # For each disease in our enrichment results, randomly select some terms as "true"
        import random
        random.seed(42)  # For reproducibility

        for result in self.enrichment_data[:100]:  # Use first 100 diseases for demo
            disease_curie = result['disease']['curie']
            enrichments = result.get('enrichment_results', {})

            # For each category, randomly mark some enriched terms as "true" in ground truth
            for category in self.categories:
                terms = enrichments.get(category, [])
                if terms:
                    # Randomly select 20-50% of enriched terms as "ground truth"
                    num_true = max(1, int(len(terms) * random.uniform(0.2, 0.5)))
                    true_terms = random.sample([t['curie'] for t in terms], min(num_true, len(terms)))
                    self.ground_truth[disease_curie][category].update(true_terms)

        # Calculate totals for mock data
        total_diseases = len(self.ground_truth)
        total_edges = sum(sum(len(self.ground_truth[disease][cat])
                             for cat in self.categories)
                         for disease in self.ground_truth)
        print(f"Created mock ground truth: {total_diseases} diseases, {total_edges} edges")

    def calculate_hits_at_k(self, disease_curie: str, predicted_terms: List[str], category: str, k: int) -> float:
        """
        Calculate hits@k for a single disease and category

        Args:
            disease_curie: Disease identifier
            predicted_terms: List of predicted term CURIEs (ranked by confidence/p-value)
            category: Category to evaluate (e.g., 'biolink:BiologicalProcess')
            k: Number of top predictions to consider

        Returns:
            Hits@k score (1.0 if any of top-k predictions are in ground truth, 0.0 otherwise)
        """
        if disease_curie not in self.ground_truth or category not in self.ground_truth[disease_curie]:
            return 0.0

        true_terms = self.ground_truth[disease_curie][category]
        if not true_terms:  # No ground truth terms for this category
            return 0.0

        top_k_predictions = set(predicted_terms[:k])

        # Return 1.0 if any top-k prediction is correct, 0.0 otherwise
        return 1.0 if len(top_k_predictions.intersection(true_terms)) > 0 else 0.0

    def evaluate_category(self, category: str, max_k: int = 50) -> Tuple[List[int], List[float]]:
        """
        Evaluate hits@k for a specific category

        Args:
            category: Category to evaluate (e.g., 'biolink:BiologicalProcess')
            max_k: Maximum k value to evaluate

        Returns:
            Tuple of (k_values, hits_at_k_scores)
        """
        print(f"Evaluating {category}...")

        k_values = list(range(1, max_k + 1))
        hits_at_k_scores = []

        for k in k_values:
            total_hits = 0
            total_diseases = 0

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

                # Sort by p-value (ascending - lower p-values are better)
                enrichments.sort(key=lambda x: x['p_value'])
                predicted_terms = [t['curie'] for t in enrichments]

                # Calculate hits@k
                hits = self.calculate_hits_at_k(disease_curie, predicted_terms, category, k)
                total_hits += hits
                total_diseases += 1

            # Calculate average hits@k across all diseases
            avg_hits_at_k = total_hits / total_diseases if total_diseases > 0 else 0.0
            hits_at_k_scores.append(avg_hits_at_k)

        return k_values, hits_at_k_scores

    def plot_hits_at_k(self, results: Dict[str, Tuple[List[int], List[float]]], output_file: str = 'hits_at_k_plot.png') -> None:
        """
        Create plots showing hits@k vs k for each category

        Args:
            results: Dictionary mapping category to (k_values, hits_at_k_scores)
            output_file: Output filename for plot
        """
        plt.figure(figsize=(12, 8))

        colors = ['blue', 'red', 'green']
        for i, category in enumerate(self.categories):
            if category in results:
                k_values, hits_scores = results[category]
                plt.plot(k_values, hits_scores,
                        color=colors[i],
                        label=category.replace('biolink:', ''),
                        marker='o',
                        markersize=3,
                        linewidth=2)

        plt.xlabel('K (Number of Top Predictions)', fontsize=12)
        plt.ylabel('Hits@K', fontsize=12)
        plt.title('Enrichment Analysis Performance: Hits@K vs K', fontsize=14)
        plt.legend(fontsize=10)
        plt.grid(True, alpha=0.3)
        plt.ylim(0, 1.05)

        # Add some styling
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_file}")
        plt.show()

    def print_summary_stats(self, results: Dict[str, Tuple[List[int], List[float]]]) -> None:
        """
        Print summary statistics for the evaluation

        Args:
            results: Dictionary mapping category to (k_values, hits_at_k_scores)
        """
        print("\\n=== Enrichment Analysis Evaluation Summary ===")
        print(f"Total diseases evaluated: {len(self.enrichment_data)}")
        print(f"Diseases with ground truth: {len(self.ground_truth)}")

        # Show ground truth statistics by category
        print("\\nGround truth statistics:")
        for category in self.categories:
            diseases_with_gt = sum(1 for disease in self.ground_truth
                                 if self.ground_truth[disease][category])
            total_terms = sum(len(self.ground_truth[disease][category])
                            for disease in self.ground_truth)
            print(f"  {category.replace('biolink:', '')}: {diseases_with_gt} diseases, {total_terms} terms")

        print("\\nHits@K Performance:")
        for category in self.categories:
            if category in results:
                k_values, hits_scores = results[category]
                hits_at_1 = hits_scores[0] if hits_scores else 0
                hits_at_5 = hits_scores[4] if len(hits_scores) > 4 else 0
                hits_at_10 = hits_scores[9] if len(hits_scores) > 9 else 0
                hits_at_20 = hits_scores[19] if len(hits_scores) > 19 else 0

                print(f"  {category.replace('biolink:', '')}:")
                print(f"    Hits@1:  {hits_at_1:.3f}")
                print(f"    Hits@5:  {hits_at_5:.3f}")
                print(f"    Hits@10: {hits_at_10:.3f}")
                print(f"    Hits@20: {hits_at_20:.3f}")

    def run_evaluation(self, max_k: int = 50) -> Dict[str, Tuple[List[int], List[float]]]:
        """
        Run complete evaluation for all categories

        Args:
            max_k: Maximum k value to evaluate

        Returns:
            Dictionary mapping category to (k_values, hits_at_k_scores)
        """
        results = {}

        for category in self.categories:
            k_values, hits_scores = self.evaluate_category(category, max_k)
            results[category] = (k_values, hits_scores)

        return results

    def calculate_precision_recall_at_k(self, category: str, k: int) -> Tuple[float, float]:
        """
        Calculate precision and recall at k for a category

        Args:
            category: Category to evaluate
            k: Number of top predictions to consider

        Returns:
            Tuple of (precision@k, recall@k)
        """
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

            true_terms = self.ground_truth[disease_curie][category]

            # Calculate precision@k and recall@k
            if predicted_terms:
                correct_predictions = set(predicted_terms).intersection(true_terms)
                precision = len(correct_predictions) / len(predicted_terms)
                recall = len(correct_predictions) / len(true_terms) if true_terms else 0

                total_precision += precision
                total_recall += recall
                num_diseases += 1

        avg_precision = total_precision / num_diseases if num_diseases > 0 else 0
        avg_recall = total_recall / num_diseases if num_diseases > 0 else 0

        return avg_precision, avg_recall

    def save_detailed_results(self, results: Dict[str, Tuple[List[int], List[float]]],
                            output_file: str = 'hits_at_k_detailed.tsv') -> None:
        """
        Save detailed hits@k results to TSV file

        Args:
            results: Dictionary mapping category to (k_values, hits_at_k_scores)
            output_file: Output TSV filename
        """
        print(f"Saving detailed results to {output_file}...")

        with open(output_file, 'w') as f:
            # Header
            f.write("K\\tBiologicalProcess\\tMolecularActivity\\tPathway\\n")

            # Get max k value
            max_k = max(len(scores[1]) for scores in results.values() if scores[1])

            for k in range(1, max_k + 1):
                row = [str(k)]
                for category in self.categories:
                    if category in results and k <= len(results[category][1]):
                        score = results[category][1][k-1]  # k-1 because list is 0-indexed
                        row.append(f"{score:.4f}")
                    else:
                        row.append("0.0000")
                f.write("\\t".join(row) + "\\n")

        print(f"Detailed results saved to {output_file}")

def main():
    """Main function to run enrichment evaluation"""

    parser = argparse.ArgumentParser(description='Evaluate enrichment results against ROBOKOP ground truth')
    parser.add_argument('--enrichment', '-e', default='fast_enrichment_results.jsonl',
                       help='Enrichment results JSONL file (default: fast_enrichment_results.jsonl)')
    parser.add_argument('--ground-truth', '-g', default='robokop_disease_term_edges_with_subclass.jsonl',
                       help='ROBOKOP ground truth JSONL file (default: robokop_disease_term_edges.jsonl)')
    parser.add_argument('--max-k', '-k', type=int, default=50,
                       help='Maximum k value for hits@k evaluation (default: 50)')
    parser.add_argument('--output-plot', '-o', default='enrichment_hits_at_k.png',
                       help='Output plot filename (default: enrichment_hits_at_k.png)')

    args = parser.parse_args()

    # Initialize evaluator
    evaluator = EnrichmentEvaluator()

    # Load data
    evaluator.load_enrichment_results(args.enrichment)
    evaluator.load_ground_truth_robokop(args.ground_truth)

    # Run evaluation
    results = evaluator.run_evaluation(args.max_k)

    # Generate plots and summary
    evaluator.plot_hits_at_k(results, args.output_plot)
    evaluator.print_summary_stats(results)

    # Save detailed TSV results
    evaluator.save_detailed_results(results, 'enrichment_hits_at_k_detailed.tsv')

    # Show some precision/recall metrics
    print("\\nPrecision/Recall at K=10:")
    for category in evaluator.categories:
        precision, recall = evaluator.calculate_precision_recall_at_k(category, 10)
        print(f"  {category.replace('biolink:', '')}: P@10={precision:.3f}, R@10={recall:.3f}")

    print(f"\\nEvaluation complete! Plot saved to {args.output_plot}")

if __name__ == "__main__":
    main()
