#!/usr/bin/env python3
"""
Create Precision@K and Recall@K Plots

Simple script to generate precision and recall curves for enrichment analysis results.
"""

import json
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from typing import Dict, List, Tuple

def load_data():
    """Load enrichment results and ground truth data"""
    print('Loading enrichment results...')
    enrichment_data = []
    with open('fast_enrichment_results.jsonl', 'r') as f:
        for line in f:
            try:
                enrichment_data.append(json.loads(line.strip()))
            except:
                pass

    print(f'Loaded {len(enrichment_data)} enrichment results')

    print('Loading ground truth...')
    ground_truth = defaultdict(lambda: defaultdict(set))
    categories = ['biolink:BiologicalProcess', 'biolink:MolecularActivity', 'biolink:Pathway']

    count = 0
    with open('robokop_disease_term_edges_with_subclass.jsonl', 'r') as f:
        for line in f:
            try:
                edge = json.loads(line.strip())
                disease_curie = edge['source_curie']
                term_curie = edge['target_curie']
                term_category = edge['target_type']

                if term_category in categories:
                    ground_truth[disease_curie][term_category].add(term_curie)
                    count += 1
                    if count % 1000000 == 0:
                        print(f'Loaded {count:,} ground truth edges...')
            except:
                pass

    print(f'Loaded ground truth for {len(ground_truth)} diseases')
    return enrichment_data, ground_truth

def calculate_precision_recall(enrichment_data, ground_truth, category: str, max_k: int = 100):
    """Calculate precision and recall at different K values for a category"""
    print(f'Calculating precision/recall for {category}...')

    k_values = list(range(1, max_k + 1))
    precision_scores = []
    recall_scores = []

    for k in k_values:
        if k % 20 == 0:
            print(f'  Processing K={k}...')

        total_precision = 0
        total_recall = 0
        num_diseases = 0

        for result in enrichment_data:
            disease_curie = result['disease']['curie']

            # Skip if no ground truth for this disease and category
            if (disease_curie not in ground_truth or
                category not in ground_truth[disease_curie] or
                not ground_truth[disease_curie][category]):
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

            true_terms = ground_truth[disease_curie][category]
            correct_predictions = set(predicted_terms).intersection(true_terms)

            # Calculate precision@k and recall@k for this disease
            precision = len(correct_predictions) / len(predicted_terms)
            recall = len(correct_predictions) / len(true_terms) if true_terms else 0

            total_precision += precision
            total_recall += recall
            num_diseases += 1

        # Average across all diseases
        avg_precision = total_precision / num_diseases if num_diseases > 0 else 0
        avg_recall = total_recall / num_diseases if num_diseases > 0 else 0

        precision_scores.append(avg_precision)
        recall_scores.append(avg_recall)

    return k_values, precision_scores, recall_scores

def create_plots(results_dict):
    """Create precision@k and recall@k plots"""

    # Create figure with subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    colors = ['blue', 'red', 'green']
    categories = ['biolink:BiologicalProcess', 'biolink:MolecularActivity', 'biolink:Pathway']

    # Plot 1: Precision@K for all categories
    for i, category in enumerate(categories):
        if category in results_dict:
            k_values, precision_scores, recall_scores = results_dict[category]
            ax1.plot(k_values, precision_scores,
                    color=colors[i],
                    label=category.replace('biolink:', ''),
                    linewidth=2)

    ax1.set_xlabel('K (Number of Top Predictions)', fontsize=12)
    ax1.set_ylabel('Precision@K', fontsize=12)
    ax1.set_title('Precision@K vs K (All Categories)', fontsize=14)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(1, 100)

    # Plot 2: Recall@K for all categories
    for i, category in enumerate(categories):
        if category in results_dict:
            k_values, precision_scores, recall_scores = results_dict[category]
            ax2.plot(k_values, recall_scores,
                    color=colors[i],
                    label=category.replace('biolink:', ''),
                    linewidth=2)

    ax2.set_xlabel('K (Number of Top Predictions)', fontsize=12)
    ax2.set_ylabel('Recall@K', fontsize=12)
    ax2.set_title('Recall@K vs K (All Categories)', fontsize=14)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(1, 100)

    # Plot 3: BiologicalProcess detailed view
    category = 'biolink:BiologicalProcess'
    if category in results_dict:
        k_values, precision_scores, recall_scores = results_dict[category]
        ax3.plot(k_values, precision_scores, color='blue', label='Precision@K', linewidth=2)
        ax3.plot(k_values, recall_scores, color='red', label='Recall@K', linewidth=2)

        ax3.set_xlabel('K', fontsize=12)
        ax3.set_ylabel('Score', fontsize=12)
        ax3.set_title('BiologicalProcess: Precision@K and Recall@K', fontsize=14)
        ax3.legend(fontsize=10)
        ax3.grid(True, alpha=0.3)
        ax3.set_xlim(1, 100)

    # Plot 4: Summary at specific K values
    k_points = [1, 5, 10, 20, 50, 100]
    precision_data = []
    recall_data = []
    category_names = []

    for category in categories:
        if category in results_dict:
            k_values, precision_scores, recall_scores = results_dict[category]
            category_names.append(category.replace('biolink:', ''))

            # Get precision/recall at K=20
            k20_precision = precision_scores[19] if len(precision_scores) > 19 else 0
            k20_recall = recall_scores[19] if len(recall_scores) > 19 else 0
            precision_data.append(k20_precision)
            recall_data.append(k20_recall)

    x = np.arange(len(category_names))
    width = 0.35

    ax4.bar(x - width/2, precision_data, width, label='Precision@20', color='blue', alpha=0.7)
    ax4.bar(x + width/2, recall_data, width, label='Recall@20', color='red', alpha=0.7)

    ax4.set_xlabel('Category', fontsize=12)
    ax4.set_ylabel('Score', fontsize=12)
    ax4.set_title('Precision@20 and Recall@20 by Category', fontsize=14)
    ax4.set_xticks(x)
    ax4.set_xticklabels(category_names)
    ax4.legend(fontsize=10)
    ax4.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig('precision_recall_at_k_plots.png', dpi=300, bbox_inches='tight')
    print("Plots saved to precision_recall_at_k_plots.png")
    plt.show()

def main():
    """Main function"""
    # Load data
    enrichment_data, ground_truth = load_data()

    # Calculate precision/recall for each category
    categories = ['biolink:BiologicalProcess', 'biolink:MolecularActivity', 'biolink:Pathway']
    results = {}

    for category in categories:
        k_values, precision_scores, recall_scores = calculate_precision_recall(
            enrichment_data, ground_truth, category, max_k=100)
        results[category] = (k_values, precision_scores, recall_scores)

        # Print some key metrics
        p1 = precision_scores[0] if precision_scores else 0
        p10 = precision_scores[9] if len(precision_scores) > 9 else 0
        p20 = precision_scores[19] if len(precision_scores) > 19 else 0
        p50 = precision_scores[49] if len(precision_scores) > 49 else 0

        r1 = recall_scores[0] if recall_scores else 0
        r10 = recall_scores[9] if len(recall_scores) > 9 else 0
        r20 = recall_scores[19] if len(recall_scores) > 19 else 0
        r50 = recall_scores[49] if len(recall_scores) > 49 else 0

        print(f"\\n{category.replace('biolink:', '')} Results:")
        print(f"  P@1={p1:.4f}, P@10={p10:.4f}, P@20={p20:.4f}, P@50={p50:.4f}")
        print(f"  R@1={r1:.4f}, R@10={r10:.4f}, R@20={r20:.4f}, R@50={r50:.4f}")

    # Create plots
    create_plots(results)

if __name__ == "__main__":
    main()