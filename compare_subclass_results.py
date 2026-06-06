#!/usr/bin/env python3
"""
Compare Original vs Subclass-Enhanced ROBOKOP Results

This script compares the original ROBOKOP disease-term edges with the
subclass-enhanced version to show how many new associations were added
through ontological reasoning.
"""

import json
from collections import defaultdict
from typing import Dict, Set, Tuple

def load_edges(jsonl_file: str) -> Dict[str, Set[str]]:
    """
    Load disease-term edges from JSONL file

    Returns:
        Dictionary mapping disease_curie -> set of term_curies
    """
    disease_terms = defaultdict(set)

    try:
        with open(jsonl_file, 'r') as f:
            for line in f:
                try:
                    edge = json.loads(line.strip())
                    disease_curie = edge['source_curie']
                    term_curie = edge['target_curie']
                    disease_terms[disease_curie].add(term_curie)
                except (json.JSONDecodeError, KeyError):
                    continue

        return disease_terms
    except FileNotFoundError:
        print(f"File {jsonl_file} not found")
        return defaultdict(set)

def load_edges_with_details(jsonl_file: str) -> Tuple[Dict[str, Set[str]], Dict[str, Dict]]:
    """
    Load edges with full details for analysis

    Returns:
        Tuple of (disease_terms_mapping, edge_details_by_key)
    """
    disease_terms = defaultdict(set)
    edge_details = {}

    try:
        with open(jsonl_file, 'r') as f:
            for line in f:
                try:
                    edge = json.loads(line.strip())
                    disease_curie = edge['source_curie']
                    term_curie = edge['target_curie']

                    disease_terms[disease_curie].add(term_curie)

                    # Store edge details with composite key
                    key = f"{disease_curie}|{term_curie}"
                    edge_details[key] = edge

                except (json.JSONDecodeError, KeyError):
                    continue

        return disease_terms, edge_details
    except FileNotFoundError:
        print(f"File {jsonl_file} not found")
        return defaultdict(set), {}

def compare_results(original_file: str, enhanced_file: str, output_file: str):
    """
    Compare original vs enhanced results and generate comparison report
    """
    print("Loading original ROBOKOP results...")
    original_edges = load_edges(original_file)

    print("Loading subclass-enhanced ROBOKOP results...")
    enhanced_edges, enhanced_details = load_edges_with_details(enhanced_file)

    if not enhanced_edges:
        print(f"Enhanced file {enhanced_file} not ready yet")
        return

    # Calculate statistics
    original_diseases = len(original_edges)
    enhanced_diseases = len(enhanced_edges)

    original_total_edges = sum(len(terms) for terms in original_edges.values())
    enhanced_total_edges = sum(len(terms) for terms in enhanced_edges.values())

    # Find new diseases and new edges
    new_diseases = set(enhanced_edges.keys()) - set(original_edges.keys())
    shared_diseases = set(enhanced_edges.keys()).intersection(set(original_edges.keys()))

    new_edges_count = 0
    new_edges_for_existing_diseases = 0
    direct_vs_inferred = {'direct': 0, 'inferred': 0}

    # Track category distributions
    category_stats = defaultdict(lambda: {'original': 0, 'enhanced': 0, 'new': 0})

    # Analyze new edges
    for disease_curie, term_curies in enhanced_edges.items():
        original_terms = original_edges.get(disease_curie, set())
        new_terms = term_curies - original_terms

        if new_terms:
            new_edges_count += len(new_terms)

            if disease_curie in original_edges:
                new_edges_for_existing_diseases += len(new_terms)

            # Analyze inference types for new edges
            for term_curie in new_terms:
                key = f"{disease_curie}|{term_curie}"
                if key in enhanced_details:
                    edge = enhanced_details[key]
                    inference_type = edge.get('inference_type', 'unknown')
                    category = edge.get('target_type', 'unknown')

                    if inference_type == 'direct':
                        direct_vs_inferred['direct'] += 1
                    elif inference_type.startswith('inferred'):
                        direct_vs_inferred['inferred'] += 1

                    category_stats[category]['new'] += 1

    # Count original and enhanced by category
    for disease_curie, term_curies in original_edges.items():
        # This is approximate since we don't have category info in original
        category_stats['total']['original'] += len(term_curies)

    for disease_curie, term_curies in enhanced_edges.items():
        for term_curie in term_curies:
            key = f"{disease_curie}|{term_curie}"
            if key in enhanced_details:
                category = enhanced_details[key].get('target_type', 'unknown')
                category_stats[category]['enhanced'] += 1

    # Generate report
    print(f"\n=== ROBOKOP Subclass Enhancement Comparison ===")
    print(f"Original results: {original_diseases:,} diseases, {original_total_edges:,} edges")
    print(f"Enhanced results: {enhanced_diseases:,} diseases, {enhanced_total_edges:,} edges")
    print(f"\nGains from subclass reasoning:")
    print(f"  New diseases with terms: {len(new_diseases):,}")
    print(f"  New edges total: {new_edges_count:,}")
    print(f"  New edges for existing diseases: {new_edges_for_existing_diseases:,}")
    print(f"  Percentage increase in edges: {((enhanced_total_edges - original_total_edges) / original_total_edges * 100):.1f}%")

    print(f"\nEdge types in new associations:")
    print(f"  Direct edges: {direct_vs_inferred['direct']:,}")
    print(f"  Inferred edges: {direct_vs_inferred['inferred']:,}")

    # Write detailed TSV report
    with open(output_file, 'w') as f:
        f.write("metric\toriginal\tenhanced\tdifference\tpercent_change\n")
        f.write(f"diseases\t{original_diseases}\t{enhanced_diseases}\t{enhanced_diseases - original_diseases}\t{((enhanced_diseases - original_diseases) / original_diseases * 100):.1f}%\n")
        f.write(f"total_edges\t{original_total_edges}\t{enhanced_total_edges}\t{enhanced_total_edges - original_total_edges}\t{((enhanced_total_edges - original_total_edges) / original_total_edges * 100):.1f}%\n")
        f.write(f"new_diseases\t0\t{len(new_diseases)}\t{len(new_diseases)}\tN/A\n")
        f.write(f"new_edges_total\t0\t{new_edges_count}\t{new_edges_count}\tN/A\n")
        f.write(f"new_edges_existing_diseases\t0\t{new_edges_for_existing_diseases}\t{new_edges_for_existing_diseases}\tN/A\n")
        f.write(f"direct_new_edges\t0\t{direct_vs_inferred['direct']}\t{direct_vs_inferred['direct']}\tN/A\n")
        f.write(f"inferred_new_edges\t0\t{direct_vs_inferred['inferred']}\t{direct_vs_inferred['inferred']}\tN/A\n")

    print(f"\nDetailed comparison saved to {output_file}")

    # Show some example inferred edges
    print(f"\nExample inferred edges (first 10):")
    inferred_count = 0
    for disease_curie, term_curies in enhanced_edges.items():
        original_terms = original_edges.get(disease_curie, set())
        new_terms = term_curies - original_terms

        for term_curie in new_terms:
            if inferred_count >= 10:
                break

            key = f"{disease_curie}|{term_curie}"
            if key in enhanced_details:
                edge = enhanced_details[key]
                if edge.get('inference_type', '').startswith('inferred'):
                    print(f"  {edge['source_name']} -> {edge['target_name']}")
                    print(f"    (inferred from {edge.get('original_disease_name', 'unknown')})")
                    inferred_count += 1

        if inferred_count >= 10:
            break

def main():
    """Main function"""
    original_file = "robokop_disease_term_edges.jsonl"
    enhanced_file = "robokop_disease_term_edges_with_subclass.jsonl"
    output_file = "subclass_enhancement_comparison.tsv"

    compare_results(original_file, enhanced_file, output_file)

if __name__ == "__main__":
    main()