#!/usr/bin/env python3
"""
Parse Disease Enrichment Analysis Results to Tab-Delimited Format

This script parses the JSON results from disease_enrichment_analysis.py
and creates tab-delimited files for easier review and analysis.

Creates multiple output files:
1. disease_summary.tsv - Disease overview with gene counts
2. enrichment_results.tsv - All enrichment results with p-values
3. error_summary.tsv - Timeout and error details
"""

import json
import csv
import sys
from typing import List, Dict, Any
import argparse

def parse_disease_summary(results: List[Dict[str, Any]], output_file: str):
    """
    Create disease summary TSV with basic information

    Args:
        results: Analysis results from JSON
        output_file: Output TSV filename
    """
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t')

        # Header
        writer.writerow([
            'Rank', 'Disease_CURIE', 'Disease_Name', 'Gene_Count',
            'BP_Enriched', 'MA_Enriched', 'Pathway_Enriched',
            'BP_Error', 'MA_Error', 'Pathway_Error', 'Has_Genes'
        ])

        # Data rows
        for i, result in enumerate(results, 1):
            disease = result['disease']
            gene_count = result['gene_count']
            enrichment_results = result.get('enrichment_results', {})
            enrichment_errors = result.get('enrichment_errors', {})

            # Count enriched terms for each category
            bp_count = len(enrichment_results.get('biolink:BiologicalProcess', []))
            ma_count = len(enrichment_results.get('biolink:MolecularActivity', []))
            pathway_count = len(enrichment_results.get('biolink:Pathway', []))

            # Check for errors
            bp_error = 'TIMEOUT' if 'timeout' in enrichment_errors.get('biolink:BiologicalProcess', '').lower() else ('ERROR' if 'biolink:BiologicalProcess' in enrichment_errors else '')
            ma_error = 'TIMEOUT' if 'timeout' in enrichment_errors.get('biolink:MolecularActivity', '').lower() else ('ERROR' if 'biolink:MolecularActivity' in enrichment_errors else '')
            pathway_error = 'TIMEOUT' if 'timeout' in enrichment_errors.get('biolink:Pathway', '').lower() else ('ERROR' if 'biolink:Pathway' in enrichment_errors else '')

            writer.writerow([
                i,
                disease['curie'],
                disease['name'],
                gene_count,
                bp_count,
                ma_count,
                pathway_count,
                bp_error,
                ma_error,
                pathway_error,
                'Yes' if gene_count > 0 else 'No'
            ])

    print(f"Disease summary written to {output_file}")

def parse_enrichment_results(results: List[Dict[str, Any]], output_file: str, min_genes: int = 1):
    """
    Create enrichment results TSV with the exact format requested:
    disease_curie, name, enrichment_type, enriched_entity, p_value, rank, number_of_enriched_entities

    Args:
        results: Analysis results from JSON
        output_file: Output TSV filename
        min_genes: Minimum gene count to include
    """
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t')

        # Header - exact format requested plus enriched entity CURIE
        writer.writerow([
            'disease_curie', 'name', 'enrichment_type', 'enriched_entity_curie',
            'enriched_entity', 'p_value', 'rank', 'number_of_enriched_entities'
        ])

        # Data rows
        for result in results:
            disease = result['disease']
            gene_count = result['gene_count']

            # Skip diseases with too few genes
            if gene_count < min_genes:
                continue

            enrichment_results = result.get('enrichment_results', {})

            for category, terms in enrichment_results.items():
                # Sort terms by p-value (ascending - most significant first)
                sorted_terms = sorted(terms, key=lambda x: x.get('p_value', float('inf')))
                total_enriched = len(sorted_terms)

                for rank, term in enumerate(sorted_terms, 1):
                    writer.writerow([
                        disease['curie'],
                        disease['name'],
                        category,
                        term.get('curie', ''),
                        term.get('name', ''),
                        term.get('p_value', ''),
                        rank,
                        total_enriched
                    ])

    print(f"Enrichment results written to {output_file}")

def parse_error_summary(results: List[Dict[str, Any]], output_file: str):
    """
    Create error summary TSV with timeout and error details

    Args:
        results: Analysis results from JSON
        output_file: Output TSV filename
    """
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t')

        # Header
        writer.writerow([
            'Disease_CURIE', 'Disease_Name', 'Gene_Count', 'Category',
            'Error_Type', 'Error_Message'
        ])

        # Data rows
        for result in results:
            disease = result['disease']
            gene_count = result['gene_count']
            enrichment_errors = result.get('enrichment_errors', {})

            if enrichment_errors:
                for category, error_msg in enrichment_errors.items():
                    error_type = 'TIMEOUT' if 'timeout' in error_msg.lower() or 'read timed out' in error_msg.lower() else 'ERROR'

                    writer.writerow([
                        disease['curie'],
                        disease['name'],
                        gene_count,
                        category,
                        error_type,
                        error_msg
                    ])

    print(f"Error summary written to {output_file}")

def parse_gene_counts_distribution(results: List[Dict[str, Any]], output_file: str):
    """
    Create gene count distribution TSV for analysis

    Args:
        results: Analysis results from JSON
        output_file: Output TSV filename
    """
    # Count diseases by gene count bins
    gene_counts = [result['gene_count'] for result in results]
    gene_counts.sort()

    # Create bins
    bins = [0, 1, 5, 10, 25, 50, 100, 200, 500, 1000, float('inf')]
    bin_labels = ['0', '1-4', '5-9', '10-24', '25-49', '50-99', '100-199', '200-499', '500-999', '1000+']
    bin_counts = [0] * len(bin_labels)

    for gene_count in gene_counts:
        for i, bin_max in enumerate(bins[1:]):
            if gene_count < bin_max:
                bin_counts[i] += 1
                break

    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t')

        # Header
        writer.writerow(['Gene_Count_Range', 'Disease_Count', 'Percentage'])

        total_diseases = len(results)

        # Data rows
        for label, count in zip(bin_labels, bin_counts):
            percentage = (count / total_diseases * 100) if total_diseases > 0 else 0
            writer.writerow([label, count, f"{percentage:.1f}%"])

    print(f"Gene count distribution written to {output_file}")

def print_summary_stats(results: List[Dict[str, Any]]):
    """Print summary statistics to console"""
    total_diseases = len(results)
    diseases_with_genes = sum(1 for r in results if r['gene_count'] > 0)
    diseases_with_errors = sum(1 for r in results if r.get('enrichment_errors'))

    gene_counts = [r['gene_count'] for r in results if r['gene_count'] > 0]
    avg_genes = sum(gene_counts) / len(gene_counts) if gene_counts else 0

    print(f"\n=== Summary Statistics ===")
    print(f"Total diseases analyzed: {total_diseases}")
    print(f"Diseases with genes: {diseases_with_genes} ({diseases_with_genes/total_diseases*100:.1f}%)")
    print(f"Diseases with no genes: {total_diseases - diseases_with_genes}")
    print(f"Diseases with enrichment errors: {diseases_with_errors}")
    if gene_counts:
        print(f"Gene count range: {min(gene_counts)} to {max(gene_counts)}")
        print(f"Average gene count: {avg_genes:.1f}")

    # Count enrichment results by category
    categories = ['biolink:BiologicalProcess', 'biolink:MolecularActivity', 'biolink:Pathway']
    for category in categories:
        success_count = sum(1 for r in results if r.get('enrichment_results', {}).get(category))
        error_count = sum(1 for r in results if r.get('enrichment_errors', {}).get(category))
        print(f"{category}: {success_count} successful, {error_count} errors")

def main():
    parser = argparse.ArgumentParser(description='Parse disease enrichment analysis results to TSV format')
    parser.add_argument('--input', '-i', default='fast_enrichment_results.jsonl',
                       help='Input JSONL file (default: fast_enrichment_results.jsonl)')
    parser.add_argument('--output-prefix', '-o', default='analysis',
                       help='Output file prefix (default: analysis)')
    parser.add_argument('--min-genes', '-m', type=int, default=1,
                       help='Minimum gene count for enrichment results (default: 1)')

    args = parser.parse_args()

    # Load results from JSONL format
    try:
        results = []
        with open(args.input, 'r') as f:
            for line in f:
                line = line.strip()
                if line:  # Skip empty lines
                    results.append(json.loads(line))
        print(f"Loaded {len(results)} disease analysis results from {args.input}")
    except FileNotFoundError:
        print(f"Error: Input file {args.input} not found")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON in {args.input}: {e}")
        sys.exit(1)

    # Generate output files
    parse_disease_summary(results, f"{args.output_prefix}_disease_summary.tsv")
    parse_enrichment_results(results, f"{args.output_prefix}_enrichment_results.tsv", args.min_genes)
    parse_error_summary(results, f"{args.output_prefix}_error_summary.tsv")
    parse_gene_counts_distribution(results, f"{args.output_prefix}_gene_distribution.tsv")

    # Print summary stats
    print_summary_stats(results)

    print(f"\nTSV files created with prefix '{args.output_prefix}_'")
    print("Files created:")
    print(f"  {args.output_prefix}_disease_summary.tsv - Disease overview")
    print(f"  {args.output_prefix}_enrichment_results.tsv - Enrichment terms with p-values")
    print(f"  {args.output_prefix}_error_summary.tsv - Timeout and error details")
    print(f"  {args.output_prefix}_gene_distribution.tsv - Gene count distribution")

if __name__ == "__main__":
    main()