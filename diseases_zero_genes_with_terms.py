#!/usr/bin/env python3
"""
Create TSV file of diseases with 0 genes but ROBOKOP terms

This script finds diseases that have 0 genes in our analysis but have
biological terms (processes, activities, pathways) in ROBOKOP knowledge graph.
Outputs a detailed TSV file showing these diseases and their associated terms.
"""

import pandas as pd
import json
from collections import defaultdict

def main():
    print('Loading disease_gene_counts.tsv...')
    df = pd.read_csv('disease_gene_counts.tsv', sep='\t')
    zero_gene_diseases = df[df['gene_count'] == 0]
    print(f'Found {len(zero_gene_diseases)} diseases with 0 genes')

    # Get the set of disease curies with 0 genes
    zero_gene_curies = set(zero_gene_diseases['disease_curie'].tolist())
    zero_gene_names = dict(zip(zero_gene_diseases['disease_curie'], zero_gene_diseases['disease_name']))

    print('Loading ROBOKOP disease-term edges...')

    # Collect all terms for each disease
    disease_terms = defaultdict(list)
    disease_names_robokop = {}

    with open('robokop_disease_term_edges.jsonl', 'r') as f:
        for line in f:
            edge = json.loads(line.strip())
            disease_curie = edge['source_curie']
            disease_name = edge['source_name']
            term_curie = edge['target_curie']
            term_name = edge['target_name']
            term_type = edge['target_type']
            predicate = edge['predicate']

            disease_names_robokop[disease_curie] = disease_name
            disease_terms[disease_curie].append({
                'term_curie': term_curie,
                'term_name': term_name,
                'term_type': term_type,
                'predicate': predicate
            })

    # Find overlap: diseases with 0 genes BUT with terms in ROBOKOP
    zero_genes_with_terms = zero_gene_curies.intersection(set(disease_terms.keys()))

    print(f'Found {len(zero_genes_with_terms)} diseases with 0 genes that have ROBOKOP terms')

    # Prepare data for TSV output
    tsv_data = []

    for disease_curie in zero_genes_with_terms:
        # Use name from gene counts file, fall back to ROBOKOP name
        disease_name = zero_gene_names.get(disease_curie, disease_names_robokop.get(disease_curie, 'Unknown'))
        terms = disease_terms[disease_curie]

        # Sort terms by type and name for consistent output
        terms.sort(key=lambda x: (x['term_type'], x['term_name']))

        # Count terms by category
        bp_count = sum(1 for t in terms if t['term_type'] == 'biolink:BiologicalProcess')
        ma_count = sum(1 for t in terms if t['term_type'] == 'biolink:MolecularActivity')
        pw_count = sum(1 for t in terms if t['term_type'] == 'biolink:Pathway')

        # Create one row per disease with all terms concatenated
        term_details = []
        for term in terms:
            term_type_short = term['term_type'].replace('biolink:', '')
            predicate_short = term['predicate'].replace('biolink:', '')
            term_details.append(f"{term['term_curie']}|{term['term_name']}|{term_type_short}|{predicate_short}")

        tsv_data.append({
            'disease_curie': disease_curie,
            'disease_name': disease_name,
            'total_terms': len(terms),
            'biological_process_count': bp_count,
            'molecular_activity_count': ma_count,
            'pathway_count': pw_count,
            'all_terms': '; '.join(term_details)
        })

    # Sort by total terms (descending) then by disease name
    tsv_data.sort(key=lambda x: (-x['total_terms'], x['disease_name']))

    # Write to TSV file
    output_file = 'diseases_zero_genes_with_robokop_terms.tsv'
    print(f'Writing results to {output_file}...')

    with open(output_file, 'w') as f:
        # Header
        f.write('disease_curie\tdisease_name\ttotal_terms\tbiological_process_count\t'
                'molecular_activity_count\tpathway_count\tall_terms\n')

        # Data rows
        for row in tsv_data:
            f.write(f"{row['disease_curie']}\t{row['disease_name']}\t{row['total_terms']}\t"
                   f"{row['biological_process_count']}\t{row['molecular_activity_count']}\t"
                   f"{row['pathway_count']}\t{row['all_terms']}\n")

    print(f'Results saved to {output_file}')

    # Print summary statistics
    print(f'\nSummary:')
    print(f'Total diseases with 0 genes but ROBOKOP terms: {len(tsv_data)}')

    if tsv_data:
        total_terms = sum(row['total_terms'] for row in tsv_data)
        avg_terms = total_terms / len(tsv_data)
        max_terms = max(row['total_terms'] for row in tsv_data)

        print(f'Average terms per disease: {avg_terms:.1f}')
        print(f'Maximum terms for one disease: {max_terms}')

        # Show top 10
        print(f'\nTop 10 diseases by term count:')
        for i, row in enumerate(tsv_data[:10]):
            print(f'{i+1:2d}. {row["total_terms"]:2d} terms: {row["disease_name"][:60]}')

if __name__ == "__main__":
    main()