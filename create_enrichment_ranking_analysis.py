#!/usr/bin/env python3
"""
Create Enrichment Ranking Analysis TSV

This script creates a detailed TSV file showing the relationship between
enrichment predictions and ROBOKOP ground truth annotations. For each disease,
it shows which ROBOKOP terms were found in the enrichment results and at what rank.

Output format: Disease | Gene_Count | ROBOKOP_Process | Enrichment_Rank
Where Enrichment_Rank shows the rank of the ROBOKOP process in our enrichment results
(or 'Not_Found' if the term wasn't predicted).
"""

import json
import pandas as pd
from collections import defaultdict
from typing import Dict, List, Set, Optional
import argparse

class EnrichmentRankingAnalyzer:
    """Analyze enrichment rankings vs ROBOKOP ground truth"""

    def __init__(self):
        """Initialize analyzer"""
        self.enrichment_data = []
        self.ground_truth = defaultdict(lambda: defaultdict(set))
        self.disease_gene_counts = {}
        self.term_names = {}  # Map from CURIE to name
        self.categories = {
            'biolink:BiologicalProcess': 'BiologicalProcess',
            'biolink:MolecularActivity': 'MolecularActivity',
            'biolink:Pathway': 'Pathway'
        }

    def load_disease_gene_counts(self, tsv_file: str):
        """Load disease gene counts"""
        print(f"Loading disease gene counts from {tsv_file}...")

        df = pd.read_csv(tsv_file, sep='\t')

        for _, row in df.iterrows():
            self.disease_gene_counts[row['disease_curie']] = row['gene_count']

        print(f"Loaded gene counts for {len(self.disease_gene_counts)} diseases")

    def load_term_names_from_robokop(self, nodes_file: str = '/Users/bizon/Projects/ROBOKOP/graph/nodes.jsonl'):
        """Load ALL term names from ROBOKOP nodes file and create complete ID->name mapping"""
        cache_file = 'robokop_term_names_cache.json'

        # Try to load from cache first
        try:
            with open(cache_file, 'r') as f:
                self.term_names = json.load(f)
                print(f"Loaded {len(self.term_names):,} term names from cache ({cache_file})")
        except FileNotFoundError:
            print(f"Cache not found. Loading complete ID->name mapping from {nodes_file}...")

            count = 0

            with open(nodes_file, 'r') as f:
                for line in f:
                    try:
                        node = json.loads(line.strip())
                        curie = node.get('id', '')
                        name = node.get('name', '')

                        # Store all nodes with valid CURIEs and names
                        if curie and name and curie.strip() and name.strip():
                            self.term_names[curie] = name
                            count += 1

                            if count % 100000 == 0:
                                print(f"  Loaded {count:,} term names...")

                    except (json.JSONDecodeError, KeyError):
                        continue

            print(f"Loaded {count:,} term names from ROBOKOP nodes")

            # Save to cache for future runs
            with open(cache_file, 'w') as f:
                json.dump(self.term_names, f)
                print(f"Saved term names cache to {cache_file}")

        # Check if we got the problematic terms
        test_terms = ['GO:0032502', 'GO:0008152']
        print("Verification of problematic terms:")
        for term in test_terms:
            if term in self.term_names:
                print(f"  ✅ {term} -> '{self.term_names[term]}'")
            else:
                print(f"  ❌ {term} -> NOT FOUND")

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
        """Load ROBOKOP ground truth"""
        print(f"Loading ground truth from {jsonl_file}...")

        count = 0
        with open(jsonl_file, 'r') as f:
            for line in f:
                try:
                    edge = json.loads(line.strip())
                    disease_curie = edge['source_curie']
                    term_curie = edge['target_curie']
                    term_category = edge['target_type']

                    if term_category in self.categories:
                        self.ground_truth[disease_curie][term_category].add(term_curie)
                        count += 1

                        if count % 1000000 == 0:
                            print(f"  Loaded {count:,} ground truth edges...")

                except (json.JSONDecodeError, KeyError):
                    continue

        print(f"Loaded {count:,} ground truth edges for {len(self.ground_truth)} diseases")

    def get_enrichment_rank(self, disease_curie: str, term_curie: str, category: str) -> Optional[int]:
        """Get the rank of a term in enrichment results for a disease"""

        # Find the enrichment result for this disease
        for result in self.enrichment_data:
            if result['disease']['curie'] == disease_curie:
                enrichments = result.get('enrichment_results', {}).get(category, [])

                if not enrichments:
                    return None

                # Sort by p-value (ascending - lower p-values are better)
                enrichments.sort(key=lambda x: x['p_value'])

                # Find the rank of this term (1-indexed)
                for rank, enrichment in enumerate(enrichments, 1):
                    if enrichment['curie'] == term_curie:
                        return rank

                return None  # Term not found in enrichment results

        return None  # Disease not found in enrichment results

    def get_term_name(self, term_curie: str) -> str:
        """Get the name of a term using the comprehensive ID->name mapping"""

        # Use the complete ROBOKOP nodes mapping
        if term_curie in self.term_names:
            return self.term_names[term_curie]

        # If not found, return the CURIE
        return term_curie

    def create_ranking_analysis(self, category: str) -> List[Dict]:
        """Create detailed ranking analysis for a category"""

        category_short = self.categories[category]
        print(f"Creating ranking analysis for {category_short}...")

        analysis_data = []

        # Process each disease that has both enrichment results and ground truth
        diseases_processed = 0

        for result in self.enrichment_data:
            disease_curie = result['disease']['curie']
            disease_name = result['disease']['name']

            # Skip if no ground truth for this disease and category
            if (disease_curie not in self.ground_truth or
                category not in self.ground_truth[disease_curie] or
                not self.ground_truth[disease_curie][category]):
                continue

            # Get gene count
            gene_count = self.disease_gene_counts.get(disease_curie, 0)

            # Get ROBOKOP terms for this disease and category
            robokop_terms = self.ground_truth[disease_curie][category]

            # For each ROBOKOP term, find its enrichment rank
            for term_curie in robokop_terms:
                enrichment_rank = self.get_enrichment_rank(disease_curie, term_curie, category)
                term_name = self.get_term_name(term_curie)

                analysis_data.append({
                    'Disease_CURIE': disease_curie,
                    'Disease_Name': disease_name,
                    'Gene_Count': gene_count,
                    'ROBOKOP_Process_CURIE': term_curie,
                    'ROBOKOP_Process_Name': term_name,
                    'Enrichment_Rank': enrichment_rank if enrichment_rank is not None else 'Not_Found'
                })

            diseases_processed += 1
            if diseases_processed % 100 == 0:
                print(f"  Processed {diseases_processed} diseases...")

        print(f"Generated {len(analysis_data)} disease-term pairs for {category_short}")
        return analysis_data

    def save_analysis_tsv(self, analysis_data: List[Dict], output_file: str):
        """Save analysis data to TSV file"""
        print(f"Saving analysis to {output_file}...")

        df = pd.DataFrame(analysis_data)

        # Sort by disease name, then by enrichment rank (with 'Not_Found' at the end)
        def sort_key(row):
            disease_name = row['Disease_Name']
            rank = row['Enrichment_Rank']
            if rank == 'Not_Found':
                return (disease_name, float('inf'))
            else:
                return (disease_name, int(rank))

        df_sorted = df.iloc[df.apply(sort_key, axis=1).argsort()]

        df_sorted.to_csv(output_file, sep='\t', index=False)
        print(f"Saved {len(df_sorted)} rows to {output_file}")

    def print_summary_stats(self, analysis_data: List[Dict], category: str):
        """Print summary statistics for the analysis"""

        category_short = self.categories[category]
        print(f"\n=== {category_short} Ranking Analysis Summary ===")

        total_pairs = len(analysis_data)
        found_pairs = [d for d in analysis_data if d['Enrichment_Rank'] != 'Not_Found']
        not_found_pairs = total_pairs - len(found_pairs)

        print(f"Total disease-term pairs: {total_pairs}")
        print(f"Found in enrichment: {len(found_pairs)} ({len(found_pairs)/total_pairs*100:.1f}%)")
        print(f"Not found in enrichment: {not_found_pairs} ({not_found_pairs/total_pairs*100:.1f}%)")

        if found_pairs:
            ranks = [int(d['Enrichment_Rank']) for d in found_pairs]
            print(f"\nRanking statistics (found terms only):")
            print(f"  Mean rank: {sum(ranks)/len(ranks):.1f}")
            print(f"  Median rank: {sorted(ranks)[len(ranks)//2]}")
            print(f"  Min rank: {min(ranks)}")
            print(f"  Max rank: {max(ranks)}")

            # Rank distribution
            top_10 = sum(1 for r in ranks if r <= 10)
            top_50 = sum(1 for r in ranks if r <= 50)
            top_100 = sum(1 for r in ranks if r <= 100)

            print(f"\nRank distribution:")
            print(f"  Top 10: {top_10} ({top_10/len(ranks)*100:.1f}%)")
            print(f"  Top 50: {top_50} ({top_50/len(ranks)*100:.1f}%)")
            print(f"  Top 100: {top_100} ({top_100/len(ranks)*100:.1f}%)")

        # Gene count analysis
        diseases = set()
        gene_counts = []
        for d in analysis_data:
            disease_key = (d['Disease_CURIE'], d['Gene_Count'])
            if disease_key not in diseases:
                diseases.add(disease_key)
                gene_counts.append(d['Gene_Count'])

        if gene_counts:
            print(f"\nGene count distribution ({len(gene_counts)} diseases):")
            print(f"  Mean gene count: {sum(gene_counts)/len(gene_counts):.1f}")
            print(f"  Median gene count: {sorted(gene_counts)[len(gene_counts)//2]}")
            print(f"  Min gene count: {min(gene_counts)}")
            print(f"  Max gene count: {max(gene_counts)}")

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Create enrichment ranking analysis TSV')
    parser.add_argument('--enrichment', '-e', default='fast_enrichment_results.jsonl',
                       help='Enrichment results JSONL file')
    parser.add_argument('--ground-truth', '-g', default='robokop_disease_term_edges_with_subclass.jsonl',
                       help='Ground truth JSONL file')
    parser.add_argument('--gene-counts', '-c', default='disease_gene_counts.tsv',
                       help='Disease gene counts TSV file')
    parser.add_argument('--category', default='biolink:BiologicalProcess',
                       choices=['biolink:BiologicalProcess', 'biolink:MolecularActivity', 'biolink:Pathway'],
                       help='Category to analyze (default: biolink:BiologicalProcess)')
    parser.add_argument('--output', '-o', default=None,
                       help='Output TSV filename (default: auto-generated)')

    args = parser.parse_args()

    # Generate output filename if not provided
    if args.output is None:
        category_short = args.category.replace('biolink:', '').lower()
        args.output = f'enrichment_ranking_analysis_{category_short}.tsv'

    # Initialize analyzer
    analyzer = EnrichmentRankingAnalyzer()
    analyzer.load_disease_gene_counts(args.gene_counts)
    analyzer.load_term_names_from_robokop()  # Load term names from ROBOKOP nodes
    analyzer.load_enrichment_results(args.enrichment)
    analyzer.load_ground_truth(args.ground_truth)

    # Create analysis
    analysis_data = analyzer.create_ranking_analysis(args.category)

    # Save results
    analyzer.save_analysis_tsv(analysis_data, args.output)
    analyzer.print_summary_stats(analysis_data, args.category)

    print(f"\nAnalysis complete! Results saved to {args.output}")

if __name__ == "__main__":
    main()