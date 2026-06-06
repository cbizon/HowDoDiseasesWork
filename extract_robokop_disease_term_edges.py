#!/usr/bin/env python3
"""
Extract Disease-Term Edges from ROBOKOP Knowledge Graph
This script processes the original ROBOKOP graph files (nodes.jsonl and edges.jsonl)
to extract ALL direct edges between diseases and BiologicalProcess/MolecularActivity/Pathway terms.

Process:
1. Parse nodes.jsonl to identify diseases and BP/MA/Pathway terms
2. Parse edges.jsonl to find direct disease-term connections (both directions)
3. Build comprehensive disease-term edge list from the knowledge graph
4. Export results in JSONL and TSV formats

This gives us the full set of disease-mechanism relationships from ROBOKOP KG.
"""

import json
import csv
import logging
import time
from typing import Dict, List, Any, Set
from collections import defaultdict
import argparse
import os

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ROBOKOPDiseaseTermExtractor:
    """Class to extract disease-term edges from ROBOKOP graph files"""

    def __init__(self, graph_dir: str = None):
        """
        Initialize the extractor

        Args:
            graph_dir: Directory containing nodes.jsonl and edges.jsonl files
        """
        self.graph_dir = graph_dir or os.environ.get("KGX_GRAPH_DIR", "/Users/bizon/Projects/ROBOKOP/graph")
        self.nodes_file = f"{self.graph_dir}/nodes.jsonl"
        self.edges_file = f"{self.graph_dir}/edges.jsonl"

        # Data storage
        self.diseases = {}  # curie -> node info
        self.terms = {}     # curie -> node info (BP/MA/Pathway)
        self.disease_subclass_edges = []  # List of disease subclass relationships
        self.direct_edges = []     # List of direct disease-term edges
        self.edges = []     # List of all disease-term edges (direct + inferred)

        # Target categories
        self.target_categories = {
            'biolink:BiologicalProcess',
            'biolink:MolecularActivity',
            'biolink:Pathway'
        }

    def parse_nodes(self) -> None:
        """
        Parse nodes.jsonl to identify diseases and target terms
        """
        logger.info("Parsing nodes.jsonl to identify diseases and terms...")

        disease_count = 0
        term_count = 0
        total_nodes = 0

        with open(self.nodes_file, 'r') as f:
            for line in f:
                total_nodes += 1

                if total_nodes % 100000 == 0:
                    logger.info(f"Processed {total_nodes:,} nodes...")

                try:
                    node = json.loads(line.strip())
                    node_id = node.get('id')
                    categories = node.get('category', [])

                    if not node_id:
                        continue

                    # Check for diseases
                    if 'biolink:Disease' in categories:
                        self.diseases[node_id] = {
                            'curie': node_id,
                            'name': node.get('name', 'Unknown'),
                            'description': node.get('description', ''),
                            'categories': categories
                        }
                        disease_count += 1

                    # Check for target terms (BP/MA/Pathway)
                    elif any(cat in self.target_categories for cat in categories):
                        # Get the primary category that matches our targets
                        primary_category = next((cat for cat in categories if cat in self.target_categories), None)
                        if primary_category:
                            self.terms[node_id] = {
                                'curie': node_id,
                                'name': node.get('name', 'Unknown'),
                                'description': node.get('description', ''),
                                'category': primary_category,
                                'all_categories': categories
                            }
                            term_count += 1

                except json.JSONDecodeError as e:
                    logger.warning(f"Error parsing node line: {e}")
                    continue

        logger.info(f"Node parsing complete!")
        logger.info(f"Total nodes processed: {total_nodes:,}")
        logger.info(f"Diseases found: {disease_count:,}")
        logger.info(f"Terms found: {term_count:,}")
        logger.info(f"  BiologicalProcess: {sum(1 for t in self.terms.values() if t['category'] == 'biolink:BiologicalProcess'):,}")
        logger.info(f"  MolecularActivity: {sum(1 for t in self.terms.values() if t['category'] == 'biolink:MolecularActivity'):,}")
        logger.info(f"  Pathway: {sum(1 for t in self.terms.values() if t['category'] == 'biolink:Pathway'):,}")

    def parse_edges(self) -> None:
        """
        Parse edges.jsonl to find disease subclass relationships and disease-term connections
        """
        logger.info("Parsing edges.jsonl for subclass relationships and disease-term connections...")

        total_edges = 0
        disease_term_edges_found = 0
        subclass_edges_found = 0

        with open(self.edges_file, 'r') as f:
            for line in f:
                total_edges += 1

                if total_edges % 500000 == 0:
                    logger.info(f"Processed {total_edges:,} edges, found {subclass_edges_found:,} subclass + {disease_term_edges_found:,} disease-term connections...")

                try:
                    edge = json.loads(line.strip())
                    subject = edge.get('subject')
                    object_id = edge.get('object')
                    predicate = edge.get('predicate', 'related_to')

                    # Check for disease -> disease subclass relationships
                    if (subject in self.diseases and object_id in self.diseases and
                        predicate in ['biolink:subclass_of', 'rdfs:subClassOf']):
                        subclass_data = {
                            'child_curie': subject,
                            'child_name': self.diseases[subject]['name'],
                            'parent_curie': object_id,
                            'parent_name': self.diseases[object_id]['name'],
                            'predicate': predicate
                        }
                        self.disease_subclass_edges.append(subclass_data)
                        subclass_edges_found += 1

                    # Check for disease -> term connections
                    elif subject in self.diseases and object_id in self.terms:
                        edge_data = {
                            'source_curie': subject,
                            'source_name': self.diseases[subject]['name'],
                            'source_type': 'biolink:Disease',
                            'target_curie': object_id,
                            'target_name': self.terms[object_id]['name'],
                            'target_type': self.terms[object_id]['category'],
                            'predicate': predicate,
                            'direction': 'disease_to_term',
                            'inference_type': 'direct'
                        }
                        self.direct_edges.append(edge_data)
                        disease_term_edges_found += 1

                    # Check for term -> disease connections (reverse direction)
                    elif subject in self.terms and object_id in self.diseases:
                        edge_data = {
                            'source_curie': object_id,  # Flip to make disease the source
                            'source_name': self.diseases[object_id]['name'],
                            'source_type': 'biolink:Disease',
                            'target_curie': subject,
                            'target_name': self.terms[subject]['name'],
                            'target_type': self.terms[subject]['category'],
                            'predicate': predicate,
                            'direction': 'term_to_disease',
                            'inference_type': 'direct'
                        }
                        self.direct_edges.append(edge_data)
                        disease_term_edges_found += 1

                except json.JSONDecodeError as e:
                    logger.warning(f"Error parsing edge line: {e}")
                    continue

        logger.info(f"Edge parsing complete!")
        logger.info(f"Total edges processed: {total_edges:,}")
        logger.info(f"Disease subclass relationships found: {subclass_edges_found:,}")
        logger.info(f"Direct disease-term connections found: {disease_term_edges_found:,}")

    def build_subclass_mapping(self) -> Dict[str, Set[str]]:
        """
        Build a mapping from each parent disease to all its subclass diseases
        Since ROBOKOP already computed transitive closure, we just need direct mappings

        Returns:
            Dictionary mapping parent_disease_curie -> set of subclass_disease_curies
        """
        logger.info("Building disease subclass mapping...")

        # Create mapping: parent -> set of all subclasses (children)
        parent_to_subclasses = defaultdict(set)
        for edge in self.disease_subclass_edges:
            child = edge['child_curie']
            parent = edge['parent_curie']
            parent_to_subclasses[parent].add(child)

        # Log statistics
        parents_with_subclasses = len(parent_to_subclasses)
        total_subclass_relationships = sum(len(subclasses) for subclasses in parent_to_subclasses.values())

        logger.info(f"Subclass mapping built: {parents_with_subclasses:,} parent diseases")
        logger.info(f"Total subclass relationships: {total_subclass_relationships:,}")

        return parent_to_subclasses

    def infer_subclass_edges(self) -> None:
        """
        Infer disease-term edges from subclass hierarchy
        If parent_disease related_to term AND child_disease subclass_of parent_disease,
        then infer child_disease related_to term

        Since ROBOKOP already computed transitive closure, we just propagate down the hierarchy
        """
        logger.info("Inferring disease-term edges from subclass hierarchy...")

        # Build subclass mapping: parent -> set of all subclasses
        parent_to_subclasses = self.build_subclass_mapping()

        # Start with direct edges
        self.edges = self.direct_edges.copy()
        inferred_edges_count = 0

        # Track existing edges to avoid duplicates
        existing_edges = set((e['source_curie'], e['target_curie']) for e in self.direct_edges)

        # For each direct disease-term edge, propagate to all subclasses
        for direct_edge in self.direct_edges:
            parent_disease = direct_edge['source_curie']

            # Get all subclasses of this disease
            subclasses = parent_to_subclasses.get(parent_disease, set())

            # For each subclass, create inferred edge
            for child_curie in subclasses:
                if child_curie in self.diseases:  # Make sure it's a valid disease
                    child_name = self.diseases[child_curie]['name']

                    # Check if this edge already exists
                    edge_key = (child_curie, direct_edge['target_curie'])
                    if edge_key not in existing_edges:
                        # Create inferred edge
                        inferred_edge = {
                            'source_curie': child_curie,
                            'source_name': child_name,
                            'source_type': 'biolink:Disease',
                            'target_curie': direct_edge['target_curie'],
                            'target_name': direct_edge['target_name'],
                            'target_type': direct_edge['target_type'],
                            'predicate': direct_edge['predicate'],
                            'direction': direct_edge['direction'],
                            'inference_type': f'inferred_from_parent',
                            'original_disease': parent_disease,
                            'original_disease_name': direct_edge['source_name']
                        }

                        self.edges.append(inferred_edge)
                        existing_edges.add(edge_key)
                        inferred_edges_count += 1

        total_edges = len(self.edges)
        direct_edges = len(self.direct_edges)

        logger.info(f"Inference complete!")
        logger.info(f"Direct edges: {direct_edges:,}")
        logger.info(f"Inferred edges: {inferred_edges_count:,}")
        logger.info(f"Total edges (direct + inferred): {total_edges:,}")

    def save_edges_jsonl(self, output_file: str) -> None:
        """
        Save edges in JSONL format

        Args:
            output_file: Output JSONL filename
        """
        logger.info(f"Saving {len(self.edges):,} edges to {output_file}...")

        with open(output_file, 'w') as f:
            for edge in self.edges:
                f.write(json.dumps(edge) + '\n')

        logger.info(f"Edges saved to {output_file}")

    def save_edges_tsv(self, output_file: str) -> None:
        """
        Save edges in TSV format for easy analysis

        Args:
            output_file: Output TSV filename
        """
        logger.info(f"Saving edges to TSV format: {output_file}...")

        with open(output_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f, delimiter='\t')

            # Header
            writer.writerow([
                'source_curie', 'source_name', 'source_type',
                'target_curie', 'target_name', 'target_type',
                'predicate', 'direction'
            ])

            # Sort edges by source disease for better organization
            sorted_edges = sorted(self.edges, key=lambda x: (x['source_name'], x['target_type'], x['target_name']))

            # Data rows
            for edge in sorted_edges:
                writer.writerow([
                    edge['source_curie'], edge['source_name'], edge['source_type'],
                    edge['target_curie'], edge['target_name'], edge['target_type'],
                    edge['predicate'], edge['direction']
                ])

        logger.info(f"Edges saved to TSV format: {output_file}")

    def generate_summary_stats(self, output_file: str) -> None:
        """
        Generate summary statistics about the extracted edges

        Args:
            output_file: Output filename for summary stats
        """
        logger.info(f"Generating summary statistics...")

        # Calculate statistics by category
        category_stats = defaultdict(lambda: {
            'edge_count': 0,
            'unique_diseases': set(),
            'unique_terms': set(),
            'predicates': defaultdict(int)
        })

        for edge in self.edges:
            category = edge['target_type']
            category_stats[category]['edge_count'] += 1
            category_stats[category]['unique_diseases'].add(edge['source_curie'])
            category_stats[category]['unique_terms'].add(edge['target_curie'])
            category_stats[category]['predicates'][edge['predicate']] += 1

        # Write summary
        with open(output_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f, delimiter='\t')

            # Header
            writer.writerow([
                'Category', 'Edge_Count', 'Unique_Diseases', 'Unique_Terms',
                'Avg_Terms_Per_Disease', 'Top_Predicate', 'Top_Predicate_Count'
            ])

            # Data rows
            for category in sorted(self.target_categories):
                if category in category_stats:
                    stats = category_stats[category]
                    unique_diseases = len(stats['unique_diseases'])
                    avg_terms = stats['edge_count'] / unique_diseases if unique_diseases > 0 else 0

                    # Find most common predicate
                    top_predicate = max(stats['predicates'].items(), key=lambda x: x[1]) if stats['predicates'] else ('N/A', 0)

                    writer.writerow([
                        category,
                        stats['edge_count'],
                        unique_diseases,
                        len(stats['unique_terms']),
                        f"{avg_terms:.1f}",
                        top_predicate[0],
                        top_predicate[1]
                    ])

        logger.info(f"Summary statistics saved to {output_file}")

    def print_summary(self) -> None:
        """Print summary statistics to console"""
        print(f"\\n=== ROBOKOP Disease-Term Edge Extraction Summary ===")
        print(f"Total edges extracted: {len(self.edges):,}")

        # Count unique diseases and terms
        unique_diseases = len(set(e['source_curie'] for e in self.edges))
        unique_terms = len(set(e['target_curie'] for e in self.edges))
        print(f"Unique diseases with term connections: {unique_diseases:,}")
        print(f"Unique terms connected to diseases: {unique_terms:,}")

        # Category breakdown
        category_counts = defaultdict(int)
        for edge in self.edges:
            category_counts[edge['target_type']] += 1

        print(f"\\nEdges by category:")
        for category in sorted(self.target_categories):
            count = category_counts[category]
            pct = (count / len(self.edges)) * 100 if self.edges else 0
            print(f"  {category}: {count:,} ({pct:.1f}%)")

        # Direction breakdown
        direction_counts = defaultdict(int)
        for edge in self.edges:
            direction_counts[edge['direction']] += 1

        print(f"\\nEdges by direction:")
        for direction, count in direction_counts.items():
            pct = (count / len(self.edges)) * 100 if self.edges else 0
            print(f"  {direction}: {count:,} ({pct:.1f}%)")

        # Top predicates
        predicate_counts = defaultdict(int)
        for edge in self.edges:
            predicate_counts[edge['predicate']] += 1

        top_predicates = sorted(predicate_counts.items(), key=lambda x: x[1], reverse=True)[:10]
        print(f"\\nTop 10 predicates:")
        for predicate, count in top_predicates:
            pct = (count / len(self.edges)) * 100 if self.edges else 0
            print(f"  {predicate}: {count:,} ({pct:.1f}%)")

        # Sample edges
        print(f"\\nSample edges:")
        for i, edge in enumerate(self.edges[:5]):
            print(f"  {edge['source_name']} --{edge['predicate']}--> {edge['target_name']} ({edge['target_type']})")

    def run_extraction(self) -> None:
        """
        Run the complete edge extraction process
        """
        start_time = time.time()
        logger.info("Starting ROBOKOP disease-term edge extraction...")

        # Step 1: Parse nodes to identify diseases and terms
        self.parse_nodes()

        # Step 2: Parse edges to find disease-term connections and subclass relationships
        self.parse_edges()

        # Step 3: Infer additional edges from subclass hierarchy
        self.infer_subclass_edges()

        elapsed = time.time() - start_time
        logger.info(f"Extraction complete in {elapsed:.1f} seconds!")

def main():
    """Main function to extract disease-term edges from ROBOKOP"""

    parser = argparse.ArgumentParser(description='Extract disease-term edges from ROBOKOP Knowledge Graph')
    parser.add_argument('--graph-dir', '-g', default=os.environ.get('KGX_GRAPH_DIR', '/Users/bizon/Projects/ROBOKOP/graph'),
                       help='Directory containing ROBOKOP graph files (default: /Users/bizon/Projects/ROBOKOP/graph)')
    parser.add_argument('--output-prefix', '-o', default='robokop_disease_term_edges',
                       help='Output file prefix (default: robokop_disease_term_edges)')

    args = parser.parse_args()

    # Initialize extractor
    extractor = ROBOKOPDiseaseTermExtractor(args.graph_dir)

    # Run extraction
    extractor.run_extraction()

    # Save results in multiple formats
    extractor.save_edges_jsonl(f"{args.output_prefix}.jsonl")
    extractor.save_edges_tsv(f"{args.output_prefix}.tsv")
    extractor.generate_summary_stats(f"{args.output_prefix}_summary.tsv")

    # Print summary
    extractor.print_summary()

    print(f"\\nROBOKOP disease-term edge extraction complete!")
    print(f"Files created:")
    print(f"  {args.output_prefix}.jsonl - All edges in JSONL format")
    print(f"  {args.output_prefix}.tsv - All edges in TSV format (sorted by disease)")
    print(f"  {args.output_prefix}_summary.tsv - Summary statistics by category")

if __name__ == "__main__":
    main()
