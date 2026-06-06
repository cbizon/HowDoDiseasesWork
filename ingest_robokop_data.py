#!/usr/bin/env python3
"""
ROBOKOP Knowledge Graph Data Ingestion Script

This script processes the local ROBOKOP graph files (nodes.jsonl and edges.jsonl)
to extract disease-gene relationships directly from the data, bypassing slow API calls.

Process:
1. Parse nodes.jsonl to identify diseases and genes
2. Parse edges.jsonl to find disease-gene connections (both directions)
3. Build disease->gene mapping
4. Sort diseases by gene count
5. Save results for fast enrichment analysis

This replaces the slow API-based gene retrieval phase.
"""

import json
import logging
from collections import defaultdict
from typing import Dict, List, Set, Any
import time
import argparse
import os

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ROBOKOPDataIngester:
    """Class to ingest ROBOKOP graph data and extract disease-gene relationships"""

    def __init__(self, graph_dir: str = None):
        """
        Initialize the data ingester

        Args:
            graph_dir: Path to directory containing nodes.jsonl and edges.jsonl
        """
        self.graph_dir = graph_dir or os.environ.get("KGX_GRAPH_DIR", "/Users/bizon/Projects/ROBOKOP/graph")
        self.nodes_file = f"{self.graph_dir}/nodes.jsonl"
        self.edges_file = f"{self.graph_dir}/edges.jsonl"

        # Data storage
        self.diseases = {}  # curie -> node info
        self.genes = {}     # curie -> node info
        self.disease_gene_edges = defaultdict(set)  # disease_curie -> set of gene_curies

    def parse_nodes(self) -> None:
        """
        Parse nodes.jsonl to identify diseases and genes
        """
        logger.info("Parsing nodes.jsonl to identify diseases and genes...")

        disease_count = 0
        gene_count = 0
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

                    # Check if this is a disease
                    if 'biolink:Disease' in categories:
                        self.diseases[node_id] = {
                            'curie': node_id,
                            'name': node.get('name', 'Unknown'),
                            'description': node.get('description', '')
                        }
                        disease_count += 1

                    # Check if this is a gene
                    elif 'biolink:Gene' in categories:
                        self.genes[node_id] = {
                            'curie': node_id,
                            'name': node.get('name', 'Unknown'),
                            'description': node.get('description', '')
                        }
                        gene_count += 1

                except json.JSONDecodeError as e:
                    logger.warning(f"Error parsing node line: {e}")
                    continue

        logger.info(f"Node parsing complete!")
        logger.info(f"Total nodes processed: {total_nodes:,}")
        logger.info(f"Diseases found: {disease_count:,}")
        logger.info(f"Genes found: {gene_count:,}")

    def parse_edges(self) -> None:
        """
        Parse edges.jsonl to find disease-gene connections in both directions
        """
        logger.info("Parsing edges.jsonl to find disease-gene connections...")

        total_edges = 0
        disease_gene_edges_found = 0

        with open(self.edges_file, 'r') as f:
            for line in f:
                total_edges += 1

                if total_edges % 500000 == 0:
                    logger.info(f"Processed {total_edges:,} edges, found {disease_gene_edges_found:,} disease-gene connections...")

                try:
                    edge = json.loads(line.strip())
                    subject = edge.get('subject')
                    object_id = edge.get('object')

                    # Check for disease -> gene connections
                    if subject in self.diseases and object_id in self.genes:
                        self.disease_gene_edges[subject].add(object_id)
                        disease_gene_edges_found += 1

                    # Check for gene -> disease connections (reverse direction)
                    elif subject in self.genes and object_id in self.diseases:
                        self.disease_gene_edges[object_id].add(subject)
                        disease_gene_edges_found += 1

                except json.JSONDecodeError as e:
                    logger.warning(f"Error parsing edge line: {e}")
                    continue

        logger.info(f"Edge parsing complete!")
        logger.info(f"Total edges processed: {total_edges:,}")
        logger.info(f"Disease-gene connections found: {disease_gene_edges_found:,}")
        logger.info(f"Diseases with genes: {len(self.disease_gene_edges):,}")

    def build_disease_gene_mapping(self) -> List[Dict[str, Any]]:
        """
        Build final disease->gene mapping sorted by gene count

        Returns:
            List of disease dictionaries with gene information, sorted by gene count
        """
        logger.info("Building disease-gene mapping...")

        diseases_with_genes = []

        for disease_curie, disease_info in self.diseases.items():
            gene_curies = list(self.disease_gene_edges.get(disease_curie, set()))
            gene_count = len(gene_curies)

            diseases_with_genes.append({
                "disease": disease_info,
                "gene_count": gene_count,
                "genes": gene_curies
            })

        # Sort by gene count (ascending - fewest genes first)
        diseases_with_genes.sort(key=lambda x: x["gene_count"])

        # Log statistics
        gene_counts = [d["gene_count"] for d in diseases_with_genes]
        diseases_with_no_genes = sum(1 for x in gene_counts if x == 0)
        diseases_with_genes_count = len(diseases_with_genes) - diseases_with_no_genes

        logger.info(f"Disease-gene mapping complete!")
        logger.info(f"Total diseases: {len(diseases_with_genes):,}")
        logger.info(f"Diseases with genes: {diseases_with_genes_count:,}")
        logger.info(f"Diseases with no genes: {diseases_with_no_genes:,}")

        if gene_counts:
            logger.info(f"Gene count range: {min(gene_counts)} to {max(gene_counts)} genes")
            avg_genes = sum(g for g in gene_counts if g > 0) / diseases_with_genes_count if diseases_with_genes_count > 0 else 0
            logger.info(f"Average gene count (for diseases with genes): {avg_genes:.1f}")

        return diseases_with_genes

    def save_results(self, diseases_with_genes: List[Dict[str, Any]], output_file: str = "robokop_disease_genes.json"):
        """
        Save disease-gene mapping to JSON file

        Args:
            diseases_with_genes: List of disease dictionaries with gene information
            output_file: Output filename
        """
        logger.info(f"Saving results to {output_file}...")

        with open(output_file, 'w') as f:
            json.dump(diseases_with_genes, f, indent=2)

        logger.info(f"Results saved successfully!")

    def print_sample_results(self, diseases_with_genes: List[Dict[str, Any]], num_samples: int = 10):
        """
        Print sample results for verification

        Args:
            diseases_with_genes: List of disease dictionaries with gene information
            num_samples: Number of samples to show
        """
        print(f"\n=== Sample Results (First {num_samples} diseases by gene count) ===")
        print(f"{'Rank':<5} {'Gene Count':<12} {'Disease Name':<50} {'CURIE':<20}")
        print("-" * 90)

        for i, result in enumerate(diseases_with_genes[:num_samples], 1):
            disease = result["disease"]
            gene_count = result["gene_count"]
            print(f"{i:<5} {gene_count:<12} {disease['name'][:47]:<50} {disease['curie']:<20}")

        # Show some with high gene counts
        high_gene_diseases = [d for d in diseases_with_genes if d["gene_count"] > 500]
        if high_gene_diseases:
            print(f"\n=== Diseases with >500 genes ===")
            for result in high_gene_diseases[:5]:
                disease = result["disease"]
                gene_count = result["gene_count"]
                print(f"{gene_count:<12} {disease['name']:<50} {disease['curie']:<20}")

    def run_ingestion(self, output_file: str = "robokop_disease_genes.json") -> List[Dict[str, Any]]:
        """
        Run the complete data ingestion process

        Returns:
            List of diseases with gene information, sorted by gene count
        """
        start_time = time.time()
        logger.info("Starting ROBOKOP data ingestion...")

        # Step 1: Parse nodes to identify diseases and genes
        self.parse_nodes()

        # Step 2: Parse edges to find disease-gene connections
        self.parse_edges()

        # Step 3: Build disease-gene mapping
        diseases_with_genes = self.build_disease_gene_mapping()

        # Step 4: Save results
        self.save_results(diseases_with_genes, output_file)

        # Print sample results
        self.print_sample_results(diseases_with_genes)

        elapsed = time.time() - start_time
        logger.info(f"Data ingestion complete in {elapsed:.1f} seconds!")

        return diseases_with_genes


def main():
    """Main function to run the data ingestion"""
    parser = argparse.ArgumentParser(description="Extract disease-gene mappings from KGX nodes/edges.")
    parser.add_argument(
        "--graph-dir",
        default=os.environ.get("KGX_GRAPH_DIR", "/Users/bizon/Projects/ROBOKOP/graph"),
        help="Directory containing nodes.jsonl and edges.jsonl",
    )
    parser.add_argument(
        "--output-file",
        default="robokop_disease_genes.json",
        help="Output JSON file for disease-gene mappings",
    )
    args = parser.parse_args()

    # Initialize ingester
    ingester = ROBOKOPDataIngester(args.graph_dir)

    # Run ingestion process
    diseases_with_genes = ingester.run_ingestion(args.output_file)

    print(f"\nData ingestion complete! {len(diseases_with_genes):,} diseases processed.")
    print(f"Results saved to '{args.output_file}'")
    print("Ready for fast enrichment analysis!")


if __name__ == "__main__":
    main()
