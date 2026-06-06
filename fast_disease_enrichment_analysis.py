#!/usr/bin/env python3
"""
Fast Disease Enrichment Analysis Pipeline

This script uses pre-processed local disease-gene data (from ingest_robokop_data.py)
to perform enrichment analysis much faster than API-based gene retrieval.

Pipeline:
1. Load pre-processed disease-gene mappings from JSON file
2. Perform enrichment analysis for BiologicalProcess, MolecularActivity, and Pathway
3. Save results with incremental saving and error tracking

This completely bypasses the slow API gene retrieval phase.
"""

import requests
import json
import time
import os
import argparse
from typing import List, Dict, Any
from collections import defaultdict
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# API URLs
ANSWER_COALESCE_URL = "https://answercoalesce.renci.org/query"

# Target categories for enrichment analysis
ENRICHMENT_CATEGORIES = [
    "biolink:BiologicalProcess",
    "biolink:MolecularActivity",
    "biolink:Pathway"
]

class FastDiseaseEnrichmentAnalyzer:
    """Fast disease enrichment analysis using pre-processed local data"""

    def __init__(
        self,
        delay_between_requests: float = 0.1,
        results_file: str = "fast_enrichment_results.jsonl",
        answer_coalesce_url: str = None,
    ):
        """
        Initialize the analyzer

        Args:
            delay_between_requests: Delay in seconds between API requests
            results_file: File path for saving results incrementally
        """
        self.delay = delay_between_requests
        self.session = requests.Session()
        self.results_file = results_file
        self.answer_coalesce_url = answer_coalesce_url or os.environ.get(
            "ANSWER_COALESCE_URL", ANSWER_COALESCE_URL
        )
        self.backup_frequency = 10  # Save backup every N diseases

    def load_disease_gene_data(self, data_file: str = "robokop_disease_genes.json") -> List[Dict[str, Any]]:
        """
        Load pre-processed disease-gene mappings

        Args:
            data_file: JSON file with disease-gene mappings from ingest_robokop_data.py

        Returns:
            List of diseases with gene information, already sorted by gene count
        """
        logger.info(f"Loading disease-gene data from {data_file}...")

        try:
            with open(data_file, 'r') as f:
                diseases_with_genes = json.load(f)

            logger.info(f"Loaded {len(diseases_with_genes):,} diseases with gene mappings")

            # Log statistics
            gene_counts = [d["gene_count"] for d in diseases_with_genes]
            diseases_with_no_genes = sum(1 for x in gene_counts if x == 0)
            diseases_with_genes_count = len(diseases_with_genes) - diseases_with_no_genes

            logger.info(f"Diseases with genes: {diseases_with_genes_count:,}")
            logger.info(f"Diseases with no genes: {diseases_with_no_genes:,}")

            if gene_counts:
                logger.info(f"Gene count range: {min(gene_counts)} to {max(gene_counts)} genes")
                avg_genes = sum(g for g in gene_counts if g > 0) / diseases_with_genes_count if diseases_with_genes_count > 0 else 0
                logger.info(f"Average gene count: {avg_genes:.1f}")

            return diseases_with_genes

        except FileNotFoundError:
            logger.error(f"Disease-gene data file {data_file} not found!")
            logger.error("Please run ingest_robokop_data.py first to generate the data file.")
            return []
        except json.JSONDecodeError as e:
            logger.error(f"Error parsing {data_file}: {e}")
            return []

    def load_existing_results(self) -> List[Dict[str, Any]]:
        """Load existing results from JSONL file if it exists"""
        if os.path.exists(self.results_file):
            try:
                results = []
                with open(self.results_file, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if line:  # Skip empty lines
                            results.append(json.loads(line))
                logger.info(f"Loaded {len(results)} existing results from {self.results_file}")
                return results
            except Exception as e:
                logger.warning(f"Error loading existing JSONL results: {e}")
                return []
        return []

    def get_processed_disease_curies(self, results: List[Dict[str, Any]]) -> set:
        """Get set of disease CURIEs that have already been processed"""
        return {result["disease"]["curie"] for result in results}

    def generate_trapi_query(self, input_node_type: str, output_node_type: str,
                           input_curies: List[str], predicate: str,
                           input_is_subject: bool = True) -> Dict[str, Any]:
        """Generate a TRAPI query for enrichment analysis"""
        envelope = {
            "message": {
                "query_graph": {
                    "nodes": {
                        "input": {},
                        "output": {}
                    },
                    "edges": {
                        "edge_0": {}
                    }
                }
            }
        }

        # Configure input node
        input_node = envelope["message"]["query_graph"]["nodes"]["input"]
        input_node["categories"] = [input_node_type]
        input_node["ids"] = ["uuid:1"]
        input_node["member_ids"] = input_curies
        input_node["set_interpretation"] = "MANY"

        # Configure output node
        output_node = envelope["message"]["query_graph"]["nodes"]["output"]
        output_node["categories"] = [output_node_type]

        # Configure edge
        query_edge = envelope["message"]["query_graph"]["edges"]["edge_0"]
        if input_is_subject:
            query_edge["subject"] = "input"
            query_edge["object"] = "output"
        else:
            query_edge["subject"] = "output"
            query_edge["object"] = "input"
        query_edge["predicates"] = [predicate]

        return envelope

    def perform_enrichment_analysis(self, genes: List[str], category: str) -> List[Dict[str, Any]]:
        """Perform enrichment analysis for a list of genes against a specific category"""
        if not genes:
            return []

        logger.debug(f"Performing enrichment analysis: {len(genes)} genes -> {category}")

        try:
            # Generate TRAPI query
            query = self.generate_trapi_query(
                input_node_type="biolink:Gene",
                output_node_type=category,
                input_curies=genes,
                predicate="biolink:related_to",
                input_is_subject=True
            )

            # Send query to AnswerCoalesce
            response = self.session.post(
                self.answer_coalesce_url,
                json=query,
                headers={"Content-Type": "application/json"},
                timeout=60  # 60 second timeout for enrichment queries
            )
            response.raise_for_status()

            data = response.json()

            # Parse results
            enriched_terms = []

            if "message" in data and "results" in data["message"]:
                knowledge_graph = data["message"].get("knowledge_graph", {})
                nodes = knowledge_graph.get("nodes", {})
                edges = knowledge_graph.get("edges", {})

                for result in data["message"]["results"]:
                    try:
                        # Get output node binding
                        output_binding = result["node_bindings"]["output"][0]["id"]
                        node_info = nodes.get(output_binding, {})

                        # Get edge information and p-value
                        edge_id = result["analyses"][0]["edge_bindings"]["edge_0"][0]["id"]
                        edge_info = edges.get(edge_id, {})

                        p_value = None
                        for attr in edge_info.get("attributes", []):
                            if attr.get("attribute_type_id") == "biolink:p_value":
                                p_value = attr.get("value")
                                break

                        enriched_terms.append({
                            "curie": output_binding,
                            "name": node_info.get("name", "Unknown"),
                            "category": category,
                            "p_value": p_value
                        })

                    except (KeyError, IndexError) as e:
                        logger.warning(f"Error parsing result: {e}")
                        continue

            logger.debug(f"Found {len(enriched_terms)} enriched terms for {category}")

            # Add delay to be respectful to the API
            time.sleep(self.delay)

            return enriched_terms

        except Exception as e:
            logger.error(f"Error in enrichment analysis for {category}: {e}")
            raise  # Re-raise to be caught by caller for error tracking

    def perform_disease_enrichment(self, disease_with_genes: Dict[str, Any]) -> Dict[str, Any]:
        """Perform enrichment analysis for a disease that already has genes"""
        disease = disease_with_genes["disease"]
        genes = disease_with_genes["genes"]
        gene_count = disease_with_genes["gene_count"]

        disease_name = disease["name"]
        disease_curie = disease["curie"]

        logger.info(f"Enrichment analysis: {disease_name} ({gene_count} genes)")

        if not genes:
            logger.warning(f"No genes found for disease {disease_name}")
            return {
                "disease": disease,
                "gene_count": 0,
                "genes": [],
                "enrichment_results": {}
            }

        # Perform enrichment analysis for each category
        enrichment_results = {}
        enrichment_errors = {}

        for category in ENRICHMENT_CATEGORIES:
            try:
                enriched_terms = self.perform_enrichment_analysis(genes, category)
                enrichment_results[category] = enriched_terms
            except Exception as e:
                logger.error(f"Error in enrichment analysis for {category}: {e}")
                enrichment_results[category] = []
                enrichment_errors[category] = str(e)

        result = {
            "disease": disease,
            "gene_count": gene_count,
            "genes": genes,
            "enrichment_results": enrichment_results
        }

        # Add enrichment errors if any occurred
        if enrichment_errors:
            result["enrichment_errors"] = enrichment_errors

        return result

    def save_incremental_result(self, result: Dict[str, Any], all_results: List[Dict[str, Any]]):
        """Save a single result incrementally to JSONL and optionally create backup"""
        # Append to results list
        all_results.append(result)

        # Append to JSONL file (much more efficient than rewriting entire file)
        try:
            with open(self.results_file, 'a') as f:
                f.write(json.dumps(result) + '\n')
            logger.debug(f"Appended result for {result['disease']['name']} to JSONL")
        except Exception as e:
            logger.error(f"Error appending to JSONL file: {e}")

        # Create periodic backup (only keep most recent)
        if len(all_results) % self.backup_frequency == 0:
            backup_file = f"{self.results_file}.backup"
            # Remove previous backup if it exists
            if os.path.exists(backup_file):
                try:
                    os.remove(backup_file)
                except Exception as e:
                    logger.warning(f"Error removing old backup: {e}")

            # Create new backup by copying current JSONL file
            try:
                import shutil
                shutil.copy2(self.results_file, backup_file)
                logger.info(f"Created backup at {backup_file} ({len(all_results)} diseases)")
            except Exception as e:
                logger.warning(f"Error creating backup: {e}")

    def run_analysis(
        self,
        max_diseases: int = None,
        resume: bool = True,
        min_genes: int = 1,
        data_file: str = "robokop_disease_genes.json",
    ) -> List[Dict[str, Any]]:
        """
        Run the fast enrichment analysis pipeline using pre-processed data

        Args:
            max_diseases: Maximum number of diseases to analyze (for testing)
            resume: Whether to resume from existing results file
            min_genes: Minimum number of genes required for enrichment analysis

        Returns:
            List of analysis results sorted by gene count (ascending)
        """
        logger.info("Starting fast disease enrichment analysis pipeline...")

        # Load pre-processed disease-gene data
        diseases_with_genes = self.load_disease_gene_data(data_file)
        if not diseases_with_genes:
            logger.error("No disease-gene data available. Exiting.")
            return []

        # Load existing results if resuming
        all_results = []
        processed_curies = set()

        if resume:
            all_results = self.load_existing_results()
            processed_curies = self.get_processed_disease_curies(all_results)
            if processed_curies:
                logger.info(f"Resuming analysis - skipping {len(processed_curies)} already processed diseases")

        # Filter diseases
        # Remove already processed diseases
        if processed_curies:
            diseases_with_genes = [d for d in diseases_with_genes if d["disease"]["curie"] not in processed_curies]
            logger.info(f"Found {len(diseases_with_genes)} new diseases to process")

        # Filter by minimum gene count
        if min_genes > 0:
            diseases_with_genes = [d for d in diseases_with_genes if d["gene_count"] >= min_genes]
            logger.info(f"After filtering by min_genes={min_genes}: {len(diseases_with_genes)} diseases remain")

        # Limit number of diseases for testing
        if max_diseases:
            total_to_process = min(max_diseases - len(all_results), len(diseases_with_genes))
            diseases_with_genes = diseases_with_genes[:total_to_process]
            logger.info(f"Limiting analysis to {total_to_process} new diseases (total target: {max_diseases})")

        if not diseases_with_genes:
            logger.info("No new diseases to process")
            return all_results

        logger.info(f"Ready to process {len(diseases_with_genes)} diseases for enrichment analysis")

        # Perform enrichment analysis for each disease
        for i, disease_with_genes in enumerate(diseases_with_genes, 1):
            current_total = len(all_results) + 1
            disease_name = disease_with_genes["disease"]["name"]
            gene_count = disease_with_genes["gene_count"]

            logger.info(f"Processing disease {i}/{len(diseases_with_genes)} (total: {current_total})")

            try:
                result = self.perform_disease_enrichment(disease_with_genes)
                self.save_incremental_result(result, all_results)

            except KeyboardInterrupt:
                logger.info("Analysis interrupted by user. Results saved up to this point.")
                break
            except Exception as e:
                logger.error(f"Error processing disease {disease_with_genes['disease']['curie']}: {e}")
                # Create a placeholder result for failed diseases
                result = {
                    "disease": disease_with_genes["disease"],
                    "gene_count": disease_with_genes["gene_count"],
                    "genes": disease_with_genes["genes"],
                    "enrichment_results": {},
                    "error": str(e)
                }
                self.save_incremental_result(result, all_results)

        # Final save not needed - results are already incrementally saved to JSONL
        logger.info(f"All results already saved incrementally to {self.results_file}")

        logger.info(f"Analysis complete! Total processed: {len(all_results)} diseases")

        return all_results

    def print_summary(self, results: List[Dict[str, Any]], top_n: int = 10):
        """Print a summary of the analysis results"""
        print(f"\n=== Fast Disease Enrichment Analysis Summary ===")
        print(f"Total diseases analyzed: {len(results)}")

        print(f"\nTop {top_n} diseases by gene count (ascending):")
        print(f"{'Rank':<5} {'Gene Count':<12} {'Disease Name':<50} {'CURIE':<20}")
        print("-" * 90)

        for i, result in enumerate(results[:top_n], 1):
            disease = result["disease"]
            gene_count = result["gene_count"]
            print(f"{i:<5} {gene_count:<12} {disease['name'][:47]:<50} {disease['curie']:<20}")

        # Show enrichment summary
        if results:
            print(f"\nEnrichment categories analyzed: {', '.join(ENRICHMENT_CATEGORIES)}")

            # Count diseases with enrichment results
            enriched_counts = {category: 0 for category in ENRICHMENT_CATEGORIES}
            error_counts = {category: 0 for category in ENRICHMENT_CATEGORIES}

            for result in results:
                for category in ENRICHMENT_CATEGORIES:
                    if result["enrichment_results"].get(category):
                        enriched_counts[category] += 1

                    # Count enrichment errors
                    if result.get("enrichment_errors", {}).get(category):
                        error_counts[category] += 1

            print("\nEnrichment results:")
            for category in ENRICHMENT_CATEGORIES:
                success_count = enriched_counts[category]
                error_count = error_counts[category]
                print(f"  {category}:")
                print(f"    Successful: {success_count}/{len(results)} diseases")
                if error_count > 0:
                    print(f"    Errors/Timeouts: {error_count} diseases")

        # Show timeout/error details
        timeout_results = [r for r in results if "enrichment_errors" in r]
        if timeout_results:
            print(f"\n=== Enrichment Errors/Timeouts ===")
            for result in timeout_results:
                disease_name = result["disease"]["name"]
                gene_count = result["gene_count"]
                print(f"\n{disease_name} ({gene_count} genes):")
                for category, error in result["enrichment_errors"].items():
                    if "timeout" in error.lower() or "read timed out" in error.lower():
                        print(f"  TIMEOUT - {category}: {error}")
                    else:
                        print(f"  ERROR - {category}: {error}")


def main():
    """Main function to run the fast analysis"""
    parser = argparse.ArgumentParser(description="Run disease enrichment through AnswerCoalesce.")
    parser.add_argument(
        "--data-file",
        default="robokop_disease_genes.json",
        help="Disease-gene JSON file created by ingest_robokop_data.py",
    )
    parser.add_argument(
        "--results-file",
        default="fast_enrichment_results.jsonl",
        help="Output JSONL file for incremental enrichment results",
    )
    parser.add_argument(
        "--answer-coalesce-url",
        default=os.environ.get("ANSWER_COALESCE_URL", ANSWER_COALESCE_URL),
        help="AnswerCoalesce /query endpoint",
    )
    parser.add_argument("--delay", type=float, default=0.1, help="Delay between API requests")
    parser.add_argument("--max-diseases", type=int, default=None, help="Maximum diseases to process")
    parser.add_argument("--min-genes", type=int, default=3, help="Minimum gene count")
    parser.add_argument("--no-resume", action="store_true", help="Do not resume from existing results")
    args = parser.parse_args()

    # Initialize analyzer
    analyzer = FastDiseaseEnrichmentAnalyzer(
        delay_between_requests=args.delay,
        results_file=args.results_file,
        answer_coalesce_url=args.answer_coalesce_url,
    )

    # Run analysis with local data
    # - Set min_genes=1 to skip diseases with 0 genes (no enrichment possible)
    # - For testing: set max_diseases=10
    # - For full analysis: set max_diseases=None
    # - Set resume=True to continue from existing file
    results = analyzer.run_analysis(
        max_diseases=args.max_diseases,
        resume=not args.no_resume,
        min_genes=args.min_genes,
        data_file=args.data_file,
    )

    # Print summary
    analyzer.print_summary(results)

    print(f"\nResults saved to '{analyzer.results_file}'")
    print("Use parse_results_to_tsv.py to convert to spreadsheet format!")


if __name__ == "__main__":
    main()
