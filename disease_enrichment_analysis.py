#!/usr/bin/env python3
"""
Disease Enrichment Analysis Pipeline

This script performs enrichment analysis to find Processes, Pathways, and Activities
that are enriched for diseases using two main APIs:
1. ROBOKOP Knowledge Graph API for disease-gene associations
2. AnswerCoalesce API for enrichment calculations

Pipeline:
1. Retrieve all diseases from ROBOKOP KG
2. For each disease, get associated genes
3. Sort diseases by number of associated genes (ascending)
4. Perform enrichment analysis for BiologicalProcess, MolecularActivity, and Pathway
"""

import requests
import json
import time
import os
from typing import List, Dict, Any, Tuple
from collections import defaultdict
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# API URLs
ROBOKOP_BASE_URL = "https://automat.renci.org/robokopkg"
ANSWER_COALESCE_URL = "https://answercoalesce.renci.org/query"

# Target categories for enrichment analysis
ENRICHMENT_CATEGORIES = [
    "biolink:BiologicalProcess",
    "biolink:MolecularActivity",
    "biolink:Pathway"
]

class DiseaseEnrichmentAnalyzer:
    """Main class for performing disease enrichment analysis"""

    def __init__(self, delay_between_requests: float = 0.1, results_file: str = "disease_enrichment_results.json"):
        """
        Initialize the analyzer

        Args:
            delay_between_requests: Delay in seconds between API requests to be respectful
            results_file: File path for saving results incrementally
        """
        self.delay = delay_between_requests
        self.session = requests.Session()
        self.results_file = results_file
        self.genes_cache_file = "disease_genes_cache.json"
        self.backup_frequency = 10  # Save backup every N diseases
        self.cache_save_frequency = 50  # Save gene cache every N diseases

    def load_existing_results(self) -> List[Dict[str, Any]]:
        """
        Load existing results from file if it exists

        Returns:
            List of existing analysis results, or empty list if file doesn't exist
        """
        if os.path.exists(self.results_file):
            try:
                with open(self.results_file, 'r') as f:
                    results = json.load(f)
                logger.info(f"Loaded {len(results)} existing results from {self.results_file}")
                return results
            except Exception as e:
                logger.warning(f"Error loading existing results: {e}")
                return []
        return []

    def get_processed_disease_curies(self, results: List[Dict[str, Any]]) -> set:
        """
        Get set of disease CURIEs that have already been processed

        Args:
            results: Existing analysis results

        Returns:
            Set of disease CURIE identifiers
        """
        return {result["disease"]["curie"] for result in results}

    def save_incremental_result(self, result: Dict[str, Any], all_results: List[Dict[str, Any]]):
        """
        Save a single result incrementally and optionally create backup

        Args:
            result: Single disease analysis result
            all_results: Complete list of results so far
        """
        # Append to results list
        all_results.append(result)

        # Save to file
        try:
            with open(self.results_file, 'w') as f:
                json.dump(all_results, f, indent=2)
            logger.debug(f"Saved incremental result for {result['disease']['name']}")
        except Exception as e:
            logger.error(f"Error saving incremental result: {e}")

        # Create periodic backup
        if len(all_results) % self.backup_frequency == 0:
            backup_file = f"{self.results_file}.backup_{len(all_results)}"
            try:
                with open(backup_file, 'w') as f:
                    json.dump(all_results, f, indent=2)
                logger.info(f"Created backup at {backup_file}")
            except Exception as e:
                logger.warning(f"Error creating backup: {e}")

    def load_genes_cache(self) -> Dict[str, Dict[str, Any]]:
        """
        Load existing gene cache from file if it exists

        Returns:
            Dictionary mapping disease CURIEs to gene information
        """
        if os.path.exists(self.genes_cache_file):
            try:
                with open(self.genes_cache_file, 'r') as f:
                    cache = json.load(f)
                logger.info(f"Loaded gene cache with {len(cache)} diseases from {self.genes_cache_file}")
                return cache
            except Exception as e:
                logger.warning(f"Error loading gene cache: {e}")
                return {}
        return {}

    def save_genes_cache(self, cache: Dict[str, Dict[str, Any]]):
        """
        Save gene cache to file

        Args:
            cache: Dictionary mapping disease CURIEs to gene information
        """
        try:
            with open(self.genes_cache_file, 'w') as f:
                json.dump(cache, f, indent=2)
            logger.info(f"Saved gene cache with {len(cache)} diseases to {self.genes_cache_file}")
        except Exception as e:
            logger.error(f"Error saving gene cache: {e}")

    def update_genes_cache(self, cache: Dict[str, Dict[str, Any]], disease_curie: str, disease_with_genes: Dict[str, Any], diseases_processed: int):
        """
        Update gene cache with new disease and optionally save to file

        Args:
            cache: Gene cache dictionary to update
            disease_curie: Disease CURIE identifier
            disease_with_genes: Disease info with genes
            diseases_processed: Number of diseases processed so far
        """
        cache[disease_curie] = disease_with_genes

        # Save cache periodically
        if diseases_processed % self.cache_save_frequency == 0:
            self.save_genes_cache(cache)

    def get_all_diseases(self) -> List[Dict[str, Any]]:
        """
        Retrieve all diseases from ROBOKOP Knowledge Graph using Cypher query

        Returns:
            List of disease nodes with their CURIE IDs and metadata
        """
        logger.info("Fetching all diseases from ROBOKOP KG...")

        cypher_query = {
            "query": "MATCH (d:`biolink:Disease`) RETURN d.id as curie, d.name as name LIMIT 1000"
        }

        try:
            response = self.session.post(
                f"{ROBOKOP_BASE_URL}/cypher",
                json=cypher_query,
                headers={
                    "accept": "application/json",
                    "Content-Type": "application/json"
                }
            )
            response.raise_for_status()

            data = response.json()
            diseases = []

            # Parse Neo4j cypher response format
            for result_set in data.get("results", []):
                for row_data in result_set.get("data", []):
                    row = row_data.get("row", [])
                    if len(row) >= 2 and row[0] and row[1]:  # curie and name
                        diseases.append({
                            "curie": row[0],
                            "name": row[1]
                        })

            logger.info(f"Found {len(diseases)} diseases")
            return diseases

        except Exception as e:
            logger.error(f"Error fetching diseases: {e}")
            return []

    def get_genes_for_disease(self, disease_curie: str) -> List[str]:
        """
        Get all genes associated with a specific disease using ROBOKOP edges API

        Args:
            disease_curie: The CURIE identifier for the disease

        Returns:
            List of gene CURIE identifiers
        """
        logger.debug(f"Fetching genes for disease: {disease_curie}")

        try:
            # Use the edges endpoint to find genes connected to this disease
            response = self.session.get(
                f"{ROBOKOP_BASE_URL}/edges/{disease_curie}",
                params={
                    "category": "biolink:Gene",  # Filter for gene connections
                    "limit": 1000  # Increase limit to get more genes
                }
            )
            response.raise_for_status()

            data = response.json()
            genes = set()  # Use set to avoid duplicates

            # Extract gene CURIEs from the edges response
            for edge_info in data.get("edges", []):
                adj_node = edge_info.get("adj_node", {})
                # Check if adjacent node is a gene
                if "biolink:Gene" in adj_node.get("category", []):
                    gene_id = adj_node.get("id")
                    if gene_id:
                        genes.add(gene_id)

            genes_list = list(genes)
            logger.debug(f"Found {len(genes_list)} genes for {disease_curie}")

            # Add delay to be respectful to the API
            time.sleep(self.delay)

            return genes_list

        except Exception as e:
            logger.error(f"Error fetching genes for disease {disease_curie}: {e}")
            return []

    def generate_trapi_query(self, input_node_type: str, output_node_type: str,
                           input_curies: List[str], predicate: str,
                           input_is_subject: bool = True) -> Dict[str, Any]:
        """
        Generate a TRAPI query for enrichment analysis

        Args:
            input_node_type: Type of input nodes (e.g., "biolink:Gene")
            output_node_type: Type of output nodes (e.g., "biolink:BiologicalProcess")
            input_curies: List of input CURIE identifiers
            predicate: Relationship predicate (e.g., "biolink:related_to")
            input_is_subject: Whether input nodes are subjects in the relationship

        Returns:
            TRAPI query dictionary
        """
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
        """
        Perform enrichment analysis for a list of genes against a specific category

        Args:
            genes: List of gene CURIE identifiers
            category: Target category for enrichment (e.g., "biolink:BiologicalProcess")

        Returns:
            List of enriched terms with metadata including p-values
        """
        if not genes:
            logger.warning("No genes provided for enrichment analysis")
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
                ANSWER_COALESCE_URL,
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
            return []

    def get_disease_genes_only(self, disease: Dict[str, Any]) -> Dict[str, Any]:
        """
        Get genes for a disease without performing enrichment analysis

        Args:
            disease: Disease dictionary with curie and name

        Returns:
            Disease info with genes and gene count
        """
        disease_curie = disease["curie"]
        disease_name = disease["name"]

        logger.debug(f"Getting genes for: {disease_name} ({disease_curie})")

        # Get genes for this disease
        genes = self.get_genes_for_disease(disease_curie)

        return {
            "disease": disease,
            "gene_count": len(genes),
            "genes": genes
        }

    def perform_disease_enrichment(self, disease_with_genes: Dict[str, Any]) -> Dict[str, Any]:
        """
        Perform enrichment analysis for a disease that already has genes retrieved

        Args:
            disease_with_genes: Disease dictionary with genes already retrieved

        Returns:
            Analysis results including enrichment results
        """
        disease = disease_with_genes["disease"]
        genes = disease_with_genes["genes"]
        gene_count = disease_with_genes["gene_count"]

        disease_name = disease["name"]
        disease_curie = disease["curie"]

        logger.info(f"Analyzing enrichment for: {disease_name} ({gene_count} genes)")

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

    def run_analysis(self, max_diseases: int = None, resume: bool = True) -> List[Dict[str, Any]]:
        """
        Run the complete enrichment analysis pipeline with two phases:
        Phase 1: Get all diseases and retrieve genes, sort by gene count
        Phase 2: Perform enrichment analysis in gene count order (fewest genes first)

        Args:
            max_diseases: Maximum number of diseases to analyze (for testing)
            resume: Whether to resume from existing results file

        Returns:
            List of analysis results sorted by gene count (ascending)
        """
        logger.info("Starting two-phase disease enrichment analysis pipeline...")

        # Load existing results if resuming
        all_results = []
        processed_curies = set()

        if resume:
            all_results = self.load_existing_results()
            processed_curies = self.get_processed_disease_curies(all_results)
            if processed_curies:
                logger.info(f"Resuming analysis - skipping {len(processed_curies)} already processed diseases")

        # === PHASE 1: GET ALL DISEASES AND THEIR GENES ===
        logger.info("=== PHASE 1: Retrieving genes for all diseases ===")

        # Load gene cache
        genes_cache = self.load_genes_cache()

        # Get all diseases
        diseases = self.get_all_diseases()
        if not diseases:
            logger.error("No diseases found. Exiting.")
            return all_results

        # Filter out already processed diseases
        if processed_curies:
            diseases = [d for d in diseases if d["curie"] not in processed_curies]
            logger.info(f"Found {len(diseases)} new diseases to process")

        if max_diseases:
            total_to_process = min(max_diseases - len(all_results), len(diseases))
            diseases = diseases[:total_to_process]
            logger.info(f"Limiting analysis to {total_to_process} new diseases (total target: {max_diseases})")

        if not diseases:
            logger.info("No new diseases to process")
            return all_results

        # Phase 1: Get genes for all diseases (using cache when available)
        cached_count = sum(1 for d in diseases if d["curie"] in genes_cache)
        need_retrieval = len(diseases) - cached_count

        if cached_count > 0:
            logger.info(f"Phase 1: Found {cached_count} diseases in gene cache, need to retrieve {need_retrieval} new diseases")
        else:
            logger.info(f"Phase 1: Retrieving genes for {len(diseases)} diseases...")

        diseases_with_genes = []

        for i, disease in enumerate(diseases, 1):
            disease_curie = disease["curie"]
            disease_name = disease["name"]

            # Check if disease genes are already cached
            if disease_curie in genes_cache:
                logger.debug(f"Using cached genes for disease {i}/{len(diseases)}: {disease_name}")
                disease_with_genes = genes_cache[disease_curie]
                diseases_with_genes.append(disease_with_genes)

                # Log gene count
                gene_count = disease_with_genes["gene_count"]
                if gene_count == 0:
                    logger.debug(f"Cached: No genes for {disease_name}")
                else:
                    logger.debug(f"Cached: {gene_count} genes for {disease_name}")
                continue

            # Need to retrieve genes for this disease
            logger.info(f"Getting genes for disease {i}/{len(diseases)}: {disease_name}")

            try:
                disease_with_genes = self.get_disease_genes_only(disease)
                diseases_with_genes.append(disease_with_genes)

                # Update cache
                self.update_genes_cache(genes_cache, disease_curie, disease_with_genes, i)

                # Log gene count
                gene_count = disease_with_genes["gene_count"]
                if gene_count == 0:
                    logger.warning(f"No genes found for {disease_name}")
                else:
                    logger.info(f"Found {gene_count} genes for {disease_name}")

            except KeyboardInterrupt:
                logger.info("Gene retrieval interrupted by user. Saving cache...")
                self.save_genes_cache(genes_cache)
                break
            except Exception as e:
                logger.error(f"Error getting genes for {disease_curie}: {e}")
                # Add disease with 0 genes and update cache
                disease_with_genes = {
                    "disease": disease,
                    "gene_count": 0,
                    "genes": [],
                    "error": str(e)
                }
                diseases_with_genes.append(disease_with_genes)
                self.update_genes_cache(genes_cache, disease_curie, disease_with_genes, i)

        # Save final gene cache
        self.save_genes_cache(genes_cache)

        # Sort diseases by gene count (ascending - fewest genes first)
        diseases_with_genes.sort(key=lambda x: x["gene_count"])
        logger.info(f"Phase 1 complete! Sorted {len(diseases_with_genes)} diseases by gene count")

        # Print gene count distribution
        gene_counts = [d["gene_count"] for d in diseases_with_genes]
        logger.info(f"Gene count range: {min(gene_counts)} to {max(gene_counts)} genes")
        logger.info(f"Diseases with 0 genes: {sum(1 for x in gene_counts if x == 0)}")

        # === PHASE 2: PERFORM ENRICHMENT ANALYSIS IN ORDER ===
        logger.info("=== PHASE 2: Performing enrichment analysis in gene count order ===")

        for i, disease_with_genes in enumerate(diseases_with_genes, 1):
            current_total = len(all_results) + 1
            disease_name = disease_with_genes["disease"]["name"]
            gene_count = disease_with_genes["gene_count"]

            logger.info(f"Phase 2: Processing disease {i}/{len(diseases_with_genes)} (total: {current_total})")
            logger.info(f"Enrichment analysis: {disease_name} ({gene_count} genes)")

            try:
                result = self.perform_disease_enrichment(disease_with_genes)
                self.save_incremental_result(result, all_results)

            except KeyboardInterrupt:
                logger.info("Enrichment analysis interrupted by user. Results saved up to this point.")
                break
            except Exception as e:
                logger.error(f"Error in enrichment analysis for {disease_with_genes['disease']['curie']}: {e}")
                # Create a placeholder result for failed enrichment
                result = {
                    "disease": disease_with_genes["disease"],
                    "gene_count": disease_with_genes["gene_count"],
                    "genes": disease_with_genes["genes"],
                    "enrichment_results": {},
                    "error": str(e)
                }
                self.save_incremental_result(result, all_results)

        # Final save (results are already sorted by gene count)
        try:
            with open(self.results_file, 'w') as f:
                json.dump(all_results, f, indent=2)
            logger.info(f"Final results saved to {self.results_file}")
        except Exception as e:
            logger.error(f"Error saving final results: {e}")

        logger.info(f"Analysis complete! Total processed: {len(all_results)} diseases")

        return all_results

    def save_results(self, results: List[Dict[str, Any]], filename: str = "disease_enrichment_results.json"):
        """
        Save analysis results to a JSON file

        Args:
            results: Analysis results from run_analysis()
            filename: Output filename
        """
        logger.info(f"Saving results to {filename}")

        with open(filename, 'w') as f:
            json.dump(results, f, indent=2)

        logger.info(f"Results saved successfully")

    def print_summary(self, results: List[Dict[str, Any]], top_n: int = 10):
        """
        Print a summary of the analysis results

        Args:
            results: Analysis results from run_analysis()
            top_n: Number of top results to display
        """
        print(f"\n=== Disease Enrichment Analysis Summary ===")
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
                total_attempted = success_count + error_count + sum(1 for r in results if r["gene_count"] > 0 and category not in r.get("enrichment_errors", {}))
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
    """Main function to run the analysis"""

    # Initialize analyzer
    analyzer = DiseaseEnrichmentAnalyzer(delay_between_requests=0.1)

    # Run analysis with incremental saving and resume capability
    # - For testing: set max_diseases=5
    # - For full analysis: remove max_diseases parameter or set to None
    # - Set resume=False to start fresh, resume=True to continue from existing file
    results = analyzer.run_analysis(max_diseases=None, resume=True)

    # Print summary
    analyzer.print_summary(results)

    # Results are already saved incrementally, no need to save again


if __name__ == "__main__":
    main()