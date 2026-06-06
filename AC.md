# AnswerCoalesce Enrichment API Documentation

## Overview

The **AnswerCoalesce API** is a statistical enrichment service provided by RENCI (Renaissance Computing Institute) that performs gene set enrichment analysis. Given a set of genes, it calculates which biological processes, molecular activities, or pathways are statistically over-represented compared to what would be expected by chance.

**Base URL:** `https://answercoalesce.renci.org/query`

## What is Enrichment Analysis?

Enrichment analysis identifies biological terms (like GO processes, molecular functions, or pathways) that appear more frequently in your gene set than would be expected by chance. This helps interpret what biological functions or processes are associated with a disease or condition.

For example, if you have 50 genes associated with Alzheimer's disease, enrichment analysis might reveal that terms like "protein folding" or "response to oxidative stress" are significantly over-represented, suggesting these processes are important in the disease mechanism.

## API Specification

### Request Format

The API uses the **TRAPI (Translator Reasoner API)** query format, a standardized JSON structure for biomedical knowledge graph queries.

#### HTTP Method
```
POST https://answercoalesce.renci.org/query
```

#### Headers
```
Content-Type: application/json
```

#### Request Body Structure

```json
{
  "message": {
    "query_graph": {
      "nodes": {
        "input": {
          "categories": ["biolink:Gene"],
          "ids": ["uuid:1"],
          "member_ids": ["NCBIGene:1234", "NCBIGene:5678", ...],
          "set_interpretation": "MANY"
        },
        "output": {
          "categories": ["biolink:BiologicalProcess"]
        }
      },
      "edges": {
        "edge_0": {
          "subject": "input",
          "object": "output",
          "predicates": ["biolink:related_to"]
        }
      }
    }
  }
}
```

### Request Parameters

#### Input Node Configuration
- **`categories`**: Always `["biolink:Gene"]` for enrichment analysis
- **`ids`**: Placeholder identifier, use `["uuid:1"]`
- **`member_ids`**: **Your list of gene CURIEs** (e.g., `["NCBIGene:1234", "NCBIGene:5678"]`)
- **`set_interpretation`**: Must be `"MANY"` to indicate set-based analysis

#### Output Node Configuration
- **`categories`**: Target enrichment category. Choose one of:
  - `"biolink:BiologicalProcess"` - GO Biological Process terms (e.g., "apoptosis", "cell division")
  - `"biolink:MolecularActivity"` - GO Molecular Function terms (e.g., "kinase activity", "DNA binding")
  - `"biolink:Pathway"` - Pathway databases (e.g., KEGG, Reactome pathways)

#### Edge Configuration
- **`subject`**: `"input"` (genes are the subject)
- **`object`**: `"output"` (enriched terms are the object)
- **`predicates`**: Use `["biolink:related_to"]` for general gene-term relationships

### Response Format

The API returns enrichment results in TRAPI message format with statistical significance values.

#### Response Structure

```json
{
  "message": {
    "results": [
      {
        "node_bindings": {
          "output": [{"id": "GO:0006915"}]
        },
        "analyses": [
          {
            "edge_bindings": {
              "edge_0": [{"id": "edge_12345"}]
            }
          }
        ]
      }
    ],
    "knowledge_graph": {
      "nodes": {
        "GO:0006915": {
          "name": "apoptotic process",
          "categories": ["biolink:BiologicalProcess"]
        }
      },
      "edges": {
        "edge_12345": {
          "subject": "uuid:1",
          "object": "GO:0006915",
          "attributes": [
            {
              "attribute_type_id": "biolink:p_value",
              "value": 0.0001234
            }
          ]
        }
      }
    }
  }
}
```

#### Parsing the Response

1. **Iterate through `message.results`**: Each result represents one enriched term
2. **Extract term CURIE**: From `result["node_bindings"]["output"][0]["id"]`
3. **Get term name**: From `knowledge_graph["nodes"][term_curie]["name"]`
4. **Extract p-value**:
   - Get edge ID from `result["analyses"][0]["edge_bindings"]["edge_0"][0]["id"]`
   - Find edge in `knowledge_graph["edges"][edge_id]`
   - Look for attribute with `"attribute_type_id": "biolink:p_value"`
   - Extract the `value` field

#### Example Parsed Result

```python
{
  "curie": "GO:0006915",
  "name": "apoptotic process",
  "category": "biolink:BiologicalProcess",
  "p_value": 0.0001234
}
```

Lower p-values indicate stronger statistical significance (more confident enrichment).

## Usage Examples

### Python Example

```python
import requests
import json

# API endpoint
ANSWER_COALESCE_URL = "https://answercoalesce.renci.org/query"

def perform_enrichment_analysis(gene_curies, category):
    """
    Perform enrichment analysis for a list of genes

    Args:
        gene_curies: List of gene CURIEs (e.g., ["NCBIGene:1234", "NCBIGene:5678"])
        category: Target category (e.g., "biolink:BiologicalProcess")

    Returns:
        List of enriched terms with p-values
    """
    # Build TRAPI query
    query = {
        "message": {
            "query_graph": {
                "nodes": {
                    "input": {
                        "categories": ["biolink:Gene"],
                        "ids": ["uuid:1"],
                        "member_ids": gene_curies,
                        "set_interpretation": "MANY"
                    },
                    "output": {
                        "categories": [category]
                    }
                },
                "edges": {
                    "edge_0": {
                        "subject": "input",
                        "object": "output",
                        "predicates": ["biolink:related_to"]
                    }
                }
            }
        }
    }

    # Send request
    response = requests.post(
        ANSWER_COALESCE_URL,
        json=query,
        headers={"Content-Type": "application/json"},
        timeout=60
    )
    response.raise_for_status()
    data = response.json()

    # Parse results
    enriched_terms = []
    knowledge_graph = data["message"].get("knowledge_graph", {})
    nodes = knowledge_graph.get("nodes", {})
    edges = knowledge_graph.get("edges", {})

    for result in data["message"].get("results", []):
        # Get output node binding (the enriched term)
        term_curie = result["node_bindings"]["output"][0]["id"]
        node_info = nodes.get(term_curie, {})

        # Get edge information and p-value
        edge_id = result["analyses"][0]["edge_bindings"]["edge_0"][0]["id"]
        edge_info = edges.get(edge_id, {})

        # Extract p-value from edge attributes
        p_value = None
        for attr in edge_info.get("attributes", []):
            if attr.get("attribute_type_id") == "biolink:p_value":
                p_value = attr.get("value")
                break

        enriched_terms.append({
            "curie": term_curie,
            "name": node_info.get("name", "Unknown"),
            "category": category,
            "p_value": p_value
        })

    return enriched_terms

# Example usage
genes = ["NCBIGene:5970", "NCBIGene:4318", "NCBIGene:596"]
results = perform_enrichment_analysis(genes, "biolink:BiologicalProcess")

for term in results[:10]:  # Show top 10 results
    print(f"{term['name']}: p-value={term['p_value']}")
```

### Analyzing Multiple Categories

```python
# Enrichment categories to analyze
ENRICHMENT_CATEGORIES = [
    "biolink:BiologicalProcess",
    "biolink:MolecularActivity",
    "biolink:Pathway"
]

# Get genes for a disease
disease_genes = ["NCBIGene:5970", "NCBIGene:4318", "NCBIGene:596"]

# Perform enrichment for each category
all_enrichment_results = {}

for category in ENRICHMENT_CATEGORIES:
    enriched_terms = perform_enrichment_analysis(disease_genes, category)
    all_enrichment_results[category] = enriched_terms
    print(f"\n{category}: {len(enriched_terms)} enriched terms")

    # Show top 5 most significant terms
    sorted_terms = sorted(enriched_terms, key=lambda x: x['p_value'] or 1.0)
    for term in sorted_terms[:5]:
        print(f"  {term['name']}: p={term['p_value']}")
```

## Best Practices

### 1. Minimum Gene Count
- **Requirement**: Use at least **3 genes** for reliable enrichment analysis
- **Why**: Fewer genes produce unstable statistics and API errors
- **Filter**: `diseases_with_genes = [d for d in diseases if d["gene_count"] >= 3]`

### 2. Rate Limiting
- Add delays between requests to be respectful to the API
- Recommended: 0.1 seconds between requests
- Example: `time.sleep(0.1)` after each API call

### 3. Timeout Handling
- Set reasonable timeout values (60 seconds recommended)
- Large gene sets may take longer to process
- Implement retry logic for transient failures

### 4. Error Handling
```python
try:
    enriched_terms = perform_enrichment_analysis(genes, category)
except requests.exceptions.Timeout:
    print(f"Timeout for {category}")
    enriched_terms = []
except requests.exceptions.HTTPError as e:
    print(f"HTTP error for {category}: {e}")
    enriched_terms = []
except Exception as e:
    print(f"Unexpected error: {e}")
    enriched_terms = []
```

### 5. Result Interpretation
- **P-values**: Lower is more significant (e.g., p < 0.05 is typically considered significant)
- **Multiple testing**: Consider Bonferroni or FDR correction when testing many terms
- **Biological relevance**: Statistical significance doesn't always mean biological importance

### 6. Gene Identifiers
- Use **NCBIGene** CURIEs (e.g., `"NCBIGene:1234"`)
- Ensure genes are properly normalized before querying
- Invalid or unrecognized genes will be silently ignored

## Common Use Cases

### 1. Disease Enrichment Analysis
Identify biological processes associated with a disease by analyzing disease-associated genes:

```python
# Get genes for Alzheimer's disease
alzheimers_genes = get_disease_genes("MONDO:0004975")

# Perform enrichment
processes = perform_enrichment_analysis(alzheimers_genes, "biolink:BiologicalProcess")
activities = perform_enrichment_analysis(alzheimers_genes, "biolink:MolecularActivity")
pathways = perform_enrichment_analysis(alzheimers_genes, "biolink:Pathway")
```

### 2. Batch Processing Multiple Diseases
Process many diseases with incremental saving:

```python
results = []
for disease in diseases:
    genes = disease["genes"]

    if len(genes) < 3:
        continue  # Skip diseases with too few genes

    enrichment_results = {}
    for category in ENRICHMENT_CATEGORIES:
        terms = perform_enrichment_analysis(genes, category)
        enrichment_results[category] = terms

    result = {
        "disease": disease,
        "gene_count": len(genes),
        "enrichment_results": enrichment_results
    }
    results.append(result)

    # Save incrementally
    with open("results.jsonl", "a") as f:
        f.write(json.dumps(result) + "\n")

    time.sleep(0.1)  # Rate limiting
```

### 3. Comparative Analysis
Compare enrichment across different gene sets:

```python
# Compare disease A vs disease B
disease_a_enrichment = perform_enrichment_analysis(disease_a_genes, "biolink:BiologicalProcess")
disease_b_enrichment = perform_enrichment_analysis(disease_b_genes, "biolink:BiologicalProcess")

# Find common enriched processes
a_terms = {term["curie"] for term in disease_a_enrichment}
b_terms = {term["curie"] for term in disease_b_enrichment}
common_terms = a_terms.intersection(b_terms)

print(f"Common enriched processes: {len(common_terms)}")
```

## Performance Considerations

### Speed
- **Typical request time**: 0.5-2 seconds per category
- **Large gene sets**: May take longer (5-10 seconds)
- **Recommended batch size**: Process 1 disease per second with 3 categories

### Scalability
- Use JSONL format for large-scale analysis (append-only, memory efficient)
- Implement incremental saving to handle interruptions
- Support resume capability for long-running analyses

### Caching
Consider caching results when:
- Running the same analysis multiple times
- Testing different evaluation methods on same enrichment data
- Comparing with ground truth or other datasets

## Troubleshooting

### Empty Results
- **Cause**: No significant enrichment found
- **Solution**: Check if genes are valid and properly formatted

### API Timeouts
- **Cause**: Large gene sets or server load
- **Solution**: Increase timeout, retry with exponential backoff

### HTTP 4xx Errors
- **Cause**: Invalid request format or parameters
- **Solution**: Validate TRAPI query structure, check gene CURIEs

### Diseases with <3 Genes
- **Cause**: Statistical tests fail with insufficient genes
- **Solution**: Filter to `min_genes >= 3` before analysis

## Related Resources

### TRAPI Specification
The AnswerCoalesce API follows the TRAPI standard:
- [TRAPI Documentation](https://github.com/NCATSTranslator/ReasonerAPI)
- Standardized biomedical knowledge graph query format
- Used across the NIH NCATS Translator project

### Biolink Model
- [Biolink Model](https://biolink.github.io/biolink-model/)
- Defines categories like `biolink:Gene`, `biolink:BiologicalProcess`
- Standard ontology for biomedical entities

### CURIE Format
- **Compact URI**: `Prefix:LocalID` (e.g., `NCBIGene:1234`, `GO:0006915`)
- Standardized identifiers for biological entities
- Enables interoperability across databases

## Integration with ROBOKOP Knowledge Graph

The AnswerCoalesce API is commonly used with ROBOKOP (Reasoning Over Biomedical Objects linked in Knowledge Oriented Pathways):

```python
# 1. Get disease-gene associations from ROBOKOP
ROBOKOP_URL = "https://automat.renci.org/robokopkg"
response = requests.get(f"{ROBOKOP_URL}/edges/{disease_curie}",
                       params={"category": "biolink:Gene"})
genes = [edge["adj_node"]["id"] for edge in response.json()["edges"]]

# 2. Perform enrichment with AnswerCoalesce
enriched_terms = perform_enrichment_analysis(genes, "biolink:BiologicalProcess")

# 3. Results provide biological interpretation of disease
for term in enriched_terms[:10]:
    print(f"{term['name']}: p={term['p_value']}")
```

## Summary

The AnswerCoalesce API provides:
- **Statistical enrichment analysis** for gene sets
- **TRAPI-compliant** query and response format
- **Multiple enrichment categories** (processes, activities, pathways)
- **P-values** for statistical significance assessment
- **Integration** with ROBOKOP and other Translator services

For production use, implement proper error handling, rate limiting, and result caching. Process diseases in order from smallest to largest gene count for optimal API stability.
