# Disease Enrichment Analysis Pipeline

Comprehensive pipeline for disease enrichment analysis using ROBOKOP Knowledge Graph. Includes enrichment prediction, ground truth comparison, hits@k evaluation, precision/recall analysis, and stratified analysis by gene count.

## Overview

This multi-phase analysis pipeline:
1. **Extracts** disease-gene associations from ROBOKOP Knowledge Graph (47,439 diseases)
2. **Predicts** enriched terms using AnswerCoalesce API for:
   - **Biological Processes** (biolink:BiologicalProcess)
   - **Molecular Activities** (biolink:MolecularActivity)
   - **Pathways** (biolink:Pathway)
3. **Extracts ground truth** disease-term edges from ROBOKOP (10.7M edges with subclass inference)
4. **Evaluates performance** using hits@k, precision@k, and recall@k metrics
5. **Provides stratified analysis** by gene count (tertile groups)

Results include statistical significance rankings, performance visualizations, and comprehensive TSV exports.

## Prerequisites

- Python 3.11+.
- `uv` for environment and dependency management.
- KGX graph files with `nodes.jsonl` and `edges.jsonl`.

Install dependencies:

```bash
uv sync
```

The original run used a local ROBOKOP graph export. New reruns should use the
Translator KGX release and CI AnswerCoalesce endpoint documented in
`docs/translator-ci-rerun.md`.

## Quick Start

### 1. Data Preparation (One-time setup)

If you haven't already processed the ROBOKOP graph data:

```bash
# Process KGX graph files to extract disease-gene mappings
uv run python ingest_robokop_data.py --graph-dir /path/to/kgx
```

This creates `robokop_disease_genes.json` containing disease-gene associations for 47,439 diseases.

### 2. Run Enrichment Analysis

```bash
# Run full analysis on all diseases with ≥3 genes (recommended)
uv run python fast_disease_enrichment_analysis.py
```

**Key parameters:**
- Minimum 3 genes required (diseases with 1-2 genes cause API errors)
- Results saved incrementally to `fast_enrichment_results.jsonl`
- Automatic resume capability if interrupted
- Backup created every 10 diseases processed

### 3. Extract Ground Truth from ROBOKOP

```bash
# Extract disease-term edges for evaluation
uv run python extract_robokop_disease_term_edges.py --graph-dir /path/to/kgx
```

This creates `robokop_disease_term_edges_with_subclass.jsonl` with 10.7M ground truth edges.

### 4. Generate TSV Reports

```bash
# Convert JSONL results to TSV format for analysis
uv run python parse_results_to_tsv.py -i fast_enrichment_results.jsonl -o results
```

This creates:
- `results_disease_summary.tsv` - Disease overview with enrichment counts
- `results_enrichment_results.tsv` - All enriched terms with p-values
- `results_error_summary.tsv` - Failed analyses (usually empty)
- `results_gene_distribution.tsv` - Gene count distribution

### 5. Run Performance Evaluation

```bash
# Hits@K analysis comparing predictions vs ground truth
uv run python compare_enrichment_vs_robokop.py --max-k 100

# Precision@K and Recall@K analysis
uv run python precision_recall_at_k_analysis.py --max-k 100

# Comprehensive precision/recall plots
uv run python create_precision_recall_plots.py
```

### 6. Stratified Analysis by Gene Count

```bash
# Gene count distribution analysis
uv run python create_gene_count_histogram.py

# Stratified precision/recall by gene count tertiles
uv run python stratified_precision_recall_analysis.py --max-k 200
```

## Translator CI Rerun

The checked-in config for the current Translator rerun is
`configs/translator_ci_latest.toml`.

```bash
uv run python scripts/download_translator_kgx.py --extract
```

Then run the pipeline with the extracted KGX files and CI AnswerCoalesce:

```bash
uv run python ingest_robokop_data.py \
  --graph-dir data/kgx/translator_kg/latest \
  --output-file artifacts/runs/translator-ci-YYYY-MM-DD/translator_disease_genes.json

uv run python fast_disease_enrichment_analysis.py \
  --data-file artifacts/runs/translator-ci-YYYY-MM-DD/translator_disease_genes.json \
  --results-file artifacts/runs/translator-ci-YYYY-MM-DD/enrichment_results.jsonl \
  --answer-coalesce-url https://answer-coalesce.ci.transltr.io/query \
  --min-genes 3
```

Large KGX downloads and generated outputs are intentionally ignored by Git.

## File Descriptions

### Core Analysis Scripts
- **`fast_disease_enrichment_analysis.py`** - Main enrichment analysis pipeline (JSONL-based)
- **`ingest_robokop_data.py`** - Extract disease-gene mappings from graph files
- **`extract_robokop_disease_term_edges.py`** - Extract ground truth disease-term edges

### Evaluation Scripts
- **`compare_enrichment_vs_robokop.py`** - Hits@k analysis vs ground truth
- **`precision_recall_at_k_analysis.py`** - Comprehensive precision/recall analysis
- **`create_precision_recall_plots.py`** - Simple precision/recall plotting
- **`stratified_precision_recall_analysis.py`** - Gene count stratified analysis
- **`create_gene_count_histogram.py`** - Gene count distribution visualization

### Utility Scripts
- **`parse_results_to_tsv.py`** - Convert JSONL results to TSV format
- **`disease_enrichment_analysis.py`** - Legacy version (slow, JSON-based)

### Key Data Files
- **`robokop_disease_genes.json`** - Disease-gene mappings (47,439 diseases)
- **`disease_gene_counts.tsv`** - Gene count summary per disease
- **`fast_enrichment_results.jsonl`** - Enrichment analysis results
- **`robokop_disease_term_edges_with_subclass.jsonl`** - Ground truth (10.7M edges)
- **`enrichment_hits_at_k_detailed.tsv`** - Hits@k evaluation results

### Generated Visualizations
- **`enrichment_hits_at_k.png`** - Hits@k performance curves
- **`precision_recall_at_k_plots.png`** - Precision/recall performance curves
- **`gene_count_histogram.png`** - Gene count distribution plots
- **`stratified_precision_recall.png`** - Stratified analysis by gene count

## Complete Analysis Workflow

### 1. Initial Setup and Enrichment
```bash
# Extract disease-gene mappings (one-time)
python ingest_robokop_data.py

# Run enrichment analysis on all eligible diseases
python fast_disease_enrichment_analysis.py

# Generate TSV reports
python parse_results_to_tsv.py -i fast_enrichment_results.jsonl -o enrichment
```

### 2. Ground Truth and Evaluation
```bash
# Extract ground truth from ROBOKOP
python extract_robokop_disease_term_edges.py

# Run hits@k evaluation
python compare_enrichment_vs_robokop.py --max-k 100 --output-plot hits_at_k.png

# Generate precision/recall analysis
python precision_recall_at_k_analysis.py --max-k 100
```

### 3. Stratified Analysis
```bash
# Analyze gene count distribution
python create_gene_count_histogram.py

# Run stratified analysis by gene count tertiles
python stratified_precision_recall_analysis.py --max-k 200
```

### 4. Quick Analysis Examples
```bash
# Test with limited diseases
python -c "
from fast_disease_enrichment_analysis import FastDiseaseEnrichmentAnalyzer
analyzer = FastDiseaseEnrichmentAnalyzer()
results = analyzer.run_analysis(max_diseases=10, min_genes=3)
print(f'Processed {len(results)} diseases')
"

# Custom evaluation parameters
python compare_enrichment_vs_robokop.py --enrichment my_results.jsonl --max-k 50

# Generate specific plots
python create_precision_recall_plots.py  # Creates 4-panel visualization
```

## Performance Summary

### Enrichment Analysis Performance
- **Speed**: ~1 disease per second (18,000x faster than original API approach)
- **Scale**: 4,221 diseases with ≥3 genes successfully analyzed
- **Memory**: JSONL format enables incremental processing
- **Reliability**: 99%+ success rate with ≥3 gene threshold

### Evaluation Results (Latest)
- **BiologicalProcess**: Hits@10=18.8%, P@10=2.5%, R@10=7.5%
- **MolecularActivity**: Hits@10=2.5%, P@10=0.4%, R@10=0.2%
- **Pathway**: Hits@10=1.4%, P@10=0.2%, R@10=0.1%
- **Ground Truth**: 10.7M disease-term edges (with subclass inference)
- **Stratification**: Tertile groups by gene count (3-5, 6-19, ≥20 genes)

## Output Formats

### Enrichment Analysis Outputs

#### Disease Summary TSV
```
Rank	Disease_CURIE	Disease_Name	Gene_Count	BP_Enriched	MA_Enriched	Pathway_Enriched	BP_Error	MA_Error	Pathway_Error	Has_Genes
1	MONDO:0006502	acute respiratory distress syndrome	3	65	20	3				Yes
```

#### Enrichment Results TSV
```
disease_curie	name	enrichment_type	enriched_entity_curie	enriched_entity	p_value	rank	number_of_enriched_entities
MONDO:0006502	acute respiratory distress syndrome	biolink:BiologicalProcess	GO:0035634	response to stilbenoid	0.000186	1	65
```

### Evaluation Analysis Outputs

#### Hits@K Results TSV
```
K	BiologicalProcess	MolecularActivity	Pathway
1	0.0348	0.0035	0.0041
10	0.1877	0.0250	0.0140
```

#### Ground Truth Edges JSONL
```json
{"source_curie": "MONDO:0007254", "target_curie": "GO:0008283", "target_type": "biolink:BiologicalProcess"}
```

### Visualization Outputs
- **Performance curves**: Hits@K, Precision@K, Recall@K vs K
- **Stratified analysis**: Performance by gene count groups
- **Distribution plots**: Gene count histograms and cumulative distributions

## Key Analysis Results

### Dataset Statistics
- **Total diseases**: 47,439 in ROBOKOP Knowledge Graph
- **Diseases with genes**: 12,744 (26.9%)
- **Eligible diseases**: 4,221 with ≥3 genes (8.9%)
- **Ground truth edges**: 10.7M disease-term associations
- **Subclass inference**: Enhances coverage via disease hierarchy

### Performance Insights
- **BiologicalProcess**: Best performing category, reasonable precision/recall
- **Gene count correlation**: Higher gene count → better enrichment performance
- **Stratification effect**: Diseases with more genes show improved metrics
- **Evaluation coverage**: Comprehensive K=1 to K=200 analysis

## Troubleshooting

### Common Issues
- **API Errors**: Diseases with <3 genes cause failures (use min_genes=3)
- **Memory**: Large ground truth files require sufficient RAM for analysis
- **Disk Space**: Monitor space for 10.7M edge ground truth file
- **Runtime**: Full stratified analysis can take hours (use background execution)

### Performance Optimization
- **Incremental saving**: All analyses support resume capability
- **Background execution**: Long-running analyses use `run_in_background=true`
- **JSONL format**: Memory-efficient for large datasets
- **Parallel processing**: Multiple analyses can run concurrently

## Data Sources and Methods

- **ROBOKOP Knowledge Graph**: 47,439 diseases, gene associations, ontology terms
- **AnswerCoalesce API**: Statistical enrichment calculations
- **Subclass Inference**: Disease hierarchy for comprehensive ground truth
- **Evaluation Metrics**: Hits@K, Precision@K, Recall@K up to K=200
- **Stratification**: Tertile-based analysis by gene count (3-5, 6-19, ≥20 genes)
