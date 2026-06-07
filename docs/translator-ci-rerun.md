# Translator CI Rerun Notes

This rerun should use the Translator Ingests merged KGX release and the CI
AnswerCoalesce service rather than the older ROBOKOP-local graph and RENCI dev
AnswerCoalesce URL.

## Verified Inputs

- Translator KGX release browser: `https://kgx-storage.rtx.ai/releases/translator_kg/latest/`
- Current `latest` release resolved during setup: `2026_04_22`
- Archive: `https://kgx-storage.rtx.ai/releases/translator_kg/latest/translator_kg.tar.zst`
- Direct nodes file: `https://kgx-storage.rtx.ai/releases/translator_kg/latest/nodes.jsonl`
- Direct edges file: `https://kgx-storage.rtx.ai/releases/translator_kg/latest/edges.jsonl`
- CI AnswerCoalesce query endpoint: `https://answer-coalesce.ci.transltr.io/query`
- CI AnswerCoalesce OpenAPI: `https://answer-coalesce.ci.transltr.io/openapi.json`

The KGX listing showed `nodes.jsonl`, `edges.jsonl`, `graph-metadata.json`,
`merge-metadata.json`, and `translator_kg.tar.zst`. The metadata reports
Biolink `4.3.6`, Babel `2025sep1`, and `29,423,079` total edges.

## Local Rerun Layout

Downloaded and derived files should stay out of Git:

```text
data/kgx/translator_kg/latest/
artifacts/runs/<run-id>/
```

Recommended run id format:

```text
translator-ci-YYYY-MM-DD
```

## Commands

Install dependencies:

```bash
uv sync
```

Download and extract the latest Translator KGX archive:

```bash
uv run python scripts/download_translator_kgx.py --extract
```

The preferred entry point is the config-driven wrapper:

```bash
uv run python run_pipeline.py --config configs/translator_ci_latest.toml --stage manifest
uv run python run_pipeline.py --config configs/translator_ci_latest.toml --stage ingest
uv run python run_pipeline.py --config configs/translator_ci_latest.toml --stage enrich
uv run python run_pipeline.py --config configs/translator_ci_latest.toml --stage ground-truth
uv run python run_pipeline.py --config configs/translator_ci_latest.toml --stage term-hierarchy
uv run python run_pipeline.py --config configs/translator_ci_latest.toml --stage tsv
uv run python run_pipeline.py --config configs/translator_ci_latest.toml --stage evaluate
uv run python run_pipeline.py --config configs/translator_ci_latest.toml --stage db
```

The wrapper writes outputs under the `run.output_dir` configured in
`configs/translator_ci_latest.toml`.

Equivalent direct command for disease-gene extraction:

```bash
uv run python ingest_robokop_data.py \
  --graph-dir data/kgx/translator_kg/latest \
  --output-file artifacts/runs/translator-ci-YYYY-MM-DD/disease_genes.json
```

Equivalent direct command for enrichment through CI AnswerCoalesce:

```bash
uv run python fast_disease_enrichment_analysis.py \
  --data-file artifacts/runs/translator-ci-YYYY-MM-DD/disease_genes.json \
  --results-file artifacts/runs/translator-ci-YYYY-MM-DD/enrichment_results.jsonl \
  --answer-coalesce-url https://answer-coalesce.ci.transltr.io/query \
  --min-genes 3
```

Extract disease-term ground truth from the same Translator KGX graph:

```bash
uv run python extract_robokop_disease_term_edges.py \
  --graph-dir data/kgx/translator_kg/latest \
  --output-prefix artifacts/runs/translator-ci-YYYY-MM-DD/disease_term_edges
```

The ground-truth extractor streams JSONL/TSV outputs so it does not have to
materialize the full subclass-inferred edge set in memory. Use `--direct-only`
for a direct-edge-only diagnostic run.

Evaluate:

```bash
uv run python compare_enrichment_vs_robokop.py \
  --enrichment artifacts/runs/translator-ci-YYYY-MM-DD/enrichment_results.jsonl \
  --ground-truth artifacts/runs/translator-ci-YYYY-MM-DD/disease_term_edges.jsonl \
  --max-k 100 \
  --output-plot artifacts/runs/translator-ci-YYYY-MM-DD/enrichment_hits_at_k.png
```

Extract generic term hierarchy where KGX provides subclass edges:

```bash
uv run python extract_term_hierarchy.py \
  --graph-dir data/kgx/translator_kg/latest \
  --output artifacts/runs/translator-ci-YYYY-MM-DD/term_hierarchy.json
```

Build the visualization database:

```bash
uv run python disease-enrichment-viz/backend/data_prep/prepare_enrichment_database.py \
  --term-hierarchy-file artifacts/runs/translator-ci-YYYY-MM-DD/term_hierarchy.json \
  --enrichment-file artifacts/runs/translator-ci-YYYY-MM-DD/enrichment_results.jsonl \
  --db-path artifacts/runs/translator-ci-YYYY-MM-DD/enrichment_database.db \
  --stats-path artifacts/runs/translator-ci-YYYY-MM-DD/database_stats.json
```
