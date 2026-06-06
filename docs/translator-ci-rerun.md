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

Build disease-gene mappings from the extracted graph:

```bash
uv run python ingest_robokop_data.py \
  --graph-dir data/kgx/translator_kg/latest \
  --output-file artifacts/runs/translator-ci-YYYY-MM-DD/translator_disease_genes.json
```

Run enrichment through CI AnswerCoalesce:

```bash
uv run python fast_disease_enrichment_analysis.py \
  --data-file artifacts/runs/translator-ci-YYYY-MM-DD/translator_disease_genes.json \
  --results-file artifacts/runs/translator-ci-YYYY-MM-DD/enrichment_results.jsonl \
  --answer-coalesce-url https://answer-coalesce.ci.transltr.io/query \
  --min-genes 3
```

Extract disease-term ground truth from the same Translator KGX graph:

```bash
uv run python extract_robokop_disease_term_edges.py \
  --graph-dir data/kgx/translator_kg/latest \
  --output-prefix artifacts/runs/translator-ci-YYYY-MM-DD/translator_disease_term_edges
```

The ground-truth extractor streams JSONL/TSV outputs so it does not have to
materialize the full subclass-inferred edge set in memory. Use `--direct-only`
for a direct-edge-only diagnostic run.

Evaluate:

```bash
uv run python compare_enrichment_vs_robokop.py \
  --enrichment artifacts/runs/translator-ci-YYYY-MM-DD/enrichment_results.jsonl \
  --ground-truth artifacts/runs/translator-ci-YYYY-MM-DD/translator_disease_term_edges.jsonl \
  --max-k 100 \
  --output-plot artifacts/runs/translator-ci-YYYY-MM-DD/enrichment_hits_at_k.png
```
