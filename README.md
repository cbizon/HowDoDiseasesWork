# How Do Diseases Work

This repository runs disease gene-set enrichment against a KGX graph and an
AnswerCoalesce endpoint built from that same graph. The same code should work
for ROBOKOP or Translator as long as the run starts with these two inputs:

1. A KGX directory containing `nodes.jsonl` and `edges.jsonl`.
2. An AnswerCoalesce `/query` URL built from that KGX.

Large KGX downloads and generated analysis outputs are intentionally ignored by
Git. Do not commit data files.

## Requirements

- Python 3.11 or newer.
- `uv` for environment and dependency management.
- A local KGX graph directory with `nodes.jsonl` and `edges.jsonl`.
- An AnswerCoalesce `/query` endpoint built from the same graph.
- Enough local disk for the KGX and derived outputs. Translator KGX-scale runs
  can require tens of GB.

Install the Python environment:

```bash
uv sync
```

All commands below assume they are run from the repository root through `uv`.

## Configure A Run

Run settings live in TOML files under `configs/`.

The checked-in Translator CI config is:

```bash
configs/translator_ci_latest.toml
```

It uses:

- KGX directory: `data/kgx/translator_kg/latest`
- AnswerCoalesce URL: `https://answer-coalesce.ci.transltr.io/query`
- Output directory: `artifacts/runs/translator-ci-latest`

For a ROBOKOP-local run, copy and edit:

```bash
cp configs/robokop_local.example.toml configs/robokop_local.toml
```

Set `kgx.graph_dir` to the local KGX directory and
`answer_coalesce.query_url` to the `/query` endpoint built from that graph.

## Get Translator KGX

For Translator runs, download and extract the latest Translator KGX release:

```bash
uv run python scripts/download_translator_kgx.py --extract
```

The downloader writes under `data/kgx/translator_kg/latest/`, which is ignored by
Git. It also stores release metadata next to the graph files so a run can be
audited later.

## Run The Pipeline

The wrapper records a manifest and runs each stage with paths from the config.
For large runs, run one stage at a time so failures are easy to inspect and
resume.

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

To see the commands without running them:

```bash
uv run python run_pipeline.py \
  --config configs/translator_ci_latest.toml \
  --stage ingest \
  --stage enrich \
  --dry-run \
  --skip-validate
```

## Stage Outputs

Given `run.output_dir = artifacts/runs/<run-id>`, the main outputs are:

- `run_manifest.json`: config, KGX paths, AnswerCoalesce URL, and output paths.
- `disease_genes.json`: disease to gene mappings extracted from KGX.
- `enrichment_results.jsonl`: one JSON object per disease, appended as results
  are derived.
- `disease_term_edges.jsonl`: ground-truth disease-term edges from the same KGX.
- `term_hierarchy.json`: term hierarchy extracted from KGX subclass edges where
  hierarchy exists.
- `enrichment_*.tsv`: spreadsheet-friendly parsed enrichment outputs.
- `enrichment_hits_at_k.png`: hits-at-k evaluation plot.
- `enrichment_database.db`: SQLite database for the visualization app.
- `database_stats.json`: database summary counts for the visualization app.

The enrichment stage writes incrementally. If it stops, rerunning the same
command resumes by reading already processed disease CURIEs from the JSONL file.

## Handling Failed Enrichment Rows

AnswerCoalesce or network failures are recorded in `enrichment_results.jsonl`
instead of being silently hidden. If a failure affected a block of diseases and
the endpoint is healthy again, create a clean resume file that drops rows with
recorded errors:

```bash
uv run python scripts/filter_clean_enrichment_results.py \
  --input artifacts/runs/<run-id>/enrichment_results.jsonl \
  --output artifacts/runs/<run-id>/enrichment_results_clean_resume.jsonl \
  --errors-output artifacts/runs/<run-id>/enrichment_results_errors.jsonl
```

Then resume enrichment using the clean file as `--results-file`:

```bash
uv run python fast_disease_enrichment_analysis.py \
  --data-file artifacts/runs/<run-id>/disease_genes.json \
  --results-file artifacts/runs/<run-id>/enrichment_results_clean_resume.jsonl \
  --answer-coalesce-url https://answer-coalesce.ci.transltr.io/query \
  --min-genes 3 \
  --delay 0.1
```

This preserves successful rows and reattempts only diseases whose previous rows
were error records.

## Hierarchy Behavior

The pipeline does not assume every enrichment category has a hierarchy.

- `extract_term_hierarchy.py` scans KGX `subclass_of` edges for configured term
  categories.
- GO-like terms keep hierarchy relationships where KGX provides them.
- Terms without hierarchy are still loaded into the database and returned by the
  API as ranked enrichment results.
- The visualization database uses generic `terms` and `term_hierarchy` tables,
  with compatibility views named `go_terms` and `go_hierarchy` for older code.

That means hierarchy-aware rollups can be added where hierarchy exists without
dropping non-hierarchical terms.

## Visualization

Build the SQLite database after enrichment and optional hierarchy extraction:

```bash
uv run python run_pipeline.py --config configs/translator_ci_latest.toml --stage db
```

Run the Flask backend against a run-specific database:

```bash
ENRICHMENT_DATABASE_PATH=artifacts/runs/<run-id>/enrichment_database.db \
ENRICHMENT_DATABASE_STATS_PATH=artifacts/runs/<run-id>/database_stats.json \
uv run python disease-enrichment-viz/backend/app.py
```

Open the frontend at:

```text
http://127.0.0.1:5000/
```

See `disease-enrichment-viz/README.md` for visualization details.

## Script Map

- `run_pipeline.py`: config-driven wrapper for the full workflow.
- `scripts/download_translator_kgx.py`: download and extract Translator KGX.
- `ingest_robokop_data.py`: historical name; extracts disease-gene mappings
  from any compatible KGX graph.
- `fast_disease_enrichment_analysis.py`: calls AnswerCoalesce and appends JSONL
  enrichment rows.
- `extract_robokop_disease_term_edges.py`: historical name; extracts
  disease-term ground truth from any compatible KGX graph.
- `extract_term_hierarchy.py`: extracts generic term hierarchy from KGX.
- `parse_results_to_tsv.py`: writes TSV summaries from enrichment JSONL.
- `compare_enrichment_vs_robokop.py`: historical name; evaluates enrichment
  results against extracted ground truth.
- `disease-enrichment-viz/backend/data_prep/prepare_enrichment_database.py`:
  builds the visualization SQLite database.

## Repository Hygiene

The `.gitignore` excludes generated JSONL, TSV, plots, SQLite databases, KGX
archives, `data/`, and `artifacts/`. Source code, config examples, and docs are
tracked; run data is not.
