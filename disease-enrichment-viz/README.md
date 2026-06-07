# Disease Enrichment Visualization

This is a local Flask and D3 web app for exploring disease enrichment results.
It visualizes hierarchy when term hierarchy is available and still shows
standalone ranked terms when no hierarchy exists.

## Inputs

The app reads a SQLite database created by:

```bash
uv run python disease-enrichment-viz/backend/data_prep/prepare_enrichment_database.py \
  --term-hierarchy-file artifacts/runs/<run-id>/term_hierarchy.json \
  --enrichment-file artifacts/runs/<run-id>/enrichment_results.jsonl \
  --db-path artifacts/runs/<run-id>/enrichment_database.db \
  --stats-path artifacts/runs/<run-id>/database_stats.json
```

The database builder accepts missing hierarchy files. In that case it loads
enrichment terms without hierarchy relationships.

## Database Schema

The current schema is term-generic:

- `diseases`: disease ID, name, description, and gene count.
- `terms`: enriched terms from AnswerCoalesce plus optional KGX term metadata.
- `term_hierarchy`: `child_id -> parent_id` subclass relationships when present.
- `enrichment_results`: disease-term p-values and ranks.

Compatibility views named `go_terms` and `go_hierarchy` are created for older
code paths that still expect those names.

## Run The App

Build the database first, then start the backend:

```bash
ENRICHMENT_DATABASE_PATH=artifacts/runs/<run-id>/enrichment_database.db \
ENRICHMENT_DATABASE_STATS_PATH=artifacts/runs/<run-id>/database_stats.json \
uv run python disease-enrichment-viz/backend/app.py
```

Open:

```text
http://127.0.0.1:5000/
```

The backend also serves the frontend, so a separate static server is not needed.

## API Endpoints

- `GET /api`: API status.
- `GET /api/disease/<disease_id>`: disease metadata.
- `GET /api/disease/<disease_id>/enrichment`: enrichment terms and optional
  hierarchy for one disease.
- `GET /api/search?q=<query>`: disease search.
- `GET /api/hierarchy?terms=<term_ids>`: hierarchy among supplied terms.
- `GET /api/stats`: database summary counts.

## Visualization Semantics

- Yellow nodes are the original enrichment results for the selected disease.
- Node size is based on `-log10(p-value)`.
- Edges are subclass relationships among returned enrichment terms.
- A virtual root is added to keep disconnected terms visible.
- If a category has no hierarchy, the network view still displays the enriched
  terms as standalone nodes connected through the virtual root.
