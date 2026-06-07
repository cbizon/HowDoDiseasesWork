#!/usr/bin/env python3
"""Run the disease enrichment workflow from a run configuration file.

The run contract is:
1. a KGX graph with `nodes.jsonl` and `edges.jsonl`
2. an AnswerCoalesce `/query` endpoint built from that KGX

Each stage writes artifacts under the configured output directory and can be
rerun independently.
"""

import argparse
import json
import subprocess
import time
import tomllib
import urllib.error
import urllib.request
from pathlib import Path


DEFAULT_STAGES = [
    "manifest",
    "ingest",
    "enrich",
    "ground-truth",
    "term-hierarchy",
    "tsv",
    "evaluate",
    "db",
]


def load_config(path):
    with open(path, "rb") as f:
        return tomllib.load(f)


def config_value(config, *paths, default=None):
    for path in paths:
        current = config
        found = True
        for part in path.split("."):
            if not isinstance(current, dict) or part not in current:
                found = False
                break
            current = current[part]
        if found:
            return current
    return default


def run_command(command, dry_run=False):
    print("\n$ " + " ".join(str(part) for part in command))
    if dry_run:
        return
    subprocess.run(command, check=True)


def build_context(config, config_path):
    run_id = config_value(config, "run.id", default=Path(config_path).stem)
    output_dir = Path(config_value(config, "run.output_dir", default=f"artifacts/runs/{run_id}"))
    graph_dir = Path(
        config_value(
            config,
            "kgx.graph_dir",
            "kgx.local_dir",
            "translator_kgx.graph_dir",
            "translator_kgx.local_dir",
            default="data/kgx/translator_kg/latest",
        )
    )

    query_url = config_value(
        config,
        "answer_coalesce.query_url",
        "answer_coalesce.url",
        default="https://answercoalesce.renci.org/query",
    )
    min_genes = int(config_value(config, "analysis.min_genes", default=3))
    delay = float(config_value(config, "analysis.delay", default=0.1))
    max_diseases = config_value(config, "analysis.max_diseases", default=None)
    categories = config_value(
        config,
        "analysis.categories",
        default=[
            "biolink:BiologicalProcess",
            "biolink:MolecularActivity",
            "biolink:Pathway",
        ],
    )

    return {
        "config_path": str(config_path),
        "run_id": run_id,
        "output_dir": output_dir,
        "graph_dir": graph_dir,
        "nodes_file": graph_dir / "nodes.jsonl",
        "edges_file": graph_dir / "edges.jsonl",
        "answer_coalesce_url": query_url,
        "min_genes": min_genes,
        "delay": delay,
        "max_diseases": max_diseases,
        "categories": categories,
        "disease_genes": output_dir / "disease_genes.json",
        "enrichment_results": output_dir / "enrichment_results.jsonl",
        "term_hierarchy": output_dir / "term_hierarchy.json",
        "ground_truth_prefix": output_dir / "disease_term_edges",
        "tsv_prefix": output_dir / "enrichment",
        "evaluation_plot": output_dir / "enrichment_hits_at_k.png",
        "visualization_db": output_dir / "enrichment_database.db",
        "visualization_stats": output_dir / "database_stats.json",
        "manifest": output_dir / "run_manifest.json",
    }


def validate_inputs(context):
    for path_key in ["nodes_file", "edges_file"]:
        path = context[path_key]
        if not path.exists():
            raise FileNotFoundError(f"Missing KGX file: {path}")

    request = urllib.request.Request(context["answer_coalesce_url"], method="HEAD")
    try:
        with urllib.request.urlopen(request, timeout=15) as response:
            status = response.status
    except urllib.error.HTTPError as e:
        # Many TRAPI endpoints do not allow HEAD but still prove DNS/TLS/routing.
        status = e.code
    except Exception as e:
        raise RuntimeError(f"AnswerCoalesce endpoint is not reachable: {e}") from e

    print(f"Validated KGX files and AnswerCoalesce endpoint ({status})")


def write_manifest(context, config):
    context["output_dir"].mkdir(parents=True, exist_ok=True)
    manifest = {
        "run_id": context["run_id"],
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "config_path": context["config_path"],
        "kgx": {
            "graph_dir": str(context["graph_dir"]),
            "nodes_file": str(context["nodes_file"]),
            "edges_file": str(context["edges_file"]),
        },
        "answer_coalesce": {
            "query_url": context["answer_coalesce_url"],
        },
        "analysis": {
            "min_genes": context["min_genes"],
            "delay": context["delay"],
            "max_diseases": context["max_diseases"],
            "categories": context["categories"],
        },
        "artifacts": {
            "disease_genes": str(context["disease_genes"]),
            "enrichment_results": str(context["enrichment_results"]),
            "term_hierarchy": str(context["term_hierarchy"]),
            "ground_truth_prefix": str(context["ground_truth_prefix"]),
            "visualization_db": str(context["visualization_db"]),
        },
        "source_config": config,
    }

    with context["manifest"].open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    print(f"Wrote run manifest to {context['manifest']}")


def run_stage(stage, context, dry_run=False):
    if stage == "manifest":
        return
    if stage == "ingest":
        run_command(
            [
                "uv",
                "run",
                "python",
                "ingest_robokop_data.py",
                "--graph-dir",
                context["graph_dir"],
                "--output-file",
                context["disease_genes"],
            ],
            dry_run,
        )
    elif stage == "enrich":
        command = [
            "uv",
            "run",
            "python",
            "fast_disease_enrichment_analysis.py",
            "--data-file",
            context["disease_genes"],
            "--results-file",
            context["enrichment_results"],
            "--answer-coalesce-url",
            context["answer_coalesce_url"],
            "--min-genes",
            str(context["min_genes"]),
            "--delay",
            str(context["delay"]),
        ]
        if context["max_diseases"] is not None:
            command.extend(["--max-diseases", str(context["max_diseases"])])
        run_command(command, dry_run)
    elif stage == "ground-truth":
        run_command(
            [
                "uv",
                "run",
                "python",
                "extract_robokop_disease_term_edges.py",
                "--graph-dir",
                context["graph_dir"],
                "--output-prefix",
                context["ground_truth_prefix"],
            ],
            dry_run,
        )
    elif stage == "term-hierarchy":
        run_command(
            [
                "uv",
                "run",
                "python",
                "extract_term_hierarchy.py",
                "--graph-dir",
                context["graph_dir"],
                "--output",
                context["term_hierarchy"],
                "--categories",
                *context["categories"],
            ],
            dry_run,
        )
    elif stage == "tsv":
        run_command(
            [
                "uv",
                "run",
                "python",
                "parse_results_to_tsv.py",
                "--input",
                context["enrichment_results"],
                "--output-prefix",
                context["tsv_prefix"],
                "--min-genes",
                str(context["min_genes"]),
            ],
            dry_run,
        )
    elif stage == "evaluate":
        run_command(
            [
                "uv",
                "run",
                "python",
                "compare_enrichment_vs_robokop.py",
                "--enrichment",
                context["enrichment_results"],
                "--ground-truth",
                f"{context['ground_truth_prefix']}.jsonl",
                "--max-k",
                "100",
                "--output-plot",
                context["evaluation_plot"],
            ],
            dry_run,
        )
    elif stage == "db":
        run_command(
            [
                "uv",
                "run",
                "python",
                "disease-enrichment-viz/backend/data_prep/prepare_enrichment_database.py",
                "--term-hierarchy-file",
                context["term_hierarchy"],
                "--enrichment-file",
                context["enrichment_results"],
                "--db-path",
                context["visualization_db"],
                "--stats-path",
                context["visualization_stats"],
            ],
            dry_run,
        )
    else:
        raise ValueError(f"Unknown stage: {stage}")


def main():
    parser = argparse.ArgumentParser(description="Run the KGX + AnswerCoalesce enrichment pipeline")
    parser.add_argument("--config", required=True, help="Run configuration TOML")
    parser.add_argument(
        "--stage",
        action="append",
        choices=DEFAULT_STAGES,
        help="Stage to run. Can be provided multiple times. Defaults to all stages.",
    )
    parser.add_argument("--skip-validate", action="store_true", help="Skip input reachability checks")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running them")
    args = parser.parse_args()

    config = load_config(args.config)
    context = build_context(config, args.config)
    stages = args.stage or DEFAULT_STAGES

    if not args.skip_validate:
        validate_inputs(context)

    if "manifest" in stages:
        write_manifest(context, config)

    for stage in stages:
        run_stage(stage, context, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
