#!/usr/bin/env python3
"""Create a clean enrichment JSONL file for resume after recorded API failures."""

import argparse
import json
from pathlib import Path


def row_has_recorded_error(row):
    return bool(row.get("error")) or bool(row.get("enrichment_errors"))


def filter_results(input_path, output_path, errors_path=None):
    kept = 0
    dropped = 0

    output_path.parent.mkdir(parents=True, exist_ok=True)
    errors_handle = None
    if errors_path:
        errors_path.parent.mkdir(parents=True, exist_ok=True)
        errors_handle = errors_path.open("w", encoding="utf-8")

    try:
        with input_path.open("r", encoding="utf-8") as input_handle, output_path.open(
            "w", encoding="utf-8"
        ) as output_handle:
            for line in input_handle:
                if not line.strip():
                    continue

                row = json.loads(line)
                if row_has_recorded_error(row):
                    dropped += 1
                    if errors_handle:
                        errors_handle.write(line)
                    continue

                kept += 1
                output_handle.write(json.dumps(row) + "\n")
    finally:
        if errors_handle:
            errors_handle.close()

    return {"kept": kept, "dropped": dropped, "input": str(input_path), "output": str(output_path)}


def main():
    parser = argparse.ArgumentParser(
        description="Filter enrichment JSONL to rows without recorded errors."
    )
    parser.add_argument("--input", required=True, type=Path, help="Original enrichment JSONL")
    parser.add_argument("--output", required=True, type=Path, help="Filtered JSONL for resume")
    parser.add_argument(
        "--errors-output",
        type=Path,
        help="Optional JSONL file containing only dropped error rows",
    )
    args = parser.parse_args()

    stats = filter_results(args.input, args.output, args.errors_output)
    print(
        "Filtered enrichment results: "
        f"kept {stats['kept']:,}, dropped {stats['dropped']:,}, wrote {stats['output']}"
    )


if __name__ == "__main__":
    main()
