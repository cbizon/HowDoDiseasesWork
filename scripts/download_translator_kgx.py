#!/usr/bin/env python3
"""Download and optionally extract the latest Translator KGX archive."""

from __future__ import annotations

import argparse
import shutil
import tarfile
from pathlib import Path
from urllib.parse import urljoin

import requests
import zstandard


DEFAULT_BASE_URL = "https://kgx-storage.rtx.ai/releases/translator_kg/latest/"
DEFAULT_OUTPUT_DIR = Path("data/kgx/translator_kg/latest")
DEFAULT_ARCHIVE_NAME = "translator_kg.tar.zst"


def download_file(url: str, output_path: Path, force: bool = False) -> None:
    if output_path.exists() and not force:
        print(f"Using existing file: {output_path}")
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output_path.with_suffix(output_path.suffix + ".part")

    print(f"Downloading {url}")
    with requests.get(url, stream=True, timeout=60) as response:
        response.raise_for_status()
        with tmp_path.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)

    tmp_path.replace(output_path)
    print(f"Saved {output_path}")


def safe_extract_tar_zst(archive_path: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    output_root = output_dir.resolve()

    print(f"Extracting {archive_path} into {output_dir}")
    with archive_path.open("rb") as compressed:
        reader = zstandard.ZstdDecompressor().stream_reader(compressed)
        with tarfile.open(fileobj=reader, mode="r|") as archive:
            for member in archive:
                target = (output_dir / member.name).resolve()
                if not str(target).startswith(str(output_root)):
                    raise ValueError(f"Refusing unsafe archive path: {member.name}")
                archive.extract(member, output_dir)

    print("Extraction complete")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-url", default=DEFAULT_BASE_URL)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--archive-name", default=DEFAULT_ARCHIVE_NAME)
    parser.add_argument("--force", action="store_true", help="Redownload existing files")
    parser.add_argument("--extract", action="store_true", help="Extract archive after download")
    args = parser.parse_args()

    archive_url = urljoin(args.base_url, args.archive_name)
    metadata_url = urljoin(args.base_url, "graph-metadata.json")
    latest_release_url = urljoin(args.base_url, "../latest-release.json")

    archive_path = args.output_dir / args.archive_name
    download_file(metadata_url, args.output_dir / "graph-metadata.json", force=args.force)
    download_file(latest_release_url, args.output_dir / "latest-release.json", force=args.force)
    download_file(archive_url, archive_path, force=args.force)

    if args.extract:
        safe_extract_tar_zst(archive_path, args.output_dir)

    nodes = args.output_dir / "nodes.jsonl"
    edges = args.output_dir / "edges.jsonl"
    if args.extract and not (nodes.exists() and edges.exists()):
        extracted = sorted(args.output_dir.rglob("*"))
        sample = "\n".join(str(path) for path in extracted[:20])
        raise FileNotFoundError(
            "Extraction finished, but nodes.jsonl and edges.jsonl were not found at "
            f"{args.output_dir}. First extracted paths:\n{sample}"
        )

    if shutil.which("du"):
        print(f"KGX directory: {args.output_dir}")


if __name__ == "__main__":
    main()
