#!/usr/bin/env python3
"""Extract term hierarchy data from KGX nodes and edges.

This is the graph-agnostic replacement for the older GO-only hierarchy extractor.
It keeps hierarchy where the KGX provides `biolink:subclass_of` edges between
terms in the requested output categories, and it is safe for categories where no
hierarchy exists.
"""

import argparse
import json
import os
import time
from collections import defaultdict
from pathlib import Path


DEFAULT_TERM_CATEGORIES = [
    "biolink:BiologicalProcess",
    "biolink:MolecularActivity",
    "biolink:Pathway",
]


def as_list(value):
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def primary_category(categories, target_categories):
    for category in categories:
        if category in target_categories:
            return category
    return categories[0] if categories else "Unknown"


def parse_nodes(nodes_file, target_categories):
    target_categories = set(target_categories)
    terms = {}
    total_nodes = 0

    print(f"Parsing term nodes from {nodes_file}")
    with open(nodes_file, "r", encoding="utf-8") as f:
        for line in f:
            total_nodes += 1
            if total_nodes % 100000 == 0:
                print(f"  Processed {total_nodes:,} nodes, found {len(terms):,} terms")

            try:
                node = json.loads(line)
            except json.JSONDecodeError:
                continue

            node_id = node.get("id")
            categories = as_list(node.get("category"))
            if not node_id or not any(category in target_categories for category in categories):
                continue

            terms[node_id] = {
                "id": node_id,
                "name": node.get("name", node_id),
                "categories": categories,
                "category": primary_category(categories, target_categories),
                "information_content": node.get("information_content"),
            }

    print(f"Node parsing complete: {total_nodes:,} nodes, {len(terms):,} target terms")
    return terms


def parse_hierarchy_edges(edges_file, terms):
    relationships = []
    total_edges = 0
    term_ids = set(terms)

    print(f"Parsing term hierarchy edges from {edges_file}")
    with open(edges_file, "r", encoding="utf-8") as f:
        for line in f:
            total_edges += 1
            if total_edges % 500000 == 0:
                print(
                    f"  Processed {total_edges:,} edges, "
                    f"found {len(relationships):,} term hierarchy edges"
                )

            try:
                edge = json.loads(line)
            except json.JSONDecodeError:
                continue

            if edge.get("predicate") not in {"biolink:subclass_of", "rdfs:subClassOf"}:
                continue

            subject = edge.get("subject")
            object_id = edge.get("object")
            if subject in term_ids and object_id in term_ids:
                relationships.append(
                    {
                        "child": subject,
                        "parent": object_id,
                        "predicate": edge.get("predicate"),
                    }
                )

    print(
        f"Edge parsing complete: {total_edges:,} edges, "
        f"{len(relationships):,} term hierarchy edges"
    )
    return relationships


def build_hierarchy_document(terms, relationships, source_nodes, source_edges):
    children_map = defaultdict(list)
    parents_map = defaultdict(list)

    for relationship in relationships:
        child = relationship["child"]
        parent = relationship["parent"]
        children_map[parent].append(child)
        parents_map[child].append(parent)

    term_ids = set(terms)
    root_nodes = sorted(term for term in term_ids if term in children_map and term not in parents_map)
    leaf_nodes = sorted(term for term in term_ids if term in parents_map and term not in children_map)

    term_records = []
    for term_id, term in sorted(terms.items()):
        term_records.append(
            {
                **term,
                "is_leaf": term_id in leaf_nodes,
                "has_hierarchy": term_id in children_map or term_id in parents_map,
            }
        )

    return {
        "metadata": {
            "source_nodes": str(source_nodes),
            "source_edges": str(source_edges),
            "total_terms": len(term_records),
            "total_relationships": len(relationships),
            "root_nodes": len(root_nodes),
            "leaf_nodes": len(leaf_nodes),
            "extraction_time": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        },
        "terms": term_records,
        "root_nodes": root_nodes,
        "leaf_nodes": leaf_nodes,
        "relationships": relationships,
        "children_map": {key: sorted(value) for key, value in children_map.items()},
        "parents_map": {key: sorted(value) for key, value in parents_map.items()},
    }


def main():
    parser = argparse.ArgumentParser(description="Extract generic term hierarchy from KGX.")
    parser.add_argument(
        "--graph-dir",
        default=os.environ.get("KGX_GRAPH_DIR", "/Users/bizon/Projects/ROBOKOP/graph"),
        help="Directory containing nodes.jsonl and edges.jsonl",
    )
    parser.add_argument("--nodes-file", help="Explicit KGX nodes.jsonl path")
    parser.add_argument("--edges-file", help="Explicit KGX edges.jsonl path")
    parser.add_argument(
        "--categories",
        nargs="+",
        default=DEFAULT_TERM_CATEGORIES,
        help="Biolink categories to include as hierarchical enrichment terms",
    )
    parser.add_argument(
        "--output",
        default="term_hierarchy.json",
        help="Output JSON file",
    )
    args = parser.parse_args()

    graph_dir = Path(args.graph_dir)
    nodes_file = Path(args.nodes_file) if args.nodes_file else graph_dir / "nodes.jsonl"
    edges_file = Path(args.edges_file) if args.edges_file else graph_dir / "edges.jsonl"

    terms = parse_nodes(nodes_file, args.categories)
    relationships = parse_hierarchy_edges(edges_file, terms)
    hierarchy = build_hierarchy_document(terms, relationships, nodes_file, edges_file)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as f:
        json.dump(hierarchy, f)

    print(f"Term hierarchy saved to {output}")
    print(
        f"Summary: {hierarchy['metadata']['total_terms']:,} terms, "
        f"{hierarchy['metadata']['total_relationships']:,} hierarchy edges"
    )


if __name__ == "__main__":
    main()
