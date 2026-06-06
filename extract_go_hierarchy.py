#!/usr/bin/env python3

import json
import argparse
from collections import defaultdict, Counter
import pandas as pd

def extract_go_hierarchy(edges_file='/Users/bizon/Projects/ROBOKOP/graph/edges.jsonl',
                        output_file='go_hierarchy.json'):
    """
    Extract GO term subclass relationships from ROBOKOP edges file
    Adapted from disease subclass extraction script
    """
    print(f"Processing {edges_file} to extract GO hierarchy...")

    # Track GO terms and their relationships
    go_terms = set()
    subclass_edges = []
    go_categories = ['biolink:BiologicalProcess', 'biolink:MolecularActivity', 'biolink:Pathway']

    # Statistics
    processed_edges = 0
    subclass_count = 0

    with open(edges_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line_num % 1000000 == 0:
                print(f"  Processed {line_num:,} edges, found {subclass_count:,} GO subclass relationships")

            try:
                edge = json.loads(line.strip())
                processed_edges += 1

                # Check if this is a subclass_of relationship
                predicate = edge.get('predicate', '')
                if predicate != 'biolink:subclass_of':
                    continue

                subject = edge.get('subject', '')
                object_term = edge.get('object', '')

                # Check if both subject and object are GO terms (start with GO:)
                if not (subject.startswith('GO:') and object_term.startswith('GO:')):
                    continue

                # Store the relationship
                subclass_edges.append({
                    'child': subject,
                    'parent': object_term,
                    'predicate': predicate
                })

                # Track all GO terms
                go_terms.add(subject)
                go_terms.add(object_term)
                subclass_count += 1

            except json.JSONDecodeError:
                continue

    print(f"\nExtraction complete!")
    print(f"  Total edges processed: {processed_edges:,}")
    print(f"  GO subclass relationships found: {subclass_count:,}")
    print(f"  Unique GO terms: {len(go_terms):,}")

    # Build hierarchy data structure
    print("Building hierarchy data structures...")

    # Parent -> children mapping
    children_map = defaultdict(list)
    # Child -> parents mapping
    parents_map = defaultdict(list)

    for edge in subclass_edges:
        child = edge['child']
        parent = edge['parent']
        children_map[parent].append(child)
        parents_map[child].append(parent)

    # Find root nodes (have children but no parents)
    root_nodes = []
    for term in go_terms:
        if term in children_map and term not in parents_map:
            root_nodes.append(term)

    # Find leaf nodes (have parents but no children)
    leaf_nodes = []
    for term in go_terms:
        if term in parents_map and term not in children_map:
            leaf_nodes.append(term)

    print(f"  Root nodes (no parents): {len(root_nodes):,}")
    print(f"  Leaf nodes (no children): {len(leaf_nodes):,}")
    print(f"  Internal nodes: {len(go_terms) - len(root_nodes) - len(leaf_nodes):,}")

    # Create comprehensive hierarchy data
    hierarchy_data = {
        'metadata': {
            'total_terms': len(go_terms),
            'total_relationships': len(subclass_edges),
            'root_nodes': len(root_nodes),
            'leaf_nodes': len(leaf_nodes),
            'extraction_date': pd.Timestamp.now().isoformat()
        },
        'terms': list(go_terms),
        'root_nodes': root_nodes,
        'leaf_nodes': leaf_nodes,
        'relationships': subclass_edges,
        'children_map': {k: v for k, v in children_map.items()},
        'parents_map': {k: v for k, v in parents_map.items()}
    }

    # Save to JSON
    print(f"Saving hierarchy data to {output_file}...")
    with open(output_file, 'w') as f:
        json.dump(hierarchy_data, f, indent=2)

    print(f"Hierarchy data saved to {output_file}")
    return hierarchy_data

def analyze_hierarchy_depth(hierarchy_data):
    """Analyze the depth distribution of the GO hierarchy"""
    print("\nAnalyzing hierarchy depth...")

    children_map = hierarchy_data['children_map']
    root_nodes = hierarchy_data['root_nodes']

    def get_max_depth(term, visited=None):
        if visited is None:
            visited = set()

        if term in visited:  # Cycle detection
            return 0

        visited.add(term)

        if term not in children_map:
            return 1

        max_child_depth = 0
        for child in children_map[term]:
            child_depth = get_max_depth(child, visited.copy())
            max_child_depth = max(max_child_depth, child_depth)

        return 1 + max_child_depth

    # Calculate depths from each root
    depths = []
    for root in root_nodes:
        depth = get_max_depth(root)
        depths.append(depth)
        print(f"  Root {root}: max depth {depth}")

    if depths:
        print(f"\nHierarchy depth statistics:")
        print(f"  Max depth: {max(depths)}")
        print(f"  Average depth: {sum(depths) / len(depths):.1f}")

def create_test_subgraph(hierarchy_data, enrichment_terms, output_file='test_subgraph.json'):
    """Create a test subgraph for a specific set of enrichment terms"""
    print(f"\nCreating test subgraph for {len(enrichment_terms)} terms...")

    # Find all ancestors and descendants of enrichment terms
    parents_map = hierarchy_data['parents_map']
    children_map = hierarchy_data['children_map']

    def get_all_ancestors(term, visited=None):
        if visited is None:
            visited = set()

        if term in visited:
            return set()

        visited.add(term)
        ancestors = set()

        if term in parents_map:
            for parent in parents_map[term]:
                ancestors.add(parent)
                ancestors.update(get_all_ancestors(parent, visited.copy()))

        return ancestors

    # Collect all relevant terms (enrichment terms + their ancestors)
    relevant_terms = set(enrichment_terms)
    for term in enrichment_terms:
        relevant_terms.update(get_all_ancestors(term))

    # Filter relationships to only include relevant terms
    relevant_relationships = [
        rel for rel in hierarchy_data['relationships']
        if rel['child'] in relevant_terms and rel['parent'] in relevant_terms
    ]

    subgraph = {
        'terms': list(relevant_terms),
        'enrichment_terms': enrichment_terms,
        'relationships': relevant_relationships,
        'metadata': {
            'original_terms': len(enrichment_terms),
            'total_terms_in_subgraph': len(relevant_terms),
            'relationships_in_subgraph': len(relevant_relationships)
        }
    }

    # Save subgraph
    with open(output_file, 'w') as f:
        json.dump(subgraph, f, indent=2)

    print(f"Test subgraph saved to {output_file}")
    print(f"  Original terms: {len(enrichment_terms)}")
    print(f"  Terms in subgraph: {len(relevant_terms)}")
    print(f"  Relationships: {len(relevant_relationships)}")

    return subgraph

def main():
    parser = argparse.ArgumentParser(description='Extract GO term hierarchy from ROBOKOP edges')
    parser.add_argument('--edges-file',
                        default='/Users/bizon/Projects/ROBOKOP/graph/edges.jsonl',
                        help='Path to ROBOKOP edges file')
    parser.add_argument('--output',
                        default='go_hierarchy.json',
                        help='Output file for hierarchy data')
    parser.add_argument('--analyze-depth', action='store_true',
                        help='Analyze hierarchy depth distribution')
    parser.add_argument('--test-terms',
                        help='Comma-separated list of GO terms for test subgraph')

    args = parser.parse_args()

    # Extract hierarchy
    hierarchy_data = extract_go_hierarchy(args.edges_file, args.output)

    # Analyze depth if requested
    if args.analyze_depth:
        analyze_hierarchy_depth(hierarchy_data)

    # Create test subgraph if terms provided
    if args.test_terms:
        test_terms = [term.strip() for term in args.test_terms.split(',')]
        create_test_subgraph(hierarchy_data, test_terms)

if __name__ == "__main__":
    main()