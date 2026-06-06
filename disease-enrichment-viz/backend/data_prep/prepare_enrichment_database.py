#!/usr/bin/env python3

import json
import sqlite3
import argparse
from pathlib import Path
import pandas as pd
from collections import defaultdict

class EnrichmentDatabaseBuilder:
    def __init__(self, db_path='../data/enrichment_database.db'):
        self.db_path = db_path
        self.conn = None

    def create_database(self):
        """Create SQLite database with optimized schema"""
        print(f"Creating database: {self.db_path}")

        # Ensure directory exists
        Path(self.db_path).parent.mkdir(parents=True, exist_ok=True)

        self.conn = sqlite3.connect(self.db_path)
        cursor = self.conn.cursor()

        # Create tables
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS diseases (
                mondo_id TEXT PRIMARY KEY,
                name TEXT NOT NULL,
                description TEXT,
                gene_count INTEGER NOT NULL
            )
        ''')

        cursor.execute('''
            CREATE TABLE IF NOT EXISTS go_terms (
                go_id TEXT PRIMARY KEY,
                name TEXT NOT NULL,
                category TEXT NOT NULL,
                is_leaf BOOLEAN NOT NULL DEFAULT 0,
                information_content REAL
            )
        ''')

        cursor.execute('''
            CREATE TABLE IF NOT EXISTS go_hierarchy (
                child_id TEXT NOT NULL,
                parent_id TEXT NOT NULL,
                PRIMARY KEY (child_id, parent_id),
                FOREIGN KEY (child_id) REFERENCES go_terms(go_id),
                FOREIGN KEY (parent_id) REFERENCES go_terms(go_id)
            )
        ''')

        cursor.execute('''
            CREATE TABLE IF NOT EXISTS enrichment_results (
                mondo_id TEXT NOT NULL,
                go_id TEXT NOT NULL,
                p_value REAL NOT NULL,
                category TEXT NOT NULL,
                rank_in_category INTEGER,
                PRIMARY KEY (mondo_id, go_id),
                FOREIGN KEY (mondo_id) REFERENCES diseases(mondo_id),
                FOREIGN KEY (go_id) REFERENCES go_terms(go_id)
            )
        ''')

        # Create indexes for fast queries
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_enrichment_disease ON enrichment_results(mondo_id)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_enrichment_pvalue ON enrichment_results(p_value)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_enrichment_category ON enrichment_results(category)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_hierarchy_child ON go_hierarchy(child_id)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_hierarchy_parent ON go_hierarchy(parent_id)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_go_category ON go_terms(category)')

        self.conn.commit()
        print("Database schema created successfully")

    def load_go_hierarchy(self, hierarchy_file='../data/go_hierarchy.json'):
        """Load GO hierarchy data into database"""
        print(f"Loading GO hierarchy from {hierarchy_file}")

        with open(hierarchy_file, 'r') as f:
            hierarchy_data = json.load(f)

        cursor = self.conn.cursor()

        # Insert GO terms
        print("Inserting GO terms...")
        leaf_nodes = set(hierarchy_data['leaf_nodes'])

        # We need term names - load from our cached names
        term_names = {}
        try:
            with open('../../robokop_term_names_cache.json', 'r') as f:
                term_names = json.load(f)
        except FileNotFoundError:
            print("Warning: Term names cache not found, using GO IDs as names")

        # Categorize GO terms based on root ancestors
        go_categories = self._categorize_go_terms(hierarchy_data)

        go_terms_to_insert = []
        for term in hierarchy_data['terms']:
            name = term_names.get(term, term)  # Use GO ID if name not found
            category = go_categories.get(term, 'Unknown')
            is_leaf = term in leaf_nodes

            go_terms_to_insert.append((term, name, category, is_leaf))

        cursor.executemany('''
            INSERT OR REPLACE INTO go_terms (go_id, name, category, is_leaf)
            VALUES (?, ?, ?, ?)
        ''', go_terms_to_insert)

        print(f"Inserted {len(go_terms_to_insert):,} GO terms")

        # Insert hierarchy relationships
        print("Inserting hierarchy relationships...")
        hierarchy_to_insert = [
            (rel['child'], rel['parent'])
            for rel in hierarchy_data['relationships']
        ]

        cursor.executemany('''
            INSERT OR REPLACE INTO go_hierarchy (child_id, parent_id)
            VALUES (?, ?)
        ''', hierarchy_to_insert)

        print(f"Inserted {len(hierarchy_to_insert):,} hierarchy relationships")
        self.conn.commit()

    def _categorize_go_terms(self, hierarchy_data):
        """Categorize GO terms based on their root ancestors"""
        print("Categorizing GO terms by root ancestors...")

        # Map root nodes to categories
        root_categories = {
            'GO:0008150': 'BiologicalProcess',
            'GO:0003674': 'MolecularFunction',
            'GO:0005575': 'CellularComponent',
            'GO:0015015': 'MolecularFunction'  # tRNA modification guide RNA binding - molecular function
        }

        children_map = hierarchy_data['children_map']
        go_categories = {}

        # BFS to assign categories
        from collections import deque

        for root, category in root_categories.items():
            if root in children_map:
                queue = deque([root])
                go_categories[root] = category

                while queue:
                    current = queue.popleft()
                    current_category = go_categories[current]

                    for child in children_map.get(current, []):
                        if child not in go_categories:
                            go_categories[child] = current_category
                            queue.append(child)

        return go_categories

    def load_enrichment_data(self, enrichment_file='../../fast_enrichment_results.jsonl'):
        """Load enrichment results into database"""
        print(f"Loading enrichment data from {enrichment_file}")

        cursor = self.conn.cursor()

        # First pass: collect diseases
        diseases_to_insert = []
        enrichment_to_insert = []

        processed = 0
        with open(enrichment_file, 'r') as f:
            for line in f:
                if processed % 500 == 0:
                    print(f"  Processed {processed:,} diseases")

                try:
                    result = json.loads(line.strip())

                    # Extract disease info
                    disease = result['disease']
                    mondo_id = disease['curie']
                    name = disease['name']
                    description = disease.get('description', '')
                    gene_count = result['gene_count']

                    diseases_to_insert.append((mondo_id, name, description, gene_count))

                    # Extract enrichment results
                    for category, terms in result['enrichment_results'].items():
                        category_clean = category.replace('biolink:', '')

                        for rank, term in enumerate(terms, 1):
                            go_id = term['curie']
                            p_value = term['p_value']

                            enrichment_to_insert.append((
                                mondo_id, go_id, p_value, category_clean, rank
                            ))

                    processed += 1

                except json.JSONDecodeError:
                    continue

        print(f"Processed {processed:,} diseases")

        # Insert diseases
        print("Inserting diseases...")
        cursor.executemany('''
            INSERT OR REPLACE INTO diseases (mondo_id, name, description, gene_count)
            VALUES (?, ?, ?, ?)
        ''', diseases_to_insert)

        print(f"Inserted {len(diseases_to_insert):,} diseases")

        # Insert enrichment results
        print("Inserting enrichment results...")
        cursor.executemany('''
            INSERT OR REPLACE INTO enrichment_results (mondo_id, go_id, p_value, category, rank_in_category)
            VALUES (?, ?, ?, ?, ?)
        ''', enrichment_to_insert)

        print(f"Inserted {len(enrichment_to_insert):,} enrichment results")
        self.conn.commit()

    def load_information_content(self, ic_file='../../robokop_information_content_cache.json'):
        """Load information content values for GO terms"""
        print(f"Loading information content from {ic_file}")

        try:
            with open(ic_file, 'r') as f:
                ic_data = json.load(f)
        except FileNotFoundError:
            print("Information content file not found, skipping...")
            return

        cursor = self.conn.cursor()

        # Update GO terms with information content
        ic_updates = []
        for go_id, ic_value in ic_data.items():
            if go_id.startswith('GO:'):
                ic_updates.append((ic_value, go_id))

        cursor.executemany('''
            UPDATE go_terms
            SET information_content = ?
            WHERE go_id = ?
        ''', ic_updates)

        print(f"Updated information content for {cursor.rowcount:,} GO terms")
        self.conn.commit()

    def create_summary_stats(self):
        """Create summary statistics and verify database"""
        print("\nCreating summary statistics...")

        cursor = self.conn.cursor()

        # Database stats
        stats = {}

        cursor.execute("SELECT COUNT(*) FROM diseases")
        stats['diseases'] = cursor.fetchone()[0]

        cursor.execute("SELECT COUNT(*) FROM go_terms")
        stats['go_terms'] = cursor.fetchone()[0]

        cursor.execute("SELECT COUNT(*) FROM go_hierarchy")
        stats['hierarchy_relationships'] = cursor.fetchone()[0]

        cursor.execute("SELECT COUNT(*) FROM enrichment_results")
        stats['enrichment_results'] = cursor.fetchone()[0]

        cursor.execute("SELECT category, COUNT(*) FROM go_terms WHERE category != 'Unknown' GROUP BY category")
        category_counts = dict(cursor.fetchall())

        cursor.execute("""
            SELECT category, COUNT(*)
            FROM enrichment_results
            GROUP BY category
        """)
        enrichment_counts = dict(cursor.fetchall())

        print(f"\n=== Database Summary ===")
        print(f"Diseases: {stats['diseases']:,}")
        print(f"GO Terms: {stats['go_terms']:,}")
        print(f"Hierarchy Relationships: {stats['hierarchy_relationships']:,}")
        print(f"Enrichment Results: {stats['enrichment_results']:,}")

        print(f"\nGO Terms by Category:")
        for category, count in category_counts.items():
            print(f"  {category}: {count:,}")

        print(f"\nEnrichment Results by Category:")
        for category, count in enrichment_counts.items():
            print(f"  {category}: {count:,}")

        # Save stats to JSON
        stats['category_counts'] = category_counts
        stats['enrichment_counts'] = enrichment_counts

        with open('../data/database_stats.json', 'w') as f:
            json.dump(stats, f, indent=2)

        print(f"\nStats saved to ../data/database_stats.json")

    def close(self):
        if self.conn:
            self.conn.close()

def main():
    parser = argparse.ArgumentParser(description='Prepare enrichment database for web app')
    parser.add_argument('--hierarchy-file', default='../data/go_hierarchy.json',
                        help='GO hierarchy JSON file')
    parser.add_argument('--enrichment-file', default='../../fast_enrichment_results.jsonl',
                        help='Enrichment results JSONL file')
    parser.add_argument('--ic-file', default='../../robokop_information_content_cache.json',
                        help='Information content cache file')
    parser.add_argument('--db-path', default='../data/enrichment_database.db',
                        help='Output database path')

    args = parser.parse_args()

    builder = EnrichmentDatabaseBuilder(args.db_path)

    try:
        # Create database schema
        builder.create_database()

        # Load GO hierarchy
        builder.load_go_hierarchy(args.hierarchy_file)

        # Load enrichment data
        builder.load_enrichment_data(args.enrichment_file)

        # Load information content
        builder.load_information_content(args.ic_file)

        # Create summary stats
        builder.create_summary_stats()

        print(f"\n✅ Database preparation complete!")
        print(f"Database saved to: {args.db_path}")

    finally:
        builder.close()

if __name__ == "__main__":
    main()