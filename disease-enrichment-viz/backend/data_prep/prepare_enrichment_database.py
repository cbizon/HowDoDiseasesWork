#!/usr/bin/env python3
"""Build the SQLite database used by the enrichment visualization app.

The database is intentionally term-generic. GO terms get hierarchy-aware
visualization because KGX contains GO subclass edges; non-hierarchical terms are
still loaded and returned as ranked enrichment results.
"""

import argparse
import json
import sqlite3
from pathlib import Path


class EnrichmentDatabaseBuilder:
    def __init__(self, db_path="../data/enrichment_database.db", stats_path="../data/database_stats.json"):
        self.db_path = db_path
        self.stats_path = stats_path
        self.conn = None

    def create_database(self):
        print(f"Creating database: {self.db_path}")
        db_path = Path(self.db_path)
        db_path.parent.mkdir(parents=True, exist_ok=True)
        if db_path.exists():
            db_path.unlink()

        self.conn = sqlite3.connect(db_path)
        cursor = self.conn.cursor()

        cursor.execute(
            """
            CREATE TABLE diseases (
                disease_id TEXT PRIMARY KEY,
                name TEXT NOT NULL,
                description TEXT,
                gene_count INTEGER NOT NULL
            )
            """
        )

        cursor.execute(
            """
            CREATE TABLE terms (
                term_id TEXT PRIMARY KEY,
                name TEXT NOT NULL,
                category TEXT NOT NULL,
                is_leaf BOOLEAN NOT NULL DEFAULT 0,
                has_hierarchy BOOLEAN NOT NULL DEFAULT 0,
                information_content REAL
            )
            """
        )

        cursor.execute(
            """
            CREATE TABLE term_hierarchy (
                child_id TEXT NOT NULL,
                parent_id TEXT NOT NULL,
                PRIMARY KEY (child_id, parent_id)
            )
            """
        )

        cursor.execute(
            """
            CREATE TABLE enrichment_results (
                disease_id TEXT NOT NULL,
                term_id TEXT NOT NULL,
                p_value REAL NOT NULL,
                category TEXT NOT NULL,
                rank_in_category INTEGER NOT NULL,
                PRIMARY KEY (disease_id, term_id, category)
            )
            """
        )

        cursor.execute("CREATE INDEX idx_enrichment_disease ON enrichment_results(disease_id)")
        cursor.execute("CREATE INDEX idx_enrichment_pvalue ON enrichment_results(p_value)")
        cursor.execute("CREATE INDEX idx_enrichment_category ON enrichment_results(category)")
        cursor.execute("CREATE INDEX idx_hierarchy_child ON term_hierarchy(child_id)")
        cursor.execute("CREATE INDEX idx_hierarchy_parent ON term_hierarchy(parent_id)")
        cursor.execute("CREATE INDEX idx_terms_category ON terms(category)")

        # Compatibility views for older code/docs that used GO-specific names.
        cursor.execute(
            """
            CREATE VIEW go_terms AS
            SELECT
                term_id AS go_id,
                name,
                category,
                is_leaf,
                information_content
            FROM terms
            """
        )
        cursor.execute(
            """
            CREATE VIEW go_hierarchy AS
            SELECT child_id, parent_id
            FROM term_hierarchy
            """
        )

        self.conn.commit()
        print("Database schema created successfully")

    def load_term_hierarchy(self, hierarchy_file=None):
        if not hierarchy_file:
            print("No term hierarchy file provided; hierarchy loading skipped")
            return

        hierarchy_path = Path(hierarchy_file)
        if not hierarchy_path.exists():
            print(f"Term hierarchy file not found: {hierarchy_file}; hierarchy loading skipped")
            return

        print(f"Loading term hierarchy from {hierarchy_file}")
        with hierarchy_path.open("r", encoding="utf-8") as f:
            hierarchy_data = json.load(f)

        raw_terms = hierarchy_data.get("terms", [])
        term_rows = []
        for term in raw_terms:
            if isinstance(term, str):
                term_id = term
                term_rows.append((term_id, term_id, "Unknown", 0, 1, None))
            else:
                term_id = term["id"]
                term_rows.append(
                    (
                        term_id,
                        term.get("name") or term_id,
                        self.clean_category(term.get("category") or "Unknown"),
                        bool(term.get("is_leaf")),
                        bool(term.get("has_hierarchy")),
                        term.get("information_content"),
                    )
                )

        relationship_rows = [
            (relationship["child"], relationship["parent"])
            for relationship in hierarchy_data.get("relationships", [])
        ]

        cursor = self.conn.cursor()
        cursor.executemany(
            """
            INSERT OR REPLACE INTO terms
                (term_id, name, category, is_leaf, has_hierarchy, information_content)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            term_rows,
        )
        cursor.executemany(
            """
            INSERT OR REPLACE INTO term_hierarchy (child_id, parent_id)
            VALUES (?, ?)
            """,
            relationship_rows,
        )
        self.conn.commit()
        print(f"Inserted {len(term_rows):,} terms")
        print(f"Inserted {len(relationship_rows):,} hierarchy relationships")

    def load_enrichment_data(self, enrichment_file):
        print(f"Loading enrichment data from {enrichment_file}")

        disease_rows = {}
        term_rows = {}
        enrichment_rows = []

        processed = 0
        skipped_errors = 0
        with open(enrichment_file, "r", encoding="utf-8") as f:
            for line in f:
                if not line.strip():
                    continue

                result = json.loads(line)
                processed += 1
                if processed % 500 == 0:
                    print(f"  Processed {processed:,} diseases")

                disease = result["disease"]
                disease_id = disease["curie"]
                disease_rows[disease_id] = (
                    disease_id,
                    disease.get("name", disease_id),
                    disease.get("description", ""),
                    result["gene_count"],
                )

                if result.get("error") or result.get("enrichment_errors"):
                    skipped_errors += 1

                for category, terms in result.get("enrichment_results", {}).items():
                    category_clean = self.clean_category(category)
                    sorted_terms = sorted(
                        terms,
                        key=lambda term: term.get("p_value") if term.get("p_value") is not None else float("inf"),
                    )

                    for rank, term in enumerate(sorted_terms, 1):
                        term_id = term["curie"]
                        p_value = term.get("p_value")
                        if p_value is None:
                            continue

                        term_rows.setdefault(
                            term_id,
                            (
                                term_id,
                                term.get("name") or term_id,
                                self.clean_category(term.get("category") or category_clean),
                                0,
                                0,
                                None,
                            ),
                        )
                        enrichment_rows.append((disease_id, term_id, p_value, category_clean, rank))

        cursor = self.conn.cursor()
        cursor.executemany(
            """
            INSERT OR REPLACE INTO diseases (disease_id, name, description, gene_count)
            VALUES (?, ?, ?, ?)
            """,
            list(disease_rows.values()),
        )
        cursor.executemany(
            """
            INSERT OR IGNORE INTO terms
                (term_id, name, category, is_leaf, has_hierarchy, information_content)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            list(term_rows.values()),
        )
        cursor.executemany(
            """
            INSERT OR REPLACE INTO enrichment_results
                (disease_id, term_id, p_value, category, rank_in_category)
            VALUES (?, ?, ?, ?, ?)
            """,
            enrichment_rows,
        )

        # Terms that participate in a hierarchy are not standalone even if they
        # were first inserted from enrichment output.
        cursor.execute(
            """
            UPDATE terms
            SET has_hierarchy = 1
            WHERE term_id IN (
                SELECT child_id FROM term_hierarchy
                UNION
                SELECT parent_id FROM term_hierarchy
            )
            """
        )

        self.conn.commit()
        print(f"Processed {processed:,} diseases")
        print(f"Inserted {len(disease_rows):,} diseases")
        print(f"Inserted or preserved {len(term_rows):,} enriched terms")
        print(f"Inserted {len(enrichment_rows):,} enrichment rows")
        if skipped_errors:
            print(f"Note: {skipped_errors:,} disease rows contained recorded errors")

    def load_information_content(self, ic_file=None):
        if not ic_file:
            print("No information content file provided; IC loading skipped")
            return

        ic_path = Path(ic_file)
        if not ic_path.exists():
            print(f"Information content file not found: {ic_file}; IC loading skipped")
            return

        print(f"Loading information content from {ic_file}")
        with ic_path.open("r", encoding="utf-8") as f:
            ic_data = json.load(f)

        cursor = self.conn.cursor()
        cursor.executemany(
            """
            UPDATE terms
            SET information_content = ?
            WHERE term_id = ?
            """,
            [(ic_value, term_id) for term_id, ic_value in ic_data.items()],
        )
        self.conn.commit()
        print(f"Updated information content for {cursor.rowcount:,} terms")

    def create_summary_stats(self):
        print("\nCreating summary statistics...")
        cursor = self.conn.cursor()
        stats = {}

        for table, key in [
            ("diseases", "diseases"),
            ("terms", "terms"),
            ("term_hierarchy", "hierarchy_relationships"),
            ("enrichment_results", "enrichment_results"),
        ]:
            cursor.execute(f"SELECT COUNT(*) FROM {table}")
            stats[key] = cursor.fetchone()[0]

        cursor.execute("SELECT category, COUNT(*) FROM terms GROUP BY category")
        stats["term_category_counts"] = dict(cursor.fetchall())

        cursor.execute("SELECT category, COUNT(*) FROM enrichment_results GROUP BY category")
        stats["enrichment_counts"] = dict(cursor.fetchall())

        cursor.execute("SELECT COUNT(*) FROM terms WHERE has_hierarchy = 1")
        stats["terms_with_hierarchy"] = cursor.fetchone()[0]

        cursor.execute("SELECT COUNT(*) FROM terms WHERE has_hierarchy = 0")
        stats["terms_without_hierarchy"] = cursor.fetchone()[0]

        print("\n=== Database Summary ===")
        print(f"Diseases: {stats['diseases']:,}")
        print(f"Terms: {stats['terms']:,}")
        print(f"Hierarchy Relationships: {stats['hierarchy_relationships']:,}")
        print(f"Enrichment Results: {stats['enrichment_results']:,}")
        print(f"Terms with hierarchy: {stats['terms_with_hierarchy']:,}")
        print(f"Terms without hierarchy: {stats['terms_without_hierarchy']:,}")

        stats_path = Path(self.stats_path)
        stats_path.parent.mkdir(parents=True, exist_ok=True)
        with stats_path.open("w", encoding="utf-8") as f:
            json.dump(stats, f, indent=2)

        print(f"\nStats saved to {stats_path}")

    @staticmethod
    def clean_category(category):
        return category.replace("biolink:", "")

    def close(self):
        if self.conn:
            self.conn.close()


def main():
    parser = argparse.ArgumentParser(description="Prepare enrichment database for web app")
    parser.add_argument(
        "--term-hierarchy-file",
        "--hierarchy-file",
        dest="term_hierarchy_file",
        default="../data/term_hierarchy.json",
        help="Generic term hierarchy JSON file",
    )
    parser.add_argument(
        "--enrichment-file",
        default="../../fast_enrichment_results.jsonl",
        help="Enrichment results JSONL file",
    )
    parser.add_argument(
        "--ic-file",
        default="../../robokop_information_content_cache.json",
        help="Optional information content cache JSON file",
    )
    parser.add_argument(
        "--db-path",
        default="../data/enrichment_database.db",
        help="Output SQLite database path",
    )
    parser.add_argument(
        "--stats-path",
        default="../data/database_stats.json",
        help="Output database statistics JSON path",
    )

    args = parser.parse_args()

    builder = EnrichmentDatabaseBuilder(args.db_path, args.stats_path)

    try:
        builder.create_database()
        builder.load_term_hierarchy(args.term_hierarchy_file)
        builder.load_enrichment_data(args.enrichment_file)
        builder.load_information_content(args.ic_file)
        builder.create_summary_stats()

        print("\nDatabase preparation complete")
        print(f"Database saved to: {args.db_path}")

    finally:
        builder.close()


if __name__ == "__main__":
    main()
