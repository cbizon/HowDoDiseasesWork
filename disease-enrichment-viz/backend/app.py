#!/usr/bin/env python3

from flask import Flask, jsonify, request, send_from_directory, send_file
from flask_cors import CORS
import sqlite3
import json
import logging
import os
from pathlib import Path

app = Flask(__name__)
CORS(app)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

APP_DIR = Path(__file__).resolve().parent
DATABASE_PATH = os.environ.get(
    "ENRICHMENT_DATABASE_PATH",
    str(APP_DIR / "data" / "enrichment_database.db"),
)
STATS_PATH = os.environ.get(
    "ENRICHMENT_DATABASE_STATS_PATH",
    str(APP_DIR / "data" / "database_stats.json"),
)

class EnrichmentAPI:
    def __init__(self, db_path=DATABASE_PATH):
        self.db_path = db_path
        self.term_names_cache = self.load_term_names_cache()
        self.schema = self.detect_schema()

    def get_connection(self):
        """Get database connection"""
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row  # Enable dict-like access
        return conn

    def detect_schema(self):
        """Detect whether the database uses generic term tables or legacy GO tables."""
        try:
            conn = self.get_connection()
            cursor = conn.cursor()
            cursor.execute(
                """
                SELECT name
                FROM sqlite_master
                WHERE type IN ('table', 'view')
                """
            )
            names = {row[0] for row in cursor.fetchall()}
            conn.close()
        except Exception as e:
            logger.warning(f"Could not inspect database schema: {e}")
            return "generic"

        if "terms" in names and "term_hierarchy" in names:
            return "generic"
        return "legacy_go"

    def get_disease_info(self, mondo_id):
        """Get basic disease information"""
        conn = self.get_connection()
        cursor = conn.cursor()

        id_column = "disease_id" if self.schema == "generic" else "mondo_id"

        cursor.execute(f'''
            SELECT
                {id_column} AS mondo_id,
                name,
                description,
                gene_count
            FROM diseases
            WHERE {id_column} = ?
        ''', (mondo_id,))

        result = cursor.fetchone()
        conn.close()

        if result:
            return dict(result)
        return None

    def get_enrichment_terms(self, mondo_id, p_threshold=1e-5, category='BiologicalProcess'):
        """Get enrichment terms for a disease below p-value threshold"""
        conn = self.get_connection()
        cursor = conn.cursor()

        if self.schema == "generic":
            cursor.execute('''
                SELECT
                    e.term_id,
                    e.term_id AS go_id,
                    e.p_value,
                    e.rank_in_category,
                    t.name,
                    t.category AS term_category,
                    t.is_leaf,
                    t.has_hierarchy,
                    t.information_content
                FROM enrichment_results e
                LEFT JOIN terms t ON e.term_id = t.term_id
                WHERE e.disease_id = ?
                AND e.category = ?
                AND e.p_value <= ?
                ORDER BY e.p_value ASC
            ''', (mondo_id, category, p_threshold))
        else:
            cursor.execute('''
                SELECT
                    e.go_id AS term_id,
                    e.go_id,
                    e.p_value,
                    e.rank_in_category,
                    g.name,
                    g.category AS term_category,
                    g.is_leaf,
                    1 AS has_hierarchy,
                    g.information_content
                FROM enrichment_results e
                JOIN go_terms g ON e.go_id = g.go_id
                WHERE e.mondo_id = ?
                AND e.category = ?
                AND e.p_value <= ?
                ORDER BY e.p_value ASC
            ''', (mondo_id, category, p_threshold))

        results = cursor.fetchall()
        conn.close()

        # Fix names using cache if database names are missing
        enrichment_terms = [dict(row) for row in results]
        for term in enrichment_terms:
            if not term.get('name'):
                term['name'] = term['term_id']
            if term['name'] == term['term_id'] and term['term_id'] in self.term_names_cache:
                term['name'] = self.term_names_cache[term['term_id']]

        return enrichment_terms

    def get_hierarchy_subgraph(self, term_ids):
        """Get hierarchy subgraph for ONLY the enrichment terms (no ancestors)"""
        if not term_ids:
            return {'terms': [], 'relationships': []}

        conn = self.get_connection()
        cursor = conn.cursor()

        # Only get information for the enrichment terms themselves
        placeholders = ','.join(['?' for _ in term_ids])
        if self.schema == "generic":
            cursor.execute(f'''
                SELECT
                    term_id,
                    term_id AS go_id,
                    name,
                    category AS term_category,
                    is_leaf,
                    has_hierarchy,
                    information_content
                FROM terms
                WHERE term_id IN ({placeholders})
            ''', term_ids)
        else:
            cursor.execute(f'''
                SELECT
                    go_id AS term_id,
                    go_id,
                    name,
                    category AS term_category,
                    is_leaf,
                    1 AS has_hierarchy,
                    information_content
                FROM go_terms
                WHERE go_id IN ({placeholders})
            ''', term_ids)

        print(f"Database query returned {cursor.rowcount} rows for {len(term_ids)} terms")

        terms_info = [dict(row) for row in cursor.fetchall()]

        # Fix names using cache if database names are missing
        for term in terms_info:
            if term['name'] == term['term_id'] and term['term_id'] in self.term_names_cache:
                term['name'] = self.term_names_cache[term['term_id']]

        # Get relationships ONLY between enrichment terms (no external ancestors)
        hierarchy_table = "term_hierarchy" if self.schema == "generic" else "go_hierarchy"
        cursor.execute(f'''
            SELECT child_id, parent_id
            FROM {hierarchy_table}
            WHERE child_id IN ({placeholders})
            AND parent_id IN ({placeholders})
        ''', term_ids + term_ids)

        relationships = [{'child': row[0], 'parent': row[1]} for row in cursor.fetchall()]

        conn.close()

        return {
            'terms': terms_info,
            'relationships': relationships,
            'original_terms': term_ids,
            'total_terms': len(terms_info),
            'hierarchy_available': bool(relationships)
        }

    def search_diseases(self, query, limit=20):
        """Search diseases by name or MONDO ID"""
        conn = self.get_connection()
        cursor = conn.cursor()

        # Search by MONDO ID or name
        id_column = "disease_id" if self.schema == "generic" else "mondo_id"

        cursor.execute(f'''
            SELECT {id_column} AS mondo_id, name, gene_count
            FROM diseases
            WHERE {id_column} LIKE ? OR LOWER(name) LIKE LOWER(?)
            ORDER BY
                CASE WHEN {id_column} = ? THEN 1 ELSE 2 END,
                CASE WHEN LOWER(name) LIKE LOWER(?) THEN 1 ELSE 2 END,
                name
            LIMIT ?
        ''', (f'%{query}%', f'%{query}%', query, f'{query}%', limit))

        results = cursor.fetchall()
        conn.close()

        return [dict(row) for row in results]

    def load_term_names_cache(self):
        """Load term names from cache file"""
        cache_files = [
            APP_DIR / 'data_prep' / 'robokop_term_names_cache.json',
            APP_DIR.parent / 'robokop_term_names_cache.json',
            APP_DIR.parent.parent / 'robokop_term_names_cache.json',
        ]

        for cache_file in cache_files:
            if cache_file.exists():
                try:
                    with open(cache_file, 'r') as f:
                        term_names = json.load(f)
                        logger.info(f"Loaded {len(term_names):,} term names from {cache_file}")
                        return term_names
                except Exception as e:
                    logger.warning(f"Failed to load {cache_file}: {e}")

        logger.warning("No term names cache file found")
        return {}

# Initialize API
api = EnrichmentAPI()

@app.route('/')
def index():
    """Serve the main web application"""
    frontend_path = os.path.join(os.path.dirname(__file__), '..', 'frontend')
    return send_file(os.path.join(frontend_path, 'index.html'))

@app.route('/src/<path:filename>')
def serve_js(filename):
    """Serve JavaScript files"""
    frontend_path = os.path.join(os.path.dirname(__file__), '..', 'frontend', 'src')
    return send_from_directory(frontend_path, filename)

@app.route('/api')
def api_status():
    """API status endpoint"""
    return jsonify({
        'status': 'Disease Enrichment Visualization API',
        'version': '1.0.0',
        'endpoints': {
            'disease_info': '/api/disease/<mondo_id>',
            'disease_enrichment': '/api/disease/<mondo_id>/enrichment',
            'search': '/api/search',
            'hierarchy': '/api/hierarchy'
        }
    })

@app.route('/api/disease/<mondo_id>')
def get_disease(mondo_id):
    """Get disease information"""
    try:
        disease_info = api.get_disease_info(mondo_id)
        if not disease_info:
            return jsonify({'error': 'Disease not found'}), 404

        return jsonify(disease_info)

    except Exception as e:
        logger.error(f"Error getting disease {mondo_id}: {e}")
        return jsonify({'error': 'Internal server error'}), 500

@app.route('/api/disease/<mondo_id>/enrichment')
def get_disease_enrichment(mondo_id):
    """Get disease enrichment data with hierarchy"""
    try:
        # Parse query parameters
        p_threshold = float(request.args.get('p_threshold', 1e-5))
        category = request.args.get('category', 'BiologicalProcess')
        include_hierarchy = request.args.get('include_hierarchy', 'true').lower() == 'true'

        # Get disease info
        disease_info = api.get_disease_info(mondo_id)
        if not disease_info:
            return jsonify({'error': 'Disease not found'}), 404

        # Get enrichment terms
        enrichment_terms = api.get_enrichment_terms(mondo_id, p_threshold, category)

        result = {
            'disease': disease_info,
            'enrichment_terms': enrichment_terms,
            'parameters': {
                'p_threshold': p_threshold,
                'category': category,
                'terms_found': len(enrichment_terms)
            }
        }

        # Add hierarchy if requested
        if include_hierarchy and enrichment_terms:
            term_ids = [term['term_id'] for term in enrichment_terms]
            hierarchy = api.get_hierarchy_subgraph(term_ids)
            result['hierarchy'] = hierarchy

        return jsonify(result)

    except Exception as e:
        logger.error(f"Error getting enrichment for {mondo_id}: {e}")
        return jsonify({'error': 'Internal server error'}), 500

@app.route('/api/search')
def search_diseases():
    """Search diseases"""
    try:
        query = request.args.get('q', '').strip()
        limit = int(request.args.get('limit', 20))

        if not query:
            return jsonify({'error': 'Query parameter q is required'}), 400

        if len(query) < 2:
            return jsonify({'error': 'Query must be at least 2 characters'}), 400

        results = api.search_diseases(query, limit)

        return jsonify({
            'query': query,
            'results': results,
            'count': len(results)
        })

    except Exception as e:
        logger.error(f"Error searching diseases: {e}")
        return jsonify({'error': 'Internal server error'}), 500

@app.route('/api/hierarchy')
def get_hierarchy():
    """Get hierarchy for specific terms"""
    try:
        term_ids = request.args.get('terms', '').split(',')
        term_ids = [term.strip() for term in term_ids if term.strip()]

        if not term_ids:
            return jsonify({'error': 'Terms parameter is required'}), 400

        hierarchy = api.get_hierarchy_subgraph(term_ids)

        return jsonify(hierarchy)

    except Exception as e:
        logger.error(f"Error getting hierarchy: {e}")
        return jsonify({'error': 'Internal server error'}), 500

@app.route('/api/stats')
def get_stats():
    """Get database statistics"""
    try:
        with open(STATS_PATH, 'r') as f:
            stats = json.load(f)
        return jsonify(stats)
    except Exception as e:
        logger.error(f"Error getting stats: {e}")
        return jsonify({'error': 'Internal server error'}), 500

if __name__ == '__main__':
    # Test database connection
    try:
        test_api = EnrichmentAPI()
        conn = test_api.get_connection()
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM diseases")
        disease_count = cursor.fetchone()[0]
        conn.close()

        logger.info(f"Database connected successfully - {disease_count} diseases loaded")

    except Exception as e:
        logger.error(f"Database connection failed: {e}")
        exit(1)

    app.run(debug=True, host='0.0.0.0', port=5000)
