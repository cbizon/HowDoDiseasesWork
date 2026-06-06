#!/usr/bin/env python3

from flask import Flask, jsonify, request, send_from_directory, send_file
from flask_cors import CORS
import sqlite3
import json
from collections import defaultdict, deque
import logging
import os

app = Flask(__name__)
CORS(app)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

DATABASE_PATH = 'data/enrichment_database.db'

class EnrichmentAPI:
    def __init__(self, db_path=DATABASE_PATH):
        self.db_path = db_path
        self.term_names_cache = self.load_term_names_cache()

    def get_connection(self):
        """Get database connection"""
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row  # Enable dict-like access
        return conn

    def get_disease_info(self, mondo_id):
        """Get basic disease information"""
        conn = self.get_connection()
        cursor = conn.cursor()

        cursor.execute('''
            SELECT * FROM diseases
            WHERE mondo_id = ?
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

        cursor.execute('''
            SELECT e.go_id, e.p_value, e.rank_in_category, g.name, g.is_leaf, g.information_content
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
            if term['name'] == term['go_id'] and term['go_id'] in self.term_names_cache:
                term['name'] = self.term_names_cache[term['go_id']]

        return enrichment_terms

    def get_hierarchy_subgraph(self, go_terms):
        """Get hierarchy subgraph for ONLY the enrichment terms (no ancestors)"""
        if not go_terms:
            return {'terms': [], 'relationships': []}

        conn = self.get_connection()
        cursor = conn.cursor()

        # Only get information for the enrichment terms themselves
        placeholders = ','.join(['?' for _ in go_terms])
        cursor.execute(f'''
            SELECT go_id, name, is_leaf, information_content
            FROM go_terms
            WHERE go_id IN ({placeholders})
        ''', go_terms)

        print(f"Database query returned {cursor.rowcount} rows for {len(go_terms)} terms")

        terms_info = [dict(row) for row in cursor.fetchall()]

        # Fix names using cache if database names are missing
        for term in terms_info:
            if term['name'] == term['go_id'] and term['go_id'] in self.term_names_cache:
                term['name'] = self.term_names_cache[term['go_id']]

        # Get relationships ONLY between enrichment terms (no external ancestors)
        cursor.execute(f'''
            SELECT child_id, parent_id
            FROM go_hierarchy
            WHERE child_id IN ({placeholders})
            AND parent_id IN ({placeholders})
        ''', go_terms + go_terms)

        relationships = [{'child': row[0], 'parent': row[1]} for row in cursor.fetchall()]

        conn.close()

        return {
            'terms': terms_info,
            'relationships': relationships,
            'original_terms': go_terms,
            'total_terms': len(terms_info)
        }

    def search_diseases(self, query, limit=20):
        """Search diseases by name or MONDO ID"""
        conn = self.get_connection()
        cursor = conn.cursor()

        # Search by MONDO ID or name
        cursor.execute('''
            SELECT mondo_id, name, gene_count
            FROM diseases
            WHERE mondo_id LIKE ? OR LOWER(name) LIKE LOWER(?)
            ORDER BY
                CASE WHEN mondo_id = ? THEN 1 ELSE 2 END,
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
            'data_prep/robokop_term_names_cache.json',
            '../robokop_term_names_cache.json',
            '../../robokop_term_names_cache.json'
        ]

        for cache_file in cache_files:
            if os.path.exists(cache_file):
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
            go_terms = [term['go_id'] for term in enrichment_terms]
            hierarchy = api.get_hierarchy_subgraph(go_terms)
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
    """Get hierarchy for specific GO terms"""
    try:
        go_terms = request.args.get('terms', '').split(',')
        go_terms = [term.strip() for term in go_terms if term.strip()]

        if not go_terms:
            return jsonify({'error': 'Terms parameter is required'}), 400

        hierarchy = api.get_hierarchy_subgraph(go_terms)

        return jsonify(hierarchy)

    except Exception as e:
        logger.error(f"Error getting hierarchy: {e}")
        return jsonify({'error': 'Internal server error'}), 500

@app.route('/api/stats')
def get_stats():
    """Get database statistics"""
    try:
        with open('data/database_stats.json', 'r') as f:
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