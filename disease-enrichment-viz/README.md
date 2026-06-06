# Disease Enrichment Visualization Web App

A web application for visualizing disease enrichment results as interactive hierarchical networks of GO terms.

## Features

- **Disease Search**: Search for diseases by MONDO ID or name
- **Interactive Network**: Hierarchical visualization of enriched GO terms
- **Customizable Parameters**: Adjustable p-value thresholds and categories
- **Rich Tooltips**: Detailed information on hover
- **Zoom & Pan**: Navigate large networks easily
- **Legend & Statistics**: Clear visual guides and network statistics

## Project Structure

```
disease-enrichment-viz/
├── backend/
│   ├── app.py                 # Flask API server
│   ├── data/
│   │   ├── enrichment_database.db    # SQLite database
│   │   ├── go_hierarchy.json         # GO hierarchy data
│   │   └── database_stats.json       # Database statistics
│   └── data_prep/
│       └── prepare_enrichment_database.py  # Database setup script
├── frontend/
│   ├── index.html            # Main web interface
│   └── src/
│       └── app.js           # JavaScript application
└── README.md
```

## Database Contents

- **4,221 diseases** with enrichment data
- **40,178 GO terms** with full hierarchy relationships
- **392,668 hierarchy relationships** (subclass_of edges)
- **5.7M enrichment results** across all categories
- **Information content values** for all GO terms

## Quick Start

### 1. Start the Backend Server

```bash
cd backend
python app.py
```

The API will be available at `http://localhost:5000`

### 2. Open the Frontend

Open `frontend/index.html` in your web browser, or serve it with a local server:

```bash
cd frontend
python -m http.server 8080
```

Then visit `http://localhost:8080`

### 3. Use the Application

1. **Search for a disease**: Type a MONDO ID (e.g., `MONDO:0005737`) or disease name
2. **Select disease**: Click on a search result to select it
3. **Adjust parameters**: Set p-value threshold and category
4. **Visualize**: Click "Visualize Network" to generate the hierarchical network
5. **Explore**: Use mouse to pan, zoom, and hover for details

## API Endpoints

- `GET /api/disease/<mondo_id>` - Get disease information
- `GET /api/disease/<mondo_id>/enrichment` - Get enrichment data with hierarchy
- `GET /api/search?q=<query>` - Search diseases
- `GET /api/hierarchy?terms=<go_ids>` - Get hierarchy for specific terms
- `GET /api/stats` - Get database statistics

## Example Usage

### Search for Huntington Disease
1. Enter "huntington" or "MONDO:0005737"
2. Select the disease from results
3. Set p-value to 1e-5
4. Choose "BiologicalProcess" category
5. Click "Visualize Network"

### Network Interpretation

- **Yellow nodes**: Original enrichment terms found for the disease
- **Node size**: Represents -log10(p-value) - larger = more significant
- **Red borders**: Leaf nodes (most specific terms)
- **Teal borders**: Internal nodes (more general terms)
- **Blue borders**: Original enrichment terms (thick border)
- **Gray lines**: Subclass relationships (child → parent)

## Technical Details

### Backend (Flask + SQLite)
- SQLite database for fast queries
- Optimized hierarchy traversal algorithms
- RESTful API with CORS support
- Error handling and logging

### Frontend (D3.js + Bootstrap)
- Force-directed network layout
- Interactive zoom and pan
- Real-time search with debouncing
- Responsive design with Bootstrap

### Data Processing
- GO hierarchy extraction from ROBOKOP edges
- Enrichment data integration
- Information content integration
- Optimized database schema with indexes

## Performance

- **Database size**: ~200MB SQLite file
- **API response times**: <1 second for most queries
- **Network rendering**: Handles 500+ nodes smoothly
- **Memory usage**: ~100MB for typical networks

## Limitations

- Currently supports BiologicalProcess, MolecularActivity, and Pathway categories
- Network layout may be cluttered for >1000 nodes
- Requires modern web browser with JavaScript enabled
- Backend must be running locally (not deployed)

## Future Enhancements

- Multi-disease comparison
- Export functionality (SVG, PNG, JSON)
- Advanced filtering options
- Pathway integration
- Deployment-ready configuration
- Real-time collaboration features

## Data Sources

- **ROBOKOP Knowledge Graph**: Disease-gene associations and GO hierarchy
- **AnswerCoalesce API**: Enrichment analysis results
- **Information Content**: From ROBOKOP node annotations