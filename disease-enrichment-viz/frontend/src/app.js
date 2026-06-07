// Disease Enrichment Visualization App
class DiseaseEnrichmentViz {
    constructor() {
        this.apiBase = 'http://localhost:5000/api';
        this.currentDisease = null;
        this.currentData = null;

        this.svg = null;
        this.width = 800;
        this.height = 600;

        this.simulation = null;

        this.initializeEventListeners();
        this.setupVisualization();
    }

    initializeEventListeners() {
        // Search input
        const searchInput = document.getElementById('diseaseSearch');
        let searchTimeout;

        searchInput.addEventListener('input', (e) => {
            clearTimeout(searchTimeout);
            searchTimeout = setTimeout(() => {
                this.searchDiseases(e.target.value);
            }, 300);
        });

        // Visualize button
        document.getElementById('visualizeBtn').addEventListener('click', () => {
            this.loadEnrichmentData();
        });

        // Parameter changes
        document.getElementById('pThreshold').addEventListener('change', () => {
            if (this.currentDisease) {
                this.loadEnrichmentData();
            }
        });

        document.getElementById('category').addEventListener('change', () => {
            if (this.currentDisease) {
                this.loadEnrichmentData();
            }
        });
    }

    setupVisualization() {
        const container = d3.select('#networkViz');
        container.selectAll('*').remove();

        // Make SVG larger to accommodate multiple components
        this.width = 1200;
        this.height = 800;

        this.svg = container.append('svg')
            .attr('width', '100%')
            .attr('height', this.height)
            .attr('viewBox', `0 0 ${this.width} ${this.height}`)
            .style('background', '#fafafa');

        this.defineArrowMarkers();

        // Add zoom behavior
        const zoom = d3.zoom()
            .scaleExtent([0.1, 4])
            .on('zoom', (event) => {
                this.svg.select('.zoom-group')
                    .attr('transform', event.transform);
            });

        this.svg.call(zoom);

        // Create zoom group for all content
        this.svg.append('g').attr('class', 'zoom-group');
    }

    defineArrowMarkers() {
        const defs = this.svg.append('defs');

        defs.append('marker')
            .attr('id', 'child-to-parent-arrow')
            .attr('viewBox', '0 -5 10 10')
            .attr('refX', 10)
            .attr('refY', 0)
            .attr('markerWidth', 7)
            .attr('markerHeight', 7)
            .attr('orient', 'auto')
            .attr('markerUnits', 'strokeWidth')
            .append('path')
            .attr('d', 'M0,-5L10,0L0,5')
            .attr('fill', '#666');
    }

    async searchDiseases(query) {
        if (query.length < 2) {
            document.getElementById('searchResults').innerHTML = '';
            return;
        }

        try {
            const response = await fetch(`${this.apiBase}/search?q=${encodeURIComponent(query)}`);
            const data = await response.json();

            this.displaySearchResults(data.results);
        } catch (error) {
            console.error('Search error:', error);
            this.showError('Search failed. Please try again.');
        }
    }

    displaySearchResults(results) {
        const container = document.getElementById('searchResults');

        if (results.length === 0) {
            container.innerHTML = '<small class=\"text-muted\">No results found</small>';
            return;
        }

        const html = results.map(disease => `
            <div class=\"border rounded p-2 mb-2 cursor-pointer\"
                 onclick=\"app.selectDisease('${disease.mondo_id}', '${disease.name}')\">
                <strong>${disease.mondo_id}</strong><br>
                <small>${disease.name} (${disease.gene_count} genes)</small>
            </div>
        `).join('');

        container.innerHTML = html;
    }

    selectDisease(mondoId, name) {
        this.currentDisease = { mondo_id: mondoId, name: name };

        document.getElementById('diseaseSearch').value = `${mondoId} - ${name}`;
        document.getElementById('searchResults').innerHTML = '';
        document.getElementById('visualizeBtn').disabled = false;

        console.log('Selected disease:', this.currentDisease);
    }

    async loadEnrichmentData() {
        if (!this.currentDisease) return;

        const pThresholdText = document.getElementById('pThreshold').value;
        const category = document.getElementById('category').value;

        // Parse p-value threshold
        let pThreshold;
        try {
            pThreshold = parseFloat(pThresholdText);
            if (isNaN(pThreshold) || pThreshold <= 0 || pThreshold > 1) {
                throw new Error('Invalid p-value');
            }
        } catch (error) {
            this.showError('Please enter a valid p-value (e.g., 1e-5, 0.001, 1e-3)');
            return;
        }

        this.showLoading(true);
        this.hideError();

        try {
            const url = `${this.apiBase}/disease/${this.currentDisease.mondo_id}/enrichment?p_threshold=${pThreshold}&category=${category}&include_hierarchy=true`;
            console.log('Fetching:', url);

            const response = await fetch(url);
            const data = await response.json();

            if (response.ok) {
                this.currentData = data;
                this.displayDiseaseInfo(data.disease, data.parameters);
                this.displayNetworkStats(data);
                this.visualizeNetwork(data);
            } else {
                this.showError(data.error || 'Failed to load enrichment data');
            }
        } catch (error) {
            console.error('Load error:', error);
            this.showError('Failed to load data. Please check if the server is running.');
        } finally {
            this.showLoading(false);
        }
    }

    displayDiseaseInfo(disease, parameters) {
        const container = document.getElementById('diseaseInfo');

        container.innerHTML = `
            <h5>${disease.name}</h5>
            <p><strong>MONDO ID:</strong> ${disease.mondo_id}</p>
            <p><strong>Genes:</strong> ${disease.gene_count}</p>
            <p><strong>Parameters:</strong> p ≤ ${parameters.p_threshold}, ${parameters.category}, ${parameters.terms_found} terms found</p>
            ${disease.description ? `<p><strong>Description:</strong> ${disease.description}</p>` : ''}
        `;

        container.style.display = 'block';
    }

    displayNetworkStats(data) {
        const container = document.getElementById('networkStats');

        if (!data.hierarchy) {
            container.style.display = 'none';
            return;
        }

        const stats = {
            'Enrichment Terms': data.enrichment_terms.length,
            'Network Nodes': data.hierarchy.total_terms,
            'Relationships': data.hierarchy.relationships.length
        };

        const html = Object.entries(stats).map(([label, value]) => `
            <div class=\"stat-item\">
                <div class=\"stat-value\">${value}</div>
                <div>${label}</div>
            </div>
        `).join('');

        container.innerHTML = html;
        container.style.display = 'flex';
    }

    visualizeNetwork(data) {
        if (!data.hierarchy || data.hierarchy.terms.length === 0) {
            this.showError('No enrichment terms available for visualization');
            return;
        }

        console.log('Visualizing network:', data.hierarchy);

        // Clear previous visualization
        this.svg.select('.zoom-group').selectAll('*').remove();

        // Prepare data
        const nodes = this.prepareNodes(data);
        const links = this.prepareLinks(data.hierarchy.relationships, nodes);

        // Add root node and connect orphan nodes
        this.addRootNode(nodes, links, data.hierarchy.relationships);

        console.log('Nodes:', nodes.length, 'Links:', links.length);

        // Create force-directed layout
        this.createForceLayout(nodes, links);
    }

    prepareNodes(data) {
        const enrichmentTerms = new Set(data.enrichment_terms.map(t => t.go_id));
        const enrichmentData = {};
        data.enrichment_terms.forEach(term => {
            enrichmentData[term.go_id] = term;
        });

        // Calculate enrichment-specific leaf nodes
        const enrichmentLeaves = this.calculateEnrichmentLeaves(data.enrichment_terms, data.hierarchy.relationships);

        return data.hierarchy.terms.map(term => {
            const isOriginal = enrichmentTerms.has(term.go_id);
            const enrichment = enrichmentData[term.go_id];
            const isEnrichmentLeaf = enrichmentLeaves.has(term.go_id);

            return {
                id: term.go_id,
                name: term.name,
                is_enrichment_leaf: isEnrichmentLeaf,
                information_content: term.information_content,
                is_original: isOriginal,
                p_value: enrichment ? enrichment.p_value : null,
                rank: enrichment ? enrichment.rank_in_category : null,
                radius: this.calculateNodeSize(enrichment)
            };
        });
    }

    calculateEnrichmentLeaves(enrichmentTerms, relationships) {
        const enrichmentIds = new Set(enrichmentTerms.map(t => t.go_id));
        const enrichmentLeaves = new Set();

        // For each enrichment term, check if any of its children are also in enrichment
        for (const term of enrichmentTerms) {
            const hasEnrichedChildren = relationships.some(rel =>
                rel.parent === term.go_id && enrichmentIds.has(rel.child)
            );

            if (!hasEnrichedChildren) {
                enrichmentLeaves.add(term.go_id);
            }
        }

        return enrichmentLeaves;
    }

    prepareLinks(relationships, nodes) {
        const nodeIds = new Set(nodes.map(n => n.id));

        return relationships.filter(rel =>
            nodeIds.has(rel.child) && nodeIds.has(rel.parent)
        ).map(rel => ({
            source: rel.child,
            target: rel.parent
        }));
    }

    addRootNode(nodes, links, relationships) {
        // Find nodes that have no parents in the enrichment results (root nodes)
        const childIds = new Set(links.map(l => l.source)); // Note: source is child, target is parent
        const orphanNodes = nodes.filter(n => !childIds.has(n.id));

        console.log(`Found ${orphanNodes.length} orphan nodes (nodes without parents in enrichment results)`);

        if (orphanNodes.length > 0) {
            // Create a virtual root node
            const rootNode = {
                id: 'VIRTUAL_ROOT',
                name: 'Enrichment Results',
                is_enrichment_leaf: false,
                information_content: null,
                is_original: false,
                p_value: null,
                rank: null,
                radius: 15, // Larger root node (reduced by half)
                is_root: true
            };

            nodes.push(rootNode);

            // Connect root to all orphan nodes (root is parent of orphans)
            orphanNodes.forEach(orphan => {
                links.push({
                    source: orphan.id,
                    target: 'VIRTUAL_ROOT'
                });
            });

            console.log(`Added virtual root node connected to ${orphanNodes.length} orphan nodes`);
        }
    }

    calculateNodeSize(enrichment) {
        if (!enrichment || !enrichment.p_value) return 12.5;

        // Size based on -log10(p_value) - reduced by half
        const logP = -Math.log10(enrichment.p_value);
        return Math.min(Math.max(logP * 4, 10), 30);
    }

    createLinks(links) {
        const g = this.svg.select('.zoom-group');

        g.append('g')
            .attr('class', 'links')
            .selectAll('line')
            .data(links)
            .enter().append('line')
            .attr('class', 'link')
            .attr('marker-end', 'url(#child-to-parent-arrow)');
    }

    createForceLayout(nodes, links) {
        console.log(`Creating force layout with ${nodes.length} nodes and ${links.length} links`);
        console.log('Sample node:', nodes[0]);

        // Create the force simulation
        this.simulation = d3.forceSimulation(nodes)
            .force('link', d3.forceLink(links).id(d => d.id).distance(100))
            .force('charge', d3.forceManyBody().strength(-300))
            .force('center', d3.forceCenter(this.width / 2, this.height / 2))
            .force('collision', d3.forceCollide().radius(d => d.radius + 5));

        // Create links
        const link = this.svg.select('.zoom-group')
            .selectAll('.link')
            .data(links)
            .enter()
            .append('line')
            .attr('class', 'link')
            .attr('marker-end', 'url(#child-to-parent-arrow)');

        console.log(`Created ${link.size()} link elements`);

        // Create nodes
        const nodeGroup = this.svg.select('.zoom-group')
            .selectAll('.node-group')
            .data(nodes)
            .enter()
            .append('g')
            .attr('class', 'node-group')
            .call(this.createDragBehavior());

        console.log(`Created ${nodeGroup.size()} node group elements`);

        // Add circles
        const circles = nodeGroup.append('circle')
            .attr('r', d => d.radius)
            .attr('fill', d => this.getNodeColor(d))
            .attr('stroke', d => this.getNodeBorderColor(d))
            .attr('class', d => this.getNodeClass(d))
            .on('mouseover', (event, d) => this.showTooltip(event, d))
            .on('mousemove', (event, d) => this.updateTooltipPosition(event))
            .on('mouseout', () => this.hideTooltip());

        console.log(`Created ${circles.size()} circle elements`);

        // Add labels
        const labels = nodeGroup.append('text')
            .attr('class', 'node-label')
            .attr('dy', '.35em')
            .text(d => this.truncateLabel(d.name, d.radius))
            .style('font-size', d => `${Math.max(8, d.radius / 3)}px`);

        console.log(`Created ${labels.size()} label elements`);

        // Set up the simulation tick handler
        this.simulation.on('tick', () => {
            link
                .attr('x1', d => this.getLinkEndpoint(d.source, d.target, d.source.radius).x)
                .attr('y1', d => this.getLinkEndpoint(d.source, d.target, d.source.radius).y)
                .attr('x2', d => this.getLinkEndpoint(d.target, d.source, d.target.radius + 4).x)
                .attr('y2', d => this.getLinkEndpoint(d.target, d.source, d.target.radius + 4).y);

            nodeGroup
                .attr('transform', d => `translate(${d.x},${d.y})`);
        });

        // Auto-freeze simulation and adjust viewport after it settles
        setTimeout(() => {
            if (this.simulation) {
                console.log('Stopping simulation after timeout, sample node positions:');
                nodes.slice(0, 3).forEach(n => console.log(`Node ${n.id}: x=${n.x}, y=${n.y}`));
                this.simulation.stop();
                this.adjustViewportToFitNodes(nodes);
                console.log('Force simulation stopped and viewport adjusted');
            }
        }, 3000); // Stop after 3 seconds

        // Also stop on low alpha (when simulation has mostly settled)
        this.simulation.on('tick', () => {
            link
                .attr('x1', d => this.getLinkEndpoint(d.source, d.target, d.source.radius).x)
                .attr('y1', d => this.getLinkEndpoint(d.source, d.target, d.source.radius).y)
                .attr('x2', d => this.getLinkEndpoint(d.target, d.source, d.target.radius + 4).x)
                .attr('y2', d => this.getLinkEndpoint(d.target, d.source, d.target.radius + 4).y);

            nodeGroup
                .attr('transform', d => `translate(${d.x},${d.y})`);

            // Stop simulation when it has settled
            if (this.simulation.alpha() < 0.01) {
                console.log('Simulation alpha dropped, stopping. Sample node positions:');
                nodes.slice(0, 3).forEach(n => console.log(`Node ${n.id}: x=${n.x}, y=${n.y}`));
                this.simulation.stop();
                this.adjustViewportToFitNodes(nodes);
                console.log('Force simulation auto-stopped at low alpha and viewport adjusted');
            }
        });
    }


    createHierarchicalLayout(nodes, links) {
        // Find disconnected components (separate trees)
        const components = this.findConnectedComponents(nodes, links);

        console.log(`Found ${components.length} disconnected components`);

        let currentX = 50;
        const componentGap = 200;
        const allPositionedNodes = [];
        const allPositionedLinks = [];

        // Layout each component separately
        components.forEach((component, index) => {
            const componentNodes = component.nodes;
            const componentLinks = component.links;

            // Find root nodes for this component (nodes with no parents in component)
            const childIds = new Set(componentLinks.map(l => l.source));
            const roots = componentNodes.filter(n => !childIds.has(n.id));

            if (roots.length === 0 && componentNodes.length > 0) {
                // No clear hierarchy, treat first node as root
                roots.push(componentNodes[0]);
            }

            // Create tree layout for this component
            const treeLayout = d3.tree()
                .size([300, this.height - 100])
                .separation((a, b) => a.parent === b.parent ? 1 : 2);

            if (roots.length > 0) {
                // Build hierarchy starting from roots
                const hierarchyRoot = this.buildComponentHierarchy(roots, componentNodes, componentLinks);
                const treeData = treeLayout(hierarchyRoot);

                // Position nodes for this component
                treeData.descendants().forEach(d => {
                    if (d.data && d.data.id !== 'component_root') {
                        const node = nodes.find(n => n.id === d.data.id);
                        if (node) {
                            node.x = d.x + currentX;
                            node.y = d.y + 50;
                            allPositionedNodes.push(node);
                        }
                    }
                });

                // Collect links for this component
                treeData.links().forEach(link => {
                    if (link.source.data.id !== 'component_root' && link.target.data.id !== 'component_root') {
                        allPositionedLinks.push({
                            source: nodes.find(n => n.id === link.source.data.id),
                            target: nodes.find(n => n.id === link.target.data.id)
                        });
                    }
                });
            } else {
                // Single node component
                componentNodes.forEach((node, i) => {
                    node.x = currentX;
                    node.y = 100 + i * 50;
                    allPositionedNodes.push(node);
                });
            }

            currentX += 250; // Space between components
        });

        // Create visualization
        this.createSimpleLinks(allPositionedLinks);
        this.createSimpleNodes(allPositionedNodes);

        // Adjust viewport to fit all nodes
        this.adjustViewportToFitNodes(allPositionedNodes);
    }

    findConnectedComponents(nodes, links) {
        const visited = new Set();
        const components = [];
        const nodeMap = new Map(nodes.map(n => [n.id, n]));

        // Build adjacency list
        const adjacency = new Map();
        nodes.forEach(n => adjacency.set(n.id, new Set()));

        links.forEach(link => {
            adjacency.get(link.source)?.add(link.target);
            adjacency.get(link.target)?.add(link.source);
        });

        // Find components using DFS
        for (const node of nodes) {
            if (!visited.has(node.id)) {
                const component = { nodes: [], links: [] };
                const stack = [node.id];
                const componentNodes = new Set();

                while (stack.length > 0) {
                    const current = stack.pop();
                    if (visited.has(current)) continue;

                    visited.add(current);
                    componentNodes.add(current);
                    component.nodes.push(nodeMap.get(current));

                    // Add neighbors to stack
                    for (const neighbor of adjacency.get(current) || []) {
                        if (!visited.has(neighbor)) {
                            stack.push(neighbor);
                        }
                    }
                }

                // Add links within this component
                component.links = links.filter(link =>
                    componentNodes.has(link.source) && componentNodes.has(link.target)
                );

                components.push(component);
            }
        }

        return components;
    }

    buildComponentHierarchy(roots, nodes, links) {
        if (roots.length === 1) {
            return d3.hierarchy(this.buildNodeTree(roots[0], nodes, links));
        } else {
            // Multiple roots - create virtual root
            const virtualRoot = {
                id: 'component_root',
                name: 'Root',
                children: roots.map(root => this.buildNodeTree(root, nodes, links))
            };
            return d3.hierarchy(virtualRoot);
        }
    }

    buildTreeStructure(nodes, links) {
        // Find root nodes (nodes with no parents in the current network)
        const children = new Set(links.map(l => l.source));
        const parents = new Set(links.map(l => l.target));
        const roots = nodes.filter(n => !children.has(n.id));

        // If multiple roots, create a virtual root
        if (roots.length > 1) {
            const virtualRoot = {
                id: 'virtual_root',
                name: 'Biological Processes',
                children: roots.map(r => this.buildNodeTree(r, nodes, links))
            };
            return d3.hierarchy(virtualRoot);
        } else if (roots.length === 1) {
            return d3.hierarchy(this.buildNodeTree(roots[0], nodes, links));
        } else {
            // Fallback: pick the most general node (highest in hierarchy)
            const rootNode = nodes[0];
            return d3.hierarchy(this.buildNodeTree(rootNode, nodes, links));
        }
    }

    buildNodeTree(node, allNodes, links) {
        const children = links
            .filter(l => l.target === node.id)
            .map(l => allNodes.find(n => n.id === l.source))
            .filter(Boolean)
            .map(child => this.buildNodeTree(child, allNodes, links));

        return {
            ...node,
            children: children.length > 0 ? children : null
        };
    }

    createTreeLinks(treeLinks) {
        const g = this.svg.select('.zoom-group');

        g.append('g')
            .attr('class', 'links')
            .selectAll('path')
            .data(treeLinks)
            .enter().append('path')
            .attr('class', 'link')
            .attr('marker-end', 'url(#child-to-parent-arrow)')
            .attr('d', d3.linkVertical()
                .x(d => d.x + 50)
                .y(d => d.y + 50)
            );
    }

    createSimpleLinks(links) {
        const g = this.svg.select('.zoom-group');

        g.append('g')
            .attr('class', 'links')
            .selectAll('line')
            .data(links)
            .enter().append('line')
            .attr('class', 'link')
            .attr('marker-end', 'url(#child-to-parent-arrow)')
            .attr('x1', d => this.getLinkEndpoint(d.source, d.target, d.source.radius).x)
            .attr('y1', d => this.getLinkEndpoint(d.source, d.target, d.source.radius).y)
            .attr('x2', d => this.getLinkEndpoint(d.target, d.source, d.target.radius + 4).x)
            .attr('y2', d => this.getLinkEndpoint(d.target, d.source, d.target.radius + 4).y);
    }

    createSimpleNodes(nodes) {
        const g = this.svg.select('.zoom-group');

        // Create node groups
        const nodeGroups = g.append('g')
            .attr('class', 'nodes')
            .selectAll('g')
            .data(nodes)
            .enter().append('g')
            .attr('class', 'node-group')
            .attr('transform', d => `translate(${d.x}, ${d.y})`);

        // Add circles
        nodeGroups.append('circle')
            .attr('class', d => this.getNodeClass(d))
            .attr('r', d => d.radius)
            .attr('fill', d => this.getNodeColor(d))
            .attr('stroke', d => this.getNodeBorderColor(d))
            .style('cursor', 'pointer')
            .on('mouseover', (event, d) => this.showTooltip(event, d))
            .on('mouseout', () => this.hideTooltip());

        // Add labels for important nodes
        nodeGroups.append('text')
            .attr('class', 'node-label')
            .attr('dy', '.35em')
            .style('display', 'block')  // Show all labels since nodes are bigger
            .style('font-size', '12px')
            .text(d => this.truncateLabel(d.name, d.radius));
    }

    createTreeNodes(descendants) {
        const g = this.svg.select('.zoom-group');

        // Create node groups
        const nodeGroups = g.append('g')
            .attr('class', 'nodes')
            .selectAll('g')
            .data(descendants)
            .enter().append('g')
            .attr('class', 'node-group')
            .attr('transform', d => `translate(${d.x + 50}, ${d.y + 50})`);

        // Add circles
        nodeGroups.append('circle')
            .attr('class', d => this.getNodeClass(d.data))
            .attr('r', d => d.data.radius)
            .attr('fill', d => this.getNodeColor(d.data))
            .attr('stroke', d => this.getNodeBorderColor(d.data))
            .on('mouseover', (event, d) => this.showTooltip(event, d.data))
            .on('mouseout', () => this.hideTooltip());

        // Add labels for important nodes
        nodeGroups.append('text')
            .attr('class', 'node-label')
            .attr('dy', '.35em')
            .style('display', d => d.data.is_original || d.data.radius > 12 ? 'block' : 'none')
            .text(d => this.truncateLabel(d.data.name, d.data.radius));
    }

    getNodeColor(d) {
        if (d.is_root) return '#e0e0e0'; // Light gray for virtual root
        if (d.is_original) return '#ffeb3b'; // Yellow for original enrichment terms
        return '#f8f8f8'; // Light gray for hierarchy nodes
    }

    getNodeBorderColor(d) {
        if (d.is_root) return '#666666'; // Dark gray for virtual root
        if (d.is_enrichment_leaf) return '#ff6b6b'; // Red for enrichment leaves
        return '#4ecdc4'; // Teal for internal nodes
    }

    getNodeClass(d) {
        const classes = ['node'];
        classes.push(d.is_enrichment_leaf ? 'leaf' : 'internal');
        if (d.is_original) classes.push('original');
        if (d.is_root) classes.push('root');
        return classes.join(' ');
    }

    getLinkEndpoint(node, otherNode, offset) {
        const dx = otherNode.x - node.x;
        const dy = otherNode.y - node.y;
        const length = Math.sqrt(dx * dx + dy * dy);

        if (!length || !Number.isFinite(length)) {
            return { x: node.x, y: node.y };
        }

        return {
            x: node.x + (dx / length) * offset,
            y: node.y + (dy / length) * offset
        };
    }

    truncateLabel(text, radius) {
        const maxLength = Math.floor(radius / 2);
        return text.length > maxLength ? text.substring(0, maxLength) + '...' : text;
    }

    createDragBehavior() {
        return d3.drag()
            .on('start', (event, d) => {
                if (!event.active) this.simulation.alphaTarget(0.3).restart();
                d.fx = d.x;
                d.fy = d.y;
            })
            .on('drag', (event, d) => {
                d.fx = event.x;
                d.fy = event.y;
            })
            .on('end', (event, d) => {
                if (!event.active) this.simulation.alphaTarget(0);
                d.fx = null;
                d.fy = null;
            });
    }

    updateVisualization() {
        this.svg.selectAll('.link')
            .attr('x1', d => this.getLinkEndpoint(d.source, d.target, d.source.radius).x)
            .attr('y1', d => this.getLinkEndpoint(d.source, d.target, d.source.radius).y)
            .attr('x2', d => this.getLinkEndpoint(d.target, d.source, d.target.radius + 4).x)
            .attr('y2', d => this.getLinkEndpoint(d.target, d.source, d.target.radius + 4).y);

        this.svg.selectAll('.node-group')
            .attr('transform', d => `translate(${d.x},${d.y})`);
    }

    showTooltip(event, d) {
        let tooltip = document.getElementById('tooltip');
        if (!tooltip) {
            // Create tooltip element if it doesn't exist
            tooltip = document.createElement('div');
            tooltip.id = 'tooltip';
            tooltip.className = 'tooltip';
            document.body.appendChild(tooltip);
        }

        // Build tooltip content
        let content = `<strong>${d.name || d.id}</strong>`;

        if (d.is_original && d.p_value) {
            content += `<br><strong>P-value:</strong> ${d.p_value.toExponential(2)}`;
            content += `<br><strong>Rank:</strong> ${d.rank}`;
        }

        if (d.information_content) {
            content += `<br><strong>Info Content:</strong> ${d.information_content.toFixed(1)}`;
        }

        const nodeType = d.is_enrichment_leaf ? 'Most Specific Finding' :
                        d.is_original ? 'Parent Enrichment Term' : 'Ancestor Term';
        content += `<br><strong>Type:</strong> ${nodeType}`;

        tooltip.innerHTML = content;

        // Position tooltip near mouse but avoid going off-screen
        const x = Math.min(event.pageX + 15, window.innerWidth - 250);
        const y = Math.max(event.pageY - 10, 10);

        // Apply styles with !important to override CSS
        tooltip.style.setProperty('position', 'absolute', 'important');
        tooltip.style.setProperty('left', x + 'px', 'important');
        tooltip.style.setProperty('top', y + 'px', 'important');
        tooltip.style.setProperty('display', 'block', 'important');
        tooltip.style.setProperty('opacity', '1', 'important');
        tooltip.style.setProperty('visibility', 'visible', 'important');
        tooltip.style.setProperty('z-index', '10000', 'important');
        tooltip.style.setProperty('background-color', 'rgba(0,0,0,0.9)', 'important');
        tooltip.style.setProperty('color', 'white', 'important');
        tooltip.style.setProperty('padding', '12px', 'important');
        tooltip.style.setProperty('border-radius', '8px', 'important');
        tooltip.style.setProperty('font-size', '12px', 'important');
        tooltip.style.setProperty('max-width', '300px', 'important');
        tooltip.style.setProperty('box-shadow', '0 4px 12px rgba(0,0,0,0.3)', 'important');
        tooltip.style.setProperty('pointer-events', 'none', 'important');
        tooltip.style.setProperty('line-height', '1.4', 'important');

        // Remove any test styling
        tooltip.style.removeProperty('width');
        tooltip.style.removeProperty('height');
        tooltip.style.removeProperty('border');
    }

    hideTooltip() {
        document.getElementById('tooltip').style.display = 'none';
    }

    updateTooltipPosition(event) {
        const tooltip = document.getElementById('tooltip');
        if (tooltip && tooltip.style.display !== 'none') {
            // Position tooltip near mouse but avoid going off-screen
            const x = Math.min(event.pageX + 15, window.innerWidth - 250);
            const y = Math.max(event.pageY - 10, 10);

            tooltip.style.setProperty('left', x + 'px', 'important');
            tooltip.style.setProperty('top', y + 'px', 'important');
        }
    }

    adjustViewportToFitNodes(nodes) {
        if (nodes.length === 0) return;

        // Calculate bounding box of all nodes (including their radius)
        const padding = 100;
        const minX = Math.min(...nodes.map(n => n.x - n.radius)) - padding;
        const maxX = Math.max(...nodes.map(n => n.x + n.radius)) + padding;
        const minY = Math.min(...nodes.map(n => n.y - n.radius)) - padding;
        const maxY = Math.max(...nodes.map(n => n.y + n.radius)) + padding;

        const width = maxX - minX;
        const height = maxY - minY;

        console.log(`Adjusting viewport: width=${width}, height=${height}, bounds=[${minX}, ${minY}, ${maxX}, ${maxY}]`);

        // Update SVG dimensions and viewBox
        this.svg
            .attr('width', '100%')
            .attr('height', Math.max(600, height))
            .attr('viewBox', `${minX} ${minY} ${width} ${height}`);
    }

    showLoading(show) {
        document.getElementById('loading').style.display = show ? 'block' : 'none';
        document.getElementById('networkViz').style.opacity = show ? '0.5' : '1';
    }

    showError(message) {
        const container = document.getElementById('errorMessage');
        container.textContent = message;
        container.style.display = 'block';
    }

    hideError() {
        document.getElementById('errorMessage').style.display = 'none';
    }
}

// Initialize the application
const app = new DiseaseEnrichmentViz();
