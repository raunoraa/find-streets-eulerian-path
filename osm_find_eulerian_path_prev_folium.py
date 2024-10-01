import osmnx as ox
import networkx as nx
import folium
from jinja2 import Template

def find_graph_boundaries(G):
    """Determine the boundaries of the graph based on node coordinates."""
    lats = [G.nodes[node]['y'] for node in G.nodes]
    lons = [G.nodes[node]['x'] for node in G.nodes]
    
    return {
        'north': max(lats),
        'south': min(lats),
        'east': max(lons),
        'west': min(lons)
    }

def is_boundary_node(node, G, boundary, margin=0.001):
    """Check if a node is close to the map boundary."""
    lat, lon = G.nodes[node]['y'], G.nodes[node]['x']
    return (lat <= boundary['south'] + margin or lat >= boundary['north'] - margin or
            lon <= boundary['west'] + margin or lon >= boundary['east'] - margin)

def handle_boundary_edges(G, boundary):
    """Adjust edges for nodes near the boundary by adding reverse edges if necessary."""
    boundary_nodes = {node for node in G.nodes if is_boundary_node(node, G, boundary)}
    H = G.copy()

    for u, v, data in G.edges(data=True):
        if u in boundary_nodes or v in boundary_nodes:
            if not H.has_edge(v, u):                
                H.add_edge(v, u, **data)
    
    return H

def filter_drivable_lanes(G):
    """Filter the graph to include only drivable lanes."""
    H = nx.MultiDiGraph()
    drivable_highways = {'motorway', 'trunk', 'primary', 'secondary', 'tertiary', 'residential'}
    
    nodes_in_filtered_edges= set()

    for u, v, key, data in G.edges(keys=True, data=True):
        if 'highway' in data and data['highway'] in drivable_highways:
            nodes_in_filtered_edges.add(u)
            nodes_in_filtered_edges.add(v)
            H.add_edge(u, v, key=key, **data)            
    
    for node in nodes_in_filtered_edges:
       H.nodes[node].update(G.nodes[node])

    default_weight = 1.0
    for u, v, key, data in H.edges(keys=True, data=True):
        if 'weight' not in data:
            H[u][v][key]['weight'] = default_weight    

    return H

def visualize_eulerian_path(G, eulerian_path, boundary):
    """Visualize the Eulerian path on the graph using Folium with JavaScript for animation control."""
    
    # Set up the initial map using Folium
    map_center = [(boundary['north'] + boundary['south']) / 2, (boundary['east'] + boundary['west']) / 2]
    m = folium.Map(location=map_center, zoom_start=15)

    polyline_ids = []  # Track polyline IDs for animation
    
    # Add each edge of the Eulerian path as a hidden polyline
    order = 0
    for u, v in eulerian_path:
        order += 1

        # Get node coordinates
        x_start, y_start = G.nodes[u]['x'], G.nodes[u]['y']
        x_end, y_end = G.nodes[v]['x'], G.nodes[v]['y']

        # Get edge data
        edge_data = G[u][v][0]
        distance = edge_data.get("length", 0)
        weight = edge_data.get("weight", 1)
        street = edge_data.get("name", "Unknown Street")

        # Add the polyline for the path segment (hidden initially)
        polyline = folium.PolyLine(
            locations=[(y_start, x_start), (y_end, x_end)],
            color='blue',
            weight=3,
            tooltip=f"Order: {order}<br>Street: {street}<br>Length: {distance*weight:.2f} m",
            opacity=0  # Initially hidden
        )
        polyline.add_to(m)
        polyline_ids.append(polyline._id)

    # Add custom JavaScript for animation control
    animation_script = Template("""
    <script>
    window.onload = function() {        

        let animation;
        let currentStep = 0;
        let prevStep = 0;
        let playing = false;
                                
        let infoText = document.getElementById('info-text');
                                
        let polylines = document.querySelectorAll(".leaflet-interactive");
        let prev_polyline = null;

        let delay = 500;  // 500 ms delay between segments
        const delaySlider = document.getElementById('delay-slider');
        const delayValueDisplay = document.getElementById('delay-value');
                                
        const stepSlider = document.getElementById('step-slider');
        stepSlider.max = polylines.length;
        const stepValueDisplay = document.getElementById('step-value');

        // Update the delay when slider is changed
        delaySlider.oninput = function() {
            delay = this.value;
            delayValueDisplay.textContent = `${delay} ms`;  // Update the display
        };               

        stepSlider.oninput = function() {
            prevStep = currentStep;
            currentStep = this.value;
            stepValueDisplay.textContent = `${currentStep}`;  // Update the display
            animateStep();
        }                  

        function animatePath() {
            if (!playing) return;  // Only animate when playing is true
            if (currentStep >= polylines.length) {
                infoText.textContent = "Finished!";
                return;
            }  // Stop when the end is reached                                                        

            if (prev_polyline) {
                prev_polyline.style.stroke = "blue";
            }

            let polyline = polylines[currentStep];
            if (polyline) {
                polyline.style["stroke-opacity"] = 1;  // Show the polyline
                polyline.style.stroke = "black";
                prev_polyline = polyline;
            }
            currentStep++;
            
            stepSlider.value = currentStep;
            stepValueDisplay.textContent = `${currentStep}`;  // Update the display                                
            
            animation = setTimeout(animatePath, delay);
        }
                                
        function animateStep(){
            if(playing) return; // Only let it animate the step, when it is not auto-playing.
            //if(currentStep >= polylines.length) {
                //infoText.textContent = "Finished!";
                //return;
            //}  // Stop here when the end is reached   
            console.log(prevStep + " " + currentStep)  ;
                                
            if (prev_polyline) {
                
                if (prevStep > currentStep) {
                    for (let i = prevStep; i > prevStep - currentStep; i--) {
                        polylines[i].style["stroke-opacity"] = 0;
                    }
                    //prev_polyline.style["stroke-opacity"] = 0;
                } else 
                {                                
                    prev_polyline.style.stroke = "blue";
                }
            }

            if (currentStep - prevStep >= 1){
                prev_polyline.style.stroke = "blue";
                for (let i = prevStep; i < currentStep; i++) {                        
                        polylines[i].style["stroke-opacity"] = 1;
                        polylines[i].style.stroke = "black";
                        prev_polyline = polylines[i];
                }
            }

            else{
                let polyline = polylines[currentStep];
                if (polyline) {
                    polyline.style["stroke-opacity"] = 1;  // Show the polyline
                    polyline.style.stroke = "black";
                    prev_polyline = polyline;
                }    
            }                
        }

        function playAnimation() {
            if (currentStep < polylines.length) {
                playing = true;
                infoText.textContent = "Playing..."            
                animatePath();
            }
        }

        function pauseAnimation() {
            playing = false;
            infoText.textContent = "Paused"
            clearTimeout(animation);
        }

        function resetAnimation() {
            playing = false;
            currentStep = 0;
            clearTimeout(animation);            
            resetAll();
        }

        function resetAll() {
            infoText.textContent = "Ready!";
            prev_polyline = null;
            polylines = document.querySelectorAll(".leaflet-interactive");
            polylines.forEach(el => {
                el.style["stroke-opacity"] = 0;
            });
        }

        document.getElementById('play-button').onclick = function() { playAnimation(); };
        document.getElementById('pause-button').onclick = function() { pauseAnimation(); };
        document.getElementById('restart-button').onclick = function() { resetAnimation(); playAnimation(); };
    };
    </script>
    """).render(polyline_ids=polyline_ids)

    # Attach the JavaScript to the map HTML
    m.get_root().html.add_child(folium.Element(animation_script))

    # Add play, pause, and restart buttons to the map
    button_html = """
        <div style="position: fixed; bottom: 60px; left: 50px; z-index: 9999;">

            <p id="info-text">Ready!</p>

            <button id="play-button" style="padding: 10px; margin-right: 10px;">Play</button>
            <button id="pause-button" style="padding: 10px; margin-right: 10px;">Pause</button>
            <button id="restart-button" style="padding: 10px;">Restart</button>

            <!-- Add slider to control delay -->
            <div style="margin-top: 20px;">
                <label for="delay-slider">Animation Delay (ms):</label>
                <input type="range" id="delay-slider" min="100" max="2000" value="500" step="100" style="margin-left: 10px;">
                <span id="delay-value">500 ms</span>
            </div>

            <!-- Add slider to control steps -->
            <div style="margin-top: 20px;">
                <label for="step-slider">Steps:</label>
                <input type="range" id="step-slider" min="0" max="1" value="0" step="1" style="margin-left: 10px;">
                <span id="step-value">0</span>
            </div>
        </div>        
    """
    m.get_root().html.add_child(folium.Element(button_html))
    
    # Return the Folium map
    return m

def main(xml_file_path, verbose=True):
    # Load the graph from the OSM XML file
    G_generated = ox.graph_from_xml(filepath=xml_file_path)

    G = filter_drivable_lanes(G_generated)
    
    # Print graph details
    print(f"Number of nodes: {len(G.nodes)}")
    print(f"Number of edges: {len(G.edges)}")
    
    # Find the boundaries of the graph
    boundary = find_graph_boundaries(G)

    if verbose:
        print("Graph boundaries:", boundary)
    
    # Handle boundary edges
    G = handle_boundary_edges(G, boundary)

    if verbose:
        print("Is G strongly connected:", nx.is_strongly_connected(G))

    eulerian_path = nx.eulerian_path(G)

    total_distance = 0

    if verbose:
        print("Eulerian path sequence:")

    eulerian_path_list = list(eulerian_path)

    for u, v in eulerian_path_list:
        edge_data = G[u][v]
        
        distance = edge_data.get(0).get("length")
        weight = edge_data.get(0).get("weight")
        street = edge_data.get(0).get("name")

        if verbose:
            print(f"Edge from {u} to {v}, length: {distance*weight} meters; street: {street}")
    
        total_distance += (distance * weight)    

    print("Total Distance:", total_distance, "meters")

    # Visualize Eulerian Path
    m = visualize_eulerian_path(G, eulerian_path_list, boundary)
    
    # Save to HTML
    m.save("eulerian_path_animation.html")

# Example call
main('map_files/osm_observable.xml', verbose=True)
