from jinja2 import Template
import networkx as nx
import geojson
from Intersection import Intersection
from shapely.geometry import shape
from shapely import Polygon
from collections import defaultdict
from itertools import combinations
from copy import deepcopy

# Graph visualization
import folium
import random

import time


# Load GeoJSON files
def load_geojson(file_path):
    with open(file_path, 'r') as f:
        return geojson.load(f)

# Function to check if lane is car-drivable
def is_car_drivable(lane):
    # Check for 'backward' access permission if present
    backward_access = lane['properties']['muv']['travel'].get('backward', {}).get('access', {}).get('LandBased', [{}])[0].get('value', 'No')
    
    # Check for 'forward' access permission if present
    forward_access = lane['properties']['muv']['travel'].get('forward', {}).get('access', {}).get('LandBased', [{}])[0].get('value', 'No')
    
    # Return True if either access is 'Yes'
    return backward_access == 'Yes' or forward_access == 'Yes'

# Load and filter car-drivable lanes
# Returns a defaultdict, where keys are road_ids and values are the lanes (and their respective data) 
#  each road has
def parse_lanes(lane_geojson):
    lanes_by_road_id = defaultdict(list)  # A dictionary where road_id is the key, and values are lists of lanes
    for lane in lane_geojson['features']:
        if is_car_drivable(lane):
            lane_data = {
                'road_id': lane['properties']['road'], 
                'lane_id': lane['properties']['index'],
                'direction': lane['properties']['direction'], 
                'src_i': None, # value will be added in parse_intersections
                'dst_i': None, # value will be added in parse_intersections
                'geometry': shape(lane['geometry'])
            }
            lanes_by_road_id[lane_data['road_id']].append(lane_data)  # Group lanes by road_id
    return lanes_by_road_id

# Load intersections and their movements
# Also assign src and destination intersections to lanes (lanes parameter will be modified)
def parse_intersections(intersection_geojson, all_lanes):    
    intersections = set()
    
    for feature in intersection_geojson['features']:
        i_type = feature['properties']['type']
        
        if i_type == 'road':
            # If the intersection type is road, then we can assign the src and dst intersections to lanes            
                        
            lanes = all_lanes[feature['properties']['id']]
                        
            for lane in lanes:
                lane['src_i'] = feature['properties']['src_i']
                lane['dst_i'] = feature['properties']['dst_i']
            
        elif i_type == 'intersection':
            if feature['properties']['intersection_kind'] != 'MapEdge':
                intersection = Intersection(feature['properties']['id'], feature['properties']['movements'], shape(feature['geometry']))
                intersections.add(intersection)
        else:
            # Let's not look at other types at the moment
            continue

    return intersections

# Create such line, which goes through the middle of the lane polygon
def centerline_from_lane_polygon(lane_polygon):
    # Shrink the polygon slightly inward to create a line (this follows the shape)
    centerline = lane_polygon.buffer(-0.000001)  # Shrink the polygon by a small amount
    
    # If the result is a polygon, extract the exterior as a line
    #if isinstance(centerline, Polygon):
        #centerline = centerline.exterior  # Take the outer boundary as a line
    
    return centerline  # This is now a LineString (the centerline)

# Create the directed multigraph
def create_graph(lanes, intersections):

    G = nx.MultiDiGraph()

    # Added nodes, which are grouped by their intersection
    # key=intersection id; value=set of nodes in an intersection
    # a node is a 3-tuple
    intersections_node_map = {}
    intersection_ids = [el.id for el in intersections]

    # We want to add intersections in a way
    # that there will be created one node per lane
    for intersection in intersections:
        id = intersection.id
        movements = intersection.movements
        unique_roads = set() # a set of unique road ids in an intersection
        for movement in movements:
            src_road, _ = movement.split(" -> ")

            # There will always be 1-1 ratio of roads in intersections
            # so we can take unique ones only from source or destination
            src_id = int(src_road.split('#')[-1])            
            unique_roads.add(src_id)        
        
        
        counter = 0

        tuples = set()
        for ur in unique_roads:
            lanes_in_ur = lanes.get(ur)            
            for lane in lanes_in_ur:
                
                # Don't add such lanes to intersection, which dont have either src or dst intersection in the intersection list
                if lane['src_i'] not in intersection_ids or lane['dst_i'] not in intersection_ids:                    
                    continue

                # Just for now, assign the same geometry to each intersection node
                # have to change later (it's just visual distortion, wont affect the graph)

                node_id = (id, lane['road_id'], lane['lane_id'])                

                # Also add the info, if the node for the current lane is for entering or leaving the intersection
                is_entering_value = False                
                if (lane['dst_i'] == id and lane['direction'] == 'Forward') or (lane['src_i'] == id and lane['direction'] == 'Backward'):
                    is_entering_value = True        
                
                # debug
                if id == 1 or id == 22:
                    print(node_id, is_entering_value)        


                tuples.add((node_id, is_entering_value, intersection.geometry, id))
                G.add_node(node_id, is_entering=is_entering_value, geometry=intersection.geometry, intersection_id=id)
                counter += 1
        
        intersections_node_map[id] = tuples

        # debug
        print(f"Nodes in intersection {id}: {counter}")
    
    # Create edges inside the intersection
    for _, nodes in intersections_node_map.items():        
        entering_nodes = []
        leaving_nodes = []
        for node in nodes:
            _, is_entering, _, _ = node
            if is_entering:
                entering_nodes.append(node)
            else:
                leaving_nodes.append(node)
            
        if len(entering_nodes) > 0 and len(leaving_nodes) > 0:
            for e_node in entering_nodes:
                e_node_id, _, e_geometry, _ = e_node                
                for l_node in leaving_nodes:
                    l_node_id, _, l_geometry, _ = l_node

                    #debug
                    if l_node_id[0] == 1 or e_node_id[0] == 1:
                        print(l_node_id, e_node_id)

                    # Need to tweak the geometry, but can be done later (TODO)
                    # after the graph seems ok.
                    # This is just placeholder geometry for now.
                    G.add_edge(e_node_id, l_node_id, geometry=e_geometry, centerline=centerline_from_lane_polygon(e_geometry))
    

    # Create edges outside the intersections (basically add the lanes)
    # We need to create an edge between all such nodes, which have the same road and lane id
    observables = {}
    for node in G.nodes(data=True):
        #print(node)
        id_tuple, node_data = node        
        #print(id_tuple)
        _, road_id, lane_id = id_tuple
        observable_tuple = (road_id, lane_id)

        # We can assume that there are only 2 nodes per each road
        if observable_tuple in observables.keys():
            popped_node = observables.pop(observable_tuple)
            node_one_id_tuple = popped_node[0]

            lane_geometry = None
            # Get the lane geometry from the lanes defaultdict
            for l in lanes[road_id]:
                if l.get('lane_id') == lane_id:
                    lane_geometry = l.get('geometry')
                    break  

            centerline_geometry = centerline_from_lane_polygon(lane_geometry)
            
            if not node_data['is_entering']:
                G.add_edge(id_tuple, node_one_id_tuple, geometry = lane_geometry, centerline=centerline_geometry)
            else:
                G.add_edge(node_one_id_tuple, id_tuple, geometry = lane_geometry, centerline=centerline_geometry)
        else:
            observables[observable_tuple] = node    

    return G

# Main function to construct the graph from geojson files
def build_city_graph(lane_geojson_file, intersection_geojson_file):
    # Load and parse the geojson files
    lanes_geojson = load_geojson(lane_geojson_file)
    intersections_geojson = load_geojson(intersection_geojson_file)

    lanes = parse_lanes(lanes_geojson)
    intersections = parse_intersections(intersections_geojson, lanes)

    # Create the directed multigraph
    G = create_graph(lanes, intersections)
    return G


def visualize_graph(G, map_boundaries, file_name="city_graph_map", visualize_nodes=True):

    # Create a Folium map centered on the middle of the boundaries
    map_center = [(map_boundaries[1] + map_boundaries[3]) / 2, 
                (map_boundaries[0] + map_boundaries[2]) / 2]

    # Initialize the Folium map
    m = folium.Map(location=map_center, zoom_start=13)

    # Add nodes to the map (assuming nodes are polygons)
    for node, data in G.nodes(data=True):
        geometry = data['geometry']
        if geometry:
            # Get the coordinates of the polygon
            coordinates = list(geometry.exterior.coords)  # Taking the first polygon's coordinates
            '''
            folium.Polygon(
                locations=[(coord[1], coord[0]) for coord in coordinates],
                color='blue',
                fill=True,
                fill_color='blue',
                fill_opacity=0.5,
                popup=f'Node {node}'
            ).add_to(m)
            '''
            if visualize_nodes:
                folium.Marker(
                    location=(coordinates[0][1] + random.uniform(-0.00003, 0.00003), coordinates[0][0] + random.uniform(-0.00003, 0.00003)),
                    popup=f"Out degree: {G.out_degree(node)}; Int.id: {data['intersection_id']}",  # Optional popup text
                    icon=folium.Icon(color='blue')  # Optional: Customize the marker color
                ).add_to(m)

    # Add edges to the map (assuming edges are also polygons)
    for u, v, data in G.edges(data=True):
        geometry = data['geometry']
        if geometry:
            coordinates = list(geometry.exterior.coords)  # Taking the first polygon's coordinates
            folium.Polygon(
                locations=[(coord[1], coord[0]) for coord in coordinates],
                color='black',
                #weight=2,
                #opacity=0.7,
                fill=True,
                fill_color='red',
                fill_opacity=0.5,
                popup=f'Edge from {u} to {v}'
            ).add_to(m)
    
    # Save the map as an html file
    m.save(file_name + '.html')    

    return m

def visualize_path(folium_map_object, G, path):   
    print() 

    for order, node in enumerate(path):
        
        if order == len(path) - 1:
            break
        
        '''
        # function to get the arithmetic mean value of the list
        get_avg = lambda numbers: sum(numbers) / len(numbers) if numbers else 0

        
        node_geometry = G.nodes[node]['geometry']

        start_coordinates = list(node_geometry.exterior.coords)  # Taking the first polygon's coordinates
        y_start_coords = []
        x_start_coords = []
        for coord in start_coordinates:
            y_start_coords.append(coord[1])
            x_start_coords.append(coord[0])        
        avg_start_coordinates = (get_avg(y_start_coords), get_avg(x_start_coords))
        '''

        next_node_id = path[order + 1]
        #next_node = G.nodes[next_node_id]        
        '''
        nn_geometry = next_node['geometry']

        end_coordinates = list(nn_geometry.exterior.coords)
        y_end_coords = []
        x_end_coords = []
        for coord in end_coordinates:
            y_end_coords.append(coord[1])
            x_end_coords.append(coord[0])
        avg_end_coordinates = (get_avg(y_end_coords), get_avg(x_end_coords))
        '''

        #print(node, next_node_id)
        # Such edge has to exist, otherwise something is wrong
        # as of now, something is wrong, it seems to be an algorithm problem as graph seems ok
        #print(node, next_node_id)
        print((node, next_node_id), (node, next_node_id) in G.edges)
        corresponding_edge = G[node][next_node_id]
        #print(corresponding_edge)
        edge_centerline_geometry = list(corresponding_edge.values())[0]['centerline']

        coords = edge_centerline_geometry.exterior.coords
        #print(coords[1], coords[0])

        polyline = folium.PolyLine(
            #locations= [avg_start_coordinates, avg_end_coordinates],
            locations = [coords[1], coords[0]],
            color = 'blue', 
            weight = 3,
            tooltip=f"Order: {order}",
            opacity = 0,
        )
        polyline.add_to(folium_map_object)


        
    # Add custom JavaScript for animation control
    animation_script = Template("""
    <script>    
    window.onload = function() {                       
                                
        let map = Object.values(window).find(v => v instanceof L.Map); // this is Folium's default map variable;        

        let animation;
        let currentStep = 0;
        let prevStep = 0;
        let playing = false;
                                
        let infoText = document.getElementById('info-text');

        /*
        // Add an SVG marker definition for an arrowhead
        let svgNS = "http://www.w3.org/2000/svg";
        let arrowDef = document.createElementNS(svgNS, "marker");
        arrowDef.setAttribute("id", "arrow");
        arrowDef.setAttribute("markerWidth", "10");
        arrowDef.setAttribute("markerHeight", "10");
        arrowDef.setAttribute("refX", "10");
        arrowDef.setAttribute("refY", "3");
        arrowDef.setAttribute("orient", "auto");
        arrowDef.setAttribute("markerUnits", "strokeWidth");

        // Create the actual arrowhead path
        let arrowPath = document.createElementNS(svgNS, "path");
        arrowPath.setAttribute("d", "M0,0 L0,6 L9,3 z"); // Simple triangle for the arrow
        arrowPath.setAttribute("fill", "yellow");  // Same color as the polyline

        // Append the arrowhead to the marker definition
        arrowDef.appendChild(arrowPath);

        // Add the marker to the map's SVG layer
        let svgDefs = document.createElementNS(svgNS, "defs");
        svgDefs.appendChild(arrowDef);
        map.getPanes().overlayPane.querySelector("svg").appendChild(svgDefs);
        */
                                
                                
        //let polylines = document.querySelectorAll(".leaflet-interactive");
        let polylines = [];
        map.eachLayer(function(layer) {
            if(layer instanceof L.Polyline) {                
                if(layer["_path"].getAttribute("stroke") === "blue")                    
                    //layer["_path"].setAttribute("marker-end", "url(#arrow)"); // Add the arrowhead to the end of the line
                    polylines.push(layer["_path"]);
            }
        });                                                                                                      
                            
        
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
                prev_polyline.setAttribute("stroke", "blue");
            }

            let polyline = polylines[currentStep];
            if (polyline) {
                polyline.setAttribute("stroke-opacity", 1);  // Show the polyline
                polyline.setAttribute("stroke", "white");
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
                        polylines[i].setAttribute("stroke-opacity", 0);
                    }
                    //prev_polyline.style["stroke-opacity"] = 0;
                } else 
                {                                
                    prev_polyline.setAttribute("stroke", "blue");
                }
            }

            if (currentStep - prevStep >= 1){
                prev_polyline.setAttribute("stroke") = "blue";
                for (let i = prevStep; i < currentStep; i++) {                        
                        polylines[i].setAttribute("stroke-opacity", 1);
                        polylines[i].setAttribute("stroke", "black");
                        prev_polyline = polylines[i];
                }
            }

            else{
                let polyline = polylines[currentStep];
                if (polyline) {
                    polyline.setAttribute("stroke-opacity", 1);  // Show the polyline
                    polyline.setAttribute("stroke", "black");
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
            polylines = [];
            map.eachLayer(function(layer) {
                if(layer instanceof L.Polyline) {                
                    if(layer["_path"].getAttribute("stroke") === "blue")
                        polylines.push(layer["_path"]);
                }
            });
            polylines.forEach(el => {
                el.setAttribute("stroke-opacity", 0);
            });
        }

        document.getElementById('play-button').onclick = function() { playAnimation(); };
        document.getElementById('pause-button').onclick = function() { pauseAnimation(); };
        document.getElementById('restart-button').onclick = function() { resetAnimation(); playAnimation(); };
    };
    </script>
    """).render()

    # Attach the JavaScript to the map HTML
    folium_map_object.get_root().html.add_child(folium.Element(animation_script))

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
    folium_map_object.get_root().html.add_child(folium.Element(button_html))

    folium_map_object.save("animated_path.html")
    
    # Return the Folium map
    return folium_map_object

    
        


folder_path = 'map_files/observable_geojson_files/'

lane_geojson_file = folder_path + 'Lane_polygons.geojson'
intersection_geojson_file = folder_path + 'Intersection_polygons.geojson'

# Build the graph
G = build_city_graph(lane_geojson_file, intersection_geojson_file)

def get_boundaries(geojson_boundaries):
    lats = []
    lons = []

    coordinates = geojson_boundaries['geometry']['coordinates'][0]    
    for coordinate in coordinates:        
        lats.append(coordinate[0])
        lons.append(coordinate[1])
    
    return (min(lats), min(lons), max(lats), max(lons))

# Define map boundaries (min_latitude, min_longitude, max_latitude, max_longitude)
# Take this data from the Boundary.geojson file
boundary_geojson_file = folder_path + 'Boundary.geojson'
loaded_boundary_geojson = load_geojson(boundary_geojson_file)
map_boundaries = get_boundaries(loaded_boundary_geojson)

# Visualize the created graph for debugging
visualize_graph(G, map_boundaries)


###
# Find the shortest path to visit each graph edge at least once
###

# Try to get the eulerian path of the graph
try:
    eulerian_path = list(nx.eulerian_path(G))    
    print(eulerian_path)

except:
# If the eulerian path doesn't exist, we need to do the following:
# 1) Filter out all such nodes from the graph, which have no outgoing edges.
# 2) Filter out all such edges from the graph, which had removed node as part of it.
# 3) Repeat steps 1-2 until there are no such nodes left in the graph, which have no outgoing edges.
# 4) Perform depth first search on the graph to find the shortest path for visiting each edge at least once.
        
    #print(G.edges)
    # Counter is just for debugging purposes
    counter = 0
    while True:

        nodes_to_remove = []

        #print()
        for node in G.nodes:            
            #print(node, G.out_degree(node))
            if G.out_degree(node) == 0:
                
                nodes_to_remove.append(node)
        #print()
        
        if len(nodes_to_remove) == 0 or counter > 0:
            break
        
        counter += 1
        G.remove_nodes_from(nodes_to_remove)       
    
    visualized_graph = visualize_graph(G, map_boundaries, "debug_map", visualize_nodes=False) 


    def calculate_surplus_and_deficit(graph, traversal_counts = None):        
        """Calculate surplus and deficit nodes based on in-degree and out-degree."""
        in_degrees = defaultdict(int)
        out_degrees = defaultdict(int)
        
        # Calculate the in-degrees and out-degrees
        for u in graph:
            out_degrees[u] += len(graph[u])  # Count outgoing edges
            for v in graph[u]:
                in_degrees[v] += 1  # Count incoming edges
        
        surplus = {}
        deficit = {}
        
        if not traversal_counts is None:
            # Now adjust the surplus and deficit based on traversal counts
            for node in set(list(in_degrees.keys()) + list(out_degrees.keys())):
                out_deg = out_degrees[node]
                in_deg = in_degrees[node]
                
                # Account for extra traversals allowed on edges
                for neighbor in graph[node]:
                    out_deg += traversal_counts[node][neighbor] - 1  # Adjust for extra traversal counts
                
                for neighbor in graph:  # Check for incoming edges
                    in_deg += traversal_counts[neighbor].get(node, 0) - 1  # Adjust for incoming traversals
                
                if out_deg > in_deg:
                    surplus[node] = out_deg - in_deg  # Node with surplus outgoing edges
                elif in_deg > out_deg:
                    deficit[node] = in_deg - out_deg  # Node with surplus incoming edges                        
        else:
            # Identify surplus and deficit nodes
            for node in set(list(in_degrees.keys()) + list(out_degrees.keys())):
                out_deg = out_degrees[node]
                in_deg = in_degrees[node]
                if out_deg > in_deg:
                    surplus[node] = out_deg - in_deg  # Node with surplus outgoing edges
                elif in_deg > out_deg:
                    deficit[node] = in_deg - out_deg  # Node with surplus incoming edges
        
        return surplus, deficit
    
    def balance_graph_with_retraced_edges(graph, surplus, deficit):
        """Balance the graph by retracing edges and adjusting traversal counts."""
        traversal_counts = defaultdict(lambda: defaultdict(int))  # traversal_counts[from][to] = count
        for u in graph:
            for v in graph[u]:
                traversal_counts[u][v] = 1  # Initialize traversal counts
        
        # Handle retracing of edges between surplus and deficit nodes
        for s_node, s_surplus in surplus.items():
            for d_node, d_deficit in deficit.items():
                if s_surplus > 0 and d_deficit > 0:
                    # Use Dijkstra to find the shortest path from s_node to d_node                    
                    shortest_path = nx.dijkstra_path(nx.DiGraph(graph), s_node, d_node)                    
                    
                    # Increment the traversal count for the edges in the path
                    for i in range(len(shortest_path) - 1):
                        u, v = shortest_path[i], shortest_path[i + 1]
                        traversal_counts[u][v] += 1  # Increment traversal count for this edge
                        s_surplus -= 1
                        d_deficit -= 1
                        if s_surplus == 0 or d_deficit == 0:
                            break

        return traversal_counts
    
    def update_graph_traversal_counts(graph, traversal_counts):
        # Handle retracing of edges between surplus and deficit nodes
        for s_node, s_surplus in surplus.items():
            for d_node, d_deficit in deficit.items():
                if s_surplus > 0 and d_deficit > 0:
                    # Use Dijkstra to find the shortest path from s_node to d_node                    
                    shortest_path = nx.dijkstra_path(nx.DiGraph(graph), s_node, d_node)                    
                    
                    # Increment the traversal count for the edges in the path
                    for i in range(len(shortest_path) - 1):
                        u, v = shortest_path[i], shortest_path[i + 1]
                        traversal_counts[u][v] += 1  # Increment traversal count for this edge
                        s_surplus -= 1
                        d_deficit -= 1
                        if s_surplus == 0 or d_deficit == 0:
                            break

    def backtrack_and_adjust(graph, traversal_counts):
        """Backtrack and adjust the graph until all edges are traversed."""       
        t_counts = traversal_counts 
        while True:
            # Recalculate surplus and deficit
            surplus, deficit = calculate_surplus_and_deficit(graph, t_counts)
            
            time.sleep(1)
            print("-"*10)
            print(surplus, deficit)
            print("-"*10)
            
            # If there's no surplus or deficit left, break the loop
            if not surplus and not deficit:
                break
            
            # Balance the graph again by adjusting traversal counts
            #t_counts = balance_graph_with_retraced_edges(graph, surplus, deficit)
            update_graph_traversal_counts(graph, t_counts)
        
        return t_counts

    def construct_weak_eulerian_path(graph, traversal_counts):
        """Construct a weak Eulerian path, visiting every edge at least once using traversal counts."""
        path = []
        current_node = next(iter(graph))  # Start at an arbitrary node (or surplus node)
        
        # Traverse the graph based on traversal counts
        while True:
            # Find the next node to visit based on the traversal count
            for next_node in graph[current_node]:
                if traversal_counts[current_node][next_node] > 0:
                    path.append((current_node, next_node))
                    traversal_counts[current_node][next_node] -= 1  # Decrement the traversal count
                    current_node = next_node
                    break
            else:
                # If no more valid edges are found, break the loop
                break
        
        return path    

    adj_table = defaultdict(set)    

    for el in G.edges:
        adj_table[el[0]].add(el[1])

    graph_list = defaultdict(list)

    for key, value in adj_table.items():
        graph_list[key] = list(value)

    debug_copy = deepcopy(graph_list)

    surplus, deficit = calculate_surplus_and_deficit(graph_list)
    traversal_counts = balance_graph_with_retraced_edges(graph_list, surplus, deficit)    
    traversal_counts = backtrack_and_adjust(graph_list, traversal_counts)
    path = construct_weak_eulerian_path(graph_list, traversal_counts)
    print(path) 
    #print(construct_eulerian_path(graph_list, traversal_counts))

    #mdg = nx.MultiDiGraph(graph_list)
    #print(nx.has_eulerian_path(mdg))
    
    #path = find_chinese_postman_tour_multidigraph(adj_table, G)
    #print(path)

    #path = dfs_multidigraph(G, adj_table, list(G.nodes)[0])    
    #print(path)
    #visualize_path(visualized_graph, G, path)