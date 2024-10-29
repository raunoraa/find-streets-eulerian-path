from jinja2 import Template
import networkx as nx
import geojson
from Intersection import Intersection
from shapely.geometry import shape
from collections import defaultdict

# Graph visualization
import folium
import random

# OSM XML parsing
import xml.etree.ElementTree as ET

# For calculating the length of a road based on its coordinates
from shapely.geometry import Polygon, MultiLineString, LineString
from shapely.ops import unary_union
from geopy.distance import geodesic
import numpy as np
from scipy.spatial import Voronoi

import time


### GLOBAL VARIABLE FOR STORING THE TOTAL LENGTH OF THE STREETS
total_streets_length = 0.0


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

def parse_osm_file(file_path, way_id):
    tree = ET.parse(file_path)
    root = tree.getroot()

    way = root.find(f"./way[@id='{way_id}']")    

    tags = {tag.attrib['k']: tag.attrib['v'] for tag in way.findall('tag')}    
    street_name = tags.get('name', None)
    return street_name


def find_lane_distance(node_coords):

    total_length = 0.0    
    for i in range(len(node_coords) - 1):
        start_node = node_coords[i]
        next_node = node_coords[i+1]
        total_length += geodesic(start_node, next_node).meters

    return total_length


# Load and filter car-drivable lanes
# Returns a defaultdict, where keys are road_ids and values are the lanes (and their respective data) 
def parse_lanes(lane_geojson, osm_xml_file_path):
    lanes_by_road_id = defaultdict(list)  # A dictionary where road_id is the key, and values are lists of lanes
    for lane in lane_geojson['features']:
        if is_car_drivable(lane):

            # Find out the lane's distance
            osm_way_id = lane['properties']['osm_way_ids'][0]
            street_name = parse_osm_file(osm_xml_file_path, osm_way_id)                        

            lane_data = {
                'road_id': lane['properties']['road'], 
                'lane_id': lane['properties']['index'],
                'direction': lane['properties']['direction'], 
                'src_i': None, # value will be added in parse_intersections
                'dst_i': None, # value will be added in parse_intersections
                'geometry': shape(lane['geometry']),
                'street_name': street_name,
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

def get_twonodes_average_coords(G, start_node, end_node):
    find_avg_coords = lambda coords : (
        sum(y for _, y in coords) / len(coords),
        sum(x for x, _ in coords) / len(coords)
    )
    start_node_data, end_node_data = G.nodes[start_node], G.nodes[end_node]
    start_node_coords, end_node_coords = start_node_data.get('geometry').exterior.coords, end_node_data.get('geometry').exterior.coords        
    avg_start, avg_end = find_avg_coords(start_node_coords), find_avg_coords(end_node_coords)
    return [avg_start, avg_end]

# Create the directed multigraph
def create_graph(lanes, intersections):
    global total_streets_length

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

                tuples.add((node_id, is_entering_value, intersection.geometry, id))
                G.add_node(node_id, is_entering=is_entering_value, geometry=intersection.geometry, intersection_id=id)
                counter += 1
        
        intersections_node_map[id] = tuples

        # debug
        #print(f"Nodes in intersection {id}: {counter}")
    
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
                    
                    # For the edges that are inside the intersection, assign the distance to be 4.0 meters and
                    # street name is None

                    total_streets_length += 4.0
                    G.add_edge(e_node_id, l_node_id, geometry=e_geometry, distance=4.0, street_name=None)
    

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
            street_name = None
            # Get the lane geometry from the lanes defaultdict
            for l in lanes[road_id]:
                if l.get('lane_id') == lane_id:
                    lane_geometry = l.get('geometry')                    
                    street_name = l.get('street_name')
                    break  

            coords = get_twonodes_average_coords(G, id_tuple, node_one_id_tuple)
            distance = find_lane_distance(coords)
            
            total_streets_length += distance

            if not node_data['is_entering']:
                G.add_edge(id_tuple, node_one_id_tuple, geometry = lane_geometry, distance=distance, street_name=street_name)
            else:
                G.add_edge(node_one_id_tuple, id_tuple, geometry = lane_geometry, distance=distance, street_name=street_name)
        else:
            observables[observable_tuple] = node    

    return G

# Main function to construct the graph from geojson files
def build_city_graph(lane_geojson_file, intersection_geojson_file, osm_xml_file_path):
    # Load and parse the geojson files
    lanes_geojson = load_geojson(lane_geojson_file)
    intersections_geojson = load_geojson(intersection_geojson_file)

    lanes = parse_lanes(lanes_geojson, osm_xml_file_path)
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
            if visualize_nodes:
                folium.Marker(
                    location=(coordinates[0][1] + random.uniform(-0.00003, 0.00003), coordinates[0][0] + random.uniform(-0.00003, 0.00003)),
                    popup=f"Out degree: {G.out_degree(node)}; Int.id: {data['intersection_id']}",  # Optional popup text
                    icon=folium.Icon(color='blue')  # Optional: Customize the marker color
                ).add_to(m)

    # Add edges to the map (assuming edges are also polygons)
    for u, v, data in G.edges(data=True):
        geometry = data['geometry']
        distance = data['distance']
        street_name = data['street_name']
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
                popup=f'Edge from {u} to {v}\nDistance: {distance} meters\nStreet name: {street_name}'
            ).add_to(m)
    
    # Save the map as an html file
    m.save(file_name + '.html')    

    return m

def visualize_path(folium_map_object, G, path):   

    total_distance = 0.0

    for i, edge in enumerate(path):
        
        corresponding_edge = G.get_edge_data(edge[0], edge[1]).get(0)
        distance = corresponding_edge.get('distance')
        street_name = corresponding_edge.get('street_name')

        total_distance += distance
        
        start_node_id, end_node_id = edge
        coords = get_twonodes_average_coords(G, start_node_id, end_node_id)

        folium.PolyLine(
            locations = coords,
            color = 'blue',
            weight = 3,
            tooltip=f"Order: {i}\nDistance: {distance}m\nStreet name: {street_name}",
            opacity = 0,
        ).add_to(folium_map_object)        

    print()
    print(f"Street lanes distance: {round(total_streets_length, 4)}m ;; Eulerian path distance: {round(total_distance, 4)}m")
    print("Coefficient:", round(total_distance/total_streets_length, 4))
    print()
        
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

osm_file_path = 'map_files/osm_observable.xml'

# Build the graph
G = build_city_graph(lane_geojson_file, intersection_geojson_file, osm_file_path)

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
    visualized_graph = visualize_graph(G, map_boundaries, "debug_map", visualize_nodes=False)   
    visualize_path(visualized_graph, G, eulerian_path)

except:
# If the eulerian path doesn't exist, we need to do the following:
# 1) Filter out all such nodes from the graph, which have no outgoing edges.
# 2) Filter out all such edges from the graph, which had removed node as part of it.
# 3) Repeat steps 1-2 until there are no such nodes left in the graph, which have no outgoing edges.
# 4) Find the shortest path for visiting each edge at least once.
        
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
    
    #print(nx.is_strongly_connected(G))
    visualized_graph = visualize_graph(G, map_boundaries, "debug_map", visualize_nodes=False) 


    def calculate_surplus_and_deficit(graph):        
        """Calculate surplus and deficit nodes based on in-degree and out-degree."""
    
        surplus = {}
        deficit = {}

        for node in graph.nodes():
            out_degree = graph.out_degree(node)
            in_degree = graph.in_degree(node)
            
            # Calculate the difference between out-degree and in-degree
            difference = out_degree - in_degree

            # Classify as surplus or deficit based on the difference
            if difference > 0:
                surplus[node] = difference  # More outgoing edges than incoming
            elif difference < 0:
                deficit[node] = abs(difference)  # More incoming edges than outgoing
        
        return surplus, deficit
    

    def balance_graph(G, surplus, deficit):
        """Balance the graph by retracing edges and adjusting traversal counts."""        

        # Handle retracing of edges between surplus and deficit nodes
        for s_node in surplus.keys():
            for d_node in deficit.keys():
                if surplus[s_node] > 0 and deficit[d_node] > 0:
                    # Use Dijkstra to find the shortest path from d_node to s_node                    
                    shortest_path = nx.dijkstra_path(G, d_node, s_node)#, G.edges[d_node, s_node]["distance"])                    
                    
                    surplus[s_node] -= 1
                    deficit[d_node] -= 1

                    for i in range(len(shortest_path) - 1):
                        u, v = shortest_path[i], shortest_path[i+1]                        

                        # Add a copy of an existing edge to the graph
                        edge_key = list(G.get_edge_data(u, v).keys())[0]
                        edge_attributes = G.get_edge_data(u, v, key=edge_key)
                        G.add_edge(u, v, **edge_attributes)
        
        for key in list(surplus.keys()):
            if surplus[key] == 0:
                del surplus[key]

        for key in list(deficit.keys()):
            if deficit[key] == 0:
                del deficit[key]  


    adj_table = defaultdict(set)    

    for el in G.edges:
        adj_table[el[0]].add(el[1])

    graph_list = defaultdict(list)

    for key, value in adj_table.items():
        graph_list[key] = list(value)

    copied_graph = nx.MultiDiGraph(G)

    surplus = {}
    deficit = {}    
    surplus, deficit = calculate_surplus_and_deficit(copied_graph)
    balance_graph(copied_graph, surplus, deficit)                
    #print(nx.is_strongly_connected(copied_graph)) # for debugging
    eulerian_path = list(nx.eulerian_path(copied_graph))        
    visualize_path(visualized_graph, G, eulerian_path)
