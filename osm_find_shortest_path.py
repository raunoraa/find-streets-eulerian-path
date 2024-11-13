from jinja2 import Template
import networkx as nx
import geojson
from shapely import Point
from Intersection import Intersection
from shapely.geometry import shape
from shapely.geometry import LineString
from shapely.geometry import Polygon
from collections import defaultdict

# Graph saving for qgis
from datetime import datetime, timedelta
import geopandas as gpd

# Graph visualization
import folium
import random

# OSM XML parsing
import xml.etree.ElementTree as ET

# For calculating the length of a road based on its coordinates
from geopy.distance import geodesic


### GLOBAL VARIABLE FOR STORING THE TOTAL LENGTH OF THE STREETS
total_streets_length = 0.0


# Load GeoJSON files
def load_geojson(file_path):
    with open(file_path, "r") as f:
        return geojson.load(f)


# Function to check if lane is car-drivable
def is_car_drivable(lane):
    # Check for 'backward' access permission if present
    backward_access = (
        lane.get("properties", {})
        .get("muv", {})
        .get("travel", {})
        .get("backward", {})
        .get("access", {})
        .get("LandBased", [{}])[0]
        .get("value", "No")
    )

    # Check for 'forward' access permission if present
    forward_access = (
        lane.get("properties", {})
        .get("muv", {})
        .get("travel", {})
        .get("forward", {})
        .get("access", {})
        .get("LandBased", [{}])[0]
        .get("value", "No")
    )

    # Return True if either access is 'Yes'
    return backward_access == "Yes" or forward_access == "Yes"


def parse_osm_file(file_path, way_id):
    tree = ET.parse(file_path)
    root = tree.getroot()

    way = root.find(f"./way[@id='{way_id}']")

    tags = {tag.attrib["k"]: tag.attrib["v"] for tag in way.findall("tag")}
    street_name = tags.get("name", None)
    
    return street_name


def find_lane_distance(node_coords):

    total_length = 0.0
    for i in range(len(node_coords) - 1):
        start_node = node_coords[i]
        next_node = node_coords[i + 1]
        total_length += geodesic(start_node, next_node).meters

    return total_length


# Load and filter car-drivable lanes
# Returns a defaultdict, where keys are road_ids and values are the lanes (and their respective data)
def parse_lanes(lane_geojson, osm_xml_file_path):
    lanes_by_road_id = defaultdict(
        list
    )  # A dictionary where road_id is the key, and values are lists of lanes
    for lane in lane_geojson["features"]:
        if is_car_drivable(lane):

            # Find out the lane's distance
            osm_way_id = lane["properties"]["osm_way_ids"][0]
            street_name = None
            #street_name = parse_osm_file(osm_xml_file_path, osm_way_id)        

            lane_data = {
                "road_id": lane["properties"]["road"],
                "lane_id": lane["properties"]["index"],
                "direction": lane["properties"]["direction"],
                "src_i": None,  # value will be added in parse_intersections
                "dst_i": None,  # value will be added in parse_intersections
                "geometry": shape(lane["geometry"]),
                "street_name": street_name,
            }
            lanes_by_road_id[lane_data["road_id"]].append(
                lane_data
            )  # Group lanes by road_id
    return lanes_by_road_id


# Load intersections and their movements
# Also assign src and destination intersections to lanes (lanes parameter will be modified)
def parse_intersections(intersection_geojson, all_lanes):
    intersections = set()

    for feature in intersection_geojson["features"]:
        i_type = feature["properties"]["type"]
        if i_type == "road":
            # If the intersection type is road, then we can assign the src and dst intersections to lanes.

            lanes = all_lanes[feature["properties"]["id"]]

            for lane in lanes:                      
                lane["src_i"] = feature["properties"]["src_i"]
                lane["dst_i"] = feature["properties"]["dst_i"]                     

        elif i_type == "intersection":            
            if feature["properties"]["intersection_kind"] != "MapEdge":
                intersection = Intersection(
                    feature["properties"]["id"],
                    feature["properties"]["movements"],
                    shape(feature["geometry"]),
                    feature["properties"]["intersection_kind"],
                )
                intersections.add(intersection)
        else:
            # Let's not look at the other types present in the file.
            continue
    return intersections


def get_twonodes_average_coords(G, start_node, end_node):
    find_avg_coords = lambda coords: (
        sum(y for _, y in coords) / len(coords),
        sum(x for x, _ in coords) / len(coords),
    )
    start_node_data, end_node_data = G.nodes[start_node], G.nodes[end_node]
    start_node_coords, end_node_coords = (
        start_node_data.get("geometry").coords,
        end_node_data.get("geometry").coords,
    )
    avg_start, avg_end = find_avg_coords(start_node_coords), find_avg_coords(
        end_node_coords
    )
    return [avg_start, avg_end]


def get_closest_edge_midpoint(polygon1, polygon2):
    # closest_point1, closest_point2 = nearest_points(polygon1, polygon2)

    closest_point = None
    min_distance = float("inf")

    # Step 1: Find the closest pair of vertices between polygon1 and polygon2 using geodesic distance
    for coord1 in polygon1.exterior.coords:
        for coord2 in polygon2.exterior.coords:
            # Calculate the geodesic distance between coord1 and coord2
            distance = geodesic(
                coord1[::-1], coord2[::-1]
            ).meters  # Reverse coord1/coord2 for (lat, lon)

            # Update the closest points if this pair is closer
            if distance < min_distance:
                min_distance = distance
                closest_point1 = Point(coord1)

    return closest_point1


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
        unique_roads = set()  # a set of unique road ids in an intersection
        for movement in movements:
            src_road, dst_road = movement.split(" -> ")
            
            src_id = int(src_road.split("#")[-1])
            dst_id = int(dst_road.split("#")[-1])
            unique_roads.add(src_id)
            unique_roads.add(dst_id)

        tuples = set()
        for ur in unique_roads:
            lanes_in_ur = lanes.get(ur)
            for lane in lanes_in_ur:

                # Don't add such lanes to intersection, which dont have either src or dst intersection in the intersection list
                if (
                    lane["src_i"] not in intersection_ids
                    or lane["dst_i"] not in intersection_ids
                ):                    
                    continue

                node_id = (id, lane["road_id"], lane["lane_id"])

                # Also add the info, if the node for the current lane is for entering or leaving the intersection
                is_entering_value = False
                if (lane["dst_i"] == id and lane["direction"] == "Forward") or (
                    lane["src_i"] == id and lane["direction"] == "Backward"
                ):
                    is_entering_value = True                
                
                lane_geometry = lane.get("geometry")
                
                int_node_geometry = get_closest_edge_midpoint(
                    lane_geometry, intersection.geometry
                )

                tuples.add(
                    (node_id, is_entering_value, int_node_geometry, id)
                )  # intersection.geometry, id))                
                
                G.add_node(
                    node_id,
                    is_entering=is_entering_value,
                    # geometry=intersection.geometry,
                    geometry=int_node_geometry,
                    intersection_id=id,
                )

        intersections_node_map[id] = tuples

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

                    # Dont add backward turns to intersection
                    if e_node_id[0] == l_node_id[0] and e_node_id[1] == l_node_id[1]:
                        continue

                    # For the edges that are inside the intersection, assign the distance to be 4.0 meters and
                    # street name is None
                    G.add_edge(
                        e_node_id,
                        l_node_id,
                        geometry=LineString([e_geometry, l_geometry]),
                        distance=4.0,
                        street_name=None,
                        edge_type="intersection",
                    )

    # Create edges outside the intersections (basically add the lanes)
    # We need to create an edge between all such nodes, which have the same road and lane id
    observables = {}
    for node in G.nodes(data=True):
        id_tuple, node_data = node
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
                if l.get("lane_id") == lane_id:
                    lane_geometry = l.get("geometry")
                    street_name = l.get("street_name")                    
                    break
            
            coords = get_twonodes_average_coords(G, id_tuple, node_one_id_tuple)
            distance = find_lane_distance(coords)

            if not node_data["is_entering"]:
                G.add_edge(
                    id_tuple,
                    node_one_id_tuple,
                    geometry=lane_geometry,
                    distance=distance,
                    street_name=street_name,
                )
            else:
                G.add_edge(
                    node_one_id_tuple,
                    id_tuple,
                    geometry=lane_geometry,
                    distance=distance,
                    street_name=street_name,
                    edge_type="road",
                )
        else:
            observables[observable_tuple] = node

    return G, intersection_ids


# Main function to construct the graph from geojson files
def build_city_graph(lane_geojson_file, intersection_geojson_file, osm_xml_file_path):
    
    # Load and parse the geojson files
    
    lanes_geojson = load_geojson(lane_geojson_file)
    print("LANE GEOJSON LOADED!")
    
    intersections_geojson = load_geojson(intersection_geojson_file)
    print("INTERSECTIONS GEOJSON LOADED!")


    lanes = parse_lanes(lanes_geojson, osm_xml_file_path)
    print("LANES PARSED!")
    
    intersections = parse_intersections(intersections_geojson, lanes)
    print("INTERSECTIONS PARSED!")


    # Create the directed multigraph
    G, int_ids = create_graph(lanes, intersections)
    
    return G, int_ids


folder_path = "map_files/observable_geojson_files/"

lane_geojson_file = folder_path + "Lane_polygons.geojson"
intersection_geojson_file = folder_path + "Intersection_polygons.geojson"

osm_file_path = "map_files/osm_observable.xml"

print()
print("INIT")
# Build the graph
G, intersection_ids = build_city_graph(lane_geojson_file, intersection_geojson_file, osm_file_path)
print("GRAPH BUILDING FINISHED!")


def get_boundaries(geojson_boundaries):
    lats = []
    lons = []

    coordinates = geojson_boundaries["geometry"]["coordinates"][0]
    for coordinate in coordinates:
        lats.append(coordinate[0])
        lons.append(coordinate[1])

    return (min(lats), min(lons), max(lats), max(lons))


# Define map boundaries (min_latitude, min_longitude, max_latitude, max_longitude)
# Take this data from the Boundary.geojson file
boundary_geojson_file = folder_path + "Boundary.geojson"
loaded_boundary_geojson = load_geojson(boundary_geojson_file)
map_boundaries = get_boundaries(loaded_boundary_geojson)
print("BOUNDARY PARSING FINISHED!")

def assign_edge_orders_with_multiple_traversals(path, G):
    """
    Assign timestamps to edges based on their traversal in the Eulerian path.
    Each edge may be traversed multiple times, so the 'order' attribute will be a list of timestamps.
    """
    current_time = datetime(1980, 1, 1, 0, 0, 0)

    for edge in path:
        # Get the corresponding edge data from the graph (assuming edge is a tuple of nodes)
        edge_data = G.get_edge_data(edge[0], edge[1]).get(0)

        # If the edge doesn't have an 'order' attribute (list), initialize it
        if 'order' not in edge_data:
            edge_data['order'] = []
        
        # Append the current timestamp to the edge's order list
        edge_data['order'].append(current_time)
        
        # Increment the datetime by one second for the next traversal
        current_time += timedelta(seconds=1)

def create_edge_features_for_qgis(path, G):
    """
    Create a list of features for edges, where each edge is duplicated for each timestamp in the 'order' list.
    Each traversal of the same edge gets its own timestamp and geometry.
    """
    edge_data = []  # List to hold the features (edges)

    # Step 1: Assign timestamps to edges based on the Eulerian path
    assign_edge_orders_with_multiple_traversals(path, G)
    
    # Step 2: For each edge, create a feature for each traversal (each timestamp)
    for u, v, edge_attrs in G.edges(data=True):
        if 'order' in edge_attrs:
            # For each timestamp in the 'order' list, duplicate the edge
            for timestamp in edge_attrs['order']:
                # Duplicate the edge with the same geometry and assign the timestamp to 'order'
                edge_attrs_copy = edge_attrs.copy()
                edge_attrs_copy['order'] = timestamp.isoformat()  # Store as string (ISO 8601)

                # Create a dictionary for the feature, including geometry and the timestamp
                edge_data.append({
                    'id': edge_attrs.get('id'),
                    'geometry': edge_attrs.get('geometry'),
                    'order': edge_attrs_copy['order'],
                    'distance': edge_attrs.get('distance'),
                    'edge_type': edge_attrs.get('edge_type'),
                })

    # Step 3: Create a GeoDataFrame for the edges
    edges_gdf = gpd.GeoDataFrame(edge_data, geometry='geometry', crs='EPSG:4326')

    return edges_gdf

def create_node_features_for_qgis(G):
    """
    Create a list of node features for the graph. Each node is represented by a point.
    """
    node_data = []
    
    for node, node_attrs in G.nodes(data=True):
        node_data.append({
            'id': node,
            'geometry': node_attrs.get('geometry')
        })
    
    # Create a GeoDataFrame for the nodes
    nodes_gdf = gpd.GeoDataFrame(node_data, geometry='geometry', crs='EPSG:4326')

    return nodes_gdf

def save_graph_to_geopackage(path, G, output_file='output.gpkg'):
    """
    Create a GeoDataFrame from the graph (nodes and edges) and save it to a GeoPackage.
    """
    # Create the GeoDataFrame for edges with timestamps
    edges_gdf = create_edge_features_for_qgis(path, G)
    
    # Create the GeoDataFrame for nodes
    nodes_gdf = create_node_features_for_qgis(G)

    # Step 4: Save both nodes and edges to a GeoPackage
    edges_gdf.to_file(output_file, layer='edges', driver="GPKG")
    nodes_gdf.to_file(output_file, layer='nodes', driver="GPKG")

    print(f"Graph saved to {output_file}")

###
# Find the shortest path to visit each graph edge at least once
###

if not nx.is_strongly_connected(G):
    components = list(nx.strongly_connected_components(G))
    largest_component = max(
        components, key=lambda c: (len(c), G.subgraph(c).size())
    )
    # print(G)
    G = G.subgraph(largest_component).copy()
    print(G)

in_counter = 0
nodes_to_remove = []
for node in G.nodes:
    # Remove the dead ends.
    if G.out_degree(node) == 0:
        if G.in_degree(node) == 1 and in_counter == 0:
            in_counter += 1
        else:
            nodes_to_remove.append(node)    
G.remove_nodes_from(nodes_to_remove)

###
# Remove edges inside the intersections
# (Sometimes makes the coefficient better, sometimes worse)
###
edges = list(G.edges(data=True))
for u, v, data in edges:
    if data.get("edge_type") == "intersection":
        if not (G.out_degree(u) <= 1 or G.in_degree(v) <= 1):
            G.remove_edge(u, v)
            if not nx.is_strongly_connected(G):
                G.add_edge(u, v, **data)

print("GRAPH MODIFYING FINISHED!")

# Calculate the total street distance
for u, v, data in G.edges(data=True):
    total_streets_length += data.get("distance")


# Try to get the eulerian path of the graph
try:        
    eulerian_path = list(nx.eulerian_path(G))
        
    print("EULERIAN PATH FOUND WITHOUT BALANCING!")
    
    save_graph_to_geopackage(eulerian_path, G)
    

except:
    ###
    # We need to balance the graph if no eulerian path was found.
    ###
    
    print("EULERIAN PATH NOT FOUND WITHOUT BALANCING!")
    print("BALANCING THE GRAPH TO FIND THE EULERIAN PATH...")   

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

        surplus_copy, deficit_copy = surplus.copy(), deficit.copy()

        # Handle retracing of edges between surplus and deficit nodes
        for s_node in surplus.keys():
            for d_node in deficit.keys():
                if surplus[s_node] > 0 and deficit[d_node] > 0:
                    # It will suffice for an eulerian path,
                    #   if we have exactly one surplus node and one deficit node left
                    #   and both of such nodes have an offset (from 0) of one.
                    if (len(surplus_copy) == 1 and len(deficit_copy) == 1) and (
                        surplus[s_node] == 1 and deficit[d_node] == 1
                    ):
                        surplus = surplus_copy
                        deficit = deficit_copy
                        return s_node

                    # Use Dijkstra algorithm to find the shortest path from d_node to s_node
                    # Calculate the shortest path based on the distance of an edge.
                    shortest_path = nx.dijkstra_path(
                        G, d_node, s_node, weight="distance"
                    )

                    surplus[s_node] -= 1
                    deficit[d_node] -= 1

                    if surplus[s_node] == 0:
                        del surplus_copy[s_node]
                    if deficit[d_node] == 0:
                        del deficit_copy[d_node]

                    for i in range(len(shortest_path) - 1):
                        u, v = shortest_path[i], shortest_path[i + 1]

                        # Add a copy of an existing edge to the graph
                        edge_key = list(G.get_edge_data(u, v).keys())[0]
                        edge_attributes = G.get_edge_data(u, v, key=edge_key)
                        G.add_edge(u, v, **edge_attributes)
        surplus = surplus_copy
        deficit = deficit_copy

        return s_node

    copied_graph = nx.MultiDiGraph(G)

    surplus = {}
    deficit = {}
    surplus, deficit = calculate_surplus_and_deficit(copied_graph)
    start_node = balance_graph(copied_graph, surplus, deficit)
    
    print("GRAPH BALANCING FINISHED!")
    
    eulerian_path = list(nx.eulerian_path(copied_graph, source=start_node))
        
    print("EULERIAN PATH FOUND!")
    
    save_graph_to_geopackage(eulerian_path, G)
