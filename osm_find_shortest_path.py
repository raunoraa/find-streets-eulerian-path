import networkx as nx
import geojson
from Intersection import Intersection
from shapely.geometry import shape, LineString, Point
from collections import defaultdict

# Graph saving for qgis
from datetime import datetime, timedelta
import geopandas as gpd

# For calculating the length of a road based on its coordinates
from geopy.distance import geodesic

# For parsing the osm.xml file
import xml.etree.ElementTree as ET

from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import time

import numpy as np
from scipy.optimize import linear_sum_assignment


def load_osm_to_dict(osm_file_path):
    # {way_id: {tag1:value1, tag2:value2, ...}, ...}
    osm_dict = {}

    # Use `iterparse` to load the XML incrementally
    context = ET.iterparse(osm_file_path, events=("start", "end"))
    for event, elem in context:
        if event == "end" and elem.tag == "way":
            osm_id = elem.get("id")
            if osm_id:
                # Collect tags into a dictionary for the current way
                tags = {tag.get("k"): tag.get("v") for tag in elem.findall("tag")}
                osm_dict[osm_id] = tags
            # Clear the element from memory
            elem.clear()

    return osm_dict


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


def find_lane_distance(node_coords):
    start_node = node_coords[0]
    next_node = node_coords[1]
    return geodesic(start_node, next_node).meters


# Load and filter car-drivable lanes
# Returns a defaultdict, where keys are road_ids and values are the lanes (and their respective data)
def parse_lanes(lane_geojson, osm_dict):
    lanes_by_road_id = defaultdict(
        list
    )  # A dictionary where road_id is the key, and values are lists of lanes

    allowed_tags = set(
        [
            "trunk",
            "trunk_link",
            "primary",
            "primary_link",
            "secondary",
            "secondary_link",
            "tertiary",
            "tertiary_link",
            "residential",
            "unclassified",
        ]
    )
    for lane in lane_geojson["features"]:
        if is_car_drivable(lane):

            osm_way_id = lane["properties"]["osm_way_ids"][0]
            highway_tag_value = osm_dict[str(osm_way_id)]["highway"]
            if highway_tag_value in allowed_tags:
                lane_data = {
                    "road_id": lane["properties"]["road"],
                    "lane_id": lane["properties"]["index"],
                    "direction": lane["properties"]["direction"],
                    "src_i": None,  # value will be added in parse_intersections
                    "dst_i": None,  # value will be added in parse_intersections
                    "geometry": shape(lane["geometry"]),
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


# Create the directed graph
def create_graph(lanes, intersections):
    G = nx.DiGraph()
    G_with_non_compulsory_edges = nx.DiGraph()

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
                    geometry=int_node_geometry,
                    intersection_id=id,
                )
                G_with_non_compulsory_edges.add_node(
                    node_id,
                    is_entering=is_entering_value,
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

                    coords = get_twonodes_average_coords(G, e_node_id, l_node_id)
                    intersection_road_distance = find_lane_distance(coords)

                    # Dont add backward turns to intersection (add to the graph, which contains non-compulsory edges)
                    if e_node_id[0] == l_node_id[0] and e_node_id[1] == l_node_id[1]:
                        G_with_non_compulsory_edges.add_edge(
                            e_node_id,
                            l_node_id,
                            geometry=LineString([e_geometry, l_geometry]),
                            distance=intersection_road_distance,  # 4.0,
                            edge_type="intersection",
                        )
                        continue

                    # For the edges that are inside the intersection, assign the distance to be 4.0 meters and
                    # street name is None
                    G.add_edge(
                        e_node_id,
                        l_node_id,
                        geometry=LineString([e_geometry, l_geometry]),
                        distance=intersection_road_distance,  # 4.0,
                        edge_type="intersection",
                    )
                    G_with_non_compulsory_edges.add_edge(
                        e_node_id,
                        l_node_id,
                        geometry=LineString([e_geometry, l_geometry]),
                        distance=intersection_road_distance,  # 4.0,
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
            # Get the lane geometry from the lanes defaultdict
            for l in lanes[road_id]:
                if l.get("lane_id") == lane_id:
                    lane_geometry = l.get("geometry")
                    break

            coords = get_twonodes_average_coords(G, id_tuple, node_one_id_tuple)
            distance = find_lane_distance(coords)

            if not node_data["is_entering"]:
                G.add_edge(
                    id_tuple,
                    node_one_id_tuple,
                    geometry=lane_geometry,
                    distance=distance,
                    edge_type="road",
                )
                G_with_non_compulsory_edges.add_edge(
                    id_tuple,
                    node_one_id_tuple,
                    geometry=lane_geometry,
                    distance=distance,
                    edge_type="road",
                )
            else:
                G.add_edge(
                    node_one_id_tuple,
                    id_tuple,
                    geometry=lane_geometry,
                    distance=distance,
                    edge_type="road",
                )
                G_with_non_compulsory_edges.add_edge(
                    node_one_id_tuple,
                    id_tuple,
                    geometry=lane_geometry,
                    distance=distance,
                    edge_type="road",
                )
        else:
            observables[observable_tuple] = node

    return G, intersection_ids, G_with_non_compulsory_edges


# Main function to construct the graph from geojson files
def build_city_graph(lane_geojson_file, intersection_geojson_file, osm_file):

    # Load and parse the geojson files and osm.xml file
    osm_dict = load_osm_to_dict(osm_file)
    print("OSM XML LOADED!")

    lanes_geojson = load_geojson(lane_geojson_file)
    print("LANE GEOJSON LOADED!")

    intersections_geojson = load_geojson(intersection_geojson_file)
    print("INTERSECTIONS GEOJSON LOADED!")

    lanes = parse_lanes(lanes_geojson, osm_dict)
    print("LANES PARSED!")

    intersections = parse_intersections(intersections_geojson, lanes)
    print("INTERSECTIONS PARSED!")

    # Create the directed multigraph
    G, int_ids, G_with_non_compulsory_edges = create_graph(lanes, intersections)

    return G, int_ids, G_with_non_compulsory_edges


def assign_edge_orders_with_multiple_traversals(path, G):
    """
    Assign timestamps to edges based on their traversal in the Eulerian path.
    Each edge may be traversed multiple times, so the 'order' attribute will be a list of timestamps.
    Order nr is just an index (starts from 1).
    """
    current_time = datetime(1980, 1, 1, 0, 0, 0)
    order_nr = 0

    for edge in path:
        # Get the corresponding edge data from the graph
        edge_data = G.get_edge_data(edge[0], edge[1])

        # If the edge doesn't have an 'order' attribute (list), initialize it
        if "order" not in edge_data:
            edge_data["order"] = []
            edge_data["order_nr"] = []

        # Append the current timestamp to the edge's order list
        edge_data["order"].append(current_time)
        edge_data["order_nr"].append(order_nr)

        # Increment the datetime by one second for the next traversal
        current_time += timedelta(seconds=1)
        order_nr += 1


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
        # For each timestamp in the 'order' list, duplicate the edge
        for i, timestamp in enumerate(edge_attrs["order"]):
            # Duplicate the edge with the same geometry and assign the timestamp to 'order'

            geometry = edge_attrs["geometry"]
            if isinstance(geometry, LineString):
                edge_attrs["geometry"] = geometry.buffer(0.000002)

            # Create a dictionary for the feature, including geometry and the timestamp
            edge_data.append(
                {
                    "node_ids": str((u, v)),
                    "geometry": edge_attrs.get("geometry"),
                    "order": timestamp,
                    "order_nr": edge_attrs["order_nr"][i],
                    "distance": edge_attrs.get("distance"),
                    "edge_type": edge_attrs.get("edge_type"),
                }
            )

    # Step 3: Create a GeoDataFrame for the edges
    edges_gdf = gpd.GeoDataFrame(edge_data, geometry="geometry", crs="EPSG:4326")

    return edges_gdf


def create_node_features_for_qgis(G):
    """
    Create a list of node features for the graph. Each node is represented by a point.
    """
    node_data = []

    for node, node_attrs in G.nodes(data=True):
        node_data.append({"id": str(node), "geometry": node_attrs.get("geometry")})

    # Create a GeoDataFrame for the nodes
    nodes_gdf = gpd.GeoDataFrame(node_data, geometry="geometry", crs="EPSG:4326")

    return nodes_gdf


def save_results_as_txt(eulerian_path_distance, coefficient):
    distances_string = f"Street lanes distance: {round(total_streets_length, 4)}m ;; Eulerian path distance: {round(eulerian_path_distance, 4)}m"
    coefficient_string = f"Coefficient: {round(coefficient, 4)}"

    # Write the results to a text file as well
    results_file = open("results.txt", "w", encoding="utf-8")
    results_file.write(distances_string)
    results_file.write("\n")
    results_file.write(coefficient_string)
    results_file.close()

    print()
    print(distances_string)
    print(coefficient_string)
    print()


def calculate_path_distance_and_coefficient(path, G):
    path_distance = 0.0

    for i, edge in enumerate(path):

        corresponding_edge = G.get_edge_data(edge[0], edge[1])
        # print(str(edge[0])+str(edge[1])+":", G.has_edge(edge[0], edge[1]), ";;", "Corresponding edge:", corresponding_edge)
        distance = corresponding_edge.get("distance")

        path_distance += distance

    coefficient = path_distance / total_streets_length

    return path_distance, coefficient


def save_graph_to_geopackage(path, G, output_file="output.gpkg"):
    """
    Create a GeoDataFrame from the graph (nodes and edges) and save it to a GeoPackage.
    """
    # Create the GeoDataFrame for edges with timestamps
    edges_gdf = create_edge_features_for_qgis(path, G)

    # Create the GeoDataFrame for nodes
    nodes_gdf = create_node_features_for_qgis(G)

    # Step 4: Save both nodes and edges to a GeoPackage
    edges_gdf.to_file(output_file, layer="edges", driver="GPKG")
    nodes_gdf.to_file(output_file, layer="nodes", driver="GPKG")

    # Save the eulerian path length and coefficient results to txt file.
    path_d, coef = calculate_path_distance_and_coefficient(path, G)
    save_results_as_txt(path_d, coef)

    print(f"Graph saved to {output_file}")


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


def compute_path(pair, G_with_non_compulsory_edges):
    deficit_node, surplus_node = pair
    distance, shortest_path = nx.single_source_dijkstra(
        G_with_non_compulsory_edges,
        source=deficit_node,
        target=surplus_node,
        weight="distance",
    )
    return distance, shortest_path

def process_chunk(chunk, G_with_non_compulsory_edges):
    chunk_distances = {}
    chunk_path_lookup = {}
    for deficit_node, surplus_node in chunk:
        # Compute shortest path for this unique pair only once
        distance, shortest_path = compute_path(
            (deficit_node, surplus_node), G_with_non_compulsory_edges
        )

        # Store the distance and path for future use
        chunk_distances[(deficit_node, surplus_node)] = distance
        chunk_path_lookup[(deficit_node, surplus_node)] = shortest_path
    return chunk_distances, chunk_path_lookup

def calculate_cost_matrix(G_with_non_compulsory_edges, deficit, surplus):
    """
    Calculate the shortest path distances between all unique deficit and surplus nodes.
    Returns a cost matrix and a mapping to the original nodes.
    """
    # Precompute the shortest paths only once per unique deficit-surplus pair
    distances = {}
    path_lookup = {}
    
    print("CALCULATING SHORTEST PATHS")    
    
    # Serial solution
    for deficit_node in deficit.keys():
        for surplus_node in surplus.keys():
            # Compute shortest path for this unique pair only once
            distance, shortest_path = compute_path(
                (deficit_node, surplus_node), G_with_non_compulsory_edges
            )

            # Store the distance and path for future use
            distances[(deficit_node, surplus_node)] = distance
            path_lookup[(deficit_node, surplus_node)] = shortest_path
    
    '''
    # Parallel solution
    print("CREATE PAIRS")
    pairs = [(deficit_node, surplus_node) for deficit_node in deficit.keys() for surplus_node in surplus.keys()]    
    print("PAIRS CREATED")
    
    num_workers = os.cpu_count()
    print("CPU count on the current system:", num_workers)
    chunk_size = len(pairs) // num_workers + 1
    chunks = [pairs[i: i + chunk_size] for i in range(0, len(pairs), chunk_size)]
    
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {}        
        for chunk in chunks:
            futures[executor.submit(process_chunk, chunk, G_with_non_compulsory_edges)] = chunk
                        
        # Process completed tasks as they finish
        for future in as_completed(futures):
            chunk_distances, chunk_path_lookup = future.result()
            # Store the distances and paths for future use
            distances.update(chunk_distances)
            path_lookup.update(chunk_path_lookup)            
    '''
    print("UNIQUE SHORTEST PATHS CALCULATED")

    # Expand the deficit and surplus nodes based on their values
    expanded_deficit_nodes = []
    expanded_surplus_nodes = []

    # Keep track of the original node for each expanded entry
    deficit_mapping = []
    surplus_mapping = []

    # Expand based on deficit values
    for node, count in deficit.items():
        for _ in range(count):
            expanded_deficit_nodes.append(node)
            deficit_mapping.append(node)

    # Expand based on surplus values
    for node, count in surplus.items():
        for _ in range(count):
            expanded_surplus_nodes.append(node)
            surplus_mapping.append(node)

    num_deficit = len(expanded_deficit_nodes)
    num_surplus = len(expanded_surplus_nodes)

    # Initialize the cost matrix
    cost_matrix = np.zeros((num_deficit, num_surplus))

    # Populate the cost matrix using pre-computed distances
    for i, deficit_node in enumerate(expanded_deficit_nodes):
        for j, surplus_node in enumerate(expanded_surplus_nodes):
            cost_matrix[i, j] = distances[(deficit_node, surplus_node)]

    print("EXPANDED COST MATRIX CREATED!")
    return (
        cost_matrix,
        path_lookup,
        expanded_deficit_nodes,
        expanded_surplus_nodes,
        deficit_mapping,
        surplus_mapping,
    )


def balance_graph(G, surplus, deficit, G_with_non_compulsory_edges, initial_graph):
    """Balance the graph by adding duplicate edges using the Hungarian algorithm."""

    surplus_node = None

    # Precompute all shortest paths and create an expanded cost matrix
    start_time = time.time()  # debugging
    (
        cost_matrix,
        path_lookup,
        expanded_deficit_nodes,
        expanded_surplus_nodes,
        deficit_mapping,
        surplus_mapping,
    ) = calculate_cost_matrix(G_with_non_compulsory_edges, deficit, surplus)

    print(
        "Time spent on calculating shortest paths:",
        round((time.time() - start_time) / 60, 1),
        "minutes",
    )  # debugging

    # Use the Hungarian algorithm to find the optimal assignment
    deficit_indices, surplus_indices = linear_sum_assignment(cost_matrix)

    # Iterate through the optimal matches
    for d_idx, s_idx in zip(deficit_indices, surplus_indices):
        deficit_node = expanded_deficit_nodes[d_idx]
        surplus_node = expanded_surplus_nodes[s_idx]

        # Retrieve the shortest path for the optimal match
        best_path = path_lookup[(deficit_node, surplus_node)]

        # Apply the selected path to balance the graph
        for i in range(len(best_path) - 1):
            u, v = best_path[i], best_path[i + 1]

            if not initial_graph.has_edge(u, v):
                edge_attributes = G_with_non_compulsory_edges.get_edge_data(u, v)
                initial_graph.add_edge(u, v, **edge_attributes)
                G.add_edge(u, v, **edge_attributes)
            else:
                # Add a copy of an existing edge to the graph.
                edge_attributes = G.get_edge_data(u, v).get(0)
                G.add_edge(u, v, **edge_attributes)

        # Reduce the surplus and deficit of the corresponding nodes.
        surplus[surplus_mapping[s_idx]] -= 1
        deficit[deficit_mapping[d_idx]] -= 1

        # Delete the dictionary entries if the node is no longer in surplus/deficit.
        if surplus[surplus_mapping[s_idx]] == 0:
            del surplus[surplus_mapping[s_idx]]
        if deficit[deficit_mapping[d_idx]] == 0:
            del deficit[deficit_mapping[d_idx]]

        # Check the early stopping condition
        if (len(surplus) == 1 and len(deficit) == 1) and (
            surplus[next(iter(surplus))] == 1 and deficit[next(iter(deficit))] == 1
        ):
            surplus_node = next(iter(surplus))
            print("Early stopping condition met!")
            break

    # Start traversing the path from the last observed surplus node.
    return surplus_node


if __name__ == "__main__":

    ### GLOBAL VARIABLE FOR STORING THE TOTAL LENGTH OF THE STREETS
    total_streets_length = 0.0

    folder_path = "map_files/observable_geojson_files/"

    lane_geojson_file = folder_path + "Lane_polygons.geojson"
    intersection_geojson_file = folder_path + "Intersection_polygons.geojson"

    osm_file = "map_files/osm_observable.xml"

    print()
    print("INIT")
    # Build the graph
    G, intersection_ids, G_with_non_compulsory_edges = build_city_graph(
        lane_geojson_file, intersection_geojson_file, osm_file
    )
    print("GRAPH BUILDING FINISHED!")

    # Find the largest strongly connected component of the graph
    if not nx.is_strongly_connected(G):
        components = list(nx.strongly_connected_components(G))
        largest_component = max(
            components, key=lambda c: (len(c), G.subgraph(c).size())
        )
        # print(G)
        G = G.subgraph(largest_component).copy()
        # print(G)
    if not nx.is_strongly_connected(G_with_non_compulsory_edges):
        components = list(nx.strongly_connected_components(G_with_non_compulsory_edges))
        largest_component = max(
            components,
            key=lambda c: (len(c), G_with_non_compulsory_edges.subgraph(c).size()),
        )
        G_with_non_compulsory_edges = G_with_non_compulsory_edges.subgraph(
            largest_component
        ).copy()

    ###
    # Remove as much edges inside the intersections as possible
    ###
    edges = list(G.edges(data=True))
    for u, v, data in edges:
        if data.get("edge_type") == "intersection":
            if not (G.out_degree(u) <= 1 or G.in_degree(v) <= 1):
                G.remove_edge(u, v)
                if not nx.is_strongly_connected(G):
                    G.add_edge(u, v, **data)
    print("GRAPH MODIFYING FINISHED!")
    print(G)

    # Calculate the total street distance
    for u, v, data in G.edges(data=True):
        total_streets_length += data["distance"]

    # Try to get the eulerian path of the graph
    try:
        eulerian_path = list(nx.eulerian_path(G))

        # print(len(list(G.edges())), len(eulerian_path))

        print("EULERIAN PATH FOUND WITHOUT BALANCING!")

        save_graph_to_geopackage(eulerian_path, G)

    except:
        ###
        # We need to balance the graph if no eulerian path was found.
        ###

        print("EULERIAN PATH NOT FOUND WITHOUT BALANCING!")
        print("BALANCING THE GRAPH TO FIND THE EULERIAN PATH...")

        copied_graph = nx.MultiDiGraph(G)

        surplus = {}
        deficit = {}
        surplus, deficit = calculate_surplus_and_deficit(copied_graph)

        start_node = balance_graph(
            copied_graph, surplus, deficit, G_with_non_compulsory_edges, G
        )

        print("GRAPH BALANCING FINISHED!")

        # Show the edge count for graph before and after the balancing.
        # print(len(list(G.edges())), len(list(copied_graph.edges())))

        eulerian_path = list(nx.eulerian_path(copied_graph, source=start_node))

        print("EULERIAN PATH FOUND!")

        save_graph_to_geopackage(eulerian_path, G)
