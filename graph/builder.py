r"""
This submodule creates the graphs.
"""

import networkx as nx
from shapely.geometry import LineString

from file_io.loader import load_geojson, load_osm_to_dict

from geo_operations.geo_utils import (
    get_closest_polygon_point,
    get_twonodes_average_coords,
    find_lane_distance,
)

from graph.data_preprocessor import preprocess_lanes, preprocess_intersections


def create_graph(lanes, intersections):
    """
    Create directed graphs representing the road network with intersections and lanes.

    Args:
        lanes (dict): Dictionary where keys are road IDs and values are lists of lane dictionaries.
        intersections (list): List of Intersection objects representing intersections.

    Returns:
        tuple: A tuple containing:
            - nx.DiGraph: The main directed graph including compulsory edges.
            - nx.DiGraph: A directed graph that also includes non-compulsory edges.
    """
    G = nx.DiGraph()
    G_with_non_compulsory_edges = nx.DiGraph()

    # Key: intersection ID, Value: set of nodes in the intersection
    intersections_node_map = {}
    intersection_ids = [el.id for el in intersections]

    # Add nodes to the graph, grouped by their intersection
    for intersection in intersections:
        id = intersection.id
        movements = intersection.movements
        unique_roads = set()  # Collect unique road IDs in the intersection

        # Extract source and destination road IDs from movements
        for movement in movements:
            src_road, dst_road = movement.split(" -> ")
            src_id = int(src_road.split("#")[-1])
            dst_id = int(dst_road.split("#")[-1])
            unique_roads.update([src_id, dst_id])

        tuples = set()
        for ur in unique_roads:
            lanes_in_ur = lanes.get(ur, [])
            for lane in lanes_in_ur:
                # Skip lanes that don't have both src and dst intersections in the intersection list
                if (
                    lane["src_i"] not in intersection_ids
                    or lane["dst_i"] not in intersection_ids
                ):
                    continue

                node_id = (id, lane["road_id"], lane["lane_id"])

                # Determine if the node represents entering or leaving the intersection
                is_entering_value = (
                    lane["dst_i"] == id and lane["direction"] == "Forward"
                ) or (lane["src_i"] == id and lane["direction"] == "Backward")

                lane_geometry = lane.get("geometry")
                int_node_geometry = get_closest_polygon_point(
                    lane_geometry, intersection.geometry
                )

                # Add the node to the graph
                tuples.add((node_id, is_entering_value, int_node_geometry))
                G.add_node(
                    node_id,
                    is_entering=is_entering_value,
                    geometry=int_node_geometry,
                )
                G_with_non_compulsory_edges.add_node(
                    node_id,
                    is_entering=is_entering_value,
                    geometry=int_node_geometry,
                )

        intersections_node_map[id] = tuples

    # Create edges within the intersection
    for _, nodes in intersections_node_map.items():
        entering_nodes = [n for n in nodes if n[1]]
        leaving_nodes = [n for n in nodes if not n[1]]

        # Add edges between entering and leaving nodes
        for e_node in entering_nodes:
            e_node_id, _, e_geometry = e_node
            for l_node in leaving_nodes:
                l_node_id, _, l_geometry = l_node

                coords = [e_geometry, l_geometry]
                # Flip latitude and longitude positions for find_lane_distance function
                # It is necessary because geodesic expects them in different order
                flipped_latlon_coords = [
                    (e_geometry.y, e_geometry.x),
                    (l_geometry.y, l_geometry.x),
                ]
                intersection_road_distance = find_lane_distance(flipped_latlon_coords)

                # Skip backward turns in the main graph but allow in non-compulsory graph
                if e_node_id[0] == l_node_id[0] and e_node_id[1] == l_node_id[1]:
                    G_with_non_compulsory_edges.add_edge(
                        e_node_id,
                        l_node_id,
                        geometry=LineString(coords),
                        distance=intersection_road_distance,
                        edge_type="intersection",
                    )
                    continue

                G.add_edge(
                    e_node_id,
                    l_node_id,
                    geometry=LineString(coords),
                    distance=intersection_road_distance,
                    edge_type="intersection",
                )
                G_with_non_compulsory_edges.add_edge(
                    e_node_id,
                    l_node_id,
                    geometry=LineString(coords),
                    distance=intersection_road_distance,
                    edge_type="intersection",
                )

    # Create edges outside intersections (based on lanes)
    # Create an edge between all nodes with the same road and lane ID
    observables = {}
    for id_tuple, node_data in G.nodes(data=True):
        _, road_id, lane_id = id_tuple
        observable_tuple = (road_id, lane_id)

        # We assume there are only two nodes per road in this context
        if observable_tuple in observables:
            node_one_id_tuple = observables.pop(observable_tuple)[0]

            # Retrieve the geometry of the lane
            lane_geometry = next(
                (l["geometry"] for l in lanes[road_id] if l["lane_id"] == lane_id), None
            )
            coords = get_twonodes_average_coords(G, id_tuple, node_one_id_tuple)
            distance = find_lane_distance(coords)

            # Add edge depending on entering status
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
            observables[observable_tuple] = (id_tuple, node_data)

    return G, G_with_non_compulsory_edges


def build_city_graph(lane_geojson_file, intersection_geojson_file, osm_file):
    """
    Build a road network graph from GeoJSON and OSM files.

    Args:
        lane_geojson_file (str): Path to the GeoJSON file containing lane data.
        intersection_geojson_file (str): Path to the GeoJSON file containing intersection data.
        osm_file (str): Path to the OSM XML file.

    Returns:
        tuple: A tuple containing:
            - nx.DiGraph: Directed graph with compulsory edges.
            - nx.DiGraph: Directed graph with non-compulsory edges.
    """
    # Load and parse the OSM XML file
    osm_dict = load_osm_to_dict(osm_file)
    print("OSM XML LOADED")

    # Load lane and intersection GeoJSON files
    lanes_geojson = load_geojson(lane_geojson_file)
    print("LANE GEOJSON LOADED")
    intersections_geojson = load_geojson(intersection_geojson_file)
    print("INTERSECTIONS GEOJSON LOADED")

    # Preprocess lanes and intersections
    lanes = preprocess_lanes(lanes_geojson, osm_dict)
    print("LANES PREPROCESSED")
    intersections = preprocess_intersections(intersections_geojson, lanes)
    print("INTERSECTIONS PREPROCESSED")

    # Create the directed graphs
    G, G_with_non_compulsory_edges = create_graph(lanes, intersections)

    return G, G_with_non_compulsory_edges
