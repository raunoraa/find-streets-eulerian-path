r"""
This submodule preprocesses the data to make it usable for further processing.
"""

import networkx as nx
from shapely.geometry import shape
from collections import defaultdict

from classes.intersection import Intersection


def is_car_drivable(lane):
    """
    Determine if a lane is drivable based on its type.

    This function checks the 'type' property within the 'properties' section of a GeoJSON lane feature.
    If the type is 'Driving', the lane is considered car-drivable.

    Args:
        lane (dict): A dictionary representing a lane feature from GeoJSON.

    Returns:
        bool: True if the lane is drivable, False otherwise.

    Example:

    ::

        drivable = is_car_drivable(lane)
    """
    return lane.get("properties", {}).get("type", "").lower() == "driving"


def preprocess_lanes(lane_geojson, osm_dict):
    """
    Load and filter car-drivable lanes based on GeoJSON features and OSM tags.

    Args:
        lane_geojson (dict): GeoJSON data containing lane features.
        osm_dict (dict): OSM data as a dictionary mapping OSM way IDs to their properties.

    Returns:
        defaultdict: A dictionary where road IDs are keys, and values are lists of filtered lanes.

    Example:

    ::

        lanes_by_road_id = preprocess_lanes(lane_geojson, osm_dict)
    """
    lanes_by_road_id = defaultdict(list)

    allowed_tags = {
        "motorway",
        "motorway_link",
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
        "living_street",
        "service",
        "track",
    }

    for lane in lane_geojson["features"]:
        if is_car_drivable(lane):
            osm_way_ids = lane["properties"]["osm_way_ids"]
            allowed_tag_present = False
            for osm_way_id in osm_way_ids:
                highway_tag_value = osm_dict[str(osm_way_id)]["highway"]
                if highway_tag_value in allowed_tags:
                    allowed_tag_present = True
                    break

            if allowed_tag_present:
                lane_data = {
                    "road_id": lane["properties"]["road"],
                    "lane_id": lane["properties"]["index"],
                    "direction": lane["properties"]["direction"],
                    "src_i": None,  # value will be added in preprocess_intersections
                    "dst_i": None,  # value will be added in preprocess_intersections
                    "geometry": shape(lane["geometry"]),
                }
                lanes_by_road_id[lane_data["road_id"]].append(lane_data)

    return lanes_by_road_id


def preprocess_intersections(intersection_geojson, all_lanes):
    """
    Load intersections and their movements, assigning source and destination intersections to lanes.

    Args:
        intersection_geojson (dict): GeoJSON data containing intersection features.
        all_lanes (defaultdict): Dictionary of lanes grouped by road ID.

    Returns:
        set: A set of Intersection objects representing valid intersections.

    Example:

    ::

        intersections = preprocess_intersections(intersection_geojson, lanes_by_road_id)
    """
    intersections = set()

    for feature in intersection_geojson["features"]:
        i_type = feature["properties"]["type"]

        if i_type == "road":
            # Assign source and destination intersections to lanes
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
            # Ignore other intersection types
            continue

    return intersections


def largest_strongly_connected_component(G):
    """
    Find and return the largest strongly connected component in the given graph.

    Args:
        G (networkx.DiGraph): The directed graph to analyze.

    Returns:
        networkx.DiGraph: A subgraph containing the largest strongly connected component.

    Example:

    ::

        largest_component = largest_strongly_connected_component(G)
    """
    return G.subgraph(max(nx.strongly_connected_components(G), key=len)).copy()


def reduce_intersection_edges(G):
    """
    Remove unnecessary edges inside intersections while maintaining strong connectivity.

    Args:
        G (networkx.DiGraph): The directed graph representing the road network.

    Returns:
        None: The function modifies the graph in place.

    Example:

    ::

        reduce_intersection_edges(G)
    """
    edges = list(G.edges(data=True))

    for u, v, data in edges:
        if data.get("edge_type") == "intersection":
            if not (G.out_degree(u) <= 1 or G.in_degree(v) <= 1):
                G.remove_edge(u, v)

                # Ensure strong connectivity is preserved
                if not nx.has_path(G, u, v):
                    G.add_edge(u, v, **data)
