r""" 
This submodule is for spatial calculations (distances, geometry handling).
"""

from shapely.geometry import Point

from geopy.distance import geodesic


def find_lane_distance(node_coords):
    """
    Calculate the geodesic distance between two nodes (representing the start and end of a road).

    Args:
        node_coords (tuple): A tuple containing two coordinates, where each coordinate is a tuple (latitude, longitude).
            Example: [(lat1, lon1), (lat2, lon2)].

    Returns:
        float: The geodesic distance between the two points in meters.

    Example:

    ::

        distance = find_lane_distance([(52.2296756, 21.0122287), (41.8919300, 12.5113300)])
    """
    start_node = node_coords[0]
    next_node = node_coords[1]
    return geodesic(start_node, next_node).meters


def get_twonodes_average_coords(G, start_node, end_node):
    """
    Calculate the average coordinates of two nodes in the graph.

    Args:
        G (networkx.Graph): The graph object containing node data.
        start_node (str): The identifier of the start node.
        end_node (str): The identifier of the end node.

    Returns:
        list: A list containing two tuples, each representing the average (latitude, longitude) coordinates of the start and end nodes.
            Example: [(lat1, lon1), (lat2, lon2)].

    Example:

    ::

        avg_coords = get_twonodes_average_coords(G, "start_node_id", "end_node_id")
    """
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
    """
    Find the closest midpoint between the edges of two polygons.

    Args:
        polygon1 (shapely.geometry.Polygon): The first polygon.
        polygon2 (shapely.geometry.Polygon): The second polygon.

    Returns:
        shapely.geometry.Point: The closest point between the two polygons.

    Example:

    ::

        closest_point = get_closest_edge_midpoint(polygon1, polygon2)
    """
    min_distance = float("inf")

    # Find the closest pair of vertices between polygon1 and polygon2 using geodesic distance
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
