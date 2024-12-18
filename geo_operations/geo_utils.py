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


def flipped_node_coords(G, node):
    """
    Flip the coordinates of the node in the graph.

    Args:
        G (networkx.Graph): The graph object containing node data.
        node (tuple): The identifier of the node.

    Returns:
        tuple: A tuple, which represents the flipped coordinates of the node (latitude, longitude) or (longitude, latitude), depending on how the input was given.        

    Example:

    ::

        flipped_coords = get_twonodes_average_coords(G, node)
    """
    node_data = G.nodes[node]
    coords = list(node_data.get("geometry").coords)[0]
    flipped_coords = [coords[1], coords[0]]
    return flipped_coords


def get_closest_polygon_point(polygon1, polygon2):
    """
    Find the closest polygon point between the edges of two polygons.

    Args:
        polygon1 (shapely.geometry.Polygon): The first polygon.
        polygon2 (shapely.geometry.Polygon): The second polygon.

    Returns:
        shapely.geometry.Point: The closest point between the two polygons. First polygon's point will be returned.

    Example:

    ::

        closest_point = get_closest_edge_midpoint(polygon1, polygon2)
    """
    min_distance = float("inf")
    closest_point = None

    # Find the closest pair of vertices between polygon1 and polygon2 using geodesic distance
    for coord1 in polygon1.exterior.coords:
        for coord2 in polygon2.exterior.coords:
            # Calculate the geodesic distance between coord1 and coord2
            distance = geodesic(
                coord1[::-1], coord2[::-1]
            ).meters  # Reverse coord1/coord2 for (lat, lon)

            # Update the closest point if this pair is closer
            if distance < min_distance:
                min_distance = distance
                closest_point = Point(coord1)

    return closest_point
