r""" 
This submodule is for QGIS export utilities.
"""

from shapely.geometry import LineString
from datetime import datetime, timedelta
import geopandas as gpd


def assign_edge_orders_with_multiple_traversals(path, G):
    """
    Assign timestamps to edges based on their traversal in the Eulerian path.
    Each edge may be traversed multiple times, so the 'order' attribute will be a list of timestamps.
    The order number is just an index (starts from 1).

    Args:
        path (list): The Eulerian path as a list of edges (tuples of nodes).
        G (networkx.Graph): The graph object containing the edges.

    Returns:
        None

    Example:

    ::

        assign_edge_orders_with_multiple_traversals(path, G)
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


def create_node_features_for_qgis(G):
    """
    Create a list of node features for the graph. Each node is represented by a point geometry.

    Args:
        G (networkx.Graph): The graph object containing the node data.

    Returns:
        geopandas.GeoDataFrame: A GeoDataFrame containing the nodes, each represented by a point geometry.

    Example:

    ::

        nodes_gdf = create_node_features_for_qgis(G)
    """
    node_data = []

    for node, node_attrs in G.nodes(data=True):
        node_data.append({"id": str(node), "geometry": node_attrs.get("geometry")})

    # Create a GeoDataFrame for the nodes
    nodes_gdf = gpd.GeoDataFrame(node_data, geometry="geometry", crs="EPSG:4326")

    return nodes_gdf


def create_edge_features_for_qgis(path, G):
    """
    Create a list of features for edges, where each edge is duplicated for each timestamp in the 'order' list.
    Each traversal of the same edge gets its own timestamp and geometry.

    Args:
        path (list): The Eulerian path as a list of edges (tuples of nodes).
        G (networkx.Graph): The graph object containing the edges and their associated attributes.

    Returns:
        geopandas.GeoDataFrame: A GeoDataFrame containing the edges, each with a timestamp and geometry.

    Example:

    ::

        edges_gdf = create_edge_features_for_qgis(path, G)
    """
    edge_data = []  # List to hold the features (edges)

    # Assign timestamps to edges based on the Eulerian path
    assign_edge_orders_with_multiple_traversals(path, G)

    # For each edge, create a feature for each traversal (each timestamp)
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

    # Create a GeoDataFrame for the edges
    edges_gdf = gpd.GeoDataFrame(edge_data, geometry="geometry", crs="EPSG:4326")

    return edges_gdf
