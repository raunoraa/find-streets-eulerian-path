r"""
This submodule deals with graph and pathfinding related calculations.
"""

import os
import numpy as np
import networkx as nx
from concurrent.futures import ProcessPoolExecutor, as_completed


def calculate_surplus_and_deficit(graph):
    """
    Calculate surplus and deficit nodes based on their in-degree and out-degree.

    Args:
        graph (networkx.DiGraph): The directed graph representing the road network.

    Returns:
        tuple: Two dictionaries:
            - surplus (dict): Nodes with surplus edges (out-degree > in-degree).
            - deficit (dict): Nodes with deficit edges (in-degree > out-degree).

    Example:

    ::

        surplus, deficit = calculate_surplus_and_deficit(G)
    """
    surplus = {}
    deficit = {}

    for node in graph.nodes():
        out_degree = graph.out_degree(node)
        in_degree = graph.in_degree(node)

        # Calculate the difference between out-degree and in-degree
        difference = out_degree - in_degree

        # Classify as surplus or deficit based on the difference
        if difference > 0:
            surplus[node] = difference
        elif difference < 0:
            deficit[node] = abs(difference)

    return surplus, deficit


def calculate_path_distance_and_coefficient(path, G, total_streets_length):
    """
    Calculate the total distance of the Eulerian path and the coefficient showing the ratio
    between the Eulerian path distance and the total length of streets.

    Args:
        path (list): List of edges in the Eulerian path.
        G (networkx.DiGraph): The graph representing the road network.
        total_streets_length (float): The total length of all streets in the graph.

    Returns:
        tuple: A tuple containing:
            - path_distance (float): Total distance of the Eulerian path.
            - coefficient (float): Ratio of path distance to total street length.

    Example:

    ::

        path_distance, coefficient = calculate_path_distance_and_coefficient(path, G, total_streets_length)
    """
    path_distance = 0.0

    for edge in path:
        corresponding_edge = G.get_edge_data(edge[0], edge[1])
        distance = corresponding_edge.get("distance")
        path_distance += distance

    coefficient = path_distance / total_streets_length

    return path_distance, coefficient


def process_chunk(deficit_node_chunk, surplus_node_set, G_with_non_compulsory_edges):
    """
    Calculate shortest paths for a chunk of deficit nodes to all surplus nodes.

    Args:
        deficit_node_chunk (list): A list of deficit nodes to process.
        surplus_node_set (set): A set of surplus nodes to calculate distances to.
        G_with_non_compulsory_edges (networkx.DiGraph): The graph containing all edges.

    Returns:
        tuple: Two dictionaries:
            - chunk_distances (dict): Shortest path distances between deficit and surplus nodes.
            - chunk_shortest_paths (dict): Shortest paths between deficit and surplus nodes.

    Example:

    ::

        chunk_distances, chunk_shortest_paths = process_chunk(deficit_chunk, surplus_nodes, G_with_non_compulsory_edges)
    """
    chunk_distances = {}
    chunk_shortest_paths = {}

    for deficit_node in deficit_node_chunk:
        # Compute shortest path for this unique pair only once
        distances_single_source, paths_single_source = nx.single_source_dijkstra(
            G_with_non_compulsory_edges,
            source=deficit_node,
            weight="distance",
        )
        for surplus_node in surplus_node_set:
            chunk_distances[(deficit_node, surplus_node)] = distances_single_source[
                surplus_node
            ]
            chunk_shortest_paths[(deficit_node, surplus_node)] = paths_single_source[
                surplus_node
            ]

    return chunk_distances, chunk_shortest_paths


def calculate_cost_matrix(G_with_non_compulsory_edges, deficit, surplus):
    """
    Calculate the shortest path distances between all unique deficit and surplus nodes.
    Returns a cost matrix and a mapping to the original nodes.

    Args:
        G_with_non_compulsory_edges (networkx.DiGraph): The graph containing all edges.
        deficit (dict): Dictionary of deficit nodes with their respective deficits.
        surplus (dict): Dictionary of surplus nodes with their respective surpluses.

    Returns:
        tuple: A tuple containing:
            - cost_matrix (numpy.ndarray): Cost matrix for matching deficit and surplus nodes.
            - shortest_paths (dict): Shortest paths between deficit and surplus nodes.
            - expanded_deficit_nodes (list): Expanded list of deficit nodes.
            - expanded_surplus_nodes (list): Expanded list of surplus nodes.
            - deficit_mapping (list): Mapping of expanded deficit nodes to original nodes.
            - surplus_mapping (list): Mapping of expanded surplus nodes to original nodes.

    Example:

    ::

        cost_matrix, shortest_paths, expanded_deficit_nodes, expanded_surplus_nodes, deficit_mapping, surplus_mapping = calculate_cost_matrix(G, deficit, surplus)
    """
    distances = {}
    shortest_paths = {}

    print("CALCULATING SHORTEST PATHS")

    num_workers = os.cpu_count()
    deficit_node_list, surplus_node_set = list(deficit.keys()), set(surplus.keys())

    # If the graph is small, avoid parallelization
    if num_workers >= len(deficit_node_list):
        print("Using serial solution")
        for deficit_node in deficit.keys():
            distances_single_source, paths_single_source = nx.single_source_dijkstra(
                G_with_non_compulsory_edges, source=deficit_node, weight="distance"
            )
            for surplus_node in surplus_node_set:
                distances[(deficit_node, surplus_node)] = distances_single_source[
                    surplus_node
                ]
                shortest_paths[(deficit_node, surplus_node)] = paths_single_source[
                    surplus_node
                ]
    else:
        print("Using parallel solution")
        print("CPU count on the current system:", num_workers)
        chunk_size = len(deficit_node_list) // num_workers + 1
        chunks = [
            deficit_node_list[i : i + chunk_size]
            for i in range(0, len(deficit_node_list), chunk_size)
        ]

        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = {
                executor.submit(
                    process_chunk,
                    chunk,
                    surplus_node_set,
                    G_with_non_compulsory_edges,
                ): chunk
                for chunk in chunks
            }

            for future in as_completed(futures):
                chunk_distances, chunk_shortest_paths = future.result()
                distances.update(chunk_distances)
                shortest_paths.update(chunk_shortest_paths)

    print("UNIQUE SHORTEST PATHS CALCULATED")

    expanded_deficit_nodes, expanded_surplus_nodes = [], []
    deficit_mapping, surplus_mapping = [], []

    # Expand based on deficit values
    for node, count in deficit.items():
        expanded_deficit_nodes.extend([node] * count)
        deficit_mapping.extend([node] * count)

    # Expand based on surplus values
    for node, count in surplus.items():
        expanded_surplus_nodes.extend([node] * count)
        surplus_mapping.extend([node] * count)

    # Initialize the cost matrix
    num_deficit = len(expanded_deficit_nodes)
    num_surplus = len(expanded_surplus_nodes)
    cost_matrix = np.zeros((num_deficit, num_surplus))

    for i, deficit_node in enumerate(expanded_deficit_nodes):
        for j, surplus_node in enumerate(expanded_surplus_nodes):
            cost_matrix[i, j] = distances[(deficit_node, surplus_node)]

    print("EXPANDED COST MATRIX CREATED")
    return (
        cost_matrix,
        shortest_paths,
        expanded_deficit_nodes,
        expanded_surplus_nodes,
        deficit_mapping,
        surplus_mapping,
    )
