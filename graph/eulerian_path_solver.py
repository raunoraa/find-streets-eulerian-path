r"""
This submodule is for balancing the graph (if necessary) and calculating the Eulerian path of the graph.
"""

import networkx as nx
from scipy.optimize import linear_sum_assignment

from graph.calculations import calculate_cost_matrix, calculate_surplus_and_deficit


def balance_graph(G, surplus, deficit, G_with_non_compulsory_edges, initial_graph):
    """
    Balance the graph by adding duplicate edges using the Hungarian algorithm.

    Args:
        G (networkx.MultiDiGraph): The graph to balance by adding duplicate edges.
        surplus (dict): Nodes with surplus edges.
        deficit (dict): Nodes with deficit edges.
        G_with_non_compulsory_edges (networkx.DiGraph): The graph containing optional edges.
        initial_graph (networkx.DiGraph): The original graph to be modified.

    Returns:
        tuple: The last observed surplus node ID.
    """
    surplus_node = None

    # Precompute all shortest paths and create an expanded cost matrix
    (
        cost_matrix,
        shortest_paths,
        expanded_deficit_nodes,
        expanded_surplus_nodes,
        deficit_mapping,
        surplus_mapping,
    ) = calculate_cost_matrix(G_with_non_compulsory_edges, deficit, surplus)

    # Use the scipy's solution to find the optimal assignment
    deficit_indices, surplus_indices = linear_sum_assignment(cost_matrix)

    # Iterate through the optimal matches
    for d_idx, s_idx in zip(deficit_indices, surplus_indices):
        deficit_node = expanded_deficit_nodes[d_idx]
        surplus_node = expanded_surplus_nodes[s_idx]

        # Retrieve the shortest path for the optimal match
        best_path = shortest_paths[(deficit_node, surplus_node)]

        # Apply the selected path to balance the graph
        for i in range(len(best_path) - 1):
            u, v = best_path[i], best_path[i + 1]

            # Add edges to the original and modified graphs if they don't exist
            if not initial_graph.has_edge(u, v):
                edge_attributes = G_with_non_compulsory_edges.get_edge_data(u, v)
                initial_graph.add_edge(u, v, **edge_attributes)
                G.add_edge(u, v, **edge_attributes)
            else:
                # Add a copy of an existing edge to the multi-directed graph.
                edge_attributes = G.get_edge_data(u, v).get(0)
                G.add_edge(u, v, **edge_attributes)

        # Update surplus and deficit nodes based on assignments
        surplus[surplus_mapping[s_idx]] -= 1
        deficit[deficit_mapping[d_idx]] -= 1

        # Remove nodes that are no longer in surplus or deficit
        if surplus[surplus_mapping[s_idx]] == 0:
            del surplus[surplus_mapping[s_idx]]
        if deficit[deficit_mapping[d_idx]] == 0:
            del deficit[deficit_mapping[d_idx]]

        # Early stopping condition when only one surplus and one deficit node remain
        if (len(surplus) == 1 and len(deficit) == 1) and (
            surplus[next(iter(surplus))] == 1 and deficit[next(iter(deficit))] == 1
        ):
            surplus_node = next(iter(surplus))
            break

    return surplus_node


def find_eulerian_path(G, G_with_non_compulsory_edges):
    """
    Find an Eulerian path in the graph. If none exists, balance the graph and retry.

    Args:
        G (networkx.Digraph): The graph to find an Eulerian path in.
        G_with_non_compulsory_edges (networkx.DiGraph): A graph with additional optional edges.

    Returns:
        tuple: A tuple containing:
            - total_streets_length (float): Total length of all edges in the graph.
            - eulerian_path (list): List of edges representing the Eulerian path.

    Example:

    ::

        total_length, eulerian_path = find_eulerian_path(G, G_with_non_compulsory_edges)
    """
    eulerian_path = None
    total_streets_length = 0.0

    try:
        # Attempt to find an Eulerian path in the current graph
        eulerian_path = list(nx.eulerian_path(G))

        # Calculate the total length of all streets in the graph
        for _, _, data in G.edges(data=True):
            total_streets_length += data["distance"]

        print("EULERIAN PATH FOUND WITHOUT BALANCING")

    except:
        # If no Eulerian path is found, the graph needs to be balanced
        print("EULERIAN PATH NOT FOUND WITHOUT BALANCING")
        print("BALANCING THE GRAPH TO FIND THE EULERIAN PATH...")

        # Create a copy of the graph as MultiDiGraph to perform modifications
        copied_graph = nx.MultiDiGraph(G)

        # Calculate nodes with surplus and deficit edges
        surplus, deficit = calculate_surplus_and_deficit(copied_graph)

        # Balance the graph and get the start node for the Eulerian path
        start_node = balance_graph(
            copied_graph, surplus, deficit, G_with_non_compulsory_edges, G
        )

        print("GRAPH BALANCING FINISHED")

        # Recalculate total streets length because balancing might add new unique edges
        for _, _, data in G.edges(data=True):
            total_streets_length += data["distance"]

        # Find the Eulerian path in the balanced graph
        eulerian_path = list(nx.eulerian_path(copied_graph, source=start_node))

        print("EULERIAN PATH FOUND")

    finally:
        return total_streets_length, eulerian_path
