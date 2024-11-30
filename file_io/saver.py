r"""
This submodule deals with saving the files.
"""

import os

from geo_operations.qgis_utils import (
    create_edge_features_for_qgis,
    create_node_features_for_qgis,
)


def save_results_as_txt(
    total_streets_length,
    eulerian_path_distance,
    coefficient,
    time_spent,
    results_txt_file_path="results.txt",
    verbose=True,
):
    """
    Save the program results to a text file and optionally print them to the console.

    This function writes the total street length, Eulerian path distance, and coefficient to a text file.
    If `verbose` is True, the results are also printed to the console.

    Args:
        total_streets_length (float): The total length of the street lanes in meters.
        eulerian_path_distance (float): The total distance of the Eulerian path in meters.
        coefficient (float): The coefficient comparing the Eulerian path distance to the total street length.
        time_spent (float): The time spent on program execution (in minutes).
        results_txt_file_path (str, optional): The path to the output text file where results will be saved.
            Default is "results.txt".
        verbose (bool, optional): If True, prints the results to the console. Default is True.

    Returns:
        None

    Example:

    ::

        save_results_as_txt(500.0, 1000.0, 2.0, 1.1, "results/output_results.txt")
    """
    time_spent_string = f"Program execution completed in {time_spent} minutes."
    distances_string = f"Street lanes distance: {round(total_streets_length, 4)}m ;; Eulerian path distance: {round(eulerian_path_distance, 4)}m"
    coefficient_string = f"Coefficient: {round(coefficient, 4)}"

    # Ensure the folders exist before writing (create them if they don't exist yet)
    os.makedirs(os.path.dirname(results_txt_file_path), exist_ok=True)

    # Write the results to a text file
    with open(results_txt_file_path, "w", encoding="utf-8") as results_file:
        results_file.write(time_spent_string + "\n")
        results_file.write(distances_string + "\n")
        results_file.write(coefficient_string + "\n")

    if verbose:
        print()
        print(time_spent_string)
        print(distances_string)
        print(coefficient_string)
        print()
        print(f"Results saved to {results_txt_file_path}")


def save_graph_to_geopackage(path, G, output_file_path="output.gpkg", verbose=True):
    """
    Save the graph's nodes and edges to a GeoPackage file.

    This function generates GeoDataFrames for both the nodes and edges of the graph. These are then
    stored as separate layers in a GeoPackage file. The edges are ordered based on the Eulerian path
    provided in the `path` argument.

    Args:
        path (list of tuple): The Eulerian path as a list of edges (tuples of nodes) that represents
            the order in which the edges are traversed.
        G (networkx.Graph): The graph object containing the nodes and edges.
        output_file_path (str, optional): The path to the output GeoPackage file. Default is "output.gpkg".
        verbose (bool, optional): If True, prints the status of the saving process to the console. Default is True.

    Returns:
        None

    Example:

    ::

        save_graph_to_geopackage(path, G, "output/graph_output.gpkg")
    """
    # Create the GeoDataFrame for edges with timestamps
    edges_gdf = create_edge_features_for_qgis(path, G)

    # Create the GeoDataFrame for nodes
    nodes_gdf = create_node_features_for_qgis(G)

    # Ensure the folders exist before writing (create them if they don't exist yet)
    os.makedirs(os.path.dirname(output_file_path), exist_ok=True)

    # Save both nodes and edges to a GeoPackage
    edges_gdf.to_file(output_file_path, layer="edges", driver="GPKG")
    nodes_gdf.to_file(output_file_path, layer="nodes", driver="GPKG")

    if verbose:
        print(f"Graph saved to {output_file_path}")
        print()
