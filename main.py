import time

from config import (
    LANES_GEOJSON_FILE_PATH,
    INTERSECTIONS_GEOJSON_FILE_PATH,
    OSM_XML_FILE_PATH,
    OUTPUT_GPKG_FILE_PATH,
    OUTPUT_TXT_FILE_PATH,
)

from file_io.loader import check_file_exists
from file_io.saver import save_graph_to_geopackage, save_results_as_txt

from graph.calculations import calculate_path_distance_and_coefficient
from graph.builder import build_city_graph
from graph.data_preprocessor import (
    largest_strongly_connected_component,
    reduce_intersection_edges,
)

from graph.eulerian_path_solver import find_eulerian_path


def main():
    """
    Main entry point for the program. This function orchestrates the entire process of loading data,
    building the graphs, finding the Eulerian path, and saving the results to files.

    Steps:
    1. Check if the required files exist.
    2. Build two graphs from the input files.
    3. Reduce the graphs to their largest strongly connected components.
    4. Remove as many edges inside the intersections as possible for one graph.
    5. Find the Eulerian path and calculate related metrics.
    6. Save results (distances, coefficient, time spent) to a text file.
    7. Save the graph with the Eulerian path to a GeoPackage file for QGIS.
    """
    # Start timer to measure execution time
    timer = time.time()

    print("\nINIT\n")

    # Check if required files exist (OSM XML, Lane GeoJSON, and Intersection GeoJSON)
    if not (
        check_file_exists(OSM_XML_FILE_PATH)
        and check_file_exists(LANES_GEOJSON_FILE_PATH)
        and check_file_exists(INTERSECTIONS_GEOJSON_FILE_PATH)
    ):
        # If any file is missing, print an error message and exit
        print("Error: One or more required files are missing.")
        return

    # Build the city graph from input files (lanes, intersections, and OSM data)
    G, G_with_non_compulsory_edges = build_city_graph(
        LANES_GEOJSON_FILE_PATH, INTERSECTIONS_GEOJSON_FILE_PATH, OSM_XML_FILE_PATH
    )
    print("GRAPH BUILDING FINISHED")

    # Reduce the graph to its largest strongly connected component (only the main connected part of the graph)
    G = largest_strongly_connected_component(G)
    G_with_non_compulsory_edges = largest_strongly_connected_component(
        G_with_non_compulsory_edges
    )

    # Remove as many edges inside the intersections as possible for one graph.
    reduce_intersection_edges(G)
    print("GRAPH MODIFYING FINISHED")

    # Find the Eulerian path in the graph
    # G will get new edges from the graph with optional edges, if it makes
    # the Eulerian path shorter.
    total_streets_length, eulerian_path = find_eulerian_path(
        G, G_with_non_compulsory_edges
    )

    # Calculate execution time (in minutes)
    time_spent = round((time.time() - timer) / 60, 4)

    # Calculate the Eulerian path distance and coefficient (comparison between path length and total street length)
    path_d, coef = calculate_path_distance_and_coefficient(
        eulerian_path, G, total_streets_length
    )

    # Save the Eulerian path distance, coefficient, and time spent to a text file
    save_results_as_txt(
        total_streets_length, path_d, coef, time_spent, OUTPUT_TXT_FILE_PATH
    )

    # Save the graph with the Eulerian path to a GeoPackage file, which can be viewed and animated in QGIS
    save_graph_to_geopackage(
        eulerian_path,
        G,
        OUTPUT_GPKG_FILE_PATH,
    )

    print("PROCESS COMPLETED\n")


if __name__ == "__main__":
    main()
