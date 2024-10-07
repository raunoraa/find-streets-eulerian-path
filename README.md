# Find all street lanes' eulerian path

This repository contains different unfinished solutions for finding all street lanes' eulerian path:
    - **`osm_find_eulerian_path_prev_folium.py`** contains a solution, where all street lanes' eulerian path is found with using only the osm.xml file as input (right now, using osm_observable.xml file as input, graph is constructed with the osmnx library, eulerian path is found with the networkx library) and animating the eulerian path with the folium library. This is not a functional solution mainly because it doesn't take into account different drivable lanes.
    - **`osm_find_eulerian_path.py`** contains a solution, where all street lanes' graph is constructed based on the Intersection_polygons.geojson file and Lane_polygons.geojson file. Graph is created with the networkx library as a networkx object. Graph creation and visualization is seems to be working correctly. However, finding the eulerian path (or the shortest path to visit all the lanes at least once) is not yet implemented, because something needs to be done about the situations where the eulerian path does not exist in the given graph.


## How to run

Here are the instructions, how to run and test the current solutions:
    - For **both** solutions, you need to clone this repository first, then open the project with visual studio code, then download the extension called **Live Server**.
    - **`osm_find_eulerian_path_prev_folium.py`**: Open the file in vsc editor called `eulerian_path_animation.html` and then run live server, which will open the html page on your web browser.
    - **`osm_find_eulerian_path.py`**: Open the file in vsc editor called `city_graph_map.html` and then run live server, which will open the html page on your web browser.