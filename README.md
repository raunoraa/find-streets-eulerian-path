# Find all streets' lanes' shortest path

This repository contains unfinished solution for finding all streets' lanes' shortest path:

- **`osm_find_shortest_path.py`** contains a solution, where all street lanes' graph is constructed based on the Intersection_polygons.geojson file and Lane_polygons.geojson file. Graph is created with the networkx library as a networkx object. Graph creation and visualization (an html file is created with folium) seems to be working correctly. It will also now generate the solution with the shortest path and visualize it. 

- **TODO** The visualization could be improved a lot (there are some bugs currently with the sliders and it isn't totally clear from the visualization if the shortest path finding algorithm is working correctly, also some of the geometry stuff is unimplemented and can be improved as well). The data for the total path length in kilometres should also be collected and shown.


## How to run

Here are the instructions, how to run and test the current solution:

- You need to clone this repository first, then open the project with visual studio code, then download the extension called **Live Server**.

- **`osm_find_shortest_path.py`**: If you plan to change anything in the code, you need to run this file before you do any of the following to generate new (updated) html files.
    - When you want to see just the graph, open the file in vsc editor called `city_graph_map.html` and then run live server.
    - When you want to see the shortest path through all of the streets' lanes', open the file in vsc editor called  `animated_path.html` and then run live server.