# Configuration file for file paths and constants used throughout the application.

# Paths, from where to load the input files (recommended to leave the folder names as they are).
# Folder "map_files" already exists, only the folder itself and the md file inside it are tracked by Git.
LANES_GEOJSON_FILE_PATH = "map_files/Lane polygons.geojson"
INTERSECTIONS_GEOJSON_FILE_PATH = "map_files/Intersection polygons.geojson"
OSM_XML_FILE_PATH = "map_files/osm.xml"

# Paths, where to save the output files (recommended to leave the folder names as they are).
# Folders are created automatically if they don't exist yet. The folder itself and the files inside are not tracked by Git.
OUTPUT_GPKG_FILE_PATH = "output/output.gpkg"
OUTPUT_TXT_FILE_PATH = "output/results.txt"
