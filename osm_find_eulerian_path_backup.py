import networkx as nx
from shapely.geometry import shape, Polygon, Point
import geojson

# Doesnt do anything useful yet 

# Task list:
# 1. Build a multidigraph graph of all drivable lanes. TODO
# 2. Find eulerian path in the graph. TODO

geojson_folder_location = "map_files/observable_geojson_files/"

# Load the geojson files
with open(geojson_folder_location + 'Boundary.geojson') as f:
    boundary_data = geojson.load(f)

with open(geojson_folder_location + 'Intersection_polygons.geojson') as f:
    intersection_data = geojson.load(f)

with open(geojson_folder_location + 'Lane_polygons.geojson') as f:
    lane_data = geojson.load(f)

# Parse boundary polygon
boundary_coords = boundary_data['geometry']['coordinates'][0]
boundary_polygon = Polygon(boundary_coords)

# Create the graph
G = nx.MultiDiGraph()

# Helper function to check if a lane is car-drivable
def is_drivable(lane_properties):
    try:
        muv = lane_properties['muv']
        forward_drive = muv['travel']['forward']['access']['LandBased'][0]['value'] == 'Yes'
        backward_drive = muv['travel']['backward']['access']['LandBased'][0]['value'] == 'Yes'
        return forward_drive or backward_drive
    except KeyError:
        return False

# Step 1: Add intersections as nodes (only features with type "intersection")
for feature in intersection_data['features']:
    if feature['properties']['type'] == 'intersection':
        intersection_id = feature['properties']['id']
        geom = shape(feature['geometry'])
        centroid = geom.centroid
        
        # Adding the intersection node
        G.add_node(intersection_id, 
                   pos=(centroid.x, centroid.y), 
                   geometry=geom, 
                   movements=feature['properties']['movements'],
                   osm_node_ids=feature['properties']['osm_node_ids'])
        print(f"Added intersection node: {intersection_id} with movements: {feature['properties']['movements']}")

# Step 2: Add lanes as edges
for feature in lane_data['features']:
    lane_properties = feature['properties']
    
    # Only include drivable lanes
    if not is_drivable(lane_properties):
        print()
        print(f"Skipped non-drivable lane with properties: {lane_properties}")
        print()
        continue
    
    # Extract lane geometry
    geom = shape(feature['geometry'])
    
    # Determine the intersections associated with this lane
    from_intersection = lane_properties['road']  # Placeholder, adjust based on your data
    to_intersection = lane_properties['road'] + 1  # Placeholder, adjust based on your data
    
    print(f"Processing lane from {from_intersection} to {to_intersection} with properties: {lane_properties}")

    # Check if the lane can be added as an edge based on allowed movements
    if from_intersection in G and to_intersection in G:
        allowed_movements = G.nodes[from_intersection]['movements']
        if f"Road #{from_intersection} -> Road #{to_intersection}" in allowed_movements:
            G.add_edge(from_intersection, to_intersection, key=lane_properties['index'], geometry=geom, properties=lane_properties)
            print(f"Added edge from {from_intersection} to {to_intersection}")
    
    # Check if the lane ends at the boundary and add a reverse lane for one-ways
    lane_end_point = Point(geom.boundary[1])  # Adjust if necessary
    if boundary_polygon.contains(lane_end_point):
        if lane_properties['muv']['travel']['forward']['access']['LandBased'][0]['value'] == 'Yes' and \
           lane_properties['muv']['travel']['backward']['access']['LandBased'][0]['value'] == 'No':
            reverse_lane_id = f"{lane_properties['index']}_reverse"
            G.add_edge(to_intersection, from_intersection, key=reverse_lane_id, geometry=geom, properties={'artificial_reverse': True})
            print(f"Added reverse edge from {to_intersection} to {from_intersection}")

# Output the final graph details
print(f"Graph created with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")
