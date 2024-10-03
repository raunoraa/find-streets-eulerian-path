import networkx as nx
import geojson
from shapely.geometry import shape
from collections import defaultdict

# Load GeoJSON files
def load_geojson(file_path):
    with open(file_path, 'r') as f:
        return geojson.load(f)

# Function to check if lane is car-drivable
def is_car_drivable(lane):
    # Check for 'backward' access permission if present
    backward_access = lane['properties']['muv']['travel'].get('backward', {}).get('access', {}).get('LandBased', [{}])[0].get('value', 'No')
    
    # Check for 'forward' access permission if present
    forward_access = lane['properties']['muv']['travel'].get('forward', {}).get('access', {}).get('LandBased', [{}])[0].get('value', 'No')
    
    # Return True if either access is 'Yes'
    return backward_access == 'Yes' or forward_access == 'Yes'

# Load and filter car-drivable lanes
# Returns a defaultdict, where keys are road_ids and values are the lanes (and their respective data) 
#  each road has
def parse_lanes(lane_geojson):
    lanes_by_road_id = defaultdict(list)  # A dictionary where road_id is the key, and values are lists of lanes
    for lane in lane_geojson['features']:
        if is_car_drivable(lane):
            lane_data = {
                'road_id': lane['properties']['road'], 
                'lane_id': lane['properties']['index'],
                'direction': lane['properties']['direction'], 
                'src_i': None, # value will be added in parse_intersections
                'dst_i': None, # value will be added in parse_intersections
                'geometry': shape(lane['geometry'])
            }
            lanes_by_road_id[lane_data['road_id']].append(lane_data)  # Group lanes by road_id
    return lanes_by_road_id

# Load intersections and their movements
# Also assign src and destination intersections to lanes (lanes parameter will be modified)
def parse_intersections(intersection_geojson, all_lanes):

    intersections = set()

    for feature in intersection_geojson['features']:
        i_type = feature['properties']['type']
        if i_type == 'road':
            # If the intersection type is road, then we can assign the src and dst intersections to lanes
            lanes = all_lanes.get(feature['properties']['id'])
            for lane in lanes:
                lane['src_i'] = feature['properties']['src_i']
                lane['dst_i'] = feature['properties']['dst_i']
            
        elif i_type == 'intersection':
            intersections.add(
                {
                    'id': feature['properties']['id'],
                    'movements': feature['properties']['movements'],
                    'geometry': shape(feature['geometry'])
                }
            )
        else:
            # Let's not look at other types
            continue

    return intersections

# Create the directed multigraph
def create_graph(lanes, intersections):

    # TODO

    G = nx.MultiDiGraph()

    # Add nodes (intersections)
    for intersection_id, data in intersections.items():
        G.add_node(intersection_id, geometry=data['geometry'])

    # Add edges (lanes)
    for lane in lanes:
        src = lane['src_intersection']
        dst = lane['dst_intersection']
        G.add_edge(src, dst, lane_id=lane['lane_id'], geometry=lane['geometry'])

    # Add movements (intersection-to-intersection connections)
    for intersection_id, data in intersections.items():
        movements = data['movements']
        for movement in movements:
            src_road, dst_road = movement.split(" -> ")
            src_id = int(src_road.split('#')[-1])
            dst_id = int(dst_road.split('#')[-1])
            if src_id in intersections and dst_id in intersections:
                G.add_edge(src_id, dst_id)

    return G

# Main function to construct the graph from geojson files
def build_city_graph(lane_geojson_file, intersection_geojson_file):
    # Load and parse the geojson files
    lanes_geojson = load_geojson(lane_geojson_file)
    intersections_geojson = load_geojson(intersection_geojson_file)

    lanes = parse_lanes(lanes_geojson)
    intersections = parse_intersections(intersections_geojson, lanes)



    # Create the directed multigraph
    G = create_graph(lanes, intersections)
    return G

# Example usage
lane_geojson_file = 'Lane_polygons.geojson'
intersection_geojson_file = 'Intersection_polygons.geojson'
city_graph = build_city_graph(lane_geojson_file, intersection_geojson_file)

# Now you have a graph `city_graph` with lanes as edges and intersections as nodes
