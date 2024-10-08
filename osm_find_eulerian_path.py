import networkx as nx
import geojson
from Intersection import Intersection
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
                #'direction': lane['properties']['direction'], 
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
                        
            lanes = all_lanes[feature['properties']['id']]
                        
            for lane in lanes:
                lane['src_i'] = feature['properties']['src_i']
                lane['dst_i'] = feature['properties']['dst_i']
            
        elif i_type == 'intersection':
            intersection = Intersection(feature['properties']['id'], feature['properties']['movements'], shape(feature['geometry']))
            intersections.add(intersection)
        else:
            # Let's not look at other types at the moment
            continue

    return intersections

# Create the directed multigraph
def create_graph(lanes, intersections):

    G = nx.MultiDiGraph()

    # Added nodes, which are grouped by their intersection
    # key=intersection id; value=set of nodes in an intersection
    # a node is a 3-tuple
    intersections_node_map = {}

    # We want to add intersections in a way
    # that there will be created one node per lane
    for intersection in intersections:
        id = intersection.id
        movements = intersection.movements
        unique_roads = set() # a set of unique road ids in an intersection
        for movement in movements:
            src_road, _ = movement.split(" -> ")

            # There will always be 1-1 ratio of roads in intersections
            # so we can take unique ones only from source or destination
            src_id = int(src_road.split('#')[-1])            
            unique_roads.add(src_id)        
        
        tuples = set()
        for ur in unique_roads:
            lanes_in_ur = lanes.get(ur)
            for lane in lanes_in_ur:
                # Just for now, assign the same geometry to each intersection node
                # have to change later (it's just visual distortion, wont affect the graph)

                node_id = (id, lane['road_id'], lane['lane_id'])                

                # Also add the info, if the node for the current lane is for entering or leaving the intersection
                is_entering_value = False                
                if lane['dst_i'] == id:
                    is_entering_value = True                

                tuples.add((node_id, is_entering_value, intersection.geometry))
                G.add_node(node_id, is_entering=is_entering_value, geometry=intersection.geometry)
        
        intersections_node_map[id] = tuples
    
    # Create edges inside the intersection
    for _, nodes in intersections_node_map.items():        
        entering_nodes = []
        leaving_nodes = []
        for node in nodes:
            _, is_entering, _ = node
            if is_entering:
                entering_nodes.append(node)
            else:
                leaving_nodes.append(node)
            
        if len(entering_nodes) > 0 and len(leaving_nodes) > 0:
            for e_node in entering_nodes:
                e_node_id, _, e_geometry = e_node                
                for l_node in leaving_nodes:
                    l_node_id, _, l_geometry = l_node

                    # Need to tweak the geometry, but can be done later (TODO)
                    # after the graph seems ok.
                    # This is just placeholder geometry for now.
                    G.add_edge(e_node_id, l_node_id, geometry=e_geometry)
    

    # Create edges outside the intersections (basically add the lanes)
    # We need to create an edge between all such nodes, which have the same road and lane id
    observables = {}
    for node in G.nodes(data=True):
        #print(node)
        id_tuple, _ = node        
        #print(id_tuple)
        _, road_id, lane_id = id_tuple
        observable_tuple = (road_id, lane_id)

        # We can assume that there are only 2 nodes per each road
        if observable_tuple in observables.keys():
            node_one_id_tuple = observables.pop(observable_tuple)[0]

            lane_geometry = None
            # Get the lane geometry from the lanes defaultdict
            for l in lanes[road_id]:
                if l.get('lane_id') == lane_id:
                    lane_geometry = l.get('geometry')
                    break
                        
            G.add_edge(id_tuple, node_one_id_tuple, geometry = lane_geometry)
        else:
            observables[observable_tuple] = node

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

folder_path = 'map_files/observable_geojson_files/'

lane_geojson_file = folder_path + 'Lane_polygons.geojson'
intersection_geojson_file = folder_path + 'Intersection_polygons.geojson'

G = build_city_graph(lane_geojson_file, intersection_geojson_file)
#eulerian_path = nx.has_eulerian_path(G)
#print(nx.is_weakly_connected(G))

#for edge, data in G.nodes(data=True):
    #print(data)



import folium

def get_boundaries(geojson_boundaries):
    lats = []
    lons = []

    coordinates = geojson_boundaries['geometry']['coordinates'][0]    
    for coordinate in coordinates:        
        lats.append(coordinate[0])
        lons.append(coordinate[1])
    
    return (min(lats), min(lons), max(lats), max(lons))

# Define map boundaries (min_latitude, min_longitude, max_latitude, max_longitude)
# Take this data from the Boundary.geojson file
boundary_geojson_file = folder_path + 'Boundary.geojson'
loaded_boundary_geojson = load_geojson(boundary_geojson_file)
map_boundaries = get_boundaries(loaded_boundary_geojson)

# Create a Folium map centered on the middle of the boundaries
map_center = [(map_boundaries[1] + map_boundaries[3]) / 2, 
              (map_boundaries[0] + map_boundaries[2]) / 2]

# Initialize the Folium map
m = folium.Map(location=map_center, zoom_start=13)

# Add nodes to the map (assuming nodes are polygons)
for node, data in G.nodes(data=True):
    geometry = data['geometry']
    if geometry:
        # Get the coordinates of the polygon
        coordinates = list(geometry.exterior.coords)  # Taking the first polygon's coordinates
        folium.Polygon(
            locations=[(coord[1], coord[0]) for coord in coordinates],
            color='blue',
            fill=True,
            fill_color='blue',
            fill_opacity=0.5,
            popup=f'Node {node}'
        ).add_to(m)

# Add edges to the map (assuming edges are also polygons)
for u, v, data in G.edges(data=True):
    geometry = data['geometry']
    if geometry:
        coordinates = list(geometry.exterior.coords)  # Taking the first polygon's coordinates
        folium.Polygon(
            locations=[(coord[1], coord[0]) for coord in coordinates],
            color='black',
            #weight=2,
            #opacity=0.7,
            fill=True,
            fill_color='red',
            fill_opacity=0.5,
            popup=f'Edge from {u} to {v}'
        ).add_to(m)

# Save or display the map
m.save('city_graph_map.html')
