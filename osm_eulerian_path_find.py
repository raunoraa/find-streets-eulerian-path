import osmnx as ox
import networkx as nx
import plotly.graph_objs as go
from copy import deepcopy

def find_graph_boundaries(G):
    """Determine the boundaries of the graph based on node coordinates."""
    lats = [G.nodes[node]['y'] for node in G.nodes]
    lons = [G.nodes[node]['x'] for node in G.nodes]
    
    return {
        'north': max(lats),
        'south': min(lats),
        'east': max(lons),
        'west': min(lons)
    }

def is_boundary_node(node, G, boundary, margin=0.001):
    """Check if a node is close to the map boundary."""
    lat, lon = G.nodes[node]['y'], G.nodes[node]['x']
    return (lat <= boundary['south'] + margin or lat >= boundary['north'] - margin or
            lon <= boundary['west'] + margin or lon >= boundary['east'] - margin)

def handle_boundary_edges(G, boundary):
    """Adjust edges for nodes near the boundary by adding reverse edges if necessary."""
    # Identify nodes close to the boundary
    boundary_nodes = {node for node in G.nodes if is_boundary_node(node, G, boundary)}
    
    # Create a copy of the graph to modify
    H = G.copy()
    
    counter = 0

    # Add reverse edges for one-way streets if they connect to boundary nodes
    for u, v, data in G.edges(data=True):
        #print(counter)
        counter += 1
        if u in boundary_nodes or v in boundary_nodes:
            # Check if the edge is one-way
            if not H.has_edge(v, u):                
                # Add the reverse edge with the same attributes
                H.add_edge(v, u, **data)
    
    return H

def filter_drivable_lanes(G):
    """Filter the graph to include only drivable lanes."""
    # Create a new graph to hold only the drivable lanes
    H = nx.MultiDiGraph()  # or MultiDiGraph if you need to preserve multiple edges
    
    # Define acceptable highway types for drivable lanes
    drivable_highways = {'motorway', 'trunk', 'primary', 'secondary', 'tertiary', 'residential'}
    
    nodes_in_filtered_edges= set()

    # Iterate over the edges and add only the drivable ones
    for u, v, key, data in G.edges(keys=True, data=True):
        # Check if the edge has a 'highway' attribute and if it's in the list of drivable highways        
        if 'highway' in list(data.keys()) and data['highway'] in list(drivable_highways):
            # Add the edge to the new graph with its attributes
            nodes_in_filtered_edges.add(u)
            nodes_in_filtered_edges.add(v)
            H.add_edge(u, v, key=key, **data)            
    
    # Add the nodes to the new graph with their attributes
    for node in nodes_in_filtered_edges:
       H.nodes[node].update(G.nodes[node])

    default_weight = 1.0
    for u, v, key, data in H.edges(keys=True, data=True):
        if 'weight' not in data:
            H[u][v][key]['weight'] = default_weight    

    return H



def visualize_eulerian_path(G, eulerian_path, boundary):
    """Visualize the Eulerian path on the graph using Plotly."""

    frames=[]
    edge_trace = []    

    order = 0

    for u, v in eulerian_path:
        
        order += 1

        # Get node coordinates
        x_start, y_start = G.nodes[u]['x'], G.nodes[u]['y']
        x_end, y_end = G.nodes[v]['x'], G.nodes[v]['y']

        # Get edge data
        edge_data = G[u][v][0]

        #print(order, edge_data)

        distance = edge_data.get("length", 0)
        weight = edge_data.get("weight", 1)
        street = edge_data.get("name", "Unknown Street")

        current_edge_trace = go.Scattermapbox(
            lon=[x_start, x_end],
            lat=[y_start, y_end],
            mode='lines',
            line=dict(width=2, color='blue'),
            hoverinfo='text',
            text=f"Order: {order}<br>Street: {street}<br>Length: {distance*weight:.2f} m"
        )

        # Create edge line trace
        edge_trace.append(current_edge_trace)

        current_edge_trace_list = deepcopy(edge_trace)

        frame = (go.Frame(
            data = current_edge_trace_list, name = str(order)
        ))

    
    # Layout for the plot
    layout = go.Layout(
        title="Eulerian Path Visualization",
        mapbox=dict(
            style="open-street-map",
            zoom=16,
            center=dict(lon=(boundary.get('west') + boundary.get('east'))/2, lat=(boundary.get('north') + boundary.get('south'))/2)
        ),
        showlegend=False,
        updatemenus=[{
            'buttons': [
                {
                    'args': [None, {'frame': {'duration': 500, 'redraw': True}, 'fromcurrent': True}],
                    'label': 'Play',
                    'method': 'animate'
                },
                {
                    'args': [[None], {'frame': {'duration': 0, 'redraw': True}, 'mode': 'immediate', 'transition': {'duration': 0}}],
                    'label': 'Pause',
                    'method': 'animate'
                }
            ],
            'direction': 'left',
            'pad': {'r': 10, 't': 87},
            'showactive': True,
            'type': 'buttons',
            'x': 0.1,
            'xanchor': 'right',
            'y': 0,
            'yanchor': 'top'
        }]
    )

    #fig = go.Figure(data=edge_trace, layout=layout)
    fig = go.Figure(data=edge_trace, layout=layout, frames=frames)
    fig.show()



def main(xml_file_path, verbose=True):
    # Load the graph from the OSM XML file
    G_generated = ox.graph_from_xml(filepath=xml_file_path)

    # Set CRS if not present and reproject the graph
    #if not G_generated.graph.get('crs'):
     #   G_generated.graph['crs'] = 'epsg:4326'  # WGS84 projection

    G = filter_drivable_lanes(G_generated)
    
    # Print graph details
    print(f"Number of nodes: {len(G.nodes)}")
    print(f"Number of edges: {len(G.edges)}")
    

    #for u, v, data in G.edges( data=True):
        #print(f"{u}, {v}: {data}")
    
    # Find the boundaries of the graph
    boundary = find_graph_boundaries(G)

    if verbose:
        print("Graph boundaries:", boundary)
    
    # Handle boundary edges
    G = handle_boundary_edges(G, boundary)

    if verbose:
        print("Is G strongly connected:", nx.is_strongly_connected(G))

    eulerian_path = nx.eulerian_path(G)

    total_distance = 0

    if verbose:
        print("Eulerian path sequence:")

    eulerian_path_list = list(eulerian_path)

    for u, v in eulerian_path_list:
        edge_data = G[u][v]
        
        distance = edge_data.get(0).get("length")
        weight = edge_data.get(0).get("weight")
        street = edge_data.get(0).get("name")

        if verbose:
            print(f"Edge from {u} to {v}, length: {distance*weight} meters; street: {street}")
    
        total_distance += (distance*weight)    

    # Print results
    print("Total Distance:", total_distance, "meters")

    visualize_eulerian_path(G, eulerian_path_list, boundary)

# Call the main function with the path to your XML file
main('osm_smaller.xml', verbose=False)
