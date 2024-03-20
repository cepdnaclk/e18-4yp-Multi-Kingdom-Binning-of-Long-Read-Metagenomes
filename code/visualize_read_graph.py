import numpy as np
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

def reduce_vector_sizes():
    all_edges = np.load('edges.npy')
    filtered_edges = [edges for edges in all_edges if edges[0] <= 800 and edges[1] <= 800]
    filtered_edges = np.array(filtered_edges)
    np.save('edges_modified.npy', filtered_edges)

# Load the data from the .npy files
edges = np.load('edges_100.npy')
classes = np.load('classes.npz')['classes'][:800]

# Create a graph object
G = nx.Graph()

# Add nodes to the graph
num_nodes = len(classes)
G.add_nodes_from(range(num_nodes))

# Add edges to the graph
for edge in edges:
    print(edge)

    G.add_edge(edge[0], edge[1])

# Set node colors based on the classes
node_colors = [classes[node] for node in G.nodes()]

# Draw the graph
pos = nx.spring_layout(G, scale=20)  # Adjust the scale parameter to increase the layout size
plt.figure(figsize=(30, 30))  # Adjust the figure size as needed
nx.draw(G, pos, node_color=node_colors, with_labels=True)

# Save the plot to a file
plt.savefig('network_plot.png') 

# Close the plot
plt.close()