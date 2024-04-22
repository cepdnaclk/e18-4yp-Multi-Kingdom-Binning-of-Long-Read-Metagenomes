import numpy as np
import sys
from tqdm import tqdm

# ----------------------------------------------------
# these should be taken as inputs

# Check if the correct number of arguments are provided
if len(sys.argv) != 4:
    print("Usage: python script.py initial_tool_results_folder output_folder initial_tool")
    sys.exit(1)

# Take paths from command line arguments
initial_tool_results = sys.argv[1]
output = sys.argv[2]
initial_binning_tool = sys.argv[3]

# ----------------------------------------------------

misbinned_nodes_file = output + "/mis_binned_reads.npy"
marker_scores_file = output + "/marker_scores.txt"
edges_file = initial_tool_results + "/edges.npy"
updated_classes_file = output + "/new_classes_step3_1.npz"

# Load the data
marker_scores = {}
with open(marker_scores_file) as f:
    for line in f:
        node_id_str, x1, marker_gene, x2 = line.strip().split("\t")
        node_id = int(node_id_str.split('_')[1]) - 1
        marker_scores[int(node_id)] = marker_gene

misbinned_nodes = np.load(misbinned_nodes_file)
classes_file = initial_tool_results + '/classes.npz'
classes = np.load(classes_file)["classes"]
edges = np.load(edges_file)

# Function to update the bins of misbinned nodes
def update_bins(misbinned_nodes):
    # Get the marker genes and current bins for the misbinned nodes
    marker_genes = np.array([marker_scores.get(node_id, None) for node_id in misbinned_nodes])
    curr_bins = classes[misbinned_nodes]

    # Get the connected nodes for the misbinned nodes
    connected_nodes_list = [edges[node_id] for node_id in misbinned_nodes]
    connected_nodes = np.concatenate(connected_nodes_list)

    # Get the bins and marker genes of the connected nodes
    connected_bins = classes[connected_nodes]
    connected_marker_genes = np.array([marker_scores.get(node_id, None) for node_id in connected_nodes])

    # Compute the new bins for the misbinned nodes
    new_bins = np.zeros_like(curr_bins)
    # Compute the new bins for the misbinned nodes
    # new_bins = curr_bins.copy()  # Initialize with the current bins


    # Case 1: Nodes with marker genes
    has_marker = marker_genes != None
    for node_idx in tqdm(np.where(has_marker)[0], desc="Processing nodes with marker genes", total=np.sum(has_marker), unit="node"):
    # for node_idx in np.where(has_marker)[0]:
        marker_gene = marker_genes[node_idx]
        node_connected_bins = connected_bins[np.isin(connected_nodes, connected_nodes_list[node_idx])]
        node_connected_marker_genes = connected_marker_genes[np.isin(connected_nodes, connected_nodes_list[node_idx])]
        same_marker_connected_bins = node_connected_bins[node_connected_marker_genes == marker_gene]
        if len(np.unique(same_marker_connected_bins)) == 1:
            new_bins[node_idx] = same_marker_connected_bins[0]

    # Case 2: Nodes without marker genes
    no_marker = marker_genes == None
    for node_idx in tqdm(np.where(no_marker)[0], desc="Processing nodes without marker genes", total=np.sum(no_marker), unit="node"):
        node_connected_bins = connected_bins[np.isin(connected_nodes, connected_nodes_list[node_idx])]
        bin_counts = np.bincount(node_connected_bins[node_connected_bins >= 0])
        if bin_counts.size > 0:
            new_bins[node_idx] = np.argmax(bin_counts)

    # Update the classes array with the new bins
    classes[misbinned_nodes] = np.where(new_bins == 0, -1, new_bins)

# Update the bins of misbinned nodes
update_bins(misbinned_nodes)

# Save the updated classes
np.savez(updated_classes_file, classes=classes)