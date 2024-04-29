import numpy as np
import sys
from tqdm import tqdm
from functools import lru_cache
from collections import defaultdict

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

vertices_file = output + "/mis_binned_reads.npy"
marker_scores = output + "/marker_scores.txt"
edges_file = initial_tool_results + "/edges.npy"
updated_classes_file = output + "/new_classes.npz"

# Load marker scores into a structured array
marker_scores_dict = {(int(line.split()[0].split('_')[1])-1): line.split()[1] for line in open(marker_scores)}

# Load edges into a structured array
edges_data = np.load(edges_file)

classes_file = initial_tool_results + '/classes.npz'

# Load the .npz file
data = np.load(classes_file)
classes = data['classes']

@lru_cache(maxsize=None)
def get_bin(ref_vertex):
    return classes[ref_vertex]

def get_connected_vertices(vertex, graph):
    return graph[vertex]

@lru_cache(maxsize=None)
def find_marker_gene(read_id):
    return marker_scores_dict.get(read_id, None)

def get_marker_genes_for_connected_vertices(connected_vertices):
    connected_vertices_info = {}
    marker_genes = [find_marker_gene(vertex) for vertex in connected_vertices]
    connected_vertices_info = {vertex: {'marker_gene': marker_gene} for vertex, marker_gene in zip(connected_vertices, marker_genes)}
    return connected_vertices_info

def get_bins_for_connected_vertices(connected_vertices_info):
    bins = [get_bin(vertex) for vertex in connected_vertices_info]
    for vertex, info in connected_vertices_info.items():
        info['bin'] = bins[list(connected_vertices_info).index(vertex)]
    return connected_vertices_info

def update_bin(new_classes_file, vertices_info):
    old_classes_file = classes_file
    data = np.load(old_classes_file)
    classes = data['classes']
    classified = data['classified']

    for ref_vertex, info in vertices_info.items():
        vertex = int(ref_vertex)
        mask = classified == vertex
        classes[mask] = info['bin']

    np.savez(new_classes_file, classes=classes, classified=classified)

def annotate_bins():

    graph = defaultdict(list)
    for u, v in tqdm(edges_data, desc="Adding edges"):
        graph[u].append(v)
        graph[v].append(u)

    ambigous_vertices = np.load(vertices_file)
    actually_updated_count = 0
    replaced_as_unlabelled = 0
    relabel_with_different_bin = 0
    mis_binned_without_markers = 0
    vertices_info = {}

    for ambigous_vertex in tqdm(ambigous_vertices, desc="Analyzing misbinned reads"):
        if (ambigous_vertex == 74670):
            print(get_connected_vertices(74670, graph))
        vertex = int(ambigous_vertex)
        old_bin = get_bin(vertex)
        marker_gene = find_marker_gene(vertex)

        if marker_gene:
            vertex_bin = -1
            connected_vertices = get_connected_vertices(vertex, graph)
            connected_read_marker_info = get_marker_genes_for_connected_vertices(connected_vertices)
            connected_read_marker_info = get_bins_for_connected_vertices(connected_read_marker_info)

            marker_genes = [info['marker_gene'] for info in connected_read_marker_info.values()]
            bins = [info['bin'] for info in connected_read_marker_info.values()]

            if marker_gene in marker_genes:
                same_marker_gene_bins = [bin for bin, mg in zip(bins, marker_genes) if mg == marker_gene]
                if len(set(same_marker_gene_bins)) == 1:
                    vertex_bin = same_marker_gene_bins[0]

        else:
            vertex_bin = old_bin
            mis_binned_without_markers += 1

        if vertex_bin != -1 and vertex_bin != old_bin:
            relabel_with_different_bin += 1

        if vertex_bin != old_bin:
            actually_updated_count += 1

        if vertex_bin == -1:
            replaced_as_unlabelled += 1

        vertices_info[vertex] = {'bin': vertex_bin}

    update_bin(updated_classes_file, vertices_info)

    # Write the results to a text file
    with open(output + 'output_details.txt', 'a') as f:
        f.write(f"Mis binned reads with marker genes: {len(ambigous_vertices) - mis_binned_without_markers}\n")
        f.write(f"Number of reads replaced with a bin different to its old bin (except minus 1): {relabel_with_different_bin}\n")
        f.write(f"Reads replaced with (-1): {replaced_as_unlabelled}\n")
        f.write(f"From {len(ambigous_vertices)} reads, {actually_updated_count} reads were updated either with -1 or with a different bin\n")

try:
    annotate_bins()
    print("Updated with new bins successfully!")
except Exception as e:
    print(f"Error: {e}")