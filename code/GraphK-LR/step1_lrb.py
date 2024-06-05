# lrbinner --> assuming that there are no unlabelled nodes
# the inconsistent labels are relabelled

import numpy as np
from collections import defaultdict
from tqdm import tqdm 
import argparse

def create_edges(initial_tool_results):
    # first, the fastq file should be converted to a fasta file (external tool)
    # second, read_id conversion (external tool)
    # third, chunking and creating reads.alns and degree files (kbm2 + build_graph.sh code from oblr)

    read_id_idx = get_idx_maps(initial_tool_results + 'read_ids')

    # create and save the edge.npy file
    edges_nparr = alignments_to_edges(initial_tool_results + 'reads.alns', read_id_idx)

    np.save(initial_tool_results + 'edges.npy', edges_nparr)

    return edges_nparr


def get_idx_maps(read_ids_file_path):
    read_id_idx = {}
    
    with open(read_ids_file_path) as read_ids_file:
        for rid in tqdm( read_ids_file):
            rid = rid.split(" ")[0].strip()  #modified: rid = rid.strip()[1:] 
            read_id_idx[rid] = len(read_id_idx)
    return read_id_idx


def alignments_to_edges(alignments_file_path, read_id_idx):
    edges = []
    with open(alignments_file_path, "r") as af:
        for line in tqdm(af):
            u, v = line.strip().split('\t')
            if u != v:
                edges.append((read_id_idx[u], read_id_idx[v]))
    
    edges_nparr = np.array(edges, dtype=np.int32)
    return edges_nparr


def get_misbinned(initial_results_folder, output_folder):
    mis_binned = []

    # Load the edges
    print("Determining edges...")
    edges = create_edges(initial_results_folder)

    print("Loading initial tool results...")
    # Load data from bins text file
    clusters = np.loadtxt(initial_tool_results + 'bins.txt', dtype=int) 
    np.save(initial_results_folder + 'classes.npy' , clusters)

    # Create the graph
    print("Creating graph...")
    graph = defaultdict(list)
    for u, v in tqdm(edges, desc="Adding edges"):
        graph[u].append(v)
        graph[v].append(u)
    print("Graph created.")

    new_clusters = clusters.copy()
    
    # Find the mislabeled nodes and compute label scores    
    for node, node_label in tqdm(enumerate(clusters), total=len(clusters), desc="Checking nodes"):
        label_scores = defaultdict(float)

        possible_labels =  set(clusters[neighbor] for neighbor in graph[node])
        if (node_label not in possible_labels and len(possible_labels)>0):
            mis_binned.append(node)
            new_clusters[node] = -1
        
        elif (len(possible_labels)>1):

            # Compute label scores for each possible label
            for label in possible_labels:
                for neighbor in set(graph[node]):
                    if clusters[neighbor] == label:
                        score = (2 ** (-1))
                        label_scores[label] += score
        
            # Get the key of the maximum value
            current_label_score  = label_scores[node_label]*1.5
            keys_greater_than_threshold = [key for key, value in label_scores.items() if value > current_label_score]
            if (len(keys_greater_than_threshold) > 0):
                new_clusters[node] = keys_greater_than_threshold[0] # relabel

    np.save(output_folder + 'misbinned_reads.npy', list(mis_binned))
    np.save(output_folder + 'relabelled_clusters.npy', list(new_clusters))
    return new_clusters

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Identify misbinned nodes""")

    parser.add_argument('--data', '-d',
                        help="folder where initial tool resulted files are stored",
                        type=str,
                        required=True)
    parser.add_argument(
        '--output', '-o', help="Output directory where the outputs from refining tool are stored", type=str, required=True)

    args = parser.parse_args()

    initial_tool_results = args.data
    output_folder = args.output

    get_misbinned(initial_tool_results, output_folder)

    