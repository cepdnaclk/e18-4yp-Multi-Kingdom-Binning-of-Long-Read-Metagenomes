# SemiBin2 --> there can be unlabelled nodes
# the inconsistent labels are relabelled

import numpy as np
from collections import defaultdict, deque
from tqdm import tqdm 
import argparse

def read_tsv_to_clusters(tsv_file, total_reads):
    
    # Initialize the clusters array with -1
    clusters = np.full(total_reads, -1, dtype=int)

    with open(tsv_file, 'r') as file:
        file.readline()
        for line in file:
            read_id, bin = line.strip().split('\t')
            read_index = int(read_id.split('_')[1]) - 1
            clusters[read_index] = int(bin)

    return clusters


def bfs_label_scores(graph, clusters, start_node):

    queue = deque([(start_node, 1)])
    visited = set()
    label_scores = defaultdict(float)
    
    while queue:
        node, depth = queue.popleft()
        visited.add(node)
        
        for neighbor in set(graph[node]):
            if neighbor not in visited:
                if clusters[neighbor] != -1:
                    score = (2 ** (-depth))
                    label_scores[clusters[neighbor]] += score
                else:
                    queue.append((neighbor, depth + 1))
                visited.add(neighbor)
    
    return label_scores

def get_misbinned(initial_tool_results, output_folder, total_reads):

    mis_binned = []

    # load the edges
    edges = np.load(initial_tool_results + 'edges.npy')
    # load the bins
    clusters = read_tsv_to_clusters(initial_tool_results + 'contig_bins.tsv', total_reads)

    # creating a dictionary with nodes as keys and neighbors as values
    graph = defaultdict(list)
    for u, v in tqdm(edges, desc="Adding edges"):
        graph[u].append(v)
        graph[v].append(u)
    print("Graph created.")

    # Find the mislabeled nodes and compute label scores
    for node, node_label in tqdm(enumerate(clusters), total=len(clusters), desc="Checking nodes"):

        # possible labels from directly connected nodes
        possible_labels =  set(clusters[neighbor] for neighbor in graph[node])

        # obviously mis binned nodes
        if len(possible_labels)>0 and -1 not in possible_labels and node_label not in possible_labels :
            mis_binned.append(node)
            clusters[node] = -1

        # nodes with the same label might not be directly connected, but connected through a path with unlabelled nodes
        elif (node_label != -1 and node_label not in possible_labels): 

            label_scores = bfs_label_scores(graph, clusters, node)

            # not even the indirectly connected nodes are having the same label
            if node_label not in label_scores and len(label_scores) > 0:
                mis_binned.append(node)
                clusters[node] = -1
            
            # nodes connected through unlabelled path is having the same label
            # elif (label_scores[node_label] < 0.5): # by this point, there are no directly connected ones with same label
            current_label_score  = label_scores[node_label]*1.5
            keys_greater_than_threshold = [key for key, value in label_scores.items() if value > current_label_score]
            if (len(keys_greater_than_threshold) > 0):
                clusters[node] = keys_greater_than_threshold[0] # relabelling


    np.save(output_folder + 'relabelled_clusters.npy', clusters)
    np.save(output_folder + 'misbinned_reads.npy', list(mis_binned))

    return mis_binned

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Identify misbinned nodes""")

    parser.add_argument('--data', '-d',
                        help="folder where initial tool resulted files are stored",
                        type=str,
                        required=True)
    parser.add_argument(
        '--output', '-o', help="Output directory where the outputs from refining tool are stored", type=str, required=True)
    
    parser.add_argument(
        '--reads', '-n', help="total no of reads", type=int, required=True)

    args = parser.parse_args()

    initial_tool_results = args.data
    output_folder = args.output
    total_reads = args.reads

    get_misbinned(initial_tool_results, output_folder, total_reads)