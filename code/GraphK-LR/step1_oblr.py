# OBLR, LRBINNER, METABCC
# for lrbinner and metabcc-lr, first the read overlap graph should be created and the edges.npy file should be already prepared
# the inconsistent labels are removed and kept unlabelled

import numpy as np
from collections import defaultdict
from tqdm import tqdm 
import argparse

from collections import defaultdict

def get_misbinned(initial_tool_results, output_folder):
    mis_binned = []
    
    # Load the edges
    print("Loading files...")
    edges = np.load(initial_tool_results + 'edges.npy')
    
    # Load the classes
    data = np.load(initial_tool_results + 'classes.npz')
    clusters = data['classes']
    print("Files loaded.")
    
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