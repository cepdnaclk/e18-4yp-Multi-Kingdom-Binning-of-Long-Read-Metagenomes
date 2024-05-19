# the inconsistent labels are removed and kept unlabelled

import numpy as np
from collections import defaultdict
from tqdm import tqdm 
import argparse

from collections import defaultdict

def get_misbinned(initial_tool_results, output_folder):
    mis_binned = []
    isolated = []
    
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
    print("Checking nodes for mislabeling...")
    
    for node, node_label in tqdm(enumerate(clusters), total=len(clusters), desc="Checking nodes"):
        label_scores = defaultdict(float)

        possible_labels =  set(clusters[neighbor] for neighbor in graph[node])
        if (len(possible_labels)==0):
            isolated.append(node)
        if (node_label not in possible_labels and len(possible_labels)>0):
            mis_binned.append(node)
            new_clusters[node] = -1
        
        elif (len(possible_labels)>1):
            # Compute label scores for each possible label
            for label in possible_labels:
                for neighbor in set(graph[node]):
                    # print(neighbor, "label: ", neighbor_label, "depth: ", depth)
                    if clusters[neighbor] == label:
                        label_scores[label] += 1
        
            # Get the key of the maximum value
            new_label = max(label_scores, key=lambda k: label_scores[k])

            if (node_label != new_label and label_scores[node_label] != label_scores[new_label]):
                mis_binned.append(node)
                new_clusters[node] = -1
        
    
    print("Isolated Nodes: ", len(isolated))

    np.save(output_folder + 'mis_binned_UPDATED_imm.npy', list(mis_binned))
    np.save(output_folder + 'new_clusters_UPDATED_imm.npy', list(new_clusters))
    return new_clusters

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""OBLR GraphSAGE Routine.""")

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