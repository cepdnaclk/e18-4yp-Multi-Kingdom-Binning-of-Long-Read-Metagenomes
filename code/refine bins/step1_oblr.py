import numpy as np
from collections import defaultdict
from tqdm import tqdm 
import argparse

def getDiff():
    # Load the two NumPy files
    arr1 = np.load('mislabeled_nodes.npy', allow_pickle=True)
    arr2 = np.load('mis_binned_reads.npy', allow_pickle=True)

    # Find the elements that are in arr2 but not in arr1
    unique_elements = np.setdiff1d(arr2, arr1)
    print(len(unique_elements))

def npy_to_txt():
    # Load the NumPy array from the .npy file
    arr = np.load('mislabeled_nodes.npy', allow_pickle=True)

    # Save the array to a text file
    np.savetxt('myarray.txt', arr, fmt='%d')

def get_misbinned(initial_tool_results, output_folder):
    isolated_nodes_count = 0

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

    # Find the mislabeled nodes
    print("Checking nodes for mislabeling...")
    mislabeled_nodes = set()
    for node, label in tqdm(enumerate(clusters), total=len(clusters), desc="Checking nodes"):
        # Count isolated nodes
        if not graph[node]:
            isolated_nodes_count += 1
        else:
            for neighbor in graph[node]:
                if clusters[neighbor] != label:
                    mislabeled_nodes.add(node)
                    break
    print("Mislabeled nodes found.")

    # Save the mislabeled nodes to a file
    print("Saving mislabeled nodes...")
    np.save(output_folder + 'mis_binned_reads.npy', list(mislabeled_nodes))
    print(f"{len(mislabeled_nodes)} mislabeled nodes found.")
    print(f"{isolated_nodes_count} isolated nodes found.")

    # Write the results to a text file
    with open(output_folder + 'output_details.txt', 'w') as f:
        f.write(f"Number of misbinned nodes: {len(mislabeled_nodes)}\n")
        f.write(f"Number of isolated nodes: {isolated_nodes_count}\n")
        f.write(f"Total number of nodes: {len(graph)}\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""OBLR GraphSAGE Routine.""")

    # data --> results from oblr
    # output --> results from refiner tool

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