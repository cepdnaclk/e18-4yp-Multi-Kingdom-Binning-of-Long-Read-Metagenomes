import numpy as np
from collections import defaultdict
from tqdm import tqdm 

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

def get_misbinned():
    # Load the edges
    print("Loading files...")
    edges = np.load('edges.npy')

    # Load the classes
    data = np.load('classes.npz')
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
        for neighbor in graph[node]:
            if clusters[neighbor] != label:
                mislabeled_nodes.add(node)
                break
    print("Mislabeled nodes found.")

    # Save the mislabeled nodes to a file
    print("Saving mislabeled nodes...")
    np.save('mis_binned_reads.npy', list(mislabeled_nodes))
    print(f"{len(mislabeled_nodes)} mislabeled nodes found.")


# get_misbinned()