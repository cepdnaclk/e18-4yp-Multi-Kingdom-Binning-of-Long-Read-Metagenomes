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

def get_misbinned(initial_results_folder, output_folder):
    # Load the edges
    print("Determining edges...")
    edges = create_edges(initial_results_folder)

    print("Loading initial tool results...")
    # Load data from bins text file
    clusters = np.loadtxt(initial_tool_results + 'bins.txt', dtype=int) # TODO: file name might be different for metabcc
    np.save(initial_results_folder + 'classes.npy' , clusters)

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
    np.save(output_folder + 'mislabeled_nodes.npy', list(mislabeled_nodes))
    print(f"{len(mislabeled_nodes)} mislabeled nodes found.")

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