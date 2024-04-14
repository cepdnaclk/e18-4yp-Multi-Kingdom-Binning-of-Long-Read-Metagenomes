import numpy as np
from Bio import SeqIO
import sys
from tqdm import tqdm

# TODO: if the bin count is 0, there will be no misbinned reads and wont be able to refine

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

# files
fasta_file = initial_tool_results + '/reads.fasta'


def get_idx_maps(read_ids_file_path):
    read_id_idx = {}
    
    with open(read_ids_file_path) as read_ids_file:
        for rid in tqdm( read_ids_file):
            rid = rid.split(" ")[0].strip()[1:]  #modified: rid = rid.strip()[1:] 
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
    np.save(initial_tool_results + 'edges.npy', edges_nparr)
    return edges_nparr


# read edges.npy to get the connecting vertices when a ref vertex is given
def get_connected_vertices(ref_vertex, edges_nparray):

    # Load the .npy file
    # file_content = np.load(edges_nparray)

    # Create a boolean mask
    mask = edges_nparray[:, 0] == ref_vertex

    # Use the mask to retrieve the second element
    result = edges_nparray[mask, 1]

    return result # result comes as a numpy array
# read clssified array in classes.npz to get the bins

classes_file = initial_tool_results + '/classes.npz'

# Load the .npz file
data = np.load(classes_file)
classes = data['classes']
classified = data['classified']
from functools import lru_cache

@lru_cache(maxsize=None)
def get_bin(ref_vertex):
    global classes, classified
    if initial_binning_tool == 'lrbinner' or initial_binning_tool == 'metabcc':
        read_cluster = np.loadtxt(initial_tool_results + 'bins.txt', dtype=int)
        bin = read_cluster[ref_vertex]

    else:
        # classes_file = initial_tool_results + '/classes.npz'

        # # Load the .npz file
        # data = np.load(classes_file)

        # Extract the arrays from the loaded data
        # classes = data['classes']
        # classified = data['classified']

        # Use NumPy advanced indexing to find the class value
        mask = classified == ref_vertex
        class_value = classes[mask]

        if (len(class_value) > 0):
            bin = class_value[0]
        else:
            bin = -1

    # result is always received as an array
    return bin


# when creating the graph, the read_id is converted to an index
# example: read_1 = 0, read_2 = 1
def get_read_id(ref_vertex):
    return f"read_{ref_vertex + 1}"

# check whether a vertex is ambigous and return the corresponding ambigous vertex too
def is_ambigous(ref_vertex, edges):
    # print(1)
    ref_vertex_bin = get_bin(ref_vertex)
    connected_vertices = get_connected_vertices(ref_vertex, edges)
    # print(2)


    for connected_vertex in connected_vertices:
        connected_vertex_bin = get_bin(connected_vertex)

        if (connected_vertex_bin != ref_vertex_bin):
            # print(3)

            return ref_vertex, connected_vertex 
    # print(4)

    return []       


# find all ambigous vertices and return an array(can be used as a list)
# also write the outputs to a file
def find_all_ambigous_vertices(edges):

    # this fasta file is assumed to be created from oblr
    all_reads_from_fasta = list(SeqIO.parse(fasta_file, "fasta"))
    ambigous_vertices_file = output + '/mis_binned_reads.npy'
    isolated_vertices_file = output + '/isolated_reads.txt'
    
    ambigous_vertices = set()
    isolated_vertices = set()
    for read in tqdm(all_reads_from_fasta):

        ref_vertex = int(read.id.split('_')[-1]) - 1
        # print("ref: ", ref_vertex)

        connected_vertices = get_connected_vertices(ref_vertex, edges)
        if len(connected_vertices) > 0:
         
            if ref_vertex not in ambigous_vertices:
                results = is_ambigous(ref_vertex, edges)
                for result in results:
                    ambigous_vertices.add(result)

        else:
           isolated_vertices.add(ref_vertex)

    # write to a file
    ambigous_vertices_nparr = np.array(list(ambigous_vertices)) 
    np.save(ambigous_vertices_file, ambigous_vertices_nparr)                   
    # with open(ambigous_vertices_file, 'w') as afile:
    #   for x in ambigous_vertices:
    #     afile.write(f"{x}\t{get_bin(x)}\n")

    with open(isolated_vertices_file, 'w') as ifile:
      for y in isolated_vertices:
        ifile.write(f"{y}\t{get_bin(y)}\n")
        
    return list(ambigous_vertices), list(isolated_vertices) # an array of read_id in numerical form [1, 2, ..]


# ----------------------------------------------------------------------------------------------------------------------------


try:

    if initial_binning_tool == 'lrbinner' or initial_binning_tool == 'metabcc':
        read_id_idx = get_idx_maps(initial_tool_results + 'read_ids')
        edges = alignments_to_edges(initial_tool_results + 'reads.alns', read_id_idx)
    else:
        edges = np.load(initial_tool_results + 'edges.npy')

    misbinned_nodes, isolated_nodes = find_all_ambigous_vertices(edges)
    print(f"Successfully completed!")
    print(f"-- {len(misbinned_nodes)} mis-binned reads")
    print(f"-- {len(isolated_nodes)} isolated reads")

except Exception as e:
    print(f"Error: {e}")
