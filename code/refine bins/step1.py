import numpy as np
from Bio import SeqIO

# ----------------------------------------------------
# these should be taken as inputs

# path to the folder where oblr results are stored
oblr_results = '../Results/zymo10_1000kto1100k/oblr'

# path to the folder where output should be stored
output = '../Results/zymo10_1000kto1100k/output'

# marker file
input_marker_file = '../e18-4yp-Multi-Kingdom-Binning-of-Long-Read-Metagenomes/data/marker_genes/bacteria_archaea.hmm'

# ----------------------------------------------------

# files
fasta_file = oblr_results + '/reads.fasta'
edges_file = oblr_results + '/edges.npy'
classes_file = oblr_results + '/classes.npz'


# read edges.npy to get the connecting vertices when a ref vertex is given
def get_connected_vertices(ref_vertex):

    # Load the .npy file
    file_content = np.load(edges_file)

    # Create a boolean mask
    mask = file_content[:, 0] == ref_vertex

    # Use the mask to retrieve the second element
    result = file_content[mask, 1]

    return result # result comes as a numpy array

# read clssified array in classes.npz to get the bins
def get_bin(ref_vertex):

    # Load the .npz file
    data = np.load(classes_file)

    # Extract the arrays from the loaded data
    classes = data['classes']
    classified = data['classified']

    # Use NumPy advanced indexing to find the class value
    mask = classified == ref_vertex
    class_value = classes[mask]

    # result is always received as an array
    return class_value[0]

# when creating the graph, the read_id is converted to an index
# example: read_1 = 0, read_2 = 1
def get_read_id(ref_vertex):
    return f"read_{ref_vertex + 1}"

# check whether a vertex is ambigous and return the corresponding ambigous vertex too
def is_ambigous(ref_vertex):

    ref_vertex_bin = get_bin(ref_vertex)
    connected_vertices = get_connected_vertices(ref_vertex)

    for connected_vertex in connected_vertices:
        connected_vertex_bin = get_bin(connected_vertex)

        if (connected_vertex_bin != ref_vertex_bin):
            return ref_vertex, connected_vertex 

    return []       


# return the bins of all connected vertices
# def bins_of_connected_vertices(ref):
#     print("ref bin: ", get_bin(ref))
#     connected = get_connected_vertices(ref)
#     for vertex in connected:
#         print(vertex, get_bin(vertex))


# find all ambigous vertices and return an array(can be used as a list)
# also write the outputs to a file
def find_all_ambigous_vertices():

    # this fasta file is assumed to be created from oblr
    all_reads_from_fasta = list(SeqIO.parse(fasta_file, "fasta"))
    ambigous_vertices_file = output + '/mis_binned_reads.txt'
    
    ambigous_vertices = set()
    for read in all_reads_from_fasta:

      ref_vertex = int(read.id.split('_')[-1]) - 1
      print("ref: ", ref_vertex)
            
      if ref_vertex not in ambigous_vertices:
        results = is_ambigous(ref_vertex)
        for result in results:
          ambigous_vertices.add(result)

    # write to a file                    
    with open(ambigous_vertices_file, 'w') as file:
      for x in ambigous_vertices:
        file.write(f"{x}\t{get_bin(x)}\n")
        
    return list(ambigous_vertices) # an array of read_id in numerical form [1, 2, ..]


# ----------------------------------------------------------------------------------------------------------------------------


find_all_ambigous_vertices()