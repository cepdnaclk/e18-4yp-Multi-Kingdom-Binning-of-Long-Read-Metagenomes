import numpy as np
import sys
from tqdm import tqdm

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


# Find marker genes of connected vertices
def get_marker_genes_for_connected_vertices(marker_scores, connected_vertices):

    # Create a dictionary to store the vertex , its marker gene and its bin
    connected_vertices_info = {}

    # Iterate through the connected vertices
    for connected_vertex in connected_vertices:

        # Find the marker gene
        marker_gene = find_marker_gene(marker_scores, connected_vertex)

        # Store the marker gene for respective vertex
        # TODO
        # connected_vertices_info = {connected_vertex: {'marker_gene': marker_gene}}  
        connected_vertices_info[connected_vertex] = {}
        connected_vertices_info[connected_vertex]['marker_gene'] = marker_gene

    # Return the dictionary which has stored the marker gene for each vertex
    return connected_vertices_info


# Get bins for each connected vertices
def get_bins_for_connected_vertices(connected_vertices_info):

    # Iterate through connected vertices
    for connected_vertex in connected_vertices_info:

        # Find the bin 
        bin = get_bin(connected_vertex)

        # Assign the bin to each vertex
        connected_vertices_info[connected_vertex]['bin'] = bin

    # Return the dictionary which has both marker genes and bins found
    return connected_vertices_info


classes_file = initial_tool_results + '/classes.npz'

# Load the .npz file
data = np.load(classes_file)
classes = data['classes']
classified = data['classified']
from functools import lru_cache

@lru_cache(maxsize=None)
def get_bin(ref_vertex):
    global classes, classified
    # if initial_binning_tool == 'lrbinner' or initial_binning_tool == 'metabcc':
    #     read_cluster = np.loadtxt(initial_tool_results + 'bins.txt', dtype=int)
    #     bin = read_cluster[ref_vertex]

    # else:

    # Use NumPy advanced indexing to find the class value
    mask = classified == ref_vertex
    class_value = classes[mask]

    if (len(class_value) > 0):
        bin = class_value[0]
    else:
        bin = -1

    # result is always received as an array
    return bin
    

# Extract the ambiguous vertices
def read_vertex(filename):

    # Take vertices as an array
    vertices = []

    # Open the ambigous_vertices.txt and extract the first entry of each row which is the vertex
    with open(filename, "r") as file:
        for line in file:
            vertex = line.strip().split("\t")[0]
            vertices.append(vertex)

    # Return the array with extracted vertices
    return vertices

# Get read_id using vertex
# when creating the graph, the read_id is converted to an index
# example: read_1 = 0, read_2 = 1
def get_read_id(ref_vertex):

    # ref_vertex should be an integer
    vertex_number = int(ref_vertex) + 1
    return f"read_{vertex_number}"

# Get the marker gene using read_id
def find_marker_gene(filename, read_id):
    # numerical_vertex = int(read_id.split("_")[1]) - 1

    # Search the read_id and find the related marker gene
    with open(filename, "r") as file:
        for line in file:
            columns = line.split()
            if len(columns) >= 3:
                numerical_read_id = int(columns[0].split("_")[1]) - 1
                if numerical_read_id == read_id:
                    return columns[1]
    return None

# Get connected vertices to a given vertex
def get_connected_vertices(edges_file, vertex):
    
    # vertex should be an integer
    # numerical_vertex = int(vertex.split("_")[1]) - 1

    # Load the .npy file
    file_content = np.load(edges_file)

    # Create a boolean mask
    mask = file_content[:, 0] == vertex

    # Use the mask to retrieve the second element
    result = file_content[mask, 1]
    
    return result

# Update the bin when the vertices is given with its bin
def update_bin(new_classes_file, vertices_info):
    classes_file = initial_tool_results + '/classes.npz'
    old_classes_file = classes_file

    # Load the classes.npz file
    data = np.load(old_classes_file)

    # Extract the arrays from the loaded data
    classes = data['classes']
    classified = data['classified']

    # Iterate through each vertex
    for ref_vertex in vertices_info:

        # vertex should be an integer to use as an index
        vertex = int(ref_vertex)

        # Use NumPy advanced indexing to find the class value
        mask = classified == vertex
        classes[mask] = vertices_info[ref_vertex]['bin']

        #print(f"Bin - {classes[mask]}")

    # Save the updated classified array back to the .npz file
    np.savez(new_classes_file, classes=classes, classified=classified)


def annotate_bins():

    # Extract the ambiguous vertices
    # ambigous_vertices = read_vertex(vertices_file)
    ambigous_vertices = np.load(vertices_file)

    # No of vertices which are changed its bin
    actually_updated_count = 0
    replaced_as_unlabelled = 0
    relabel_with_different_bin = 0
    mis_binned_without_markers = 0

    # Store the bin assigned to vertices
    vertices_info = {}

    # Iterate through ambiguous vertices
    for ambigous_vertex in tqdm(ambigous_vertices):
        vertex = int(ambigous_vertex) # since these are read from txt file, type is str

        old_bin = get_bin(vertex)

        # Find the marker gene for the ambiguous vertex 
        marker_gene = find_marker_gene(marker_scores,vertex)
        
        # If there exist a marker gene for the ambiguous vertex
        if marker_gene:

            # Mark the bin as ambiguous by assigning -1
            vertex_bin = -1

           # print(f"Marker gene found = {marker_gene} for ref_vertex {vertex}")

            # Get connected vertices for the ambigous vertex
            connected_vertices = get_connected_vertices(edges_file, vertex)

            # Prepare a dictionary with both marker gene and bin for each connected_vertex 
            connected_read_marker_info = {} 

            # Get marker genes
            connected_read_marker_info = get_marker_genes_for_connected_vertices(marker_scores, connected_vertices)
            # Get bins     
            connected_read_marker_info = get_bins_for_connected_vertices(connected_read_marker_info)

            # Iterate through connected vertices
            for connected_vertex, info in connected_read_marker_info.items():
                
                # Get the marker gene to the connected vertex
                connected_marker_gene = info['marker_gene']
                                
                # Compare the marker gene of connected vertex with the original marker gene
                if connected_marker_gene == marker_gene:
                    
                    #print(f"Same marker gene found - {marker_gene}, bin - {info['bin']}")

                    # Same marker gene found at first, so assign its bin to the ambiguous vertex
                    if(vertex_bin == -1):
                        vertex_bin = info['bin']
                        #print(f"Bin got assigned, new bin - {vertex_bin}")

                    # Same marker gene found , but the bin is different, so keep the vertex as ambiguous
                    elif(vertex_bin != info['bin']):
                        vertex_bin = -1
                        ambigous = True
                        #print(f"Vertex {vertex} is ambiguous")
                        break
        
        # if no marker gene, dont change the bin
        else:
            vertex_bin = old_bin
            mis_binned_without_markers += 1

        #print("------------------------------------------------------------------------------------")

        # Count the no of vertices assigned by bins
        if vertex_bin != -1 and vertex_bin != old_bin:
            relabel_with_different_bin += 1
        
        if vertex_bin != old_bin:
            actually_updated_count += 1

        if vertex_bin == -1:
            replaced_as_unlabelled += 1

        # Store the bin's info for each vertex
        vertices_info[vertex] = {'bin': vertex_bin}

    # Update the bin of the vertex
    update_bin(updated_classes_file, vertices_info)
    print(f"Mis binned reads without any marker gene = {mis_binned_without_markers}")
    print(f"Number of reads replaced with a bin different to its old bin (except minus 1)= {relabel_with_different_bin}") 
    print(f"Reads replaced with (-1) = {replaced_as_unlabelled}") 
    print(f"From {len(ambigous_vertices)} reads, {actually_updated_count} of reads were updated either by -1 or different bin.") 

# ----------------------------------------------------------------------------------------------------------------------------

try:
    annotate_bins()
    print("Updated with new bins successfully!")
except Exception as e:
    print(f"Error: {e}")