import numpy as np

vertices_file = "ambigous_vertices.txt"
marker_scores = "zymo10_1000kto1100k_genes_bacteria_test.hmm.hmmout_marker_scores.txt"
edges_file = "edges.npy"
classes_file = "classes.npz"

# Annotate bins according to the assignment of marker genes
def annotate_bin(vertices_file, marker_scores, edges_file, classes_file):
    # Extract the ambiguous vertices
    vertices = read_vertex(vertices_file)

    # Update the class.npz by updating bins of ambiguous vertices with -1
    update_bins_of_ambiguous_vertices(vertices, classes_file)

    # No of vertices which are changed its bin
    count = 0

    # Store the bin assigned to vertices
    vertices_info = {}

    # Iterate through ambiguous vertices
    for vertex in vertices:
        
        # Find the marker gene for the ambiguous vertex 
        marker_gene = find_marker_gene(marker_scores,get_read_id(vertex))

        # Mark the bin as ambiguous by assigning -1
        vertex_bin = -1
        
        # If there exist a marker gene for the ambiguous vertex
        if marker_gene:

            print(f"Marker gene found = {marker_gene} for ref_vertex {vertex}")

            # Get connected vertices for the ambigous vertex
            connected_vertices = get_connected_vertices(edges_file, vertex)

            # Prepare a dictionary with both marker gene and bin for each connected_vertex 
            connected_read_marker_info = {} 
            # Get marker genes
            connected_read_marker_info = get_marker_genes_for_connected_vertices(marker_scores, connected_vertices)
            # Get bins     
            connected_read_marker_info = get_bins_for_connected_vertices(connected_read_marker_info, classes_file)      

            # Iterate through connected vertices
            for connected_vertex, info in connected_read_marker_info.items():
                
                # Get the marker gene to the connected vertex
                connected_marker_gene = info['marker_gene']
                                
                # Compare the marker gene of connected vertex with the original marker gene
                if connected_marker_gene == marker_gene:
                    
                    print(f"Same marker gene found - {marker_gene}, bin - {info['bin']}")

                    # Same marker gene found at first, so assign its bin to the ambiguous vertex
                    if(vertex_bin == -1):
                        vertex_bin = info['bin']
                        print(f"Bin got assigned, new bin - {vertex_bin}")

                    # Same marker gene found , but the bin is different, so keep the vertex as ambiguous
                    elif(vertex_bin != info['bin']):
                        vertex_bin = -1
                        ambigous = True
                        print(f"Vertex {vertex} is ambiguous")
                        break

        # Count the no of vertices assigned by bins
        if vertex_bin != -1:
            count = count +1
        
        # Store the bin's info for each vertex
        vertices_info[vertex] = {'bin': vertex_bin}

    # Update the bin of the vertex
    update_bin(classes_file, vertices_info)

    print(f"Count = {count}")    



# Update bins of ambiguous vertices in to -1
def update_bins_of_ambiguous_vertices(vertices, classes_file):
    # Load the classes.npz file
    data = np.load(classes_file)

    # Extract the arrays from the loaded data
    classes = data['classes']
    classified = data['classified']

    for ref_vertex in vertices:
        # vertex should be an integer to use as an index
        vertex = int(ref_vertex)

        # Use NumPy advanced indexing to find the class value
        mask = classified == vertex
        classes[mask] = -1
    
    # Save the updated classified array back to the .npz file
    np.savez(classes_file, classes=classes, classified=classified)


# Update the bin when the vertices is given with its bin
def update_bin(classes_file, vertices_info):

    # Load the classes.npz file
    data = np.load(classes_file)

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

        print(f"Bin - {classes[mask]}")

    # Save the updated classified array back to the .npz file
    np.savez(classes_file, classes=classes, classified=classified)




# Find marker genes of connected vertices
def get_marker_genes_for_connected_vertices(marker_scores, connected_vertices):

    # Create a dictionary to store the vertex , its marker gene and its bin
    connected_vertices_info = {}

    # Iterate through the connected vertices
    for connected_vertex in connected_vertices:

        # Find the marker gene
        marker_gene = find_marker_gene(marker_scores, get_read_id(connected_vertex))

        # Store the marker gene for respective vertex
        # connected_vertices_info = {connected_vertex: {'marker_gene': find_marker_gene(marker_scores, get_read_id(connected_vertex))}}
        connected_vertices_info[connected_vertex] = {}
        connected_vertices_info[connected_vertex]['marker_gene'] = marker_gene

    # Return the dictionary which has stored the marker gene for each vertex
    return connected_vertices_info




# Get bins for each connected vertices
def get_bins_for_connected_vertices(connected_vertices_info, classes_file):

    # Iterate through connected vertices
    for connected_vertex in connected_vertices_info:

        # Find the bin 
        bin = get_bin(connected_vertex, classes_file)

        # Assign the bin to each vertex
        connected_vertices_info[connected_vertex]['bin'] = get_bin(connected_vertex, classes_file)

    # Return the dictionary which has both marker genes and bins found
    return connected_vertices_info




# Get the bin
def get_bin(ref_vertex, class_file):

    # ref_vertex should be an integer to use as an index
    ref_vertex = int(ref_vertex)

    # Load the .npz file
    data = np.load(class_file)

    # Extract the arrays from the loaded data
    classes = data['classes']
    classified = data['classified']

    # Use NumPy advanced indexing to find the class value
    mask = classified == ref_vertex
    class_value = classes[mask]

    # result is always received as an array
    return class_value[0]
                


# Extract the ambiguous vertices
def read_vertex(filename):

    # Take vertices as an array
    vertices = []

    # Open the ambigous_vertices.txt and extract the first entry of each row which is the vertex
    with open(filename, "r") as file:
        for line in file:
            values = line.split()
            vertex = values[0]
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

    # Search the read_id and find the related marker gene
    with open(filename, "r") as file:
        for line in file:
            columns = line.split()
            if len(columns) >= 3:
                if columns[1] == read_id:
                    return columns[2]
    return None

# Get connected vertices to a given vertex
def get_connected_vertices(edges_file, vertex):
    
    # vertex should be an integer
    vertex = int(vertex)

    # Load the .npy file
    file_content = np.load(edges_file)

    # Create a boolean mask
    mask = file_content[:, 0] == vertex

    # Use the mask to retrieve the second element
    result = file_content[mask, 1]
    
    return result


# Call the function - This will update the classes.npz with the new bin assignment
annotate_bin(vertices_file, marker_scores, edges_file, classes_file)
