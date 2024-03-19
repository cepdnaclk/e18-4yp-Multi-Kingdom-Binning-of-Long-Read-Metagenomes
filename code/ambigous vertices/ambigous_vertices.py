import numpy as np

ambigous_vertices_file = "ambigous_vertices.txt"

# read data.npz to get the connecting vertices when a ref vertex is given
def read_data_npz(ref_vertex):

    # Load the .npz file
    file_content = np.load('data.npz')

    # Get the list of array names
    array_names = file_content.files

    # get the edges
    data = file_content[array_names[0]]

    # Create a boolean mask
    mask = data[:, 0] == ref_vertex

    # Use the mask to retrieve the second element
    result = data[mask, 1]

    return result # result comes as a numpy array

# read edges.npy to get the connecting vertices when a ref vertex is given
def get_connected_vertices(ref_vertex):

    # Load the .npy file
    file_content = np.load('edges.npy')

    # Create a boolean mask
    mask = file_content[:, 0] == ref_vertex

    # Use the mask to retrieve the second element
    result = file_content[mask, 1]

    return result # result comes as a numpy array

# read edges.npy to get the connecting vertices when a ref vertex is given
def get_bin(ref_vertex):

    # Load the .npz file
    data = np.load('classes.npz')

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
    # print(ref_vertex_bin)

    for connected_vertex in connected_vertices:
        connected_vertex_bin = get_bin(connected_vertex)
        # print(connected_vertex_bin)

        if (connected_vertex_bin != ref_vertex_bin):
            return ref_vertex, connected_vertex 

    return []       

# find all ambigous vertices and return a numpy array(can be used as a list)
# also write the outputs to a file
def find_all_ambigous_vertices():

    with open(ambigous_vertices_file, 'w') as file:

        ambigous_vertices = []
        for ref_vertex in range (0, 100000):
            print("ref: ", ref_vertex)
            if ref_vertex not in ambigous_vertices:
                results = is_ambigous(ref_vertex)
                for result in results:
                    # ambigous_vertices.append(result)
                    ambigous_vertices.append(get_read_id(result))
                    file.write(f"{result}\t{get_bin(result)}\n")
        
    return ambigous_vertices

# return the bins of all connected vertices
def bins_of_connected_vertices(ref):
    print("ref bin: ", get_bin(ref))
    connected = get_connected_vertices(ref)
    for vertex in connected:
        print(vertex, get_bin(vertex))


print(len(find_all_ambigous_vertices()))


# [22, 56238, 225, 79252, 304, 45639, 307, 56270, 308, 99766, 392, 88937, 395, 88937, 498, 35904, 605, 23215, 639, 68700, 640, 90285, 754, 65502]
''' ['read_23', 'read_56239', 'read_226', 'read_79253', 'read_305', 'read_45640', 'read_308', 'read_56271', 'read_309', 'read_99767', 'read_393', 
  'read_88938', 'read_396', 'read_88938', 'read_499', 'read_35905', 'read_606', 'read_23216', 'read_640', 'read_68701', 'read_641', 'read_90286',
  'read_755', 'read_65503'] '''