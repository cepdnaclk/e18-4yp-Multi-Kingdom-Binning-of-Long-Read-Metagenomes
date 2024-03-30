import numpy as np

# Read the file
data = {}
MG_dict = {}

with open("reads_genes_bacteria.hmmout_marker_scores.txt", "r") as file:
    marker_count = 0
    for line in file:
        line = line.strip().split()
        read_index = int(line[0].split("_")[1]) - 1
        marker_gene = line[2]
        if marker_gene not in MG_dict:
            MG_dict[line[2]] = marker_count
            marker_count += 1
        # data.append([int(read_index), float(line[3])])
        if MG_dict[marker_gene] in data:
            data[MG_dict[marker_gene]].append(read_index)
        else:
            data[MG_dict[marker_gene]] = []
            data[MG_dict[marker_gene]].append(read_index)
        

# Get the number of rows and columns
num_rows = 100000
num_cols = len(MG_dict)

# Create an empty NumPy array
matrix = np.zeros((num_rows, num_cols), dtype=int)

# Fill the matrix based on the dictionary values
row_index = 0
for col, values in data.items():
    for value in values:
        matrix[value, col] = 1
        # row_index += 1

# Save the matrix to a .npy file
np.save("presence_of_markers.npy", matrix)

print(matrix[1])