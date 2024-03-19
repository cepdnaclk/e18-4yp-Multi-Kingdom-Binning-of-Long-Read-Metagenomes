import numpy as np
import sys

if len(sys.argv) != 3:
    print("Usage: python script.py classes.npz read_ids")
    sys.exit(1)

classes_file = sys.argv[1]
read_ids_file = sys.argv[2]

# Unzipping and reading the classes.npz file
with np.load(classes_file) as data:
    class_array = data['classes']

# Reading read_ids file and extracting read names
with open(read_ids_file, 'r') as file:
    read_names = [line.split()[0][1:] for line in file if line.strip()]

# Mapping values from class_array to read_names
mapped_values = dict(zip(read_names, class_array))

# Writing mapped values to bins.txt
with open('bins.txt', 'w') as output_file:
    for read_name in read_names:
        output_file.write(f"{read_name} {mapped_values.get(read_name)}\n")