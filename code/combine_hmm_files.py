import os

# Function to read contents of a file
def read_file(file_path):
    with open(file_path, 'r') as file:
        return file.read()

# Get folder path from the user
folder_path = input("Enter the folder path containing the HMM files: ")

# Ensure the folder path is valid
if not os.path.isdir(folder_path):
    print("Invalid folder path!")
    exit()

# Output file name
output_file = "combined_hmm.txt"

# List to store contents of all HMM files
all_contents = []

# Get the total number of HMM files in the folder
total_files = sum(1 for filename in os.listdir(folder_path) if filename.endswith(".hmm"))

# Initialize counter for the current file
current_file_count = 0

# Iterate through each file in the folder
for filename in os.listdir(folder_path):
    if filename.endswith(".hmm"):
        file_path = os.path.join(folder_path, filename)
        # Read content of the file
        content = read_file(file_path)
        # Append content to the list
        all_contents.append(content)
        
        # Increment current file count
        current_file_count += 1
        
        # Calculate percentage completion
        completion_percentage = (current_file_count / total_files) * 100
        
        # Print progress
        print(f"Progress: {completion_percentage:.2f}% completed", end='\r')

# Combine contents
combined_contents = "\n//\n".join(all_contents)

# Write combined contents to the output file
with open(output_file, 'w') as out_file:
    out_file.write(combined_contents)

print("\nCombined HMM file created successfully!")
