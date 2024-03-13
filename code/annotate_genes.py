import subprocess
import sys

# Check if the required arguments are provided
if len(sys.argv) != 5:
    print("Incorrect number of agruments provided")
    print("Usage: python script.py [in_fasta_file] [out_faa_file] [in_hmm_file] [out_hmmout_file]")
    sys.exit(1)

# Get the input file path, output file path, and the starting read index
fasta_file = sys.argv[1]
faa_file = sys.argv[2]
marker_hmm_file = sys.argv[3]
hmmout_file = sys.argv[4]

# Define the command arguments
prodigal_cmd = ["prodigal", "-i", fasta_file, "-a", faa_file]
hmmsearch_cmd = ["hmmsearch", marker_hmm_file, faa_file]

# Run the Prodigal command
try:
    subprocess.run(prodigal_cmd, check=True)
except subprocess.CalledProcessError as e:
    print(f"Prodigal command failed with error code {e.returncode}")
    exit(1)

# Run the hmmsearch command and capture the output
try:
    hmmsearch_output = subprocess.run(hmmsearch_cmd, check=True, stdout=subprocess.PIPE, universal_newlines=True)
    with open(hmmout_file, "w") as f:
        f.write(hmmsearch_output.stdout)
except subprocess.CalledProcessError as e:
    print(f"hmmsearch command failed with error code {e.returncode}")
    exit(1)

print("Commands executed successfully!")