import sys
import random
from Bio import SeqIO

# Check if the required arguments are provided
if len(sys.argv) != 3:
    print("Usage: python script.py input.fastq output.fastq")
    sys.exit(1)

# Get the input file path and output file path
input_file = sys.argv[1]
output_file = sys.argv[2]

# Parse the input FASTQ file
records = list(SeqIO.parse(input_file, "fasta"))

# Randomly select 10000 sequences
random_records = random.sample(records, 10000)

# Write the selected sequences to the output file
with open(output_file, "w") as out_handle:
    SeqIO.write(random_records, out_handle, "fasta")

print(f"Wrote 10000 random sequences to {output_file}")