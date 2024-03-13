import sys
from Bio import SeqIO

# Check if the required arguments are provided
if len(sys.argv) != 4:
    print("Usage: python script.py input.fastq output.fastq start_read")
    sys.exit(1)

# Get the input file path, output file path, and the starting read index
input_file = sys.argv[1]
output_file = sys.argv[2]
start_read = int(sys.argv[3])
read_count = int(sys.argv[4])
end_read = start_read + read_count

# Parse the input FASTQ file
records = SeqIO.parse(input_file, "fastq")

# Write the desired range of reads to the output file
with open(output_file, "w") as out_handle:
    count = 0
    for record in records:
        if count >= start_read and count < end_read:
            SeqIO.write(record, out_handle, "fastq")
        count += 1

print(f"Wrote reads {start_read} to {end_read - 1} to {output_file}")
