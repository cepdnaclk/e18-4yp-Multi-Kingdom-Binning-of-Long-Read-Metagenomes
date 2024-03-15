from Bio import SeqIO

"""
Convert a fastq file to fasta format
"""

with open("../../../Results/datasets/zymo10_1000kto102k.fastq", "r") as input_handle:
    with open("../../../Results/datasets/zymo10_1000kto102k.fasta", "w") as output_handle:
        sequences = SeqIO.parse(input_handle, "fastq")
        count = SeqIO.write(sequences, output_handle, "fasta")
        
print("Converted %i records" % count)