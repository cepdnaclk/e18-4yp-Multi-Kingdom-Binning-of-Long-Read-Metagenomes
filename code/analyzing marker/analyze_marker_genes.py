import logging
import os
import sys
from Bio import SeqIO


# create logger
logger = logging.getLogger("== analyze_marker_genes.py ==")

def convert_fastq_fasta(input_file):
    split_file_name = input_file.rsplit('.', 1)
    output_fasta_file_name = split_file_name[0]
    output_file = f"{output_fasta_file_name}.fasta"
    
    with open(input_file, "r") as input_handle:
        with open(output_file, "w") as output_handle:
            sequences = SeqIO.parse(input_handle, "fastq")
            count = SeqIO.write(sequences, output_handle, "fasta")
    print("Converted %i records" % count)
    
    return output_fasta_file_name

def generate_gene_faa_file(fasta_file):
    # Define the command as a list (safer than a string)
    gene_file = f"{fasta_file}_genes"
    fragCmd = (
                "prodigal"
                + " -i "
                + fasta_file
                + ".fasta"
                + " -a "
                + gene_file
                + ".faa"
            )
    logger.debug(f"exec cmd: {fragCmd}")
    os.system(fragCmd)
    
    return gene_file #returning the file name: string
    
# 1e1137d3-1334-4bcd-be69-ef08a8665b18_1 # 3 # 203 # 1 # ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.3>
# CAAFRSVTYLLQLITQLLPGKYLLNIQVFNVLRVVFDEFTTRFHAVTHQHFKSFIGINSV
# FDIDLK*

def get_id_map(gene_file):
    # Open the input file
    output_file_name = f"{gene_file}_id_mapper.txt"
    read_id_data = {}
    
    with open(f"{gene_file}.faa", "r") as input_file:
        # Read each line in the input file
        for line in input_file:
            # Check if the line starts with "<"
            if line.startswith(">"):
                # Split the line by spaces
                parts = line.strip().split()
                target_id = parts[0].lstrip(">")
                target_id = target_id.split("_")
                target_id = target_id[0]
                # Iterate through each part to find the one containing "ID="
                for part in parts:
                    if "ID=" in part:
                        # Extract the desired information from the part
                        id_part = part.split(";")[0]  # Extracting ID part before ";"
                        id_value = id_part.split("=")[1]  # Extracting value after "="
                        # Extracting the required parts before and after "_"
                        read_id, number = id_value.split("_")
                        # Write the extracted information to the output file
                        read_id_data[target_id] = f"read_{read_id}"
                        break  # Exit loop once ID information is found

    with open(output_file_name, 'w') as output_file:
        for key,value in read_id_data.items():
            output_file.write(f"{key}\t{value}\n")

    print(f"{output_file_name} is generated with all the unique read ids.")
    return read_id_data


def generate_hmmout_file(marker_file, gene_file):
    #hmmsearch --domtblout output_filename.hmmout marker_genes/bacteria.hmm reads.faa
    hmmout_file = f"{gene_file}_{marker_file}"
    hmmCmd = (
                "hmmsearch "
                + " --domtblout "
		+ hmmout_file
                + ".hmmout "
		+ marker_file
                + " "
                + gene_file
                + ".faa"
                
            )

    logger.debug(f"exec cmd: {hmmCmd}")
    os.system(hmmCmd)
    
    return f"{hmmout_file}.hmmout"

def process_hmmout_file(file_path, threshold, read_id_dict):
    """
    Process the .hmmout file and get unique reads, their most matching marker gene, and score
    :param file_path: Path to the .hmmout file
    :param threshold: Threshold value for score
    :return: Dictionary containing unique reads, their most matching marker gene, and score
    """
    # Dictionary to store data for each target read
    score_data = {}

    try:
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('#'):
                    continue
                
                data = line.split()
                target_read = data[0].split("_")
                
                # Extracting target ID, calculating score, and retrieving marker gene
                target_id = target_read[0]
                temp_score = (int(data[18]) - int(data[17]) + 1 )/int(data[5])*100
                temp_marker_gene = data[3]
                
                read_id = read_id_dict[target_id]
                
                # Update data if target_id exists, otherwise add new entry if score is above threshold
                if target_id in score_data:
                    if score_data[target_id]['score'] < temp_score:
                        score_data[target_id]['score'] = temp_score
                        score_data[target_id]['marker_gene'] = temp_marker_gene
                        score_data[target_id]['read_id'] = read_id
                else:
                    if temp_score >= threshold:
                        score_data[target_id] = {'read_id': read_id,'marker_gene': temp_marker_gene, 'score': temp_score}
                        
    except FileNotFoundError:
        print("File not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
    
    # Write the processed data to a text file
    output_file_path = f"{file_path}_marker_scores.txt"
    with open(output_file_path, 'w') as output_file:
        for key,value in score_data.items():
            output_file.write(f"{key}\t{value['read_id'].ljust(20)}\t{value['marker_gene'].ljust(20)}\t{value['score']}\n")

    
    return f"Processed data has been written to {output_file_path}"




# Define the main function to call other functions
def main():
    if len(sys.argv) < 3:
        print("Usage: python script.py <reads.fastq> <marker.hmm>")
        sys.exit(1)

    # Get the input from the command line argument
    input_reads_file = sys.argv[1]
    input_marker_file = sys.argv[2]
    
    fasta_file = convert_fastq_fasta(input_reads_file)
    gene_file = generate_gene_faa_file(fasta_file)
    read_id_mapping_dict = get_id_map(gene_file)
    final_hmmout_file = generate_hmmout_file(input_marker_file, gene_file)
    process_hmmout_file(final_hmmout_file,50, read_id_mapping_dict)
    
# Call the main function if this script is run directly
if __name__ == "__main__":
    main()
