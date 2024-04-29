import numpy as np
import os
import sys


# ----------------------------------------------------
# these should be taken as inputs

# Check if the correct number of arguments are provided
if len(sys.argv) != 4:
    print("Usage: python script.py oblr_results_folder output_folder marker_file_path")
    sys.exit(1)

# Take paths from command line arguments
oblr_results = sys.argv[1]
output = sys.argv[2]
input_marker_file = sys.argv[3]

# ----------------------------------------------------


# files
fasta_file = oblr_results + '/reads.fasta'

def generate_gene_faa_file(fasta_file):
    # here, fasta_file means the path to the file   
    faa_file = f"{output}/reads_genes"

    fragCmd = (
                "prodigal"
                + " -i "
                + fasta_file
                + " -a "
                + faa_file
                + ".faa"
            )
    
    # logger.debug(f"exec cmd: {fragCmd}")
    return_code = os.system(fragCmd)
    
    return return_code, faa_file # returning the file name: string

def generate_hmmout_file(marker_file, faa_file):
    # marker_file is the file paths
    
    #hmmsearch --domtblout output_filename.hmmout marker_genes/bacteria.hmm reads.faa
    hmmout_file = f"{output}/reads_markers"
    hmmCmd = (
                "hmmsearch "
                + " --domtblout "
		+ hmmout_file
                + ".hmmout "
		+ marker_file
                + " "
                + faa_file
  
                
            )

    # logger.debug(f"exec cmd: {hmmCmd}")
    return_code = os.system(hmmCmd)
    return return_code, f"{hmmout_file}.hmmout"

def process_hmmout_file(hmmout_filepath, threshold):
    """
    Process the .hmmout file and get unique reads, their most matching marker gene, and score
    :param file_path: Path to the .hmmout file
    :param threshold: Threshold value for score
    :return: Dictionary containing unique reads, their most matching marker gene, and score
    """
    # Dictionary to store data for each target read
    score_data = {}

    try:
        with open(hmmout_filepath, 'r') as file:
            for line in file:
                if line.startswith('#'):
                    continue
                
                data = line.split()
                target_read = data[0].rsplit("_", 1)
                
                # Extracting target ID, calculating score, and retrieving marker gene
                target_id = target_read[0]
                temp_score = (int(data[16]) - int(data[15]) + 1 )/int(data[5])*100
                temp_marker_gene = data[3]
                
                # read_id = read_id_dict[target_id]
                
                # Update data if target_id exists, otherwise add new entry if score is above threshold
                if target_id in score_data:
                    if score_data[target_id]['score'] < temp_score:
                        score_data[target_id]['score'] = temp_score
                        score_data[target_id]['marker_gene'] = temp_marker_gene
                        score_data[target_id]['read_id'] = target_id
                else:
                    if temp_score >= threshold:
                        score_data[target_id] = {'read_id': target_id,'marker_gene': temp_marker_gene, 'score': temp_score}
                        
    except FileNotFoundError:
        print("HMMOUT file not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
    
    # Write the processed data to a text file
    output_file_path = f"{output}/marker_scores.txt"
    with open(output_file_path, 'w') as output_file:
        for key,value in score_data.items():
            output_file.write(f"{key}\t{value['read_id'].ljust(20)}\t{value['marker_gene'].ljust(20)}\t{value['score']}\n")

    return f"Processed data has been written to {output_file_path}"


# ----------------------------------------------------------------------------------------------------------------------------

# generate .faa file
# return_code1, faa_file = generate_gene_faa_file(fasta_file)

faa_file = output + "/reads_genes.faa"

# Check if the file exists
if os.path.exists(faa_file):
  # Check if the file is not empty
  if os.path.getsize(faa_file) > 0:
    return_code1 = 0
    print("faa file already created")
  else:
    print("creating faa file from fasta file")
    return_code1, faa_file = generate_gene_faa_file(fasta_file)
else:
    print("creating faa file from fasta file")
    return_code1, faa_file = generate_gene_faa_file(fasta_file)
        

if return_code1 == 0:
    # faa_file = output + "/reads_genes"
    # generate hmmout file
    return_code2, hmmout_file = generate_hmmout_file(input_marker_file, faa_file)

    if return_code2 == 0:
        # proecss hmmout file and get marker genes of each read
        threshold = 50
        try:
            output = process_hmmout_file(hmmout_file, threshold)
            print(output)
        except Exception as e:
            print(f"Error: {e}")
    else:
        print ("Error occured in HMMSEARCH command. Please try again")

else:
    print ("Error occured in PRODIGAL command. Please try again")