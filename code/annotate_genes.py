import subprocess
import sys
import os
from analyze_hmmout import analyze_hmmout_file

def run_prodigal(fasta_file, faa_file):
    # Define the Prodigal command arguments
    prodigal_cmd = ["prodigal", "-i", fasta_file, "-a", faa_file]

    # Run the Prodigal command
    try:
        subprocess.run(prodigal_cmd, check=True)
        print("Prodigal executed successfully!")
    except subprocess.CalledProcessError as e:
        print(f"Prodigal failed with error code {e.returncode}")
        exit(1)

def run_hmmer(faa_file, marker_hmm_file, hmmout_file):
    # Define the HMMER (hmmsearch) command arguments
    hmmsearch_cmd = ["hmmsearch", "--domtblout",  hmmout_file, marker_hmm_file, faa_file]

    # Run the HMMER (hmmsearch) command
    try:
        subprocess.run(hmmsearch_cmd, check=True)
        print("hmmsearch executed successfully!")
    except subprocess.CalledProcessError as e:
        print(f"hmmsearch failed with error code {e.returncode}")
        exit(1)

def write_dict_to_file(dictionary, file_path, column_names):
    with open(file_path, "a+") as file:
        first_line = "Read_id"
        for col_name in column_names:
            first_line += f"\t{col_name}"
        
        file.write(first_line+"\n")

        for key, value in dictionary.items():
            line = f"{key}\t{value[0]}\t{value[1]}\t{value[2]}\n"
            file.write(line)


def main():

    # Check if the required arguments are provided
    if len(sys.argv) != 4:
        print("Incorrect number of agruments provided")
        print("Usage: python script.py [input_fasta_file] [folder_wiht_markers] [results_folder]")
        sys.exit(1)

    # Get the input file path, output file path, and the starting read index
    fasta_file = sys.argv[1]
    marker_hmm_folder = sys.argv[2]
    results_folder = sys.argv[3] + "/"

    result_file_name = "final_result.txt"

    # for all the files in marker_files folder, run the below code
    file_count = 0
    reads_scores_dict = {}
    kingdoms = []

    for filename in os.listdir(marker_hmm_folder):
        hmm_file_path = os.path.join(marker_hmm_folder, filename)
        path_for_results = results_folder + filename.split(".")[0]
        kingdoms.append(filename.split(".")[0])

        run_prodigal(fasta_file, path_for_results+".faa")

        run_hmmer(path_for_results+".faa", hmm_file_path, path_for_results+".hmmout")
        # run_prodigal_and_hmmer(fasta_file, results_folder+filename+".faa", hmm_file_path, results_folder+filename+".hmmout")

        reads_scores_dict = analyze_hmmout_file(reads_scores_dict, path_for_results+".hmmout", file_count, 0.5)
        file_count += 1

    write_dict_to_file(reads_scores_dict, results_folder+result_file_name, kingdoms)
    
main()