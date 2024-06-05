import logging
import os
import sys
from Bio import SeqIO


# create logger
logger = logging.getLogger("== analyze_marker_genes.py ==")

def convert_fastq_fasta(input_file):
    
    output_fasta_file_name = "reads"
    output_file = f"{output_fasta_file_name}.fasta"
    
    with open(input_file, "r") as input_handle:
        with open(output_file, "w") as output_handle:
            sequences = SeqIO.parse(input_handle, "fastq")
            count = SeqIO.write(sequences, output_handle, "fasta")
    print("Converted %i records" % count)
    
    return output_file

def generate_gene_faa_file(fasta_file):
    # Define the command as a list (safer than a string)
    gene_file = f"fasta_reads_genes.faa"
    fragCmd = (
                "prodigal"
                + " -i "
                + fasta_file
                + " -a "
                + gene_file
            )
    logger.debug(f"exec cmd: {fragCmd}")
    os.system(fragCmd)
    
    return gene_file #returning the file name: string
    
# ------------------------------------------------HMMSEARCH-----------------------------------------

def generate_hmmout_file(marker_file, gene_file):
    #hmmsearch --domtblout output_filename.hmmout marker_genes/bacteria.hmm reads.faa
    hmmout_file = f"hmmsearch_output.hmmout"
    hmmCmd = (
                "hmmsearch "
                + " --domtblout "
		+ hmmout_file
                + " "
		+ marker_file
                + " "
                + gene_file
                
            )

    logger.debug(f"exec cmd: {hmmCmd}")
    os.system(hmmCmd)
    
    return hmmout_file

# ---------------------------------------------MMSEQS2----------------------------------------------


def generate_query_db(gene_file):
    
    query_db = f"queryDB"
    
    query_db_cmd = ("mmseqs createdb "+gene_file + " " + query_db)
    logger.debug(f"exec cmd: {query_db_cmd}")
    os.system(query_db_cmd)

    return query_db

def generate_mmseqs_result(query_db, target_db_list):
    
    result_db = "resultdb"
    
    create_output_dir = ("mkdir mmseqs2_files")
    logger.debug(f"exec cmd: {create_output_dir}")
    os.system(create_output_dir)
    count = 0
    for target_db in target_db_list:
      target_db = str(target_db)
      final_file_name = f"{result_db}_{count}"
      
      # usage of thread is limited to 8 threads
      # special note: --exhaustive-search should be used for PHROGS
      mmseq_search_cmd = ("mmseqs " + "search " + query_db + " " + target_db + " mmseqs2_files/" + final_file_name + " tmp_"+final_file_name + " --threads 8")
      logger.debug(f"exec cmd: {mmseq_search_cmd}")
      os.system(mmseq_search_cmd)
      
      final_db_file = ("ls mmseqs2_files/"+ final_file_name +"* | grep -v \'" + final_file_name +".dbtype "+"\\|"+ final_file_name +".index\'"+" > mmseqs2_files/"+final_file_name) #renaming with _
      logger.debug(f"exec cmd: {final_db_file}")
      os.system(final_db_file)    
    
      convertalist = ("mmseqs " + "convertalis " + query_db + " " + target_db +  " mmseqs2_files/" + final_file_name + " mmseqs2_files/" + final_file_name + ".tab --format-output query,target,evalue,gapopen,pident,fident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen")
      logger.debug(f"exec cmd: {convertalist}")
      os.system(convertalist)
      count = count + 1
      
    final_merge = ("cat mmseqs2_files/*.tab > mmseqs2_files/" +  result_db + ".tab")
    logger.debug(f"exec cmd: {final_merge}")
    os.system(final_merge)
      
    remove_tmp = ("rm -r "+ "tmp*")
    logger.debug(f"exec cmd: {remove_tmp}")
    os.system(remove_tmp)
    
    return f"mmseqs2_files/{result_db}.tab"




def process_hmmout_file(file_path,score_data, threshold):
    """
    Process the .hmmout file and get unique reads, their most matching marker gene, and score
    :param file_path: Path to the .hmmout file
    :param threshold: Threshold value for score
    :return: Dictionary containing unique reads, their most matching marker gene, and score
    """
    # Dictionary to store data for each target read
    # score_data = {}

    try:
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('#'):
                    continue
                
                data = line.split()
                
                read_id = "_".join(data[0].split("_")[:-1])
                
                # Extracting target ID, calculating score, and retrieving marker gene
                
                temp_score = (int(data[16]) - int(data[15]) + 1 )/int(data[5])*100
                temp_marker_gene = data[3]
                
                # read_id = read_id_dict[target_id]
                
                # Update data if target_id exists, otherwise add new entry if score is above threshold
                if read_id in score_data:
                    if score_data[read_id]['score'] < temp_score:
                        score_data[read_id]['score'] = temp_score
                        score_data[read_id]['marker_gene'] = temp_marker_gene
                        
                else:
                    if temp_score >= threshold:
                        score_data[read_id] = {'marker_gene': temp_marker_gene, 'score': temp_score}
                        
    except FileNotFoundError:
        print("File not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
    
    # # Write the processed data to a text file
    # output_file_path = f"hmmsearch_marker_scores.txt"
    # with open(output_file_path, 'w') as output_file:
    #     for key,value in score_data.items():
    #         output_file.write(f"{key.ljust(20)}\t{value['marker_gene'].ljust(20)}\t{value['score']}\n")

    
    return f"Dictionary updated with hmmout file results"


def process_mmseqs_file(filename, score_data, threshold):
    try:
        with open(filename, 'r') as f:
            for line in f:
                row = line.strip().split('\t')
                
                read_id = "_".join(row[0].split("_")[:-1])
                
                temp_score = 100 * (float(row[11]) - float(row[10])) / float(row[12])
                
                temp_marker_gene = row[1]
                
                if read_id in score_data:
                    if score_data[read_id]['score'] < temp_score:
                        score_data[read_id]['score'] = temp_score
                        score_data[read_id]['marker_gene'] = temp_marker_gene
                        
                else:
                    if temp_score >= threshold:
                        score_data[read_id] = {'marker_gene': temp_marker_gene, 'score': temp_score}

    except FileNotFoundError:
        print("File not found.")
    except Exception as e:
        print(f"An error occurred: {e}")


    #   with open(output_filename, 'w') as f:
    #     for row in filtered_data:
    #       f.write('\t'.join(map(str, row)) + '\n')

    print(f"Dictionary updated with mmseqs .tab file results")



# Define the main function to call other functions
def main():
    if len(sys.argv) < 3:
        print("Usage: python script.py <reads.fastq | .fasta | .faa> <marker.hmm> <target db list>")
        sys.exit(1)

    # Get the input from the command line argument
    read_file = sys.argv[1]
    input_marker_file = sys.argv[2]
    target_db = sys.argv[3:]
    
    gene_file = ""
    
    read_file_type = read_file.split(".")[-1]
    if read_file_type == "fastq":
      #fastq file is given and need to do the conversion
      read_file = convert_fastq_fasta(read_file)
      gene_file = generate_gene_faa_file(read_file)
    elif read_file_type == "fasta":
      # fasta file is given
      gene_file = generate_gene_faa_file(read_file)
    elif read_file_type == "faa":
      # directly .faa file is given
      gene_file = read_file 
    else:
      print("Entered file format is incorrect. Only .fasta, .fastq and .faa are valid")
      return
    
    # creating the DB for the query reads
    query_db = generate_query_db(gene_file)
    
    # get the tab version representation
    resultDBtab = generate_mmseqs_result(query_db, target_db)
    
    final_hmmout_file = generate_hmmout_file(input_marker_file, gene_file)
    
    score_data = {}
    # threshold = 50
    
    process_hmmout_file(final_hmmout_file, score_data, 50)
    
    process_mmseqs_file(resultDBtab, score_data, 0)
    
    output_file_path = f"final_marker_scores.txt"
    with open(output_file_path, 'w') as output_file:
        for key,value in score_data.items():
            output_file.write(f"{key.ljust(20)}\t{value['marker_gene'].ljust(20)}\t{value['score']}\n")
    
    
if __name__ == "__main__":
    main()
