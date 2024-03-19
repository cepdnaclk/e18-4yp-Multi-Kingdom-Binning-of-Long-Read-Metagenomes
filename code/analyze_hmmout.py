# Get contigs containing marker genes
def get_reads_with_marker_genes(
    reads_and_scores_dict, hmmout_file, position, mg_overlap_threshold
):

    with open(f"{hmmout_file}", "r") as file:
        for line in file.readlines():
            if not line.startswith("#"):
                strings = line.strip().split()
                sequence = strings[0]

                # extract the read_id
                seq_name = sequence.split("_")
                seq_name = seq_name[: len(seq_name) - 1]
                read_id = "_".join(seq_name)

                # Marker gene length
                marker_gene_length = int(strings[5])

                # Mapped marker gene length
                mapped_marker_length = int(strings[16]) - int(strings[15])

                score = mapped_marker_length/marker_gene_length

                if (score > mg_overlap_threshold):

                    if (read_id not in reads_and_scores_dict):
                        reads_and_scores_dict[read_id] = [0,0,0]
                        reads_and_scores_dict[read_id][position] = score
                    
                    else:
                        if reads_and_scores_dict[read_id][position] < score:
                            reads_and_scores_dict[read_id][position] = score

    return reads_and_scores_dict

# hmmsearch --domtblout output_filename.hmmout marker_genes/bacteria.hmm reads.faa

def analyze_hmmout_file(previous_result, hmmout_file, position):

    mg_overlap_threshold = 0.5

    reads_and_scores_dict = get_reads_with_marker_genes(previous_result, hmmout_file, position, mg_overlap_threshold)
    return reads_and_scores_dict