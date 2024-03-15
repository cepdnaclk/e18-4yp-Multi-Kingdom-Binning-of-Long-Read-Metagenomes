# Get contigs containing marker genes
def get_reads_with_marker_genes(
    reads_and_scores_dict, hmmout_file, position
):

    with open(f"{hmmout_file}", "r") as file:
        for line in file.readlines():
            if not line.startswith("#"):
                strings = line.strip().split()
                sequence = strings[0]

                seq_id = sequence.split("_")
                seq_name = seq_name[: len(seq_name) - 1]
                read_id = "_".join(seq_name)

                score = float(strings[7])

                if (read_id not in reads_and_scores_dict):
                    reads_and_scores_dict[read_id] = [0,0,0]
                    reads_and_scores_dict[read_id][position] = score
                    
                else:
                    if reads_and_scores_dict[read_id][position] < score:
                        reads_and_scores_dict[read_id][position] = score

    return reads_and_scores_dict

# hmmsearch --domtblout domtblout.hmmout marker\ genes/bacteria.hmm reads.faa 1> domtblout.out 2> domtblout.err

def analyze_hmmout_file(previous_result, hmmout_file, position):

    reads_and_scores_dict = get_reads_with_marker_genes(previous_result, hmmout_file, position)
    return reads_and_scores_dict