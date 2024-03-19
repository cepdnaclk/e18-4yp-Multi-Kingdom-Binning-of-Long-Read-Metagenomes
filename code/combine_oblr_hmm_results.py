bins_file = "oblr/bins.txt"
gene_file = "prodigal_hmmer/hmmout_result.txt"
output = "final.txt"

def write_dict_to_file(dictionary, file_path):
    with open(file_path, "w") as file:

        for key, value in dictionary.items():
            line = f"{key}"
            for val in value:
                line += f"\t{val}" 

            line += "\n"
            file.write(line)

def bin_based_files(dictionary):
    bins = {1: [], 2: [], 3: []}

    for key, value in dictionary.items():
        bin_num = float(value[0])
        if bin_num in bins:
            
            line = f"{key}"
            for val in value:
                line += f"\t{val}" 

            # line += "\n"
            bins[bin_num].append(line)

    for bin_num, lines in bins.items():
        with open(f"bin{bin_num}.txt", "w") as file:
            file.write("\n".join(lines))


final_dict = {}

bins = open(bins_file, 'r')
gene = open(gene_file, 'r')

for line in bins.readlines():

    line = line.strip().split(" ")
    final_dict[line[0]]= []
    final_dict[line[0]].append(line[1])

kingdoms = gene.readline().strip().split("\t")[1:]


for line in gene.readlines():
    splitted = line.strip().split("\t")
    scores = [float(x) for x in splitted[1:4]]  # Convert scores to float
    max_score = max(scores)
    max_scored_indices = [i + 1 for i, x in enumerate(scores) if x == max_score]
    max_scored_kingdoms = [kingdoms[i-1] for i in max_scored_indices]
    final_dict[splitted[0]].extend(max_scored_kingdoms)


# write_dict_to_file(final_dict, output)


bin_based_files(final_dict)