input_file_path = "../../../Results/zymo10_1000kto1100k/mapping/all_reads.output"
output_file_path = "../../../Results/zymo10_1000kto1100k/mapping/delete.tsv"
count = 0
tst = []
# Open input and output files
with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
    # Iterate over each line in the input file
    bins = set()
    for line in input_file:
        # Remove the identifier at the beginning and the number at the end of each line
        # if (len(line.split('\t')) > 2):
        line = line.strip()
        read_id=line.split('\t')[0]
        modified_line = line.split('\t')[1]  # Assuming tab-separated data
        if (modified_line == 'Staphylococcus_aureus_plasmid1' or modified_line == 'Staphylococcus_aureus_plasmid2' or modified_line == 'Staphylococcus_aureus_plasmid3'):
            modified_line = 'Staphylococcus_aureus_chromosome'
        elif (modified_line == 'BS.pilon.polished.v3.ST170922'):
            modified_line = 'Bacillus_subtilis_complete_genome'
        elif ('tig' in modified_line):
            modified_line = 'Saccharomyces_cerevisiae_draft_genome'
        elif (modified_line == 'Escherichia_coli_plasmid'):
            modified_line = 'Escherichia_coli_chromosome'
        
        if (modified_line != 'None' and modified_line != 'POOR MAPPING'):
            bins.add(modified_line)
            output_file.write(f"{read_id}\t{modified_line}\n")
        else:
            count += 1
            tst.append(read_id)

print("Processing completed. Modified content saved to:", output_file_path)
print(tst)
print(count)