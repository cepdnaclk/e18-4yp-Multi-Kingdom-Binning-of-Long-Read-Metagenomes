def compare_find_bact_arch(result_file, genes_file_bact, genes_file_arch, output_file):
  """
  Compares accession numbers from result.tab with genes.txt and adds a "bacteria" column to results.tab

  Args:
      result_file (str): Path to the result.tab file.
      genes_file (str): Path to the genes.txt file.
      output_file (str): Path to the output file (modified result.tab).
  """

  # Read gene data into a dictionary for efficient lookup
  genes_dict_arch = {}
  genes_dict_bact = {}
  count_arch = 0
  count_bact = 0
  with open(genes_file_bact, 'r') as f:
    for line in f:
      gene_name, accession, _ = line.strip().split('\t')
      accession = accession.split(".")[0]
      genes_dict_bact[accession] = "bacteria"

  with open(genes_file_arch, 'r') as f:
    for line in f:
      gene_name, accession, _ = line.strip().split('\t')
      accession = accession.split(".")[0]
      genes_dict_arch[accession] = "archaea"

  with open(result_file, 'r') as f_in, open(output_file, 'w') as f_out:
    for line in f_in:
      identifier, accession_1, *_ = line.strip().split('\t')
      accession_1 = accession_1.split(".")[0]
      if accession_1 in genes_dict_bact:
        source = genes_dict_bact[accession_1]
        count_bact = count_bact + 1
      elif accession_1 in genes_dict_arch:
        source = genes_dict_arch[accession_1]
        count_arch = count_arch + 1
      else:
        source = "not_found"

      f_out.write(line.strip() + '\t' + source + '\n')
  return [count_bact, count_arch]

# Specify file paths
result_file = "resultDB.tab"
genes_file_arch = "genes_arch.txt"
genes_file_bact = "genes_bac.txt"
output_file = "result_with_bacteria.tab"
genes_dict_file = "genes_dict.json"
counts = compare_find_bact_arch(result_file, genes_file_bact, genes_file_arch, output_file)

print(f"Bacteria: {counts[0]}, Archaea: {counts[1]}")
print(f"Results processed and written to: {output_file}")

