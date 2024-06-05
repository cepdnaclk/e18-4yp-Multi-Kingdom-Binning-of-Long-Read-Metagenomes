step1 --> identify misbinned nodes
python step1_oblr.py -o <path-to-output-folder>/ -d <path-to-initial-tool-results-folder>/

step2 --> scan for markers
python step2_analyze_marker_genes.py <path-to-fasta-read-file> <path-to-hmm-file> <path-to-vogdb-profile-db> <path-to-phrog-profile-db> <path-to-fungi-db>

step3 --> relabel using markers
python step3_relabelling_using_MG.py <path-to-initial-tool-results-folder>/ <path-to-output-folder>/

step4 --> label propagation
python step4_label_prop.py -o <path-to-output-folder>/ -d <path-to-initial-tool-results-folder>/ -b <name-of-initial-tool>




