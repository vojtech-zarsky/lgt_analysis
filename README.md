Prerequisities: ete2 or ete3

Syntax example:
python LGTanalysis2.py --tree_files \<file1\>,\<file2\>,... --query_sequence_ids \<id1\>,\<id2\>,... 
  --taxonomy_mappings \<table1\>,\<table2\>,... OR (--tax_id_mappings \<table1\>,\<table2\>,... --taxdump_dir \<dir\>)
  --taxonomy_of_selected_groups \<table\> --intaxon eukaryota

This script was originally designed to analyze possible LGTs from prokaryotes to eukaryotes. This new version (v2) supports multi-file input and allows to analyze LGTs on different taxonomic levels through the optional --intaxon argument.
