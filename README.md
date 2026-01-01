This script was originally designed to analyze possible LGTs from prokaryotes to eukaryotes. This new version (v2) supports multi-file input and allows to analyze LGTs on different taxonomic levels through the optional --intaxon argument.

```
usage:  LGTanalysis (https://github.com/vojtech-zarsky/lgt_analysis)
Prerequisities:
 ete2 or ete3

Syntax:
python LGTanalysis2.py --tree_files <file1>,<file2>,...
  --query_sequence_ids <id1>,<id2>,...
  --taxonomy_mappings <table1>,<table2>,... OR (--tax_id_mappings <table1>,<table2>,... --taxdump_dir <dir>)
  --taxonomy_of_selected_groups <table>
  --intaxon eukaryota

       [-h] --tree_files TREE_FILES --query_sequence_ids QUERY_SEQUENCE_IDS
       [--taxonomy_mappings TAXONOMY_MAPPINGS]
       [--tax_id_mappings TAX_ID_MAPPINGS] [--taxdump_dir TAXDUMP_DIR]
       --taxonomy_of_selected_groups TAXONOMY_OF_SELECTED_GROUPS
       [--intaxon INTAXON] [--mask_viruses] [--support_cutoff SUPPORT_CUTOFF]
       [--orthology_cutoff ORTHOLOGY_CUTOFF] [--output_suffix OUTPUT_SUFFIX]
       [--gain_weight GAIN_WEIGHT] [--loss_weight LOSS_WEIGHT]

optional arguments:
  -h, --help            show this help message and exit
  --tree_files TREE_FILES
                        A single file name or a comma-separated list of tree
                        files. (e.g. IQ-Tree *.treefile)
  --query_sequence_ids QUERY_SEQUENCE_IDS
                        A single ID or a comma-separated list of IDs of the
                        query sequence. (In the same order as the respective
                        tree_files.)
  --taxonomy_mappings TAXONOMY_MAPPINGS
                        A table or a comma-separated list of tables of
                        sequence IDs with their taxonomies. (e.g. A2FH21_TRIVA
                        <tab>Eukaryota<tab>Metamonada<tab>Parabasalia<tab>Tric
                        homonadida) (In the same order as the respective
                        tree_files.)
  --tax_id_mappings TAX_ID_MAPPINGS
                        Alternatively a table or a comma-separated list of
                        tables of sequence IDs with their NCBI taxids can be
                        provided. (e.g. A2FH21_TRIVA<tab>5722) In that case
                        the Taxonomy.py (https://github.com/vojtech-
                        zarsky/lgt_analysis-tools) file is needed and the
                        --taxdump_dir directory must be specified. You could
                        also specify taxon name, however taxids are better.
                        (In the same order as the respective tree_files.)
  --taxdump_dir TAXDUMP_DIR
                        Directory to the NCBI taxdump folder which can be
                        downloaded here:
                        https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
  --taxonomy_of_selected_groups TAXONOMY_OF_SELECTED_GROUPS
                        A tab-separated table of selected eukaryotic and
                        prokaryotic groups with their taxonomy. (e.g.
                        https://github.com/vojtech-
                        zarsky/lgt_analysis/LGTanalysis.groups.tsv) An arrow
                        indicates inclusion of the taxon on the left site to
                        the group on the right. The domain-level taxon should
                        be first.
  --intaxon INTAXON     Only consider LGTs from the outside of this taxon.
                        Default: The domain-level taxon of the query (e.g.
                        Euakryota).
  --mask_viruses        This tag causes to ignore viral sequences when
                        evaluating "out taxons".
  --support_cutoff SUPPORT_CUTOFF
                        Nodes with support bellow this will be removed (mind
                        e.g. 0.85 vs 85 support). Default: 0.9/90
  --orthology_cutoff ORTHOLOGY_CUTOFF
                        Cutoff of the orthology score. Set between 0.2 (more
                        inclusive) - 0.3 (a little less inclusive).
  --output_suffix OUTPUT_SUFFIX
                        (Default: *.lgt_analysis_<intaxon>.txt)
  --gain_weight GAIN_WEIGHT
  --loss_weight LOSS_WEIGHT
```
