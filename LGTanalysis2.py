
def parseArguments():
    parser = argparse.ArgumentParser(prog=' LGTanalysis (https://github.com/vojtech-zarsky/vojta-tools)'\
                                     '\nPrerequisities:\n ete2 or ete3\n'\
                                     '\nSyntax:\npython LGTanalysis2.py --tree_files <file1>,<file2>,... --query_sequence_ids <id1>,<id2>,... --taxonomy_mappings <table1>,<table2>,... OR (--tax_id_mappings <table1>,<table2>,... --taxdump_dir <dir>) --taxonomy_of_selected_groups <table> --intaxon eukaryota\n'\
                                    )
    parser.add_argument('--tree_files', required=True, help='A single file name or a comma-separated list of tree files. (e.g. IQ-Tree *.treefile)')
    parser.add_argument('--query_sequence_ids', required=True, help='A single ID or a comma-separated list of IDs of the query sequence. (In the same order as the respective tree_files.)')
    parser.add_argument('--taxonomy_mappings', required=False, help='A table or a comma-separated list of tables of sequence IDs with their taxonomies. (e.g. A2FH21_TRIVA<tab>Eukaryota<tab>Metamonada<tab>Parabasalia<tab>Trichomonadida) (In the same order as the respective tree_files.)')
    parser.add_argument('--tax_id_mappings', required=False, help='Alternatively a table or a comma-separated list of tables of sequence IDs with their NCBI taxids can be provided. (e.g. A2FH21_TRIVA<tab>5722) In that case the Taxonomy.py (https://github.com/vojtech-zarsky/vojta-tools) file is needed and the --taxdump_dir directory must be specified. You could also specify taxon name, however taxids are better. (In the same order as the respective tree_files.)')
    parser.add_argument('--taxdump_dir', required=False, help='Directory to the NCBI taxdump folder which can be downloaded here: https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz')
    parser.add_argument('--taxonomy_of_selected_groups', required=True, help='A tab-separated table of selected eukaryotic and prokaryotic groups with their taxonomy. (e.g. https://github.com/vojtech-zarsky/vojta-tools/blob/master/LGTanalysis.groups.tsv) An arrow indicates inclusion of the taxon on the left site to the group on the right. The domain-level taxon should be first.')
    parser.add_argument('--intaxon', required=False, help='Only consider LGTs from the outside of this taxon. Default: The domain-level taxon of the query (e.g. Euakryota).')
    parser.add_argument('--mask_viruses', action='store_true', help='This tag causes to ignore viral sequences when evaluating "out taxons".')
    parser.add_argument('--support_cutoff', type=float, required=False, help='Nodes with support bellow this will be removed (mind e.g. 0.85 vs 85 support). Default: 0.9/90')
    parser.add_argument('--orthology_cutoff', type=float, default=0.24, help='Cutoff of the orthology score. Set between 0.2 (more inclusive) - 0.3 (a little less inclusive).')
    parser.add_argument('--output_suffix', required=False, help='(Default: *.lgt_analysis_<intaxon>.txt)')
    parser.add_argument('--gain_weight', type=float, default=1, help='')
    parser.add_argument('--loss_weight', type=float, default=1, help='')
    args = parser.parse_args()
    if args.taxonomy_mappings==None and (args.tax_id_mappings==None or args.taxdump_dir==None):
        print('Error: At least one way to assign taxonomies to leaves is required!')
        raise
    args.source_cutoff = 0.5 ##
    args.tree_files = args.tree_files.split(',')
    args.query_sequence_ids = args.query_sequence_ids.split(',')
    if len(args.query_sequence_ids) != len(args.tree_files):
        print('Error: Multi-file inputs need to contain the same number of trees, query_sequence_ids and taxonomy-mapping tables!')
        raise
    if args.taxonomy_mappings != None:
        args.taxonomy_mappings = args.taxonomy_mappings.split(',')
        if len(args.taxonomy_mappings) != len(args.tree_files):
            print('Error: Multi-file inputs need to contain the same number of trees, query_sequence_ids and taxonomy-mapping tables!')
            raise
    if args.tax_id_mappings != None:
        args.tax_id_mappings = args.tax_id_mappings.split(',')
        if len(args.tax_id_mappings) != len(args.tree_files):
            print('Error: Multi-file inputs need to contain the same number of trees, query_sequence_ids and taxonomy-mapping tables!')
            raise
    return args

def removeZeroSisterNodes(tree, add_branch_lengths=True):
    for node in tree.traverse():
        ancestors = node.get_ancestors()
        if len(ancestors) > 0 and len(ancestors[0].get_children()) == 1:
            ancestor = ancestors[0]
            for child in node.get_children():
                if add_branch_lengths:
                    child.dist += node.dist
                ancestor.add_child(child)
            ancestor.remove_child(node)
    return tree

def readTaxonomy(taxonomy_table, skip_viruses=True, limit_to=None):
    groups = []
    taxon2group = {}
    group2taxonomy = {}
    tree = Tree()
    tree.name = "root"
    for line in open(taxonomy_table):
        taxonomy = list(map(lambda x:x.strip(), line.strip().split('\t')))
        if skip_viruses and 'viruses' in taxonomy:
            continue
        if limit_to != None:
            limit_to_check = False
            for taxon in taxonomy:
                if taxon.lower() == limit_to.lower():
                    limit_to_check = True
            if not limit_to_check:
                continue
        if '->' in taxonomy[-1]:
            (taxon, group) = list(map(lambda x:x.strip(), taxonomy[-1].split('->')))
            taxon2group[taxon] = group
            del taxonomy[-1]
        else:
            taxon2group[taxonomy[-1]] = taxonomy[-1]
        group = taxonomy[-1]
        if group not in group2taxonomy:
            groups.append(group)
            group2taxonomy[group] = taxonomy
            if True:
                node = tree
                for tax_level in range(len(taxonomy)):
                    child_check = False
                    for child in node.children:
                        if child.name == taxonomy[tax_level]:
                            node = child
                            child_check = True
                            break
                    if not child_check:
                        child_new = Tree()
                        child_new.name = taxonomy[tax_level]
                        node.add_child(child_new)
                        node = child_new
    tree = removeZeroSisterNodes(tree, add_branch_lengths=False)
    return (groups, group2taxonomy, taxon2group, tree)

def removeUnsupported(tree, node_support_cutoff):
    for node in tree.traverse():
        if node == tree or node.is_leaf():
            continue
        if node.support < node_support_cutoff:
            ancestor = node.get_ancestors()[0]
            for child in node.get_children():
                child.dist += node.dist
                ancestor.add_child(child)
            ancestor.remove_child(node)
    return tree

def createPseudonodes(node):
    if node.is_leaf():
        return node
    for child in node.get_children():
        createPseudonodes(child)
    if len(node.get_children()) > 2:
        dDominantTaxon2Children = {}
        for child in node.get_children():
            sDominantTaxon = 0
            if 1 in child.taxonomy and child.taxonomy[1] >= 0.5:
                sDominantTaxon = 1
            if sDominantTaxon not in dDominantTaxon2Children:
                dDominantTaxon2Children[sDominantTaxon] = []
            dDominantTaxon2Children[sDominantTaxon].append(child)
        if len(dDominantTaxon2Children) > 1:
            for (sDominantTaxon, lDominantTaxonChildren) in dDominantTaxon2Children.items():
                if len(lDominantTaxonChildren) > 1:
                    newChild = Tree()
                    newChild.dist = min(list(map(lambda x:x.dist, node.get_children())))/2.
                    for child in lDominantTaxonChildren:
                        child.dist -= newChild.dist
                        newChild.add_child(child.copy())
                        node.remove_child(child)
                    node.add_child(newChild)
    return node

def setTaxonomy(node, tax_levels):
    if node.is_leaf():
        return node
    node.taxonomy = []
    for tax_level in range(tax_levels):
        node.taxonomy.append({})
    children = node.get_children()
    children_tax_freq_sum = 0 ## typically the number of children, but some taxons may be masked
    for child in children:
        setTaxonomy(child, tax_levels)
        for (taxon, freq) in child.taxonomy[0].items():
            children_tax_freq_sum += freq
    for child in children:
        for tax_level in range(tax_levels):
            for (taxon, freq) in child.taxonomy[tax_level].items():
                if taxon not in node.taxonomy[tax_level]:
                    node.taxonomy[tax_level][taxon] = 0
                if children_tax_freq_sum > 0:
                    node.taxonomy[tax_level][taxon] += float(freq)/children_tax_freq_sum
    return node

## create pseudonodes in multifurcations connecting children with dominant 'eukaryota' or 'prokaryota' taxonomy respectively
def createPseudonodes(node, taxlvl=0):
    if node.is_leaf():
        return node
    for child in node.get_children():
        createPseudonodes(child)
    if len(node.get_children()) > 2:
        dominant_taxon2children = {}
        for child in node.get_children():
            dominant_taxon = 0
            if 1 in child.taxonomy[taxlvl] and child.taxonomy[taxlvl][1] >= 0.5:
                dominant_taxon = 1
            if dominant_taxon not in dominant_taxon2children:
                dominant_taxon2children[dominant_taxon] = []
            dominant_taxon2children[dominant_taxon].append(child)
        if len(dominant_taxon2children) > 1:
            for (dominant_taxon, dominant_taxon_children) in dominant_taxon2children.items():
                if len(dominant_taxon_children) > 1:
                    new_child = Tree()
                    new_child.dist = min(list(map(lambda x:x.dist, node.get_children())))/2.
                    new_child.support = 0
                    for child in dominant_taxon_children:
                        child.dist -= new_child.dist
                        new_child.add_child(child.copy())
                        node.remove_child(child)
                    node.add_child(new_child)
    return node

def countSumPresence(node):
    node.presence_sum = 0
    if node.is_leaf():
        if node.presence:
            node.presence_sum = 1 
    else:
        for child in node.get_children():
            countSumPresence(child)
            node.presence_sum += child.presence_sum

def findLosses(tree):
    leaves_used = set()
    for leaf in tree.get_leaves():
        if not leaf.presence and leaf.name not in leaves_used:
            last_node = leaf
            for ancestor in leaf.get_ancestors():
                if ancestor.presence_sum > 0:
                    break
                last_node = ancestor
            for leafTemp in last_node.get_leaves():
                leaves_used.add(leafTemp.name)
            last_node.loss = 1

def countSumLosses(node):
    node.losses_sum = 0
    if not hasattr(node, 'loss'):
        node.loss = 0
    node.losses_sum += node.loss
    if node.is_leaf():
        if node.loss == 1: ###
            node.losses_sum = 1 ###
        #pass
    else:
        for child in node.get_children():
            countSumLosses(child)
            node.losses_sum += child.losses_sum

def selectGains(tree):
    gains = []
    for node in tree.traverse():
        #node.gainScore = None
        if node.presence_sum > 0:
            if node.is_leaf():
                gains.append(node)
            else:
                count_children_with_presence = 0
                #iCountChildrenWithPresence = 0
                for child in node.get_children():
                    if child.presence_sum >= 1:
                        count_children_with_presence += 1
                if count_children_with_presence >= 2:
                    gains.append(node)
    return gains

def scoreGainCombos(gains, presence_sum, cutoff_score, gain_weight=1, loss_weight=1, indices=[], combo2score={}):
    presences = set()
    overlap_check = False
    losses_sum = 0
    for index in indices:
        node = gains[index]
        losses_sum += node.losses_sum
        for leaf in node.get_leaves():
            if leaf.presence:
                if leaf.name in presences:
                    overlap_check = True
                else:
                    presences.add(leaf.name)
        if overlap_check:
            break
    ## stop adding gain to the combo if the number of events get bigger than in the case of single gain at the root of eukaryotes
    if not overlap_check and losses_sum+len(indices) <= cutoff_score:
        if len(presences) == presence_sum:
            combo = tuple(sorted(list(map(lambda x:gains[x].name, indices))))
            ## the score of gain combination is the number of events - gains and losses
            ## I set them to an equal weight.
            combo2score[combo] = losses_sum*loss_weight + len(indices)*gain_weight
        for index in range(len(gains)):
            if len(indices) == 0 or index > indices[-1]:
                scoreGainCombos(gains, presence_sum, cutoff_score, gain_weight=gain_weight, loss_weight=loss_weight, indices=indices+[index], combo2score=combo2score)
    return combo2score

def find_lgts(args, tree_index):
    query_sequence_id = args.query_sequence_ids[tree_index]
    if args.output_suffix == None:
        out_file = open('{}.lgt_analysis_{}.txt'.format(args.tree_files[tree_index], args.intaxon), 'w')
    else:
        out_file = open('{}{}'.format(args.tree_files[tree_index], args.output_suffix), 'w')
    # get taxonomy for each sequence and assign it to one of the selected groups
    seqid2taxonomy = {}
    if args.taxonomy_mappings != None:
        for line in open(args.taxonomy_mappings[tree_index]):
            line = line.strip().split('\t')
            seqid = line[0]
            taxonomy = line[1:]
            seqid2taxonomy[seqid] = list(map(lambda x:x.lower(), taxonomy))
    else:
        for line in open(args.tax_id_mappings[tree_index]):
            if line.strip() == '':
                continue
            (seqid, taxid) = re.split(r'\s+', line.strip(), 1)
            taxonomy = None
            if taxid.strip().isdigit():
                taxonomy = taxonomy_parser.getTaxonomyByTaxId(taxid)
            else:
                taxonomy = taxonomy_parser.getTaxonomy(taxid)
            taxonomy = list(map(lambda x:x[1], taxonomy))[::-1]
            seqid2taxonomy[seqid] = list(map(lambda x:x.lower(), taxonomy))
    # mask viruses
    #if args.mask_viruses:
    #    for seqid, taxonomy in seqid2taxonomy.items():
    #        if 'viruses' in taxonomy[:2]:
    #            seqid2taxonomy[seqid] = []
    # define intaxon
    if args.intaxon == None:
        domains = ['archaea', 'bacteria', 'eukaryota', 'viruses']
        for taxon in seqid2taxonomy[query_sequence_id]:
            if taxon in domains:
                args.intaxon = taxon
    if args.intaxon == None:
        print('Error: Cannot identify the default intaxon!')
        raise
    args.intaxon = args.intaxon.lower()
    intaxon_level = seqid2taxonomy[query_sequence_id].index(args.intaxon)
    out_file.write('Query base group: {}\n'.format(args.intaxon))
    out_file.write('Query base group taxlevel: {}\n'.format(intaxon_level))
    seqid2group = {}
    for seqid, taxonomy in seqid2taxonomy.items():
        group = None
        for taxon in taxonomy[::-1]:
            if taxon in taxon2group:
                group = taxon2group[taxon]
                break
        if group == None:
            out_file.write('Taxon out of scope: {} {}\n'.format(seqid, ','.join(taxonomy)))
        seqid2group[seqid] = group
    # set base groups (intaxon vs outtaxon) # -1 - taxon out of the scope, 0 - non-query base group, 1 - query base group
    seqid2base_group = {}
    for seqid, taxonomy in seqid2taxonomy.items():
        if seqid2group[seqid] == None:
            seqid2base_group[seqid] = -1
        elif  len(taxonomy) > intaxon_level and taxonomy[intaxon_level] == args.intaxon:
            seqid2base_group[seqid] = 1
        else:
            seqid2base_group[seqid] = 0
    # read the gene tree
    gene_tree = Tree(args.tree_files[tree_index])
    # check the support values
    if args.support_cutoff == None:
        supports = []
        for node in gene_tree.traverse():
            supports.append(node.support)
            if max(supports) > 1.01:
                args.support_cutoff = 90
            else:
                args.support_cutoff = 0.9
    out_file.write('Support cutoff: {}\n'.format(args.support_cutoff))
    # root by the query
    try:
        gene_tree.set_outgroup(gene_tree&query_sequence_id)
    except:
        print('Error: Cannot find the query sequence id {}!'.format(query_sequence_id))
        raise
    # assign base group and taxonomy to leaves
    # also record outtaxons
    out_taxons = []
    for leaf in gene_tree.get_leaves():
        leaf.group = seqid2base_group[leaf.name]
        group_freq = 1
        if leaf.group == 0:
            out_taxons.append(seqid2group[leaf.name])
        elif leaf.group == -1:
            group_freq = 0
        leaf.taxonomy = [{leaf.group:group_freq}]
    out_taxons = set(out_taxons)
    # remove unsupported nodes
    removeUnsupported(gene_tree, node_support_cutoff=args.support_cutoff)
    # infer node groups
    setTaxonomy(gene_tree, tax_levels=1)
    # create pseudonodes solving multifurcations when appropriate
    createPseudonodes(gene_tree)
    out_file.write('Parsed gene tree: {}\n'.format(gene_tree.write()))
    setTaxonomy(gene_tree, tax_levels=1)
    # find intaxon orthologs
    seqid2orthology = {}
    in_taxons = []
    for leaf in gene_tree.get_leaves():
        if 1 in leaf.taxonomy[0]:
            min_in_freq = 1
            for ancestor in leaf.get_ancestors():
                in_freq = 0
                if 1 in ancestor.taxonomy[0]:
                    in_freq = ancestor.taxonomy[0][1]
                min_in_freq = min(min_in_freq, in_freq)
            if min_in_freq >= args.orthology_cutoff:
                in_taxons.append(seqid2group[leaf.name])
                out_file.write('Ortholog: {} {}\n'.format(leaf.name, min_in_freq))
            seqid2orthology[leaf.name] = min_in_freq
    in_taxons = set(in_taxons)
    out_file.write('In taxons: {}\n'.format(','.join(in_taxons)))
    out_file.write('Out taxons: {}\n'.format(','.join(out_taxons)))
    # now we have the number of in_taxons and out_taxons, so we can count the "directionality score"
    directionality_score = float(len(out_taxons))/(len(in_taxons)+len(out_taxons))
    out_file.write('Directionality score: {}\n'.format(directionality_score))
    out_file.write('# A low number of the Directionality score means that there is higher taxonomic diversity in the "in taxon" group and this could point to LGT in the other direction.')
    # get taxonomy tree of the intaxon
    intaxon_tree = None
    for node in taxonomy_tree.traverse():
        if node.name.lower() == args.intaxon:
            intaxon_tree = node
    if intaxon_tree == None:
        print('Error: Intaxons {} not found in the taxonomy_of_selected_groups!'.format(args.intaxon))
        raise
    # map presence of the orthologs on the taxonomy tree of the intaxon
    for leaf in intaxon_tree.get_leaves():
        leaf.presence = 0
        if leaf.name in in_taxons:
            leaf.presence = 1
    # Count number of losses in case the gene was present in the last common ancestor of the intaxon
    countSumPresence(intaxon_tree)
    findLosses(intaxon_tree)
    countSumLosses(intaxon_tree)
    # get all possible gains
    gains = selectGains(intaxon_tree)
    # find best combination of non-ancestral gains
    combo2score = scoreGainCombos(gains, intaxon_tree.presence_sum, intaxon_tree.losses_sum+1, gain_weight=args.gain_weight, loss_weight=args.loss_weight)
    if ('root',) in combo2score:
        del combo2score[('root',)]
    if (args.intaxon,) in combo2score:
        del combo2score[(args.intaxon,)]
    # 
    ancestral_losses = []
    for node in intaxon_tree.traverse():
        if node.loss:
            #out_file.write('Loss in case of ancestral presence: {}\n'.format(node.name))
            ancestral_losses.append( node.name )
    out_file.write('Losses in case of ancestral presence: {}\n'.format(','.join(ancestral_losses)))
    best_combo = None
    new_losses = []
    if len(combo2score) > 0:
        best_combo = max( combo2score.keys(), key=lambda x:(-combo2score[x], len(x)) )
        out_file.write('Gains in case of later acquisition: {}\n'.format(','.join(best_combo)))
        for gain in intaxon_tree.traverse():
            if best_combo != None and gain.name in best_combo: ####
                for node in gain.traverse():
                    if node.loss:
                        new_losses.append( node.name )
        out_file.write('Losses in case of later acquisition: {}\n'.format(','.join(new_losses)))

    nonancestral_score = 0
    if best_combo != None:
        nonancestral_score = float(len(ancestral_losses)*args.loss_weight)/( len(ancestral_losses)*args.loss_weight+len(new_losses)*args.loss_weight+(len(best_combo)-1)*args.gain_weight )
    out_file.write('Nonancestral score: {}\n'.format(nonancestral_score))
    out_file.write('# A low number of the Nonancestral score means that the presence/absence pattern of the gene suggests acquisition in the common ancestor of the "in taxon" group.\n')
    # Evaluate the donor taxon.
    longest_taxonomy = 0
    for taxonomy in seqid2taxonomy.values():
        longest_taxonomy = max(longest_taxonomy, len(taxonomy))
    for leaf in gene_tree.get_leaves():
        current_taxon = None
        leaf.taxonomy = []
        if leaf.name in seqid2orthology and seqid2orthology[leaf.name] >= args.orthology_cutoff:
            for tax_level in range(longest_taxonomy):
                leaf.taxonomy.append({})
        else:
            for tax_level in range(longest_taxonomy):
                leaf.taxonomy.append({})
                if tax_level < len(seqid2taxonomy[leaf.name]):
                    current_taxon = tuple(seqid2taxonomy[leaf.name][:tax_level+1]) ##
                leaf.taxonomy[tax_level] = {current_taxon:1}
    tree = setTaxonomy(gene_tree, longest_taxonomy)
    taxonomy_wo_inorthologs = tree.taxonomy
    donor_taxon = ['']
    donor_taxon_freq = ''
    donor_taxon_taxlevel = ''
    for tax_level in range(len(taxonomy_wo_inorthologs)):
        if len(taxonomy_wo_inorthologs[tax_level]) > 0 and max(taxonomy_wo_inorthologs[tax_level].values()) >= args.source_cutoff:
            donor_taxon, donor_taxon_freq = max(taxonomy_wo_inorthologs[tax_level].items(), key=lambda x:x[1])
            donor_taxon_taxlevel = tax_level
        if tax_level <= 5: ##
            line_out = sorted(taxonomy_wo_inorthologs[tax_level].items(), key=lambda x:(-x[1], x[0]))
            line_out = list(map(lambda x:x[0][-1]+':'+str(x[1]), line_out))
            line_out = ','.join(line_out)
            out_file.write('Sources at {} tax level: {}\n'.format(tax_level, line_out))
    out_file.write('Donor taxon: {}\n'.format(','.join(donor_taxon)))
    out_file.write('Donor taxon tax level: {}\n'.format(donor_taxon_taxlevel))
    out_file.write('Donor taxon fraction: {}\n'.format(donor_taxon_freq))
    out_file.write('Finished.\n')

if __name__=="__main__":
    import argparse
    args = parseArguments()
    import sys
    import re
    import traceback
    try:
        from ete3 import Tree
    except ImportError:
        from ete2 import Tree
    if args.taxdump_dir != None:
        from Taxonomy import Taxonomy
        taxonomy_parser = Taxonomy(args.taxdump_dir)
    # Read a system of taxonomic groups and their relationships. How we define these groups can influence the results a lot.
    # The table gives us the phylogenetic structure of those groups. The rightmost taxons, that I later on call "groups", should match an NCBI taxonomy group or our supplied taxonomy mapping.
    # Sometimes an arrow indicates merger of certain groups. E.g. I wanted to merge apicomplexa+chrompodellids into colpodellida, which however isn't annotated in the NCBI taxonomy.
    groups, group2taxonomy, taxon2group, taxonomy_tree = readTaxonomy(taxonomy_table=args.taxonomy_of_selected_groups, skip_viruses=args.mask_viruses)
    # multiple inputs
    for tree_index, tree_file in enumerate(args.tree_files):
        #find_lgts(args, tree_index)
        try:
            find_lgts(args, tree_index)
        except Exception as e:
            print(traceback.format_exc())
            print('Error running LGT analysis for: {}. Skipping...'.format(tree_file))
