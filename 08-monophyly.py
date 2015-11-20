#!/usr/bin/env python
'''Test trees for polyphily and reexport with taxa information.'''

from __future__ import print_function
import Bio.Phylo as phy
from blastplus import BlastProcessor, BLASTP_OUTPUT_FILE
import os
import re                       # For extracting taxa

TREE_DIR = '../data/fasttree-output/'
TREE_FORMAT = 'newick'
SEQ_DIR = '../data/mafft-output/'
SEQ_EXT = 'mafft'
OUTPUT_DIR = '../data/phylo-output/'

# Read all the files in as trees
filenames = os.listdir(TREE_DIR)
files = map(os.path.join, [TREE_DIR] * len(filenames), filenames)
trees = map(phy.read, files, [TREE_FORMAT] * len(files))

# Read Aeromonas names to be compared in tree terminals for monophily
bp = BlastProcessor(BLASTP_OUTPUT_FILE)
queries = bp.df['query id'].unique()

# Test for monophyly
digits = len(str(len(trees)))
is_monophyletic = []
for tree, filename in zip(trees, files):
    all_terminals = tree.get_terminals()

    # Ignore <= 3 taxa trees.
    if len(all_terminals) <= 3:
        is_monophyletic += ['Ignored due to <= 3 terminals']
        continue

    # Read taxa names of gi blast results
    tree_filename = os.path.basename(filename)
    seq_filename = os.path.splitext(tree_filename)[0] + '.' + SEQ_EXT
    seq_file = os.path.join(SEQ_DIR, seq_filename)
    with open(seq_file, 'rU') as fp:
        lines = fp.readlines()
        comments = [line[1:-1] for line in lines if line[0] == '>']
    gene_to_taxa = {}
    regex_gene = r'(gi\|[0-9]+\|[a-z]+\|[^\|]+\|)'
    regex_sep = r'[^\[]+\['
    regex_taxa = r'([^\]]+)'
    for comment in comments:
        if comment.count(' ') > 0:
            pairs = re.findall(regex_gene + regex_sep + regex_taxa, comment)
            genes = [] # DEBUG
            for pair in pairs:
                gene, taxa = pair
                gene_to_taxa[gene] = taxa
                genes += [gene] # DEBUG
            assert len(pairs) == len([gene_to_taxa[gene] for gene in genes])

    aeromonas_terminals = []
    for i in xrange(len(all_terminals)):
        terminal = all_terminals[i]
        if terminal.name in queries:
            aeromonas_terminals += [terminal]
        else:
            # Append taxa names to gi blast results
            gene = re.findall(regex_gene, all_terminals[i].name)[0]
            if gene_to_taxa.has_key(gene):
                all_terminals[i].name += gene_to_taxa[gene].replace(' ', '_')
            else:
                print('Warning: unable to find taxa name of gene {}'.format(gene))

    tree.root_at_midpoint()
    is_monophyletic += [tree.is_monophyletic(aeromonas_terminals)]

    # Print paraphyletic trees
    if is_monophyletic[-1] == False:
        num_aero = len(aeromonas_terminals)
        num_all = len(all_terminals)
        print("{diff:>{digits}} blast hits, {aero:>{digits}} Aeromonas, ".format(
            aero=num_aero, digits=digits, diff=num_all-num_aero), end='')
        print(tree_filename)

# Save all the paraphyletic trees
output_files = []
num_paraphyletic = 0
for i, tree in enumerate(trees):
    if is_monophyletic[i]:
        sub_dir = 'monophyletic'
    else:
        sub_dir = 'paraphyletic'
        num_paraphyletic += 1
    filename = filenames[i]
    output_filename = os.path.splitext(filename)[0] + '.' + TREE_FORMAT
    output_file = os.path.join(OUTPUT_DIR, sub_dir, filename)
    output_files += [output_file]

map(phy.write, trees, output_files,
    [TREE_FORMAT] * len(trees))
print('Finished writing {} trees'.format(len(output_files)))
print('{} of {} trees are paraphyletic'.format(num_paraphyletic, len(trees)))
