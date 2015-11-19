#!/usr/bin/env python
'''Test trees for monophily.'''
DEBUG = False
TREE_FORMAT='newick'

import Bio.Phylo as phy
from blastplus import BlastProcessor, BLASTP_OUTPUT_FILE
import os
TREE_DIR = '../data/fasttree-output/'

# Read all the files in as trees
filenames = os.listdir(TREE_DIR)
files = map(os.path.join, [TREE_DIR] * len(filenames), filenames)
trees = map(phy.read, files, [TREE_FORMAT] * len(files))

# Read Aeromonas names to be compared in tree terminals for monophily
bp = BlastProcessor(BLASTP_OUTPUT_FILE)
queries = bp.df['query id'].unique()

# Test for monophyly
digits = len(str(len(trees)))
is_monophylic = []
for tree in trees:
    all_terminals = tree.get_terminals()
    aeromonas_terminals = []
    for terminal in all_terminals:
        if terminal.name in queries:
            aeromonas_terminals += [terminal]
    if DEBUG:
        print("{0:>{digits}}/{1:>{digits}} terminals are Aeromonas".format(
            len(aeromonas_terminals), len(all_terminals), digits=digits))
    is_monophylic += [tree.is_monophyletic(aeromonas_terminals)]

num_monophylic = len(is_monophylic) - is_monophylic.count(False)
print('{0} of {1} trees are monophylic'.format(num_monophylic, len(is_monophylic)))
