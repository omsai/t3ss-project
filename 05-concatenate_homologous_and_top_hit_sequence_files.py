#!/usr/bin/env python
'''Create sequence group files for alignment.  Each will have the
original homologous group sequences, as well as all their matching
subject id sequences.

'''
from blastplus import BlastProcessor, BLASTP_OUTPUT_FILE
import glob
import os
import pandas as pd

QUERY_FILE_GLOB = '../data/blastp-input/*.faa'
QUERY_DIR = os.path.dirname(QUERY_FILE_GLOB)
SUBJ_DIR = '../data/analyze-output/'
ALIGN_DIR = '../data/mafft-input/'

# Create a DataFrame of the filenames with their homologous sequences
columns = {'filename', 'query id'}
df_queries = pd.DataFrame(columns=columns)
query_files = glob.glob(QUERY_FILE_GLOB)
COMMENT = '>'
for query_file in query_files:
    with open(query_file, 'rU') as f:
        query_ids = []
        for line in f.xreadlines():
            if line[0] == COMMENT:
                query_ids += [line[1:-1]] # Leave out newline character
    filenames = len(query_ids) * [os.path.basename(query_file)]
    df_query = pd.DataFrame({
        'filename': filenames,
        'query id': query_ids,
    })
    df_queries = df_queries.append(df_query)
df_queries = df_queries.reset_index()

# Create DataFrame of blast results
bp = BlastProcessor(BLASTP_OUTPUT_FILE)
df_blast_results = bp.top_hits()

# Merge the DataFrames to get the filename and blast hits
df = pd.merge(df_queries, df_blast_results)
df = df[['filename', 'subject id']]
mask_dups = df.duplicated()
df = df[~mask_dups]             # Get unique values

# Concatenate the query and subject sequences into the alignment files
for query_filename, group in df.groupby('filename'):
    query_file = os.path.join(QUERY_DIR, query_filename)
    align_filename = os.path.splitext(query_filename)[0] + '.cat.faa'
    align_file = os.path.join(ALIGN_DIR, align_filename)
    with open(query_file, 'rU') as f:
        query = f.readlines()   
    subject_filenames = group['subject id'].values
    hits = []
    for subject_filename in subject_filenames:
        subj_file = os.path.join(SUBJ_DIR, subject_filename) + '.faa'
        with open(subj_file, 'rU') as f:
            hit = f.readlines()
        hits += hit
    with open(align_file, 'w') as f:
        f.writelines(query + hits)
    print('Wrote {}'.format(align_file))
