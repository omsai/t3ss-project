'''Analyze blastp output files.'''

import pandas as pd
import os

# Disable table wrapping
pd.set_option('expand_frame_repr', False)

BLASTP_OUTPUT_FILE = '../data/blastp-output/all.fasta.blastp'

def blastp_output(blastp_file):
    '''Read output from NIH BLAST+ `blastp` command-line program run with
    the `-outfmt 7` option.

    '''
    # The column header is commented out, therefore pandas.read_csv()
    # will ignore it.  As a workaround, read in the that commented
    # header line.
    names = None
    with open(blastp_file) as f:
        for i in range(10):
            line = f.readline()
            if line[:9] == '# Fields:':
                names = line[10:].split(', ')
                names[-1] = names[-1].strip('\r\n')
                break
    # Read in the file_name.
    return pd.read_csv(blastp_file, sep='\t', comment='#', names=names)

if __name__ == '__main__':
    # Read each file, adding a homologous group column
    df = blastp_output(BLASTP_OUTPUT_FILE)
    # Get the best hit of each query
    q_grouped = df.groupby('query id')
    max_idx = q_grouped['bit score'].idxmax().values
    max_idx.sort()
    df_max = df.iloc[max_idx]

    # Write the subject ids to disk for samtools faidx
    arr_subj = df_max['subject id'].unique()
    s_subj = pd.Series(arr_subj)
    s_subj.to_csv('../data/analyze-output/subject_ids.ref',
                  header=False, index=False)
