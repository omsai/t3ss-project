'''Analyze blastp output files.'''

import pandas as pd
import os

BLASTP_EXT = '.faa.blastp'

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
    df = pd.DataFrame()
    dirname = '../data/blastp-output/'
    for blastp_file in os.listdir(dirname):
        df_homologs = blastp_output(os.path.join(dirname, blastp_file))
        # Strip file extension
        homologs = blastp_file[:-len(BLASTP_EXT)]
        df_homologs['query homologs'] = homologs
        # Append to single DataFrame
        df = pd.concat([df, df_homologs])
    # Re-enumerate the index
    df = df.reset_index()
    # Filter by evalue
    df_filtered = df[df['evalue'] < 1e-5]
    # Filter by max bit score in each group
    grouped = df_filtered.groupby('query id')
    max_idx = grouped['bit score'].idxmax().values
    df_max = df.iloc[max_idx]
