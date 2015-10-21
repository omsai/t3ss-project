'''Analyze blastp output files.'''

import pandas as pd

def blastp_output(file_name):
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
    # Read all files
    blastp_file = '../data/analyze-input/all.faa.blastp'
    df = blastp_output(blastp_file)
    # Filter by evalue
    df_filtered = df[df['evalue'] < 1e-5]
    # Filter by max bit score in each group
    grouped = df_filtered.groupby('query id')
    max_idx = grouped['bit score'].idxmax().values
    df_max = df.iloc[max_idx]
