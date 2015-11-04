'''Read and process blastp output file.'''

import pandas as pd

# Disable table wrapping
pd.set_option('expand_frame_repr', False)

BLASTP_OUTPUT_FILE = '../data/blastp-output/all.fasta.blastp'

class BlastProcessor(object):
    '''Pandas DataFrame and processing of of NIH BLAST+ `blastp`
command-line program output.

    Requires `blastp` to be run with the `-outfmt 7` option.

    '''
    
    def __init__(self, blastp_file):
        # The column header is commented out, therefore
        # pandas.read_csv() will ignore it.  As a workaround, read in
        # the that commented header line.
        names = None
        with open(blastp_file) as f:
            for i in range(10):
                line = f.readline()
                if line[:9] == '# Fields:':
                    names = line[10:].split(', ')
                    names[-1] = names[-1].strip('\r\n')
                    break
        # Read in the file_name.
        self.df = pd.read_csv(blastp_file, sep='\t', comment='#', names=names)

    def top_hits(self):
        # Get the best hit of each query
        q_grouped = self.df.groupby('query id')
        max_idx = q_grouped['bit score'].idxmax().values
        max_idx.sort()
        return self.df.iloc[max_idx]

    def uniq_hits(self):
        df_max = self.top_hits()
        arr_subj = df_max['subject id'].unique()
        return pd.Series(arr_subj)
