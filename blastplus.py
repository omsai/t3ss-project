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
        # Add query coverage column.
        self.df['query coverage'] = (self.df['s. end'] - self.df['s. start']) \
                                        / self.df['query length']

    def top_hits(self,
                 querylen_frac=0.5,
                 top_identity_cutoff=0.5,
                 top_bitscore_frac=0.85):
        '''Return best hits of each query of homologous proteins.

        Parameters
	----------
        querylen_frac : float (0 to 1)
	    Filter hits greater than ((end - start) / querylen).
        top_identity_cutoff : float (0 to 1)
	    Filter out groups with top hit bitscores below this
	    fraction.
        top_bitscore_frac : float (0 to 1)
	    Filter query group members by this fraction of their top
	    hit bitscore.

        Returns
        -------
        DataFrame

        '''
        # Filter all hits by query coverage.
        mask_coverage = self.df['query coverage'] >= querylen_frac

        # Filter group members by top hit bitscore.
        top_hit_bitscores = self.df.groupby('query id')[['query id', 'bit score']].head(1)
        bitscore_min = top_hit_bitscores.set_index('query id') * top_bitscore_frac
        df2 = self.df.set_index(['query id', 'subject id'])
        # See http://stackoverflow.com/a/12463255 and
        # http://stackoverflow.com/a/27865433
        bitscore_min = bitscore_min.reindex(df2.index, level=0)
        mask_bitscore = self.df['bit score'] >= bitscore_min['bit score']

        # Filter group members by top hit identity percent.
        df_identity = self.df.groupby('query id').filter(
            lambda x: x['% identity'].head(1) > top_identity_cutoff * 100,
            dropna=False)
        mask_identity = self.df.isin(df_identity)['% identity']

        # Apply all filters
        mask = mask_coverage & mask_bitscore & mask_identity
        top_idx = mask[mask].index.values
        top_idx.sort()
        return self.df.iloc[top_idx]

    def uniq_hits(self):
        df_max = self.top_hits()
        arr_subj = df_max['subject id'].unique()
        return pd.Series(arr_subj)


if __name__ == '__main__':
    bp = BlastProcessor(BLASTP_OUTPUT_FILE)
