'''samtools faidx [...] on our lab server crashes due to high memory
usage.  Instead of using samtools, we can directly read the nr text
file directly using `xreadlines()` to conserve memory.

This implementation assumes the sequence of interest may be anywhere
in the comment line instead of just at the beginning, so we use regex
matches and set intersection to test for matches.

'''
import os                       # To write files
import re

db_file='/Volumes/Macintosh HD 2/thiberio/blast_dbs/nr'
refs_file='../data/analyze-output/subject_ids.ref'
output_dir='../data/analyze-output'

# Read in the sequence IDs
with open(refs_file, 'rU') as f:
    sequence_ids = set(f.read().splitlines())

# Search for the sequence IDs
COMMENT = '>'
PATTERN = re.compile('gi\|[^ ]+')
active_sequence_id = None
sequence = []
with open(db_file, 'rU') as f:
    for line in f.xreadlines():
        if line[0] == COMMENT:
            line_sequence_ids = set(re.findall(PATTERN, line))
            # Write any previous sequence to disk.
            if active_sequence_id != None:
                output_file = os.path.join(
                    output_dir,
                    '{0}.faa'.format(active_sequence_id)
                )
                with open(output_file, 'w') as fo:
                    fo.writelines(sequence)
            active_sequence_id = None
            sequence = []
            matching_sequence_ids = sequence_ids & line_sequence_ids
            if len(matching_sequence_ids) > 0:
                active_sequence_list = list(matching_sequence_ids)
                active_sequence_list.sort()
                active_sequence_id = '_'.join(active_sequence_list)
        elif active_sequence_id != None:
            sequence += [line]
