#!/usr/bin/env python
'''Read and process blastp output file.'''

from blastplus import BlastProcessor, BLASTP_OUTPUT_FILE

bp = BlastProcessor(BLASTP_OUTPUT_FILE)

# Write the subject ids to disk for sequence retrieval
bp.uniq_hits().to_csv('../data/analyze-output/subject_ids.ref',
                      header=False, index=False)
