#!/bin/bash
# Create blastp output files using NIH BLAST+ `blastp` command-line
# program.

BLAST_INPUT_DIR="/Volumes/Macintosh HD 2/thiberio/t3ss_fasta_files/"
BLAST_DB="../blast_dbs/nr"
BLAST_NEGATIVE_GILIST="~/data/blastp-input/sequence.gi.txt"
PARALLEL_LOG_DIR="~/code/log"

( cd "${BLAST_INPUT_DIR}"
  find . -iname '*.faa' |
      LC_ALL=C \
	    parallel \
	    --jobs 12 \
	    blastp \
	    -query "{}" \
	    -db "${BLAST_DB}" \
	    -outfmt 7 \
	    -negative_gilist "${BLAST_NEGATIVE_GILIST}" \
	    "> ~/data/blastp-output/{/}.blastp"
  cd - > /dev/null		# Return to original directory.
)
