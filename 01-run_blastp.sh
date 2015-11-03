#!/bin/bash
# Create blastp output files using NIH BLAST+ `blastp` command-line
# program.

BLAST_DIR="/Volumes/Macintosh HD 2/thiberio/blast_dbs/"
BLAST_DB="nr"
BLAST_INPUT_FILE="/Users/omsai/data/blastp-input/all.fasta"
BLAST_OUTPUT_FILE="/Users/omsai/data/blastp-output/$(basename ${BLAST_INPUT_FILE}).blastp"
BLAST_NEGATIVE_GILIST="/Users/omsai/data/blastp-input/sequence.gi"

# Use a subshell (enclosing the commands in parantheses) so that
# failure of any of these commands does not change the working
# directory.
(
    cd "${BLAST_DIR}"
    pwd
    blastp -num_threads 12 -db ${BLAST_DB} -query ${BLAST_INPUT_FILE} -negative_gilist ${BLAST_NEGATIVE_GILIST} -out ${BLAST_OUTPUT_FILE} -outfmt 7\ qseqid\ sseqid\ qlen\ pident\ length\ mismatch\ gapopen\ qstart\ qend\ sstart\ send\ evalue\ bitscore
    cd - > /dev/null		# Return to original directory.
)
