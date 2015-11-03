#!/bin/bash
# Merge blastp output files from NIH BLAST+ `blastp` command-line
# program and create an .xz archive.

INPUT_DIR="../data/blastp-output"
OUTPUT_DIR="../data/analyze-input"
OUTPUT_FILE="all.faa.blastp"

find $INPUT_DIR -iname '*.faa.blastp' -exec cat {} \; > $OUTPUT_DIR/$OUTPUT_FILE
(
    cd $OUTPUT_DIR
    xz --keep $OUTPUT_FILE
    cd -
)
