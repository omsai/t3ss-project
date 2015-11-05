#!/usr/bin/env bash

# mafft does not perform well multi threaded, so we will use it
# serially and parallelize it using GNU Parallel.

in_dir='../data/mafft-input'
out_dir='../data/mafft-output'
log='log/parallel-mafft.log'

rm -f ${log}

ls ${in_dir}/*.faa | \
    parallel \
	--jobs 50% \
	--joblog ${log} \
	--resume-failed \
	mafft --auto --reorder "{} > ${out_dir}/{/.}.mafft"
