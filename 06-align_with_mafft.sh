#!/usr/bin/env bash

# mafft does not perform well multi threaded, so we will use it
# serially and parallelize it using GNU Parallel.

in_dir='../data/mafft-input'
out_dir='../data/mafft-output'
log_dir='log'

ls ${in_dir}/*.faa | \
    parallel \
	--jobs 50% \
	--joblog ${log_dir}/parallel-mafft.log \
	--resume-failed \
	mafft --auto --reorder "{} > ${out_dir}/{/.}.mafft"
