#!/usr/bin/env bash

# FastTree is not compiled with openmp, so we will run it serially and
# parallelize it using GNU Parallel.

in_dir='../data/mafft-output'
out_dir='../data/fasttree-output'
log='log/parallel-fasttree.log'

rm -f ${log}

ls ${in_dir}/*.mafft | \
    parallel \
	--jobs 50% \
	--joblog ${log} \
	--resume-failed \
	../apps/fasttree/FastTree "{} > ${out_dir}/{/.}.fasttree"
