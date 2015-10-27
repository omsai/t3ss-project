refs_file=/Users/omsai/data/analyze-output/subject_ids.ref
output_file=/Users/omsai/data/analyze-output/samtools-faidx.faa

ids_inline=$(cat ${refs_file} | tr '\n' ' ')
(
    cd '/Volumes/Macintosh HD 2/thiberio/blast_dbs'
    time samtools faidx nr ${ids_inline} | tee ${output_file}
    cd - > /dev/null
)
