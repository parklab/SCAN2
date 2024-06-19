#!/bin/bash

if [ $# -ne 4 ]; then
    echo "usage: $0 joint_depth_matrix.tab.gz reference_genome.fasta.fai column_number_for_sample output.bw"
    exit 1
fi

# Fledgling support for the snakemake 'script:' directive
#joint_depth_matrix="${snakemake_input[joint_matrix]}"
#fasta_index="${snakemake_input[fai]}"
#sample_name="${snakemake_wildcards[sample]}"
#output="${snakemake_output[0]}"

joint_depth_matrix=$1
fasta_index=$2
sample_column=$3
output=$4
echo $joint_depth_matrix
echo $fasta_index
echo $sample_column
echo $output

if [ -f $output ]; then
    echo "output file $output already exists, please delete it first"
    exit 1
fi

# For supporting per-chrom bigwigs
#tabix $joint_depth_matrix $chrom \

gunzip -c $joint_depth_matrix \
    | tail -n +2 \
    | cut -f 1,2,$sample_column \
    | awk -v OFS='\t' '{ print $1, $2-1, $2, $3; }' \
    | bigtools bedgraphtobigwig --sorted start --nzooms 0 - $fasta_index $output
