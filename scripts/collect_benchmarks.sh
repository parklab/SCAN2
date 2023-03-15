#!/bin/bash

if [ $# -ne 2 ]; then
    echo "usage: $0 scan2_dir output.txt"
    exit 1
fi

sd=$1
out=$2

if [ -f $out ]; then
    echo "output file $out already exists, please delete it first"
    echo "NOTE: benchmark collecting is only supported for SCAN2 --analysis=call_mutations and --analysis=makepanel"
    exit 1
fi

# Which analysis was run?
analysis=$(cat $sd/scan.yaml | shyaml get-value analysis)

# Bulk sample ID
bulk=$(cat $sd/scan.yaml | shyaml get-value bulk_sample)

# Single cells analyzed
scs=$(cat $sd/scan.yaml | shyaml keys sc_bams)

# Chromosomes to analyze
chroms=$(cat $sd/scan.yaml | shyaml get-values chrs)

# Number of windows the genome is divided into for GATK HaplotypeCaller and DepthOfCoverage
gatk_chunks=$(cat $sd/scan.yaml | shyaml get-value gatk_chunks)

# all benchmark files have identical headers. However, we need to choose a
# benchmark file shared across all analysis types (or do some special handling).
echo "job_category job subjob $(head -1 $sd/gatk/gather_benchmark.mmq60.tsv | tr '\t' ' ')" >> $out




if [ "x$analysis" == "xmakepanel" ]; then
    echo "GATK (MMQ >= 60)"
    echo "panel gather_mmq60 1 $(tail -n +2 $sd/gatk/gather_benchmark.mmq60.tsv | tr '\t' ' ')" >> $out
    for i in `seq 1 $gatk_chunks`; do
        echo "panel haplotypecaller_mmq60 $i $(tail -n +2 $sd/gatk/scatter_benchmark.mmq60_chunk${i}.tsv | tr '\t' ' ')" >> $out
    done

    echo "Make panel"
    echo "panel makepanel 1 $(tail -n +2 $sd/panel/benchmark_make_panel.txt | tr '\t' ' ')" >> $out
elif [ "x$analysis" == "xcall_mutations" ]; then
    echo "GATK (MMQ >= 60)"
    echo "gatk gather_mmq60 1 $(tail -n +2 $sd/gatk/gather_benchmark.mmq60.tsv | tr '\t' ' ')" >> $out
    for i in `seq 1 $gatk_chunks`; do
        echo "gatk haplotypecaller_mmq60 $i $(tail -n +2 $sd/gatk/scatter_benchmark.mmq60_chunk${i}.tsv | tr '\t' ' ')" >> $out
    done

    echo "GATK (MMQ >= 1)"
    echo "gatk gather_mmq1 1 $(tail -n +2 $sd/gatk/gather_benchmark.mmq1.tsv | tr '\t' ' ')" >> $out
    for i in `seq 1 $gatk_chunks`; do
        echo "gatk haplotypecaller_mmq1 $i $(tail -n +2 $sd/gatk/scatter_benchmark.mmq1_chunk${i}.tsv | tr '\t' ' ')" >> $out
    done

    echo "Depth profiling"
    echo "depth_profile joint_depth_matrix 1 $(tail -n +2 $sd/depth_profile/benchmark_joint_depth_matrix.txt | tr '\t' ' ')" >> $out
    for sc in $scs; do
        echo "depth_profile summarize $sc $(tail -n +2 $sd/depth_profile/benchmark_${sc}_region_summarize.txt | tr '\t' ' ')" >> $out
    done
    for i in `seq 1 $gatk_chunks`; do
        echo "depth_profile gatk_depthofcoverage $i $(tail -n +2 $sd/depth_profile/gatk_depthofcoverage/benchmark_depthofcoverage_chunk${i}.txt | tr '\t' ' ')" >> $out
    done


    echo "SHAPEIT"
    for chr in $chroms; do
        echo "phasing shapeit $chr $(tail -n +2 $sd/shapeit/$chr/benchmark.tsv | tr '\t' ' ')" >> $out
    done


    echo "AB model fitting"
    for sc in $scs; do
        echo "ab_model gather ${sc} $(tail -n +2 $sd/ab_model/$sc/benchmark_abmodel_gather.txt | tr '\t' ' ')" >> $out
        for chr in $chroms; do
            echo "ab_model scatter ${sc}_${chr} $(tail -n +2 $sd/ab_model/$sc/benchmark_abmodel_scatter_${chr}.tsv | tr '\t' ' ')" >> $out
        done
    done


    echo "Summarizing GATK tables"
    echo "call_mutations vcf_to_table mmq1 $(tail -n +2 $sd/call_mutations/benchmark_gatkvcf_to_tab_mmq1.txt | tr '\t' ' ')" >> $out
    echo "call_mutations vcf_to_table mmq60 $(tail -n +2 $sd/call_mutations/benchmark_gatkvcf_to_tab_mmq60.txt | tr '\t' ' ')" >> $out
    echo "call_mutations integrate_tables 1 $(tail -n +2 $sd/call_mutations/benchmark_integrate_tables.txt | tr '\t' ' ')" >> $out


    echo "CIGAR analysis"
    for id in $bulk $scs; do
        echo "call_mutations gather_cigars ${sc} $(tail -n +2 $sd/call_mutations/$id/benchmark_cigar_gather.txt | tr '\t' ' ')" >> $out
        for chr in $chroms; do
            echo "call_mutations scatter_cigars ${sc}_${chr} $(tail -n +2 $sd/call_mutations/$id/benchmark_cigars.${chr}.txt | tr '\t' ' ')" >> $out
        done
    done


    echo "Call mutations"
    for sc in $scs; do
        echo "call_mutations pregenotyping $sc $(tail -n +2 $sd/call_mutations/$sc/benchmark_pregenotyping.txt | tr '\t' ' ')" >> $out
        echo "call_mutations genotyping $sc $(tail -n +2 $sd/call_mutations/$sc/benchmark_genotype.txt | tr '\t' ' ')" >> $out
    done
else
    echo "ERROR: unsupported analysis type in $sd/scan.yaml ($analysis). Only analysis=call_mutations and analysis=makepanel are supported"
    exit 1
fi
