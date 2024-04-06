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

# Was SHAPEIT or Eagle used to phase?
phaser=$(cat $sd/scan.yaml | shyaml get-value phaser)

# Which method was used to fit the AB model? grid_search or gradient_descent?
abmodel_method=$(cat $sd/scan.yaml | shyaml get-value abmodel_method)

# The use_fit option allows users to skip the expensive parameter fit, instead
# using previously calculated values. There are no benchmark files in that case.
# This is done on a per-single cell basis, so there can be a mixture of pre-
# computed AB fits and AB fits using abmodel_method above.
abmodel_use_fit_samples=$(cat $sd/scan.yaml | shyaml keys abmodel_use_fit)

# Was basepair resolution depth measured by GATK DepthOfCoverage or samtools depth?
depth_method=$(cat $sd/scan.yaml | shyaml get-value depth_method)

# Which GATK version (or Sentieon) and strategy (joint vs. GVCF) was used?
gatk=$(cat $sd/scan.yaml | shyaml get-value gatk)

# Bulk sample ID
bulk=$(cat $sd/scan.yaml | shyaml get-value bulk_sample)

# All single cells analyzed
scs=$(cat $sd/scan.yaml | shyaml keys sc_bams)

# Chromosomes to analyze
chroms=$(cat $sd/scan.yaml | shyaml get-values chrs)

# Some tools don't run on chrY. Get the correct name of the Y chrom from the config
chrY=$(cat $sd/scan.yaml | shyaml get-value chrY)

# Number of windows the genome is divided into for GATK HaplotypeCaller and DepthOfCoverage
analysis_regions=$(cat $sd/scan.yaml | shyaml get-values analysis_regions)

# Was spatial sensitivity calculated?
spatial_sens=$(cat $sd/scan.yaml | shyaml get-value compute_sensitivity)

# all benchmark files have identical headers. However, we need to choose a
# benchmark file shared across all analysis types (or do some special handling).
echo -e "category\tjob\tsubjob\tfile\t$(head -1 $sd/gatk/gatk_gather.mmq60.benchmark.tsv)" >> $out

# writes to the file defined by the global "out" variable
writeline () {
    category=$1
    job=$2
    subjob=$3
    benchfile=$4

    # -e: interpret escape sequences like '\t'
    echo -e "$category\t$job\t$subjob\t$benchfile\t$(tail -n +2 $benchfile)" >> $out
}


if [ "x$analysis" == "xmakepanel" ]; then
    echo "GATK (MMQ >= 60)"
    writeline "panel" "gather_mmq60" "1" "$sd/gatk/concat_benchmark.mmq60.tsv"
    for region in $analysis_regions; do
        writeline "panel" "haplotypecaller_mmq60" "$region" "$sd/gatk/scatter_benchmark.mmq60.region_${region}.tsv"
    done

    echo "Make panel"
    writeline "panel" "bcftotab" "1" "$sd/panel/benchmark_bcf_to_tab.txt"
    writeline "panel" "makepanel" "1" "$sd/panel/benchmark_make_panel.txt"

elif [ "x$analysis" == "xcall_mutations" ]; then
    for mmq in 1 60; do
        echo "GATK (MMQ >= ${mmq})"
        for region in $analysis_regions; do
            if [ "x$gatk" == "xgatk3_joint" ] || [ "x$gatk" == "xgatk4_joint" ] || [ "x$gatk" == "xsentieon_joint" ]; then
                writeline "gatk" "haplotypecaller_mmq${mmq}" "$region" "$sd/gatk/scatter_benchmark.mmq${mmq}.region_${region}.tsv"
            fi
            if [ "x$gatk" == "xgatk4_gvcf" ] || [ "x$gatk" == "xsentieon_gvcf" ]; then
                writeline "gatk" "genotypegvcfs_mmq${mmq}" "$region" "$sd/gatk/genotypegvcfs_benchmark.mmq${mmq}.region_${region}.tsv"
                for id in $bulk $scs; do
                    writeline "gatk" "haplotypecaller_mmq${mmq}" "${id}_${region}" "$sd/gatk/$id/scatter_benchmark.mmq${mmq}.region_${region}.tsv"
                done
            fi
            if [ "x$gatk" == "xgatk4_gvcf" ]; then
                writeline "gatk" "combinegvcfs_mmq${mmq}" "$region" "$sd/gatk/combinegvcfs_benchmark.mmq${mmq}.region_${region}.tsv"
            fi
        done

        writeline "gatk" "gather_mmq${mmq}" "1" "$sd/gatk/gatk_gather.mmq${mmq}.benchmark.tsv"
    done

    echo "Depth profiling (depth_method=$depth_method)"
    writeline "depth_profile" "joint_depth_matrix" "1" "$sd/depth_profile/benchmark_joint_depth_matrix.txt"
    for sc in $scs; do
        writeline "depth_profile" "summarize" "$sc" "$sd/depth_profile/benchmark_${sc}_region_summarize.txt"
    done

    # GATK uses analysis regions because it is very slow, samtools depth uses chromosomes
    if [ "x$depth_method" == "xgatkdocov" ]; then
        for region in $analysis_regions; do
            writeline "depth_profile" "gatk_depthofcoverage" "$region" "$sd/depth_profile/gatk_depthofcoverage/benchmark_depthofcoverage.region_${region}.txt"
        done
    elif [ "x$depth_method" == "xsamtoolsdepth" ]; then
        for chr in $chroms; do
            writeline "depth_profile" "samtoolsdepth" "$chr" "$sd/depth_profile/samtoolsdepth/benchmark_samtoolsdepth.region_${chr}.txt"
        done
    else
        echo "Unrecognized depth_method"
    fi

    echo 'Binned counts'
    for id in $bulk $scs; do
        writeline "depth_profile" "binned_counts" "$id" "$sd/depth_profile/binned_counts/${id}.benchmark.txt"
    done


    # I think this will all be applicable to eagle as well, but not sure yet.
    echo "Phasing (phaser=$phaser)"
    if [ "x$phaser" == "xshapeit" ]; then
        writeline "phasing" "phasing_prepare_bulk" "snv" "$sd/$phaser/analyzable_sites.mmq60.filtered_bulk_snv_only.benchmark.txt"
        writeline "phasing" "phasing_prepare_bulk" "indel" "$sd/$phaser/analyzable_sites.mmq60.filtered_bulk_indel_only.benchmark.txt"

        for chr in $chroms; do
            writeline "phasing" "prepare" "$chr" "$sd/$phaser/$chr/final_for_phasing.analyzable_sites.mmq60.benchmark.txt"
            writeline "phasing" "qual_cutoff_snv" "$chr" "$sd/$phaser/$chr/analyzable_sites.mmq60.filtered_bulk_snv_only.qual_cutoff.benchmark.txt"
            writeline "phasing" "qual_cutoff_indel" "$chr" "$sd/$phaser/$chr/analyzable_sites.mmq60.filtered_bulk_indel_only.qual_cutoff.benchmark.txt"
            writeline "phasing" "$phaser" "$chr" "$sd/$phaser/$chr/benchmark_phaser.tsv"
            if [ "x$chr" != "x$chrY" ]; then
                writeline "phasing" "postprocess" "$chr" "$sd/$phaser/$chr/postprocess.benchmark.txt"
            fi
        done
        writeline "phasing" "gather_final" "1" "$sd/$phaser/phased_filtered.benchmark.txt"
    fi


    echo "AB model fitting (method=$abmodel_method)"
    ab_suffix=""
    if [ "x$abmodel_method" == "xgradient_descent" ]; then
        ab_suffix="_gradient"
    fi
    for sc in $scs; do
        if [ $(echo "$abmodel_use_fit_samples" | tr ' ' '\n' | grep -c "^$sc\$") == "1" ]; then
            echo "skipping AB model for sample $sc (use_fit)"
        else
            writeline "ab_model" "gather$ab_suffix" "${sc}" "$sd/ab_model$ab_suffix/$sc/benchmark_abmodel_gather.txt"
            for chr in $chroms; do
                writeline "ab_model" "scatter$ab_suffix" "${sc}_${chr}" "$sd/ab_model$ab_suffix/$sc/benchmark_abmodel_scatter_${chr}.tsv"
            done
        fi
    done


    echo "Call mutations"
    writeline "call_mutations" "phase_info" "1" "$sd/call_mutations/phase_info.benchmark.txt"
    writeline "call_mutations" "analyzable_sites" "mmq60" "$sd/call_mutations/analyzable_sites.mmq60.benchmark.txt"
    writeline "call_mutations" "analyzable_sites" "mmq1" "$sd/call_mutations/analyzable_sites.mmq1.benchmark.txt"
    writeline "call_mutations" "tablefy" "mmq1" "$sd/call_mutations/benchmark_tablefy_sites_mmq1.txt"
    writeline "call_mutations" "tablefy" "mmq60" "$sd/call_mutations/benchmark_tablefy_sites_mmq60.txt"
    writeline "call_mutations" "integrate_tables" "1" "$sd/call_mutations/benchmark_integrate_tables.txt"


    echo "CIGAR analysis (per-sample)"
    for id in $bulk $scs; do
        writeline "call_mutations" "gather_cigars" "${sc}" "$sd/call_mutations/$id/benchmark_cigar_gather.txt"
        for chr in $chroms; do
            writeline "call_mutations" "scatter_cigars" "${sc}_${chr}" "$sd/call_mutations/$id/benchmark_cigars.${chr}.txt"
        done
    done


    echo "Genotyping (per-sample)"
    for sc in $scs; do
        writeline "call_mutations_precompute" "precompute_ab_ests_and_models_gather" "${sc}" "$sd/call_mutations/$sc/ab_ests_and_models.benchmark.txt"
        writeline "call_mutations_precompute" "precompute_excess_cigar_scores_gather" "${sc}" "$sd/call_mutations/$sc/excess_cigar_scores.benchmark.txt"
        if [ "x$spatial_sens" == "xTrue" ]; then
            writeline "call_mutations_sensitivity" "sensitivity_abmodel_gather" "${sc}" "$sd/call_mutations/$sc/abmodel_covariates.benchmark.txt"
            writeline "call_mutations_sensitivity" "sensitivity_depth_gather" "${sc}" "$sd/call_mutations/$sc/depth_covariates.benchmark.txt"
        fi

        for chr in $chroms; do
            writeline "call_mutations_precompute" "precompute_ab_ests_and_models_scatter" "${sc}_${chr}" "$sd/call_mutations/$sc/ab_ests_and_models.${chr}.benchmark.txt"
            writeline "call_mutations_precompute" "precompute_excess_cigar_scores_scatter" "${sc}_${chr}" "$sd/call_mutations/$sc/excess_cigar_scores.${chr}.benchmark.txt"
            if [ "x$spatial_sens" == "xTrue" ]; then
                writeline "call_mutations_sensitivity" "sensitivity_abmodel_scatter" "${sc}_${chr}" "$sd/call_mutations/$sc/abmodel_covariates.${chr}.benchmark.txt"
                writeline "call_mutations_sensitivity" "sensitivity_depth_scatter" "${sc}_${chr}" "$sd/call_mutations/$sc/depth_covariates.${chr}.benchmark.txt"
            fi
        done
        writeline "call_mutations" "pregenotyping" "$sc" "$sd/call_mutations/$sc/benchmark_pregenotyping.txt"
        writeline "call_mutations" "genotyping" "$sc" "$sd/call_mutations/$sc/benchmark_call_mutations.txt"
    done

else
    echo "ERROR: unsupported analysis type in $sd/scan.yaml ($analysis). Only analysis=call_mutations and analysis=makepanel are supported"
    exit 1
fi
