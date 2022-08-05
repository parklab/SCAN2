#!/bin/bash

if [ $# -ne 2 ]; then
    echo "usage: $0 path_to_reference chromosome"
    echo "Do NOT include the 'chr' prefix in chromosome"
    echo "This script is intended for downloading hg38 phasing panel files ONLY."
    exit 1
fi

chr=$2
hg38_ref_path=$1

# chrX was phased by Eagle apparently and has a different file name
downloadpath=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased
if [ "x$chr" != "xX" ]; then
    downloadfile=${downloadpath}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz
else
    downloadfile=${downloadpath}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.eagle2-phased.v2.vcf.gz
fi


echo "Downloading for chromosome $chr.."
wget -O chr${chr}_GRCh38.genotypes.vcf.gz $downloadfile
wget -O chr${chr}_GRCh38.genotypes.vcf.gz.tbi ${downloadfile}.tbi

echo "Converting and normalizing variants for chromosome $chr.."
(bcftools view --no-version -h chr${chr}_GRCh38.genotypes.vcf.gz \
    | grep -v "^##contig=<ID=[GNh]" \
    | sed 's/^##contig=<ID=MT/##contig=<ID=chrM/;s/^##contig=<ID=\([0-9XY]\)/##contig=<ID=chr\1/' ; \
 bcftools view --no-version -H -c 2 chr${chr}_GRCh38.genotypes.vcf.gz ) \
| bcftools norm --no-version -Ou -m -any \
| bcftools norm --no-version -Ob -o chr${chr}_GRCh38.genotypes.bcf -d none -f $hg38_ref_path

echo "Indexing bcf.."
bcftools index -f chr${chr}_GRCh38.genotypes.bcf
