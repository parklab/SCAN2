#!/bin/bash
#SBATCH -c 1
#SBATCH -t 6:00:00
#SBATCH -A park_contrib
#SBATCH -p park
#SBATCH --mem-per-cpu=12G

if [ $# -ne 6 ]; then
    echo "usage: $0 mmq60_table mmq1_table sc_dir bulk_sample output_dir fdr"
    echo "note: sc_dir must match the sample name used in the gatk tables"
    exit 1
fi

mmq60=$1
mmq1=$2
scdir=$3
bulk=$4
outdir=$5
fdr=$6

mkdir -p $outdir

module load gcc
module load R/3.3.3

scansnv_dir=/home/ljl11/balance/scan-snv

echo "Step 1: finding somatic sites"
$scansnv_dir/get_somatic_positions.sh $mmq60 $mmq1 $scdir $bulk $outdir

# creates outdir/somatic_positions.txt
posfile=$outdir/somatic_positions.txt
if [ ! -f $posfile ]; then
    echo "ERROR: expected positions file not created: $posfile"
    exit 1
fi



echo "Step 2: counting CIGAR ops at each somatic candidate"
cigartxt=$outdir/somatic_cigars.txt
cigarfile=$outdir/somatic_cigars.tab
$scansnv_dir/get_cigars.sh $posfile $scdir/wg.bam $cigartxt
if [ ! -f $cigartxt ]; then
    echo "ERROR: expected raw CIGAR op file not created: $cigartxt"
    exit 1
fi
$scansnv_dir/count_cigars.py $cigartxt > $cigarfile



echo "Step 3: sampling hSNPs for CIGAR op distribution comparison"
echo "This can take a while.."
$scansnv_dir/get_hsnp_positions.sh $scdir $outdir
hsnp_posfile=$outdir/hsnp_positions.txt
if [ ! -f $hsnp_posfile ]; then
    echo "ERROR: expected positions file not created: $hsnp_posfile"
    exit 1
fi



echo "Step 4: counting CIGAR ops at each hSNP"
hsnp_cigartxt=$outdir/hsnp_cigars.txt
hsnp_cigarfile=$outdir/hsnp_cigars.tab
$scansnv_dir/get_cigars.sh $hsnp_posfile $scdir/wg.bam $hsnp_cigartxt
if [ ! -f $hsnp_cigartxt ]; then
    echo "ERROR: expected raw CIGAR op file not created: $hsnp_cigartxt"
    exit 1
fi
$scansnv_dir/count_cigars.py $hsnp_cigartxt > $hsnp_cigarfile



echo "Step 5: estimating AB at each somatic candidate"
# creates outdir/somatic_ab.rda
$scansnv_dir/get_somatic_ab.sh $scdir $outdir
abfile=$outdir/somatic_ab.rda
if [ ! -f $abfile ]; then
    echo "ERROR: expected AB estimate file not created: $abfile"
    exit 1
fi



echo "Step 6: genotyping"
# creates outdir/somatic_gt.rda
$scansnv_dir/get_somatic_gt.sh $mmq60 $mmq1 $scdir $bulk $outdir $fdr
gtfile=$outdir/somatic_gt.rda
if [ ! -f $gtfile ]; then
    echo "ERROR: expected genotype file not created: $gtfile"
    exit 1
fi
