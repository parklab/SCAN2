#!/bin/bash
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH -A park_contrib
#SBATCH -p park
#SBATCH --mem-per-cpu=12G

if [ $# -ne 8 ]; then
    echo "usage: $0 mmq60_vcf mmq1_vcf sc_dir sc_sample bulk_sample output_dir fdr nhsnps"
    echo "nhsnps - the number of hSNPs to sample to get CIGAR indel and clipping distributions"
    exit 1
fi

mmq60_vcf=$1
mmq1_vcf=$2
scdir=$3
sc=$4
bulk=$5
outdir=$6
fdr=$7
nhsnps=$8

mkdir -p $outdir

tmpfile=`mktemp -p $outdir`

echo "Step 0: converting VCFs to tables"
mmq60="$outdir/$(basename $mmq60_vcf .vcf).tab"
java -jar $GATK_PATH/gatk.jar -R $GATK_PATH/human_g1k_v37_decoy.fasta \
    -T SelectVariants \
    -V $mmq60_vcf \
    -selectType SNP -restrictAllelesTo BIALLELIC \
    -env -trimAlternates \
    -select 'vc.getGenotype("'$bulk'").isCalled()' \
    -o $tmpfile
totab.sh $tmpfile $mmq60

mmq1="$outdir/$(basename $mmq1_vcf .vcf).tab"
java -jar $GATK_PATH/gatk.jar -R $GATK_PATH/human_g1k_v37_decoy.fasta \
    -T SelectVariants \
    -V $mmq1_vcf \
    -selectType SNP -restrictAllelesTo BIALLELIC \
    -env -trimAlternates \
    -select 'vc.getGenotype("'$bulk'").isCalled()' \
    -o $tmpfile
totab.sh $tmpfile $mmq1
rm $tmpfile

echo "Step 1: finding somatic sites"
get_somatic_positions.sh $mmq60 $mmq1 $sc $bulk $outdir

# creates outdir/somatic_positions.txt
posfile=$outdir/somatic_positions.txt
if [ ! -f $posfile ]; then
    echo "ERROR: expected positions file not created: $posfile"
    exit 1
fi



echo "Step 2: counting CIGAR ops at each somatic candidate"
cigartxt=$outdir/somatic_cigars.txt
cigarfile=$outdir/somatic_cigars.tab
get_cigars.sh $posfile $scdir/wg.bam $cigartxt
if [ ! -f $cigartxt ]; then
    echo "ERROR: expected raw CIGAR op file not created: $cigartxt"
    exit 1
fi
count_cigars.py $cigartxt > $cigarfile



echo "Step 3: sampling $nhsnps hSNPs for CIGAR op distribution comparison"
echo "This can take a while.."
get_hsnp_positions.sh $scdir $outdir $nhsnps
hsnp_posfile=$outdir/hsnp_positions.txt
if [ ! -f $hsnp_posfile ]; then
    echo "ERROR: expected positions file not created: $hsnp_posfile"
    exit 1
fi



echo "Step 4: counting CIGAR ops at each hSNP"
hsnp_cigartxt=$outdir/hsnp_cigars.txt
hsnp_cigarfile=$outdir/hsnp_cigars.tab
get_cigars.sh $hsnp_posfile $scdir/wg.bam $hsnp_cigartxt
if [ ! -f $hsnp_cigartxt ]; then
    echo "ERROR: expected raw CIGAR op file not created: $hsnp_cigartxt"
    exit 1
fi
count_cigars.py $hsnp_cigartxt > $hsnp_cigarfile



echo "Step 5: estimating AB at each somatic candidate"
# creates outdir/somatic_ab.rda
get_somatic_ab.sh $scdir $outdir
abfile=$outdir/somatic_ab.rda
if [ ! -f $abfile ]; then
    echo "ERROR: expected AB estimate file not created: $abfile"
    exit 1
fi



echo "Step 6: genotyping"
# creates outdir/somatic_gt.rda
get_somatic_gt.sh $mmq60 $mmq1 $sc $bulk $outdir $fdr
gtfile=$outdir/somatic_gt.rda
if [ ! -f $gtfile ]; then
    echo "ERROR: expected genotype file not created: $gtfile"
    exit 1
fi
