# scan-snv
Single cell somatic genotyper

# Running the demo

## STEP 0. Download demo BAM files

    BAM and index files for the kindred system single cells and matched bulk
    can be found at http://compbio.med.harvard.edu/scan-snv/. To make the demo
    run in reasonable time, only reads from chromosome 22 are provided.

    At least one single cell and the unamplified bulk must be used. We recommend
    downloading hunamp and il-12. The two other kindred system samples are
    provided.

## STEP 1. Compile `laplace_cpu`

**Dependencies**: LAPACKE (v3.6.1), OpenBLAS (v0.2.19)
**Optional dependencies**: Intel C compiler

    1. Edit mkl-gridfit-gauss/Makefile (or Makefile.gcc if using gcc) to point to
       install paths for LAPACKE and OpenBLAS:

```
OPENBLAS=/n/app/openblas/0.2.19   # Set these 2 lines appropriately
LAPACKE=/n/app/lapacke/3.6.1
```
    2. Compile
```
make -f Makefile.gcc    # or Makefile.icc if using Intel C compiler
```

    3. **IMPORTANT** Add OpenBLAS to the linker path.

```
# Should be the same as $OPENBLAS with /lib appended
export LD_LIBRARY_PATH=/n/app/openblas/0.2.19/lib  
```


## STEP 2: Run GATK HaplotypeCaller on single cell and matched bulk data

**Dependencies**: Java (v1.8), GATK (v3.8-0-ge9d806836)
**Data dependencies**: human reference genome (GRCh37 with decoy), dbSNP (v147, b37)

    1. Edit `gatk/run_gatk.sh`. You will need to supply paths to the dependencies
       described above. This can be accomplished by editing the following lines:
```
RESOURCES=~/balance/resources
GATK=$RESOURCES/GATK3.8.jar
HG19=$RESOURCES/human_g1k_v37_decoy.fasta
DBSNP=$RESOURCES/dbsnp_147_b37_common_all_20160601.vcf
```
    2. Run the script two times, once using MMQ=60 and once using MMQ=1. The script
       will automatically create the directory `output_dir` to hold the generated
       VCF files. **Approximate run time**: 30 minutes each (~1 hour total) using
       a single core.
```
run_gatk.sh 60 output_dir hunamp.chr22.bam il-12.chr22.bam
run_gatk.sh 1 output_dir hunamp.chr22.bam il-12.chr22.bam
```
