# scan-snv
Single cell somatic genotyper

# Running the demo
The provided demo shows how the pipeline is run in practice. However,
because parameter fitting (step 3) requires significant compute time,
the demo only analyzes reads from chromosome 22. Using a server with 8
cores, the demo took approximately **X HOURS** to run.
Because some properties of hSNPs and sSNVs (e.g., VAF distributions) are
measured from sites across all chromosome, the demo's output will not exactly
reconstitute the findings reported in the manuscript.

In each step, the dependency versions refer to the **tested** versions.
Other versions may work as well.

## STEP 0. Download demo BAM files
At least one single cell and the unamplified bulk must be used. We recommend
downloading hunamp and il-12. The two other kindred system samples are also
provided. Both the BAM and index file must be downloaded.

* **[REQUIRED]** Unamplified cell line bulk\
    BAM: http://compbio.med.harvard.edu/scan-snv/hunamp.chr22.bam \
    Index: http://compbio.med.harvard.edu/scan-snv/hunamp.chr22.bam.bai
* **[RECOMMENDED]** Kindred single cell IL-12\
    BAM: http://compbio.med.harvard.edu/scan-snv/il-12.chr22.bam \
    Index: http://compbio.med.harvard.edu/scan-snv/il-12.chr22.bam.bai
* **[OPTIONAL]** Kindred single cell IL-11\
    BAM: http://compbio.med.harvard.edu/scan-snv/il-11.chr22.bam \
    Index: http://compbio.med.harvard.edu/scan-snv/il-11.chr22.bam.bai
* **[OPTIONAL]** Kindred single cell-derived clone IL-1c\
    BAM: http://compbio.med.harvard.edu/scan-snv/il-1c.chr22.bam \
    Index: http://compbio.med.harvard.edu/scan-snv/il-1c.chr22.bam.bai

## STEP 1. Compile `laplace_cpu`

**Dependencies**: LAPACKE (v3.6.1), OpenBLAS (v0.2.19)\
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

**Dependencies**: Java (v1.8), GATK (v3.8-0-ge9d806836)\
**Data dependencies**: human reference genome (GRCh37 with decoy; e.g. `ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.fasta.gz`), dbSNP (v147, b37: e.g. `ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/common_all_20170710.vcf.gz`)

1. Edit `gatk/run_gatk_demo.sh`. You will need to supply paths to the dependencies
   described above. This can be accomplished by editing the following lines:
```
RESOURCES=~/balance/resources
GATK=$RESOURCES/GATK3.8.jar
HG19=$RESOURCES/human_g1k_v37_decoy.fasta
DBSNP=$RESOURCES/dbsnp_147_b37_common_all_20160601.vcf
```
2. Run the script two times, once using MMQ=60 and once using MMQ=1. The script
   will automatically create the directory `output_dir` to hold the generated
   VCF files. **Approximate run time**: 1 hour each using 4 cores.
```
run_gatk_demo.sh 60 output_dir hunamp.chr22.bam il-12.chr22.bam
run_gatk_demo.sh 1 output_dir hunamp.chr22.bam il-12.chr22.bam
```


## STEP 3: Define and phase hSNPs

**Dependencies**: SHAPEIT2 (v2, r837)
**Data dependencies**: SHAPEIT2 formatted 1000 genomes reference haplotype panel

1. Install and configure SHAPEIT2. 
    * Download and unzip SHAPEIT (e.g., https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r904.glibcv2.12.linux.tar.gz).
      The unzipped path is `SHAPEIT_ROOT`.
    * Download and unzip the 1000 genomes haplotype panel. (e.g., https://mathgen.stats.ox.ac.uk/impute/data_download_1000G_phase1_integrated_SHAPEIT2_16-06-14.html)
      The unzipped path is `REFPANEL_ROOT`.
    * Edit `shapeit/phase_chr.sh` and `shapeit/extract_chr.sh` and set
      `SHAPEIT_ROOT` and `REFPANEL_ROOT`:
```
SHAPEIT_ROOT=/path/to/shapeit
REFPANEL_ROOT=/path/to/1000g_reference_panel
```
2. Run SHAPEIT2 on heterozygous SNVs found in the bulk sample. Note that hSNPs are
   taken only from the MMQ=60 GATK output.
```
# Select only biallelic SNVs for phasing
java -jar /path/to/GATK -T SelectVariants -R /path/to/human_g1k_v37_decoy.fasta \
    -T SelectVariants -V output_dir/hc_raw.mmq60.vcf \
    -selectType SNP -restrictAllelesTo BIALLELIC -env -trimAlternates \
    -select 'vc.getGenotype("hunamp").isCalled()' \
    -o output_dir/hc_raw.mmq60.snp.vcf

# Only phase sites reported in the bulk
echo "hunamp" > samples_to_phase.txt

shapeit/phase_chr.sh output_dir/hc_raw.mm60.snp.vcf output_dir 22
shapeit/extract_chr.sh output_dir 22
```
   If successful, a file named `phased_hsnps.chr22.vcf` should exist in
   `output_dir`.


## STEP 4: Fit covariance function parameters via grid search
**Dependencies**: Python (v2.7), R (v3.3.3)

1. Create training data files.
```
# Creates a VCF with only 
# NOTE: the sample name "h25" corresponds to il-12
shapeit/get_hsnps_singlecell.sh h25 output_dir/hc_raw.mmq60.snp.vcf output_dir/phased_hsnps.chr22.vcf output_dir

awk -f gridfit_slurm/totab.awk output_dir/h25_hsnps.vcf > output_dir/h25_hsnps.tab

Rscript gridfit_slurm/torda.R output_dir/h25_hsnps.tab output_dir/training_chr%d
```

2. Run the grid searching algorithm to find parameter MLEs. Using the parameters
   below (ngrids=8, points-per-grid=400), the run took ~20 minutes. The parameters
   used in the manuscript were ngrids=20, points-per-grid=1000.
```
mkdir -p output_dir/gridfit/chr22

./gridfit_slurm/runchr.py \
    --bindata=output_dir/training_chr22.bin \
    --local \
    --ngrids 8 \          # Set this to the number of cores available
    --points-per-grid 400 \
    --resume \
    --combine=../gridfit_slurm/combine.R \
    --mkl-laplace=../mkl-gridfit-gauss/laplace_cpu_gcc \
    --outprefix=output_dir/gridfit/chr22
```
    If the run is successful, the file `output_dir/gridfit/chr22/fit.rda` will be created.

3. Make the final fit file.
```
Rscript ./gridfit_slurm/make_fits.R
```


## STEP 5: Run SCAN-SNV
