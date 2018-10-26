# scan-snv
Single cell somatic genotyper


## Installation
Version numbers in parentheses denote the versions used in the manuscript. They
are not necessarily required to run.

**Dependencies**: Python (v2.7), R (v3.3.3), LAPACKE (v3.6.1), OpenBLAS (v0.2.19),
    Java (v1.8), GATK (v3.8-0-ge9d806836), SHAPEIT2 (v2-r837)

**Optional dependencies**: Intel C compiler (2016)

**Data dependencies**:
    * Human reference genome version GRCh37 with decoy
        e.g. ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.fasta.gz
    * dbSNP (v147, b37)
        e.g. ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/common_all_20170710.vcf.gz
    * SHAPEIT2 formatted 1000 genomes reference haplotype panel
        e.g., https://mathgen.stats.ox.ac.uk/impute/data_download_1000G_phase1_integrated_SHAPEIT2_16-06-14.html

> **IMPORTANT!** The environment variables listed below (`LD_LIBRARY_PATH`,
> `PATH`, `GATK_PATH`, `SHAPEIT_ROOT`, `REFPANEL_ROOT`) must always be set
> appropriately before running SCAN-SNV.

1. Build and install the SCAN-SNV R package
```
$ cd /root/of/git/repo
$ R CMD build rpkg
$ R CMD INSTALL scansnv-0.1.tar.gz
```
2. Install the Laplace approximator. First edit `gridfit-gauss/Makefile.gcc`
   if you are compiling with gcc (`gridfit-gauss/Makefile.icc if you are
   using Intel's C compiler) and set the paths to OpenBLAS and LAPACKE by
   modifying these two lines:
```
OPENBLAS=/n/app/openblas/0.2.19
LAPACKE=/n/app/lapacke/3.6.1
```
   Add the OpenBLAS library path to the linker path.
```
# Should be the same as $OPENBLAS with /lib appended
$ export LD_LIBRARY_PATH=/path/to/openblas/lib  
```
   Compile the approximator program.
```
$ cd gridfit-gauss
$ make -f Makefile.gcc   # or -f Makefile.icc if Intel C compiler is available
```
3. Add the helper scripts and laplace approximator binary to the global path
```
# From the git repo root
$ cd ..
$ export PATH=$PATH:`realpath scripts`:`realpath bin`
```
4. Make a GATK .jar file and relevant databases available in a single path.
   You may place this directory anywhere you'd like. It is important,
   however, that the files in this directory are named as stated below
   (i.e., gatk.jar, dbsnp.vcf, human_g1k_v37_decoy.fasta).
```
$ mkdir gatkpath
$ cd gatkpath
# Put GATK (preferably version 3.8) here
$ cp /path/to/gatk.jar gatk.jar
# Put a copy of the human reference genome, GRCh37 with decoy here
# All 3 files are necessary: the fasta, fasta.fai and dict.
$ cp /path/to/ref.fasta human_g1k_v37_decoy.fasta
$ cp /path/to/ref.fasta.fai human_g1k_v37_decoy.fasta.fai
$ cp /path/to/ref.dict human_g1k_v37_decoy.dict
# Put a copy of dbSNP (in VCF format) here
$ cp /path/to/dbsnp dbsnp.vcf
$ export GATK_PATH=`pwd`
```
5. Install SHAPEIT2 and the 1000 genomes haplotype panel.
    * Download and unzip SHAPEIT.  The path to the top level of the unzipped archive
      is `SHAPEIT_ROOT`.
        e.g., https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r904.glibcv2.12.linux.tar.gz
    * Download and unzip the 1000 genomes haplotype panel. The path to the top
      level of the unzipped archive is `REFPANEL_ROOT`.
        e.g., https://mathgen.stats.ox.ac.uk/impute/ALL.integrated_phase1_SHAPEIT_16-06-14.nosing.tgz
    * Set the following two environment variables.
```
$ export SHAPEIT_ROOT=/path/to/shapeit
$ export REFPANEL_ROOT=/path/to/refpanel
```



## Running the demo
The provided demo shows how the pipeline is run in practice. However,
because parameter fitting (step 3) requires significant compute time,
the demo only analyzes reads from chromosome 22. Using a server with 8
cores, the demo took approximately **X HOURS** to run.
Because some properties of hSNPs and sSNVs (e.g., VAF distributions) are
measured from sites across all chromosome, the demo's output will not exactly
reconstitute the findings reported in the manuscript.

In each step, the dependency versions refer to the **tested** versions.
Other versions may work as well.


### Step 0: Download data
1. Create a directory to contain the demo files and outputs.
```
$ cd /root/of/git/repo
$ mkdir demo
```
2. Download demo BAM files. Save the downloaded files to the `demo` directory.
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


### STEP 1: Run GATK HaplotypeCaller on single cell and matched bulk data
1. Configure the number of threads you wish to use for GATK by editing
   `scripts/run_gatk_demo.sh` and replacing the indicated line:
```
ncores=8    # Replace this line with the desired number of cores
mem=22G     # If using ncores > 1, increase ~linearly up to ~24G
```
**8 cores can process the chr22 data in approximately 1 hour per run.**

2. Run GATK two times, once using MMQ=60 and once using MMQ=1. The script
   will automatically create the directory `output_dir` to hold the generated
   VCF files.
```
$ cd /path/to/git/repo
$ run_gatk_demo.sh 60 demo demo/hunamp.chr22.bam demo/il-12.chr22.bam
$ run_gatk_demo.sh 1 demo demo/hunamp.chr22.bam demo/il-12.chr22.bam
```


### STEP 2: Define and phase hSNPs
1. Run SHAPEIT2 on potential heterozygous SNVs found in the bulk sample. Note
   that potential hSNPs are taken only from the MMQ=60 GATK output.
```
run_shapeit.sh demo/hc_raw.mm60.vcf demo hunamp 22
```
   If successful, a file named `phased_hsnps.chr22.vcf` should exist in `demo`.


## STEP 3: Fit covariance function parameters via grid search
1. Create training data files.
```
# NOTE: the sample name "h25" corresponds to il-12
get_hsnps_singlecell.sh h25 demo/hc_raw.mmq60.snp.vcf demo/phased_hsnps.chr22.vcf demo
```

2. Run the grid searching algorithm to find parameter MLEs. Using the parameters
   below (ngrids=8, points-per-grid=400), the run took ~30 minutes. The parameters
   used in the manuscript were ngrids=20, points-per-grid=1000.
```
mkdir -p output_dir/gridfit/chr22

gridfit_chr.py \
    --bindata=demo/training_chr22.bin \
    --local \
    --ngrids 8 \                 # Set this to the number of cores available
    --points-per-grid 400 \
    --resume \
    --laplace=laplace_cpu_gcc \  # Or laplace_cpu_icc if Intel C compiler was used
    --outprefix=demo/gridfit/chr22
```
   If the run is successful, the file `output_dir/gridfit/chr22/fit.rda` will be created.

3. Make the final fit file.
```
make_fits.R demo/gridfit demo/fits.rda
```


## STEP 4: Run SCAN-SNV
1. Convert GATK VCFs into table format.
```
# Make expected symlinks to the single cell BAMs
$ cd demo
$ ln -s il-12.chr22.bam wg.bam
$ ln -s il-12.chr22.bam.bai wg.bam.bai
$ cd ..
```
2. Run SCAN-SNV
```
# Run SCAN-SNV
$ scan_snv.sh demo/hc_raw.mmq60.vcf demo/hc_raw.mmq1.vcf demo h25 hunamp demo 0.1
```
