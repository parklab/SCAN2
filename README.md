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

## STEP 0. Installation and demo data
**Dependencies**: LAPACKE (v3.6.1), OpenBLAS (v0.2.19)\
**Optional dependencies**: Intel C compiler

1. Build and install the SCAN-SNV R package
```
$ cd /root/of/git/repo
$ R CMD build rpkg
$ R CMD INSTALL scansnv-0.1.tar.gz
```
2. Install the Laplace approximator. First edit gridfit-gauss/Makefile.gcc (if you are
   compiling with gcc) or gridfit-gauss/Makefile.icc (if you are using Intel's C
   compiler) and set the paths to OpenBLAS and LAPACKE by modifying these two lines
   to point to your local installation of OpenBLAS and LAPACKE.
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
$ export PATH=$PATH:`realpath scripts`:`realpath bin`
```
4. Create a directory to contain the demo files and outputs.
```
$ cd /root/of/git/repo
$ mkdir demo
```
5. Make a GATK .jar file and relevant databases available in a single path.
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
6. Install SHAPEIT2 and the 1000 genomes haplotype panel.
    * Download and unzip SHAPEIT (e.g., https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r904.glibcv2.12.linux.tar.gz).
      The unzipped path is `SHAPEIT_ROOT`.
    * Download and unzip the 1000 genomes haplotype panel (e.g., https://mathgen.stats.ox.ac.uk/impute/data_download_1000G_phase1_integrated_SHAPEIT2_16-06-14.html).
      The unzipped path is `REFPANEL_ROOT`.
```
$ export SHAPEIT_ROOT=/path/to/shapeit
$ export REFPANEL_ROOT=/path/to/refpanel
```
7. Download demo BAM files. Save the downloaded files to the `demo` directory.
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



## STEP 1: Run GATK HaplotypeCaller on single cell and matched bulk data
**Dependencies**: Java (v1.8), GATK (v3.8-0-ge9d806836)\
**Data dependencies**: human reference genome (GRCh37 with decoy; e.g. `ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.fasta.gz`), dbSNP (v147, b37: e.g. `ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/common_all_20170710.vcf.gz`)

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


## STEP 2: Define and phase hSNPs
**Dependencies**: SHAPEIT2 (v2, r837)
**Data dependencies**: SHAPEIT2 formatted 1000 genomes reference haplotype panel
1. Run SHAPEIT2 on potential heterozygous SNVs found in the bulk sample. Note
   that potential hSNPs are taken only from the MMQ=60 GATK output.
```
phase_chr.sh demo/hc_raw.mm60.vcf hunamp demo 22
```
   If successful, a file named `phased_hsnps.chr22.vcf` should exist in
   `demo`.


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
   below (ngrids=8, points-per-grid=400), the run took ~30 minutes. The parameters
   used in the manuscript were ngrids=20, points-per-grid=1000.
```
mkdir -p output_dir/gridfit/chr22

./gridfit_slurm/runchr.py \
    --bindata=output_dir/training_chr22.bin \
    --local \
    --ngrids 8 \          # Set this to the number of cores available
    --points-per-grid 400 \
    --resume \
    --combine=./gridfit_slurm/combine.R \
    --mkl-laplace=./mkl-gridfit-gauss/laplace_cpu_gcc \
    --outprefix=output_dir/gridfit/chr22
```
   If the run is successful, the file `output_dir/gridfit/chr22/fit.rda` will be created.

3. Make the final fit file.
```
Rscript ./gridfit_slurm/make_fits.R
```


## STEP 5: Run SCAN-SNV

1. Convert GATK VCFs into table format.
```
# A similar process has already been run for mmq60
java -jar /path/to/GATK -R /path/to/human_g1k_v37_decoy.fasta \
    -T SelectVariants -V output_dir/hc_raw.mmq1.vcf \
    -selectType SNP -restrictAllelesTo BIALLELIC -env -trimAlternates \
    -select 'vc.getGenotype("hunamp").isCalled()' \
    -o output_dir/hc_raw.mmq1.snp.vcf

totab.sh output_dir/hc_raw.mmq1.snp.vcf output_dir/mmq1.tab
totab.sh output_dir/hc_raw.mmq60.snp.vcf output_dir/mmq60.tab

# Make expected symlinks to the single cell BAMs
ln -s ../il-12.chr22.bam wg.bam
ln -s ../il-12.chr22.bam.bai wg.bam.bai
```
