# scan-snv
Single cell somatic genotyper


## Installation
**Operating systems tested**: GNU/Linux, kernel version 3.10.0, CentOS 7. Note that pre-compiled SHAPEIT2 binaries are only made available for Linux systems, although in principle other phasing algorithms can be used instead.

**Installation time**: installation should take only a few minutes on a modern
computer.

Version numbers in parentheses are the tested versions and were used to produce
the results in the manuscript. They are not necessarily required to run.

**Dependencies**: Python (v2.7), R (v3.3.3), LAPACKE (v3.6.1), OpenBLAS (v0.2.19),
    Java (v1.8), GATK (v3.8-0-ge9d806836), SHAPEIT2 (v2-r837), samtools (v1.3)

**Optional dependencies**: Intel C compiler (2016)

**Data dependencies**:

* Human reference genome version GRCh37 with decoy
  e.g. ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.fasta.gz and associated \*.fai and \*.dict files.
* dbSNP (v147, b37)
  e.g. ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/common_all_20170710.vcf.gz
* SHAPEIT2 formatted 1000 genomes reference haplotype panel
  e.g., https://mathgen.stats.ox.ac.uk/impute/ALL.integrated_phase1_SHAPEIT_16-06-14.nosing.tgz

> **IMPORTANT!** The environment variables listed below (`LD_LIBRARY_PATH`,
> `PATH`, `GATK_PATH`, `SHAPEIT_ROOT`, `REFPANEL_ROOT`) must always be set
> appropriately before running SCAN-SNV. The `samtools`, `Rscript`, `python`
> and `java` binaries must also be in your $PATH.

1. Build and install the SCAN-SNV R package
```
$ cd /root/of/git/repo
$ R CMD build rpkg
$ R CMD INSTALL scansnv_0.1.tar.gz
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
$ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/openblas/lib  
```
   Compile the approximator program.
```
$ cd gridfit-gauss
$ make -f Makefile.gcc   # or -f Makefile.icc if Intel C compiler is available
```
3. Add the helper scripts and laplace approximator binary to the global path
```
# Return to the git repo root
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
$ cp /path/to/gatk.jar gatk.jar
$ cp /path/to/ref.fasta human_g1k_v37_decoy.fasta
$ cp /path/to/ref.fasta.fai human_g1k_v37_decoy.fasta.fai
$ cp /path/to/ref.dict human_g1k_v37_decoy.dict
$ cp /path/to/dbsnp dbsnp.vcf
$ cp /path/to/dbsnp/index dbsnp.vcf.idx
$ export GATK_PATH=`pwd`
```
5. Install SHAPEIT2 and the 1000 genomes haplotype panel.
    * Download and unzip SHAPEIT.  The path to the top level of the unzipped archive
      is `SHAPEIT_ROOT`.
    * Download and unzip the 1000 genomes haplotype panel. The path to the top
      level of the unzipped archive is `REFPANEL_ROOT`.
    * Set the environment variables `SHAPEIT_ROOT` and `REFPANEL_ROOT`.
```
# For example, using the version of SHAPEIT2 available on 10/22/2018
$ cd /path/to/git/repo
$ mkdir shapeit
$ cd shapeit
$ wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r904.glibcv2.12.linux.tar.gz
$ tar xzvf shapeit.v2.r904.glibcv2.12.linux.tar.gz
$ wget https://mathgen.stats.ox.ac.uk/impute/ALL.integrated_phase1_SHAPEIT_16-06-14.nosing.tgz
$ tar xzvf ALL.integrated_phase1_SHAPEIT_16-06-14.nosing.tgz
$ export SHAPEIT_ROOT=`realpath shapeit.v2.904.2.6.32-696.18.7.el6.x86_64`
$ export REFPANEL_ROOT=`realpath ALL.integrated_phase1_SHAPEIT_16-06-14.nosing`
```



# Running the demo
**Approximate run time**: 5 hours with 8 threads.

The provided demo shows how the pipeline is run in practice. However,
because parameter fitting (step 3) requires significant compute time,
the demo only analyzes reads from chromosome 22.
Because some properties of hSNPs and sSNVs (e.g., VAF distributions) are
normally measured from sites across all chromosomes, the demo's output will
not exactly match the findings reported in the manuscript. In addition, a
smaller number of grid points are used in the parameter MLE estimation,
meaning the correlation function will also differ from the manuscript.
For reference, the demo script computes 250 points times the number of
requested threads, while 20,000 points per chromosome are used for normal
analyses.

Before running `demo.sh`, ensure that environment variables are set as
described in **Installation**. You should also edit the `demo.sh` file to
reflect the number of CPU cores you wish to use. Simply modify the line
```
ncores=8
```

To run:
```
$ cd /path/to/scan-snv
$ demo.sh
```

The final output is stored in Rdata format in `demo/scan-snv/somatic_gt.rda`.
The set of final, called sSNVs can be printed using
```
$ Rscript -e 'load("demo/scan-snv/somatic_gt.rda"); gt$somatic[gt$somatic$pass,]'
```
The expected output is a matrix of ~5-10 PASSed variants with several
covariates.



# Step by step usage
### STEP 1: Run GATK HaplotypeCaller on single cell and matched bulk data
1. Configure parallelism for GATK. Region-based parallelization was used for the
   manuscript's analysis, although an alternative and simpler method to
   parallelize is to increase thread count.
   * Region-based parallelization is only supported for clusters with
     SLURM installed.
   * The first parameter of the `run_gatk.sh` script sets the number of
     compute threads for GATK HaplotypeCaller.
2. Run GATK two times: once using MMQ=60 and once using MMQ=1.
```
$ run_gatk.sh n_threads 60 output_directory input_bam1 input_bam2 [ input_bam3 ... ]
$ run_gatk.sh n_threads 1 output_directory input_bam1 input_bam2 [ input_bam3 ... ]
```


### STEP 2: Define and phase hSNPs
1. Run SHAPEIT2 on potential heterozygous SNPs found in the bulk sample. The
   `hc_raw.mmq60.vcf` below refers to the VCF produced by `run_gatk.sh 60 ...`.
   `bulk_sample_name` should be the sample string in the VCF corresponding to
   the bulk sample.
```
$ for chrom in `seq 1 22`; do
>   run_shapeit.sh hc_raw.mmq60.vcf output_directory bulk_sample_name $chrom
> done
```
   If successful, files named `phased_hsnps.chr[1-22].vcf` should exist in
   the specified output directory.
2. Combine the chromosome-specific VCFs into a single VCF. For example, using
   GATK:
```
$ java -cp /path/to/gatk.jar org.broadinstitute.gatk.tools.CatVariants \
      -R /path/to/human/reference.fasta  \
      $(for chrom in `seq 1 22`; do echo -V phased_hsnps.chr$chrom.vcf; done) \
      -out phased_hsnps.vcf \
      -assumeSorted
```


### STEP 3: Fit covariance function parameters via grid search
1. Create training data files. `single_cell_sample` is the string found in
   the mmq60 VCF corresponding to the single cell sample.
```
$ get_hsnps_singlecell.sh single_cell_sample hc_raw.mmq60.snp.vcf phased_hsnps.vcf output_directory
```

2. Run the grid searching algorithm to find parameter MLEs. In a typical
   analysis, we use 20,000 points per grid level (ngrids=20,
   points-per-grid=1000). For each chromosome, run:
```
mkdir -p output_directory/gridfit/chr22

gridfit_chr.py \
    --bindata=training_chr22.bin \
    --local \
    --ngrids 8 \                 # Set this to the number of cores available
    --points-per-grid 400 \
    --resume \
    --laplace=laplace_cpu_gcc \  # Or laplace_cpu_icc if Intel C compiler was used
    --outprefix=output_directory/gridfit/chr22
```
   If the run is successful, the files
   `output_directory/gridfit/chr[1-22]/fit.rda` will be created.

3. Make the final fit file.
```
make_fits.R demo/gridfit demo/fits.rda
```


### STEP 4: Run SCAN-SNV
1. Convert GATK VCFs into table format.
```
# Make expected symlinks to the single cell BAMs
$ cd output_directory
$ ln -s /path/to/single/cell/BAM wg.bam
$ ln -s /path/to/single/cell/BAI wg.bam.bai
$ cd ..
```
2. Run SCAN-SNV. The final parameter is the target FDR.
```
# Run SCAN-SNV
$ scan_snv.sh hc_raw.mmq60.vcf hc_raw.mmq1.vcf output_directory \
    single_cell_sample bulk_sample output_directory 0.1
```
