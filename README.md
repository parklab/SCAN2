# SCAN2
Genotyper for **somatic** SNV and indel discovery in PTA-amplified single cells.

SCAN2 should not be used for genotyping germline mutations, as it excludes any
mutation with any read support in matched bulk samples.

**IMPORTANT**: SCAN2 should only be applied to diploid autosomes (i.e., human
autosomes 1-22). Sex chromosomes are currently not supported, but in principle
SCAN2 is applicable to female X chromosomes.

> **Note**
>
> The instructions below guide the user through a minimal demo of the
> SCAN2 tool with the primary goal of checking the installation. For a more
> detailed walkthrough including cross-sample panel building and mutation
> signature-based rescue, please [see the wiki](../../wiki).

## Clonal vs. non-clonal somatic mutation detection
SCAN2 currently will not call any mutation candidate with **any** supporting
reads in the matched bulk sample. At standard WGS sequencing depths of 30x,
this excludes nearly all germline variants, the majority of early
post-zygotic mutations and, if the matched bulk is from a closely related
tissue, many lower frequency clonal somatic mutations as well.
We therefore often refer
to SCAN2 calls as "non-clonal," although this is only approximately the
case since clonal mutations with low frequency may also have 0 reads in
the matched bulk.

It is not necessarily a good
idea to allow mutations with >0 reads in bulk since it may allow some
germline heterozygous mutations or relatively rare artifacts (of
sequencing, alignment, etc.) to be miscalled as somatic mutations. This
is particularly true in single cells with very low (i.e., <100) somatic
mutation burden.

This restriction may be removed in the future after further testing.


# License
SCAN2 is freely available for non-commercial use.



# Installation
SCAN2 is distributed as a conda package. Installation requires the conda
package management tool and a Linux-flavored OS.

**Operating systems tested**
* GNU/Linux, kernel version 3.10.0, CentOS 7.


## Install miniconda
**IMPORTANT** Ensure that your chosen install prefix has sufficient
disk space: SCAN2 requires ~15G of disk space to analyze GRCh37 samples and an
additional ~10G is required for each extra genome.

```
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
# Accept the license by typing "yes"
# Choose to run conda init (enter yes a second time during script)
```
Log-out and back in to source .bashrc and put conda on $PATH


## Install SCAN2
Create a conda environment for SCAN2 and install necessary packages
```
conda create -n scan2 -c conda-forge -c bioconda -c jluquette -c dranew -c soil scan2

# Once the environment is created, it must be activated every time
# you wish to use SCAN2.
conda activate scan2
```

## Install reference genome files for SigProfilerMatrixGenerator
SigProfilerMatrixGenerator is used to classify indels into the 83-class indel mutation
signature format ID83. This classification requires a reference genome to determine
the sequence context around each indel. The following command shows how to install a
SigProfilerMatrixGenerator reference genome for GRCh37; GRCh38 is also available. See
https://github.com/AlexandrovLab/SigProfilerMatrixGenerator for details on other
genome assemblies and custom assemblies.

**N.B.** this installation step is required even if indels are not analyzed.

**REMEMBER** the scan2 conda environment must be activated (see above) when
running these commands!
```
$ python
from SigProfilerMatrixGenerator import install as genInstall
genInstall.install('GRCh37', rsync=False, bash=True)
quit()
```


# Download external data dependencies
SCAN2 has been tested on the NCBI human reference build 37 and hg38.

## Human genome version GRCh37
The SCAN2 demo requires GRCh37.

### Human reference version GRCh37 with decoy
Download the human reference genome if needed. The reference genome should match
the genome used for read alignment.
```
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.fasta.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.fasta.fai.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.dict.gz
gunzip *.gz
```

### dbSNP **common variants**
dbSNP build 147 (common variants only) was used in the publication. It is critical to use the database of **common variants**; the full dbSNP often intersects with somatic mutations and will cause reduced sensitivity. 
Instructions for generating the required tribble index are provided at the bottom of this page.
```
wget https://ftp.ncbi.nih.gov/snp/pre_build152/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz
gunzip common_all_20180423.vcf.gz
# Create Tribble index (see bottom of this page)
```
This VCF can be filtered to match build 147 (used in Luquette et al. 2022).
Each line contains a dbSNPBuildID=XXX tag; simply filter for XXX <= 147.

### SHAPEIT's haplotype reference panel
```
wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz
wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3_chrX.tgz
```

Unzip everything and move the chrX SHAPEIT files into the main SHAPEIT
directory.
```
tar xzvf 1000GP_Phase3.tgz
tar xzvf 1000GP_Phase3_chrX.tgz
mv genetic_map_chrX_* 1000GP_Phase3_chrX* 1000GP_Phase3
```

## Human reference version hg38
To run the SCAN2 demo, you will need the **GRCh37** genome, not GRCh38.

NOTE: GRCh38 alignments use the 'chr' prefix for chromosome names (i.e., chr1,
chr2, ..., rather than 1, 2, ...).

### Internal R reference genome
```
conda install -c conda-forge -c bioconda bioconductor-bsgenome.hsapiens.ucsc.hg38
```

### SigProfilerMatrixGenerator reference genome
```
$ python
from SigProfilerMatrixGenerator import install as genInstall
genInstall.install('GRCh38', rsync=False, bash=True)
quit()
```

### Human reference version GRCh38
Download the human reference genome if needed. The reference genome should match
the genome used for read alignment.
```
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict
```

### dbSNP **common variants**
Instructions for generating the required tribble index are provided at the bottom of this page.
```
wget https://ftp.ncbi.nih.gov/snp/pre_build152/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz
gunzip common_all_20180418.vcf.gz
# IMPORTANT: add a 'chr' prefix to sites. The dbSNP VCF currently
# does not contain ##contig=<ID=chrXXX... header lines, so we only need
# to update sites. This may change in the future.
sed  -e 's/\(^[^#]\)/chr\1/' common_all_20180418.vcf > common_all_20180418.chrprefix.vcf
# Create Tribble index (see bottom of this page)
```

### Eagle2 phasing panel
The user must create a reference panel manually for Eagle phasing. SCAN2 provides
a script (`scan2_download_eagle_refpanel.sh`) to do this, which essentially
implements the strategy recommended by Eagle2's authors (see
https://alkesgroup.broadinstitute.org/Eagle/#x1-320005.3.2). The script is located
in the same directory as the scan2 binary (`which scan2`), but is
runnable by simply typing the script's name when the scan2 conda environment is
activated.
```
mkdir eagle_1000g_panel
cd eagle_1000g_panel

for chr in {1..22} X; do  
    scan2_download_eagle_refpanel.sh path/to/Homo_sapiens_assembly38.fasta $chr
done
```

Download Eagle's genmap file:
```
wget https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz
```

Finally, you must supply the following arguments to `scan2 config` to use Eagle instead of SHAPEIT:
* `--phaser eagle`
* `--eagle-refpanel path/to/eagle_1000g_panel` - this should be the directory created by the commands above.
* `--eagle-genmap /path/to/genetic_map_hg38_withX.txt.gz`



# Running the SCAN2 demo
## Download the demo chr22 BAMs
These files are aligned to GRCh37, so one
must downloaded the GRCh37 external data files as directed above to run
the demo.

We provide two MDA-amplified (not PTA-amplified!) single cell BAMs and
one matched bulk from Dong et al. (*Nature Methods* 2017) for the demo.
These are chosen (rather than our PTA data) because they are publicly
available on SRA while our PTA cells are only available only through
protected access at dbGaP.

```
wget http://compbio.med.harvard.edu/scan-snv/hunamp.chr22.bam
wget http://compbio.med.harvard.edu/scan-snv/hunamp.chr22.bam.bai
wget http://compbio.med.harvard.edu/scan-snv/il-11.chr22.bam
wget http://compbio.med.harvard.edu/scan-snv/il-11.chr22.bam.bai
wget http://compbio.med.harvard.edu/scan-snv/il-12.chr22.bam
wget http://compbio.med.harvard.edu/scan-snv/il-12.chr22.bam.bai
```

## Run SCAN2
Replace instances of /path/to/... with the paths
downloaded above. This demo runs in ~15 minutes on a 20 core
machine by restricting analysis to a 2 MB segment of chr22.
```
# Creates a directory for the analysis
scan2 -d demo init
cd demo

# Configure analysis parameters. Can be run multiple times to change
# parameters.
scan2 config \
    --verbose \
    --ref /path/to/human_g1k_v37_decoy.fasta \
    --dbsnp /path/to/dbsnp_138_b37.vcf \
    --shapeit-refpanel  /path/to/1000GP_Phase3 \
    --regions 22:30000000-30999999,22:31000000-31999999 \
    --bulk-bam /path/to/hunamp.chr22.bam \
    --sc-bam /path/to/il-11.chr22.bam 
    --sc-bam /path/to/il-12.chr22.bam 

# Attempts to validate the chosen parameters and input files. Not
# exhaustive, but helpful to prevent errors.
scan2 validate

# Run the analysis. Set the number of cores you wish you use; for
# more advanced parallelization on a cluster or the cloud, see the
# sections below.
scan2 run --joblimit <N cores>
```
See `scan2 -h` for more details on configuration options and runtime arguments.

After SCAN2 completes, single sample results are stored in 
Rdata format `demo/call_mutations/[single_cell_sample_name]/scan2_object.rda`.

Run R and load the RData file. IMPORTANT! Use the R installation in the SCAN2 conda
environment for the SCAN2 R library.
```
$ R   
R> library(scan2)
# The two single cell samples are 'h24' (il-11) and 'h25' (il-12)
R> load('demo/call_mutations/h25/scan2_object.rda')
```

We now provide an S4 class called 'SCAN2' that handles both the
genotyping logic and allows convenient user interaction. Each
single cell in a SCAN2 run will have a 'results' object saved in the
corresponding scan2_object.rda file.


The S4 class provides a summary of several features of the data.
```
R> results
# SCAN2 
#   Region:
GRanges object with 22 ranges and 0 metadata columns:
       seqnames      ranges strand
          <Rle>   <IRanges>  <Rle>
   [1]        1 1-249250621      *
   [2]        2 1-243199373      *
   [3]        3 1-198022430      *
   [4]        4 1-191154276      *
   [5]        5 1-180915260      *
   ...      ...         ...    ...
  [18]       18  1-78077248      *
  [19]       19  1-59128983      *
  [20]       20  1-63025520      *
  [21]       21  1-48129895      *
  [22]       22  1-51304566      *
  -------
  seqinfo: 22 sequences from GRCh37.p13 genome
#   Genome: hs37d5 
#   Single cell ID: h25 
#   Bulk ID: hunamp 
#   GATK: 2517 raw sites
#   AB model training hSNPs: 753 phasing: 0|1=411, 1|0=342
#       VAF correlation between neighboring hSNPs:
#           <10 bp 1 <100 bp 0.983 <1000 bp 0.897 <10 kbp 0.806 <100 kbp 0.203 
#         36 resampled hSNPs
#         7 resampled hIndels
#   AB model parameters:
#       average (over chromosomes): a=-2.374, b=669.314, c=1.125, d=10193.908
#   Allele balance:
#       mean (0 is neutral): -0.073 
#       uncertainty (Q25, median, Q75): 0.372 0.558 1.162 
#       mean at training hSNPs: 0.058 
#   Mutation models: computed
#   CIGAR data: single cell: 2517 sites, bulk: 2517 sites
#   Static filters: 
#          snv:       12 retained     1873 removed        0 NA
#        indel:        7 retained      625 removed        0 NA
#   FDR prior data: 
#          snv:       25 candidates      657 germline hets        4 max burden
#        indel:       13 candidates      108 germline hets        0 max burden
#   Depth profile: 
#          snv: genome <:   min. bulk DP 1.3%,   min. sc DP 8.6%,   either 9.4%
#        indel: genome <:   min. bulk DP 1.3%,   min. sc DP 15.4%,   either 15.9%
#   Somatic mutation calls: mode=new, target.fdr=0.010
#          snv:        0 called       11 resampled training calls
#        indel:        0 called        1 resampled training calls
#       ALL indel calls have been suppressed (insufficient number of single cells in cross-sample filter
#   Somatic mutation burden: 
#          snv:     NA somatic,   NA% sens,   NA callable Gbp,   NA muts/haploid Gbp,   NA muts per genome
#        indel:     NA somatic,   NA% sens,   NA callable Gbp,   NA muts/haploid Gbp,   NA muts per genome
#       indel mutation burden estimates ARE NOT VALID (cross-sample panel insufficiency)!
#   Mutation rescue by signature: (not computed)
```

This SCAN2 object contains only VAF-based somatic mutation calls. For
signature-based rescue, see `scan2 rescue -h`. Further, this SCAN2 object
will not contain indel calls because no cross sample filter was provided
(see `scan2 config --help`, `--cross-sample-panel`).

There should be no called somatic SNVs in this demo dataset.
```
R> snvs <- results@gatk[pass == TRUE & muttype == 'snv']
R> snvs
Empty data.table (0 rows and 73 cols): chr,pos,dbsnp,refnt,altnt,mq...
```


# Parallelization
In a practical setting, parallelization will be required. SCAN2 leverages
Snakemake to offer parallelization via any of the following:
* A single machine with multiple cores. To do this, increase the
  `--joblimit` parameter but do not invoke cluster or cloud
  arguments.  A total memory limit can also be set via
  `--memlimit`, which should be supplied in megabytes.
* Clusters with distributed resource management software (e.g., SLURM,
  LSF, GridEngine). If your resource management software supports
  DRMAA (https://www.drmaa.org), the authors recommend using the
  `--drmaa`. **However, additional libraries may be necessary to
  interface with the DRMAA wrapper.**
    * E.g., in the SCAN2 publication, a SLURM cluster was accessed via
      the slurm-drmaa 1.1.1 package (https://github.com/natefoo/slurm-drmaa). This
      package was maintained by local cluster system admins.
      An example command line for running SCAN2 with DRMAA follows. Note the `--mem`
      argument, which uses Snakemake's `{resources.mem}` placeholder to allow each
      individual job in the SCAN2 pipeline to specify the amount of memory needed and
      the `{threads}` placeholder.
        ```
        scan2 run --joblimit 100 --drmaa ' -p YOUR_QUEUE_NAME -A YOUR_ACCOUNT -t 10:00 -c {threads} --mem={resources.mem}'
        ```
    * If your scheduler is not DRMAA-compatible (or if the appropriate
      DRMAA interface is unavailable), Snakemake's `--cluster` option
      offers similar functionality to `--drmaa`, but with fewer features.
* Cloud environments.

See Snakemake's documentation for more details on cluster and cloud execution: https://snakemake.readthedocs.io/en/stable/.


# Generating a Tribble index for dbSNP
dbSNP VCFs must be indexed by Tribble (*not* tabix) for GATK. The dbSNP
found in the GATK's resource bundle is already indexed. If you wish to use
a different dbSNP version (as in the above section on downloading external
data dependencies), use IGVtools as detailed below:

```
# IMPORTANT! IMPORTANT! Do not install igvtools into your SCAN2 conda
# environment. Create a new environment for igvtools! At the time of this writing,
# igvtools depends on an incompatible version of the JDK and will break the
# GATK installation.
conda deactivate                               # exit the SCAN2 environment
conda create -n igvtools -c bioconda igvtools  # create a new environment with only igvtools installed
conda activate igvtools
igvtools index /path/to/your/dbsnp.vcf
conda deactivate                               # exit the igvtools environment
conda activate scan2                           # reactivate the SCAN2 environment
```
