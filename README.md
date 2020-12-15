# SCAN2
Genotyper for **somatic** SNV and indel discovery in PTA-amplified single cells.

SCAN2 should not be used for genotyping germline mutations, as it excludes any
mutation with any read support in matched bulk samples.

**IMPORTANT**: SCAN2 should only be applied to diploid chromosomes (i.e., human
autosomes 1-22 and the X chromosome in female samples).



# Installation
SCAN2 is distributed as a conda package. Installation requires the conda
package management tool and a Linux-flavored OS.

**Operating systems tested**
* GNU/Linux, kernel version 3.10.0, CentOS 7. Note that pre-compiled SHAPEIT2 binaries are only made available for Linux systems, although in principle other phasing algorithms can be used instead.
* Ubuntu 16.04.4 LTS AWS instance.


## Installing miniconda
```
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
# Accept the license by typing "yes"
# Choose an install prefix (the default is often fine)
# Choose to run conda init (enter yes a second time during script)
# Log-out and back in to source .bashrc and put conda on $PATH
```

## Installing SCAN2
Create a conda environment for SCAN2 and install necessary packages
```
$ conda create -n scan2 -c conda-forge -c bioconda -c jluquette -c dranew scan2 shapeit

# Activate the scan2 conda environment
$ conda activate scan2
```

## Downloading external data dependencies
SCAN2 has been tested on the NCBI human reference build 37, **support is being added for hg38**.

Download reference genome.
```
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.fasta.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.fasta.fai.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.dict.gz
```

Download dbSNP **common variants**.
Note that dbSNP build 147 (common variants only) was used in
the publication. However, NCBI
does not guarantee long term hosting of dbSNP builds, so we recommend
downloading the version of dbSNP included in the Broad's GATK resource
bundle. To use other builds of dbSNP, you will need to generate a tribble
index (see below).
```
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.idx.gz
```

Download SHAPEIT's haplotype reference panel.
```
$ wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz
$ wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3_chrX.tgz
```

Unzip everything and move the chrX SHAPEIT files into the main SHAPEIT
directory.
```
$ gunzip *.gz
$ tar xzvf 1000GP_Phase3.tgz
$ tar xzvf 1000GP_Phase3_chrX.tgz
$ mv genetic_map_chrX_* 1000GP_Phase3_chrX* 1000GP_Phase3
```

# Running the SCAN2 demo
Download the demo chr22 BAMs. Our PTA data is only available through
protected access at dbGaP. However, any sequencing data (whether
single cell or not) can be used to test the pipeline installation.

We provide the following publically available data for MDA-amplified
single cells (Dong et al *Nature Methods* 2017) for the demo:
```
$ wget http://compbio.med.harvard.edu/scan-snv/hunamp.chr22.bam
$ wget http://compbio.med.harvard.edu/scan-snv/hunamp.chr22.bam.bai
$ wget http://compbio.med.harvard.edu/scan-snv/il-12.chr22.bam
$ wget http://compbio.med.harvard.edu/scan-snv/il-12.chr22.bam.bai
```

Run SCAN2. Replace instances of /path/to/... with the paths
downloaded above. This demo runs in less than 5 minutes on an 8 core
machine by restricting analysis to a 500 KB segment of chr22 and by
using an impractically coarse grid for covariance function fitting.
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
    --abmodel-chunks=1 \
    --abmodel-samples-per-chunk=100 \
    --abmodel-steps=1 \
    --callable-regions True \
    --score-all-sites \
    --regions 22:30000001-30500000 \
    --bulk-bam /path/to/hunamp.chr22.bam \
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
Rdata format `demo/snv/[single_cell_sample_name]/somatic_genotypes.rda`.

Run R and load the RData file. IMPORTANT! Use the R installation in the SCAN2 conda
environment for the SCAN2 R library.
```
$ R   
library(scan2)
#
# Attaching package: 'scan2'
#
# The following object is masked from 'package:stats':
#
#    df

# 'h25' is the name of the single cell sample
load('demo/snv/h25/somatic_genotypes.rda')
```

We now provide an S4 class called 'SCAN2' that handles both the
genotyping logic and allows convenient user interaction. Each
single cell in a SCAN2 run will have a 'gt' object saved in the
corresponding somatic_genotypes.rda file.


The S4 class below provides a summary of several features of the data.
```
gt
# SCAN2 
#   Single cell ID: h25 
#   Bulk ID: hunamp 
#   GATK: 429 raw sites
#   GATK with low mapping quality: 429 raw sites
#   AB model training hSNPs: 92 phased sites (hap1=57, hap2=35) 
#       VAF correlation between neighboring hSNPs:
#           <100 bp 0.993 <1000 bp 0.842 <10 kbp 0.587 <100 kbp 0.572 
#         49 resampled hSNPs
#   Allele balance:
#       mean (0 is neutral): 0.212 
#       uncertainty (Q25, median, Q75): 0.289 0.801 0.895 
#       mean at training hSNPs: -0.028 
#       correlation with VAF at training hSNPs 0.815 
#   Mutation models: computed
#   CIGAR data: 429 sites
#   Static filters: 0 retained 429 removed 0 NA
```

 Now it is possible to extract the genotyped sites called using the
static filters (such
 as minimum depths, exclusion of dbSNP sites, etc.) and the two
 single cell artifact models. lysis.fdr and mda.fdr <= 0.01 correspond
 to a target false discovery rate (FDR) of 1%. Note that SCAN2 does
 **NOT** provide formal FDR control.


 There should be no called somatic SNVs in this demo dataset.
```
snvs <- df(gt)
snvs[snvs$static.filter & snvs$lysis.fdr <= 0.01 & snvs$mda.fdr <= 0.01,]
 [1] chr               pos               dbsnp             refnt            
 [5] altnt             mq                mqrs              h25              
 [9] scref             scalt             hunamp            bref             
[13] balt              dp                af                bulk.dp          
[17] bulk.af           training.site     scref             scalt            
[21] bref              balt              ab                gp.mu            
[25] gp.sd             abc.pv            lysis.pv          lysis.beta       
[29] mda.pv            mda.beta          M.cigars          ID.cigars        
[33] HS.cigars         other.cigars      dp.cigars         M.cigars.bulk    
[37] ID.cigars.bulk    HS.cigars.bulk    other.cigars.bulk dp.cigars.bulk   
[41] id.score.y        id.score.x        hs.score.y        hs.score.x       
[45] id.score          hs.score          static.filter     nt               
[49] na                lysis.fdr         mda.fdr          
<0 rows> (or 0-length row.names)
```

**NOTE**: VCF output is forthcoming.
```
# Called sSNVs can be extracted from the data frame via
R> load('demo/snv/[single_cell_sample_name]/somatic_genotypes.rda')
R> somatic[somatic$pass,]
# The demo should not produce any passing variants.
```

### WARNING!
* The conda environment (named scan2 in these instructions) must always
  be `conda activate`d before running SCAN2.

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
      the slurm-drmaa 1.1.1 package.
    * If your scheduler is not DRMAA-compatible (or if the appropriate
      DRMAA interface is unavailable), Snakemake's `--cluster` option
      offers similar functionality to `--drmaa`, but with fewer features.
* Cloud environments.

See Snakemake's documentation for more details on cluster and cloud execution: https://snakemake.readthedocs.io/en/stable/.


# Using a custom dbSNP version
## Generating a Tribble index for dbSNP
dbSNP VCFs must be indexed by Tribble (*not* tabix) for GATK. The dbSNP
found in the GATK's resource bundle is already indexed. If you wish to use
a different dbSNP version, the file can be indexed by `igvtools`.

```
$ mamba install -c bioconda igvtools
$ igvtools index /path/to/your/dbsnp.vcf
```
