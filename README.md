# SCAN2
Genotyper for _somatic_ SNV and indel discovery in PTA-amplified single cells.

SCAN2 should not be used for genotyping germline mutations, as it excludes any
mutation with any read support in matched bulk samples.



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
Create a conda environment for SCAN2
```
$ conda deactivate   # The "base" environment will be active after login
$ conda create -n scan2
$ conda activate scan2
```
Install the SCAN2 package
```
$ conda install -c bioconda -c conda-forge/label/cf201901 -c jluquette scan2
```
Register your GATK installation
```
$ wget 'https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-1-0-gf15c1c3ef' -O GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
$ tar xjvf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
$ gatk-register GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
# Test the install
$ gatk --version
# Above should print 3.8-1-0-gf15c1c3ef
```

## Downloading external data dependencies
SCAN2 has been tested on the NCBI human reference build 37.

Download reference genome.
```
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.fasta.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.fasta.fai.gz
$ wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.dict.gz
```

Download dbSNP *common variants*.
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
Download the demo chr22 BAMs.
```
$ wget http://compbio.med.harvard.edu/scan-snv/hunamp.chr22.bam
$ wget http://compbio.med.harvard.edu/scan-snv/hunamp.chr22.bam.bai
$ wget http://compbio.med.harvard.edu/scan-snv/il-12.chr22.bam
$ wget http://compbio.med.harvard.edu/scan-snv/il-12.chr22.bam.bai
```

Run SCAN2. Replace instances of /path/to/... with the paths
downloaded above. This demo runs in about 5 minutes on a single core
machine by restricting analysis to a 1 MB segment of chr22 and by
using an impractically coarse grid for covariance function fitting.
```
scan2 --analysis-dir demo init
scan2 --analysis-dir demo config \
    --verbose \
    --ref /path/to/human_g1k_v37_decoy.fasta \
    --dbsnp /path/to/dbsnp_138_b37.vcf \
    --shapeit-refpanel  /path/to/1000GP_Phase3 \
    --sc-bam il-12.chr22_30M.bam \
    --bulk-bam hunamp.chr22_30M.bam \
    --regions 22:30000001-30100000
scan2 --analysis-dir demo validate
scan2 --analysis-dir demo run
```

See `scansnv -h` for more details on arguments.

After SCAN-SNV completes, single sample results are available in the
Rdata file `demo/scansnv/[single_cell_sample_name]/somatic_genotypes.rda`.
SNVs that pass SCAN-SNV's calling thresholds will have `pass=TRUE` in the
`somatic` data frame (see below).

**NOTE**: a VCF output option is forthcoming.
```
# Called sSNVs can be extracted from the data frame via
R> load('demo/scansnv/[single_cell_sample_name]/somatic_genotypes.rda')
R> somatic[somatic$pass,]
# The demo should not produce any passing variants.
```

# WARNING!
* The conda environment (named scan2 in these instructions) must always
  be `conda activate`d before running SCAN2.
* Real world analyses will require parallelization.
    * On a machine with multiple cores, increasing the `--joblimit` parameter
      will run multiple parts of the analysis in parallel.
    * For clusters with distributed resource management software (e.g., SLURM),
      SCAN-SNV exposes Snakemake's parallelization options
      `--cluster` and `--drmaa`.


# Using a custom dbSNP version
## Generating a Tribble index for dbSNP
dbSNP VCFs must be indexed by Tribble (*not* tabix) for GATK. The dbSNP
found in the GATK's resource bundle is already indexed. If you wish to use
a different dbSNP version, the file can be indexed by `igvtools`.

```
$ conda install -c bioconda igvtools
$ igvtools index /path/to/your/dbsnp.vcf
```
