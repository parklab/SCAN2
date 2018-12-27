module load gcc
module load samtools
module load R/3.3.3
module load java
export PATH=$PATH:`realpath bin`:`realpath scripts`
export GATK_PATH=`realpath gatkpath`
export SHAPEIT_ROOT=`realpath shapeit/shapeit.v2.904.2.6.32-696.18.7.el6.x86_64`
export REFPANEL_ROOT=`realpath shapeit/1000GP_Phase3`
export DRMAA_LIBRARY_PATH=/n/app/slurm-drmaa/1.0.7/lib/libdrmaa.so
check_env.sh
