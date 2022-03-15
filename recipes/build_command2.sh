#!/bin/bash
#SBATCH -p priopark
#SBATCH -A park_contrib
#SBATCH -t 124:00:00
#SBATCH --mem=32G

/usr/bin/time conda-build --override-channels --R 3.5.1 -c conda-forge -c bioconda -c jluquette -c dranew -c soil scan2
