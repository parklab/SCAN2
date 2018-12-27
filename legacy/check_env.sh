#!/bin/bash

retval=0

echo "Checking environment for necessary variables, PATH augmentations and libraries."

echo "-------------------------------------------------------------------------------"

test_path() {
    x=$(which $1 2>/dev/null)
    if [ $? -ne 0 ]; then
        echo "ERROR: $1 not found in PATH"
        retval=1
    else
        echo "Using $x"
    fi
}

which laplace_cpu_gcc >& /dev/null
gcc=$?
which laplace_cpu_icc >& /dev/null
icc=$?
if [ $gcc -eq 1 ] && [ $icc -eq 1 ]; then
    echo "ERROR: could not find laplace_cpu_gcc or laplace_cpu_icc in PATH"
    retval=1
else
    if [ $gcc -eq 0 ]; then
        echo "Using laplace_cpu_gcc"
    else
        echo "Using laplace_cpu_icc"
    fi
fi

if [ -z ${LD_LIBRARY_PATH+x} ]; then
    echo "ERROR: the \$LD_LIBRARY_PATH environment variable is unset"
    retval=1
else
    x=$(/sbin/ldconfig -N -v $(sed 's/:/ /g' <<< $LD_LIBRARY_PATH) 2>/dev/null |grep libopenblas)
    if [ "x$x" == "x" ]; then
        echo "ERROR: libopenblas.so not found (usually this means LD_LIBRARY_PATH is incorrect)"
        retval=1
    else
        echo "Using OpenBLAS: $x"
    fi
fi



if [ -z ${GATK_PATH+x} ]; then
    echo "ERROR: the \$GATK_PATH environment variable is unset"
    retval=1
else
    echo "Using GATK_PATH=$GATK_PATH"
    if [ ! -f "$GATK_PATH/gatk.jar" ]; then
        echo "ERROR: gatk.jar not found in GATK_PATH"
        retval=1
    fi
    if [ ! -f "$GATK_PATH/human_g1k_v37_decoy.fasta" ]; then
        echo "ERROR: human_g1k_v37_decoy.fasta not found in GATK_PATH"
        retval=1
    fi
    if [ ! -f "$GATK_PATH/human_g1k_v37_decoy.fasta.fai" ]; then
        echo "ERROR: human_g1k_v37_decoy.fasta.fai not found in GATK_PATH"
        retval=1
    fi
    if [ ! -f "$GATK_PATH/human_g1k_v37_decoy.dict" ]; then
        echo "ERROR: human_g1k_v37_decoy.dict not found in GATK_PATH"
        retval=1
    fi
    if [ ! -f "$GATK_PATH/dbsnp.vcf" ]; then
        echo "ERROR: dbsnp.vcf not found in GATK_PATH"
        retval=1
    fi
    if [ ! -f "$GATK_PATH/dbsnp.vcf.idx" ]; then
        echo "ERROR: dbsnp.vcf.idx not found in GATK_PATH"
        retval=1
    fi
fi



if [ -z ${SHAPEIT_ROOT+x} ]; then
    echo "ERROR: the \$SHAPEIT_ROOT environment variable is not set"
    retval=1
else
    echo "Using SHAPEIT_ROOT=$SHAPEIT_ROOT"
    if [ ! -f "$SHAPEIT_ROOT/bin/shapeit" ]; then
        echo "ERROR: shapeit binary not found in SHAPEIT_ROOT"
        retval=1
    fi
fi



if [ -z ${REFPANEL_ROOT+x} ]; then
    echo "ERROR: the \$REFPANEL_ROOT environment variable is unset"
    retval=1
else
    echo "Using REFPANEL_ROOT=$REFPANEL_ROOT"
    if [ ! -f "$REFPANEL_ROOT/genetic_map_chr1_combined_b37.txt" ]; then
        echo "ERROR: did not find expected files in REFPANEL_ROOT"
        retval=1
    fi
fi


if [ -z ${DRMAA_LIBRARY_PATH+x} ]; then
    echo "WARNING: the \$DRMAA_LIBRARY_PATH environment variable is unset"
    echo "WARNING: the cluster scripts for grid fitting will not work!"
else
    echo "Using DRMAA_LIBRARY_PATH=$DRMAA_LIBRARY_PATH"
    if [ ! -f "$DRMAA_LIBRARY_PATH" ]; then
        echo "ERROR: DRMAA_LIBRARY_PATH specified but file does not exist"
        retval=1
    else
        echo "CLUSTER GRIDFITTING ENABLED"
    fi
fi


##################################################################
# PATH must contain quite a few things..
##################################################################
test_path "samtools"
test_path "Rscript"
test_path "python"
test_path "java"
test_path "demo.sh"
test_path "run_gatk.sh"
test_path "run_shapeit.sh"
test_path "scan_snv.sh"

which Rscript >& /dev/null
if [ $? -ne 0 ]; then
    Rscript -e "library('scansnv')" >& /dev/null
    if [ $? -ne 0 ]; then
        echo "ERROR: could not load SCAN-SNV R package"
        retval=1
    else
        echo "Found R library scansnv"
    fi
fi

exit $retval
