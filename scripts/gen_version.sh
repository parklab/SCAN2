#!/bin/bash

if [ $# -ne 1 ]; then
    echo "usage: $0 path_to_SCAN2_gitrepo"
    exit 1
fi

dir=$1

metafile=$dir/recipes/scan2/meta.yaml

echo "version='$(cat $metafile | grep -v GIT_FULL_HASH | shyaml get-value package.version)'"
echo "buildnum='$(cat $metafile | grep -v GIT_FULL_HASH | shyaml get-value build.number)'"
echo "githash='$(GIT_DIR=`realpath $dir/.git` git rev-parse HEAD)'"
