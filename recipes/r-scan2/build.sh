#!/bin/bash

if [[ $target_platform =~ linux.* ]]; then
  export DISABLE_AUTOBREW=1
  mv DESCRIPTION DESCRIPTION.old
  grep -v '^Priority: ' DESCRIPTION.old > DESCRIPTION
  $R CMD INSTALL --build .
else
  mkdir -p $PREFIX/lib/R/library/scan2
  mv * $PREFIX/lib/R/library/scan2
fi
