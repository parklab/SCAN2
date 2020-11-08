mkdir -p $PREFIX/bin
cp bin/scan2 $PREFIX/bin
mkdir -p $PREFIX/lib/scansnv
cp scripts/*.{R,sh,py} $PREFIX/lib/scansnv/
cp snakemake/Snakefile $PREFIX/lib/scansnv/
cp snakemake/cluster.yaml $PREFIX/lib/scansnv/
