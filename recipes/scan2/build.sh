mkdir -p $PREFIX/bin
cp bin/scan2 $PREFIX/bin
mkdir -p $PREFIX/lib/scan2
cp scripts/*.{R,sh,py} $PREFIX/lib/scan2/
cp snakemake/Snakefile $PREFIX/lib/scan2/
mkdir -p $PREFIX/lib/scan2/scripts
cp snakemake/scripts/*.R $PREFIX/lib/scan2/scripts
