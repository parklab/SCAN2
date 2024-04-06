mkdir -p $PREFIX/bin
cp bin/scan2 $PREFIX/bin
cp bin/digest_calls.R $PREFIX/bin
cp bin/scan2_download_eagle_refpanel.sh $PREFIX/bin
# Creates a file with nothing but a few version strings
echo "version='$PKG_VERSION'" > bin/version.py
echo "buildnum='$PKG_BUILDNUM'" >> bin/version.py
echo "githash='$GIT_FULL_HASH'" >> bin/version.py
mkdir -p $PREFIX/lib/scan2
cp scripts/*.{R,sh,py} $PREFIX/lib/scan2/
cp snakemake/Snakefile $PREFIX/lib/scan2/
cp snakemake/snakefile.* $PREFIX/lib/scan2/
mkdir -p $PREFIX/lib/scan2/resources
cp -r resources/analysis_regions $PREFIX/lib/scan2/resources/
cp -r resources/binned_counts $PREFIX/lib/scan2/resources/
mkdir -p $PREFIX/lib/scan2/scripts
cp snakemake/scripts/*.R $PREFIX/lib/scan2/scripts
