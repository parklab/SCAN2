./configure --prefix=$PREFIX
make
make install

mkdir -p "${PREFIX}/etc/conda/activate.d"
echo "export DRMAA_LIBRARY_PATH=$PREFIX/lib/libdrmaa.so" > "$PREFIX/etc/conda/activate.d/${PKG_NAME}_activate.sh"
mkdir -p "${PREFIX}/etc/conda/deactivate.d"
echo "unset DRMAA_LIBRARY_PATH" > "$PREFIX/etc/conda/deactivate.d/${PKG_NAME}_deactivate.sh"
