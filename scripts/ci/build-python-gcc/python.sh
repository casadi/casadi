#!/usr/bin/env bash

# Script to download and build Python 3 from source

set -ex

version="${1:-3.9}"
builddir="/tmp"
prefix="$HOME/.local"

case $version in 
  3.7)
    full_version=3.7.12
    python="Python-$full_version";;
  3.8)
    full_version=3.8.12
    python="Python-$full_version";;
  3.9)
    full_version=3.9.7
    python="Python-$full_version";;
  3.10)
    full_version=3.10.0
    python="Python-${full_version}";;
esac

# Download and extract the Python source code
mkdir -p "$builddir"
cd $builddir
if [ ! -d "$python" ]; then
    wget -O- "https://www.python.org/ftp/python/$full_version/$python.tgz" | tar -xz
fi

cd "$python"
./configure --prefix="$prefix" \
    --enable-ipv6 \
    --enable-shared \
    --with-lto --enable-optimizations \
    'LDFLAGS=-Wl,-rpath,\$$ORIGIN/../lib'

make -j$(($(nproc) * 2))
# make altinstall
make install
