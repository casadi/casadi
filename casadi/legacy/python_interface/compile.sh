#!/bin/bash
rm -rf build
rm -rf casadi_python.cpp
rm -rf casadi_python.so

CFLAGS="-I$CASADI"  \
LDFLAGS="-L$CASADI/build/lib"     \
    python setup.py build_ext -i

#python setup.py build_ext --inplace