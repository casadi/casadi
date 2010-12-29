#!/bin/bash

if [ ! -f doxy2swig.py ]; then
  wget http://www.aero.iitb.ac.in/~prabhu/software/code/python/doxy2swig.py
fi

python doxy2swig.py XML/index.xml doc.i
