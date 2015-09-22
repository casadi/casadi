#!/bin/bash

#if [ ! -f extra/doxy2swig.py ]; then
#  cd extra && wget http://www.aero.iitb.ac.in/~prabhu/software/code/python/doxy2swig.py && cd ..
#fi

cd extra && python doxy2swigX.py ../XML_internal/index.xml ../../../swig/doc.i ../../../swig/internal.i ../../../swig/deprecated.i && cd ..
cd extra && python doxy2swigX.py --merge ../XML_internal/index.xml ../../../swig/doc_merged.i ../../../swig/internal.i ../../../swig/deprecated.i && cd ..
