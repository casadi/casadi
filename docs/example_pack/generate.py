#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             K.U. Leuven. All rights reserved.
#     Copyright (C) 2011-2014 Greg Horn
#
#     CasADi is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 3 of the License, or (at your option) any later version.
#
#     CasADi is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public
#     License along with CasADi; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
#
import zipfile
import os
import casadi
release = casadi.__version__

import zlib
compression = zipfile.ZIP_DEFLATED
zf = zipfile.ZipFile('example_pack.zip', mode='w')


if "+" in release:
  releasedocurl = "http://docs.casadi.org/api/html/annotated.html"
else:
  releasedocurl = "http://docs.casadi.org/v%s/api/html/annotated.html" % release

zf.writestr('README.txt',"""

This is a collection of python examples for CasADi %s.
Consult the documentation online at %s.

To run an example vdp_single_shooting.py, you have several options:

a) fire up a terminal, navigate to this directory, and type "python vdp_single_shooting.py"
b) fire up a terminal, navigate to this directory, and type "ipython" and then "run vdp_single_shooting.py"
c) start spyder, open vdp_single_shooting.py and press "run".

    
""" % (release,releasedocurl))

base = "../examples/python"
for root, dirs, files in os.walk(base): # Walk directory tree
  for f in files:
    if f.endswith(".py"):
       zf.write(os.path.join(root,f),os.path.join("python",f),compress_type=compression)
base = "../examples/matlab"
for root, dirs, files in os.walk(base): # Walk directory tree
  for f in files:
    if f.endswith(".m"):
       zf.write(os.path.join(root,f),os.path.join("matlab",f),compress_type=compression)
       
base = "../documents"
for root, dirs, files in os.walk(base): # Walk directory tree
  for f in files:
    if f.endswith(".pdf"):
       zf.write(os.path.join(root,f),f,compress_type=compression)

#zf.write("../cheatsheet/python.pdf","cheatsheet.pdf",compress_type=compression)
zf.write("../users_guide/users_guide.pdf","user_guide.pdf",compress_type=compression)

zf.close()
