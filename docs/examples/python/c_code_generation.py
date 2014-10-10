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
from casadi import *
from numpy import *
from os import system
import time
import sys


compileme = True  # Flag wether compiling should take place or not

if len(sys.argv)>1 and sys.arg[0]=='nc':
  compileme = False
else:
  print "Info: Use 'python c_code_generation.py nc' to omit compiling"
  
x = SX.sym("x",7,7)
f = det(x)
x = vec(x)
x0 = [random.rand() for xi in x.data()]

fcn = SXFunction([x],[f])
fcn.init()

# adjoint
gf = fcn.grad()

gfcn = SXFunction([x],[gf])
gfcn.init()

srcname = "grad_det.c"
gfcn.generateCode(srcname)

objname_no_opt = "grad_det_no_opt.so"
print "Compiling without optimization: ", objname_no_opt
t1 = time.time()
if compileme:
  system("gcc -fPIC -shared " + srcname + " -o " + objname_no_opt)
t2 = time.time()
print "time = ", (t2-t1)*1e3, " ms"

objname_O3_opt = "grad_det_O3_opt.so"
print "Compiling with O3 optimization: ", objname_O3_opt
t1 = time.time()
if compileme:
  system("gcc -fPIC -shared -O3 " + srcname + " -o " + objname_O3_opt)
t2 = time.time()
print "time = ", (t2-t1)*1e3, " ms"

objname_Os_opt = "grad_det_Os_opt.so"
print "Compiling with Os optimization: ", objname_Os_opt
t1 = time.time()
if compileme:
  system("gcc -fPIC -shared -Os " + srcname + " -o " + objname_Os_opt)
t2 = time.time()
print "time = ", (t2-t1)*1e3, " ms"

# Read function
efcn_no_opt = ExternalFunction("./"+objname_no_opt)
efcn_O3_opt = ExternalFunction("./"+objname_O3_opt)
efcn_Os_opt = ExternalFunction("./"+objname_O3_opt)
efcn_no_opt.init()
efcn_O3_opt.init()
efcn_Os_opt.init()
f_test = [gfcn,efcn_no_opt,efcn_O3_opt,efcn_Os_opt]

# Just-in-time compilation with OpenCL
if False:
  print "Just-in-time compilation with OpenCL"
  t1 = time.time()
  gfcn_opencl = SXFunction([x],[gf])
  gfcn_opencl.setOption("just_in_time_opencl",True)
  gfcn_opencl.init()
  t2 = time.time()
  print "time = ", (t2-t1)*1e3, " ms"
  f_test.append(gfcn_opencl)

for f in f_test:
  f.setInput(x0)
  t1 = time.time()
  nrep = 10000
  for r in range(nrep):
    f.evaluate()
  t2 = time.time()
  print "result = ", f.output().data()
  dt = (t2-t1)/nrep
  print "time = ", dt*1e3, " ms"
  
  num_op = gfcn.getAlgorithmSize()
  print "number of elementary operations: ", num_op
  print "time per elementary operations: ", dt/num_op*1e9, " ns"
  
  
