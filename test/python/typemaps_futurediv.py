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
import casadi as c
import numpy
import numpy as n
from numpy import array, matrix
import unittest
from types import *
from helpers import *

scipy_available = True
try:
	from scipy.sparse import *
	from scipy import *
except:
	scipy_available = False
import warnings

class typemaptests(casadiTestCase):

  def setUp(self):
    pass

  def test_floordiv(self):
    self.message("make sure that floor_div raises errors")
    for x in [SX.sym("x"),MX.sym("x"),SX.sym("x")]:
      for y in [2,2.0,DM(3),numpy.array([2.0])]:
        print((x,y))
        self.assertRaises(Exception,lambda : x//y)
        self.assertRaises(Exception,lambda : y//x)

  def test_autoconversion(self):
    self.message("Auto conversion DM")
    x=np.array([2.3])
    s = DM([[1,2],[3,4]])
    n = np.array(s)

    self.checkarray(x[0]*s,s*x[0],"")
    self.checkarray(x[0]*s,n*x[0],"")

    self.checkarray(x[0]/s,1/(s/x[0]),"")
    self.checkarray(x[0]/s,x[0]/n,"")

    self.checkarray(x[0]-s,-(s-x[0]),"")
    self.checkarray(x[0]-s,x[0]-n,"")

    self.checkarray(x[0]+s,s+x[0],"")
    self.checkarray(x[0]+s,x[0]+n,"")

    w=np.array([2.3])[0]
    w+=s
    self.checkarray(w,2.3+n,"")

    w=np.array([2.3])[0]
    w-=s
    self.checkarray(w,2.3-n,"")

    w=np.array([2.3])[0]
    w*=s
    self.checkarray(w,2.3*n,"")

    w=np.array([2.3])[0]
    w/=s
    self.checkarray(w,2.3/n,"")

    x=[2.3]

    self.checkarray(x[0]*s,s*x[0],"")
    self.checkarray(x[0]*s,n*x[0],"")

    self.checkarray(x[0]/s,1/(s/x[0]),"")
    self.checkarray(x[0]/s,x[0]/n,"")

    self.checkarray(x[0]-s,-(s-x[0]),"")
    self.checkarray(x[0]-s,x[0]-n,"")

    self.checkarray(x[0]+s,s+x[0],"")
    self.checkarray(x[0]+s,x[0]+n,"")


    w=2.3
    w+=s
    self.checkarray(w,2.3+n,"")

    w=2.3
    w-=s
    self.checkarray(w,2.3-n,"")

    w=2.3
    w*=s
    self.checkarray(w,2.3*n,"")

    w=2.3
    w/=s
    self.checkarray(w,2.3/n,"")

  def test_autoconversionMX(self):
    self.message("Auto conversion MX")
    s = DM([[1,2],[3,4]])
    x = SX(3)
    y = MX(3)

    def doit(z,s,fun):
      if type(z) in [type(SX()),type(SX())]:
        ztype = [type(SX()),type(SX())]

      if type(z) in [type(MX())]:
        ztype = [type(MX())]

      r = fun(z,s)

      if type(z) is type(SX()) and type(s) is type(SX()):
        self.assertTrue(type(r) is type(SX()))


      self.assertTrue(type(r) in ztype)

      hasNum = True
      if type(s) in [type(SX()),type(MX()),type(SX())]:
        hasNum = False

      if hasNum:
        dummy = [1.3,2.7,9.4,1.0]

        f=Function('f', [z],[r])
        f_in = DM(f.sparsity_in(0),dummy[0:f.nnz_in(0)])
        f_out = f(f_in)

        f_=Function('f', [z],[z])
        f__in = DM(f_.sparsity_in(0),dummy[0:f.nnz_in(0)])
        f__out = f_(f__in)


        self.checkarray(fun(f__out,DM(s)),f_out,"operation")
      else:
        dummy = [1.3,2.7,9.4,1.0]
        dummy2 = [0.3,2.4,1.4,1.7]

        f=Function('f',[z,s],[r])
        f_in = [DM(f.sparsity_in(0),dummy[0:f.nnz_in(0)]), DM(f.sparsity_in(1),dummy2[0:f.nnz_in(1)])]
        f_out = f(*f_in)

        f_=Function('f', [z,s],[z,s])
        f__in= [DM(f_.sparsity_in(0),dummy[0:f.nnz_in(0)]), DM(f_.sparsity_in(1),dummy2[0:f.nnz_in(1)])]
        f__out = f_(*f__in)

        self.checkarray(fun(f__out[0],f__out[1]),f_out,"operation")


    def tests(z,s):
      doit(z,s,lambda z,s: -z)
      doit(z,s,lambda z,s: z+s)
      doit(z,s,lambda z,s: s+z)
      doit(z,s,lambda z,s: s*z)
      doit(z,s,lambda z,s: z*s)
      doit(z,s,lambda z,s: z-s)
      doit(z,s,lambda z,s: s-z)
      doit(z,s,lambda z,s: z/s)
      doit(z,s,lambda z,s: s/z)
      doit(z,s,lambda z,s: z**s)
      doit(z,s,lambda z,s: s**z)
      doit(z,s,lambda z,s: np.fmin(s,z))
      doit(z,s,lambda z,s: np.fmax(s,z))
      doit(z,s,lambda z,s: constpow(s,z))
      doit(z,s,lambda z,s: constpow(z,s))

    nums = [np.array([[1,2],[3,4]]),DM([[1,2],[3,4]]), DM(4), np.array(4),4.0,4]

    ## numeric & SX
    for s in nums:
      for z in [SX.sym("x"), SX.sym("x"), SX.sym("x",2,2)]:
        print("z = %s, s = %s" % (str(z),str(s)))
        print("  z = %s, s = %s" % (type(z),type(s)))
        tests(z,s)

    # numeric & MX
    for s in nums:
      for z in [MX.sym("x",2,2)]:
        print("z = %s, s = %s" % (str(z),str(s)))
        print("  z = %s, s = %s" % (type(z),type(s)))
        tests(z,s)

    # SX & SX
    for s in [SX.sym("x"), SX.sym("x"), SX.sym("x",2,2)]:
      for z in [SX.sym("x"),SX.sym("x"), SX.sym("x",2,2)]:
        print("z = %s, s = %s" % (str(z),str(s)))
        print("  z = %s, s = %s" % (type(z),type(s)))
        tests(z,s)

    ## MX & MX
    for s in [MX.sym("x"),MX.sym("x",2,2)]:
      for z in [MX.sym("x"),MX.sym("x",2,2)]:
        print("z = %s, s = %s" % (str(z),str(s)))
        print("  z = %s, s = %s" % (type(z),type(s)))
        tests(z,s)

    for (s,x,y) in [
                  (np.array([[1,2],[3,4]]),SX.sym("x",2,2),MX.sym("x",2,2))
                  ]:
      for z,ztype in zip([x,y],[[type(SX()),type(SX())],[type(MX())]]):
        print("z = %s, s = %s" % (str(z),str(s)))
        print("  z = %s, s = %s" % (type(z),type(s)))
        doit(z,s,lambda z,s: -z)
        -s
        doit(z,s,lambda z,s: z+s)
        doit(z,s,lambda s,z: s+z)
        doit(z,s,lambda z,s: s*z)
        doit(z,s,lambda s,z: z*s)
        doit(z,s,lambda z,s: z-s)
        doit(z,s,lambda s,z: s-z)
        doit(z,s,lambda z,s: z/s)
        doit(z,s,lambda s,z: s/z)

if __name__ == '__main__':
    unittest.main()
