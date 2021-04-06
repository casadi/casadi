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
import unittest
from types import *
from helpers import *
import pickle
import os
scipy_interpolate = False
try:
  import scipy.interpolate
  scipy.interpolate.RectBivariateSpline
  scipy_interpolate = True
except:
  pass

def vector_widths(li):
  for e in li:
    GlobalOptions.setVectorWidthReal(e)
    yield e
    GlobalOptions.setVectorWidthReal(1)

class Vectorizetests(casadiTestCase):

  def test_call_empty(self):
    


    if 1:
      for ab in AutoBrancher():

        n = 2#ab.branch([1, 2])
        nmap = ab.branch([3, 4])

        x = MX.sym("x",n)

        DM.rng(1)
        A = DM.rand(n,n)

        expand = ab.branch()
        mode = ab.branch(range(2))
        

        for vw in vector_widths([1, 4]):
          f = Function("f",[x],[A @ x**2,sqrt(x+sign(x)+tan(x))])
          if expand:
            f = f.expand()
          fm = f.map(nmap)
          F = [fm, fm.forward(5)][mode]#, fm.forward(5), fm.reverse(5)][mode]#, fm.forward(5).map(7).reverse(2)][mode]

          DM.rng(0)
          args = [DM.rand(F.sparsity_in(i)) for i in range(F.n_in())]
          
          print("name",F.name(),vw)
          if vw==1:
            ref = F.call(args)
          else:
            # -DSLEEF_ENABLE_OMP_SIMD -Dtan=Sleef_tan_u10 -I/usr/local/include -lsleef
            res = self.do_codegen(F,inputs=args,extra_options=["-Wno-maybe-uninitialized"])
            res_unsafe = self.do_codegen(F,inputs=args,extra_options=["-Wno-maybe-uninitialized -fno-math-errno -funsafe-math-optimizations -ffinite-math-only"])
            res_unsafe_vec = self.do_codegen(F,inputs=args,extra_options=["-Wno-maybe-uninitialized -DSLEEF_ENABLE_OMP_SIMD -Dtan=Sleef_tan_u10 -I/usr/local/include -lsleef -Wl,-rpath=/usr/local/lib/"],vectorize_check=4,compiler="gcc-8",static=True,includes="sleef.h",std="gnu18")
        self.check_outputs(res,ref)
        self.check_outputs(res_unsafe_vec,ref,digits=12)
        self.check_outputs(res_unsafe_vec,res_unsafe,digits=12)
    n = 2

    x = MX.sym("x",n)
    y = MX.sym("y",2*n)
    ys = vertsplit(y,n)

    DM.rng(1)
    A = DM.rand(n,n)

    for ab in AutoBrancher():

      nmap = 3#ab.branch([3,5])

      expand = True#ab.branch()
      mode = 5#ab.branch(range(5))
      #mode = 2
      reduce_in = False#ab.branch()
      reduce_out = True#ab.branch()

      for vw in vector_widths([1, 4]):
        print("cache reset")
        f = Function("f",[x, y],[A @ (x-ys[0])**2, cos(ys[1]*x+sign(x))])
        if expand:
          f = f.expand()
        
        if reduce_in:
          fm = f.map(nmap,[False,True],[False,reduce_out])
          Y = MX.sym("y",2*n)
        else:
          fm = f.map(nmap)
          Y = MX.sym("y",2*n,nmap)

        X = MX.sym("x",n,nmap)


        w = fm(X,Y)
        z = fm(w[0]/w[1],X[0,0]*Y)
        g = Function('g',[X,Y],z)


        if mode==0:
          F = fm
        elif mode==1:
          F = fm.forward(4)
        elif mode==2:
          F = g.forward(4)
        elif mode==3:
          F = fm.reverse(4)
        elif mode==4:
          F = g.reverse(4)
        elif mode==5:
          F = g.forward(4).map(4).reverse(4)

        DM.rng(0)
        args = [DM.rand(F.sparsity_in(i)) for i in range(F.n_in())]
        if vw==1:
          ref = F.call(args)
        else:
          res = self.do_codegen(F,inputs=args,extra_options=["-Wno-maybe-uninitialized"])
          res_unsafe = self.do_codegen(F,inputs=args,extra_options=["-Wno-maybe-uninitialized -fno-math-errno -funsafe-math-optimizations -ffinite-math-only"])
          res_unsafe_vec = self.do_codegen(F,inputs=args,extra_options=["-Wno-maybe-uninitialized"],vectorize_check=4,compiler="gcc-8",static=True,std="gnu18")
      self.check_outputs(res,ref,digits=12)
      self.check_outputs(res_unsafe_vec,ref,digits=12)
      self.check_outputs(res_unsafe_vec,res_unsafe,digits=12)

if __name__ == '__main__':
    unittest.main()
