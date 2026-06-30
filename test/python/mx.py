#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             KU Leuven. All rights reserved.
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
import casadi as ca
import numpy as np
from numpy import inf, pi
import casadi as c
import numpy
from numpy import random, array, linalg, matrix, zeros, ones, ndarray, eye
import unittest
import warnings
from types import *
from helpers import *
from copy import deepcopy

import sys

if sys.version_info >= (3, 0):
  TupleType = tuple


warnings.filterwarnings("ignore",category=DeprecationWarning)

scipy_available = True
try:
	from scipy.sparse import csr_matrix
except:
	scipy_available = False

def checkarray(self,zr,zt,name):
    if len(zr.shape)==1 and (zt.shape[0]==1 or zt.shape[1]==1) and zr.shape[0]==zt.shape[1]*zt.shape[0]:
      zr=ca.reshape(zr,(zt.shape));
    self.assertEqual(zt.shape[0],zr.shape[0],"%s dimension error. Got %s, expected %s" % (name,str(zt.shape),str(zr.shape)))
    self.assertEqual(len(zt.shape),len(zr.shape),"%s dimension error. Got %s, expected %s" % (name,str(zt.shape),str(zr.shape)))
    self.assertEqual(zt.shape[1],zr.shape[1],"%s dimension error. Got %s, expected %s" % (name,str(zt.shape),str(zr.shape)))
    for i in range(zr.shape[0]):
      for j in range(zr.shape[1]):
        self.assertAlmostEqual(zt[i,j],zr[i,j],10,"%s evaluation error. %s <-> %s" % (name, str(zt),str(zr)))

def checkMXoperations(self,ztf,zrf,name):
    x = ca.MX.sym("x",1,3)
    z=ca.vertcat(*[x*(i+1) for i in range(8)])
    f = ca.Function("f", [x],[ztf(z)])
    L=[1,2,3]
    f_out = f(L)
    zt = f_out.full()
    zr = array([[L[0]*(i+1),L[1]*(i+1),L[2]*(i+1)] for i in range(8)])
    checkarray(self,zrf(zr),zt,name)
    return (zt,zrf(zr))

def checkMXoperations2(self,ztf,zrf,name):
    x = ca.MX.sym("x",3,1)
    z = ca.horzcat(*[x*i for i in range(8)])
    f = ca.Function("f", [x],[ztf(z)])
    L=[1,2,3]
    f_out = f(L)
    zt = f_out.full()
    zr = array([[L[0]*i,L[1]*i,L[2]*i] for i in range(8)]).T
    checkarray(self,zrf(zr),zt,name)
    return zt

def checkMXoperations3(self,ztf,zrf,name):
    x = ca.MX.sym("x",3,1)
    p = ca.horzcat(*[x[0,0],x[1,0],x[2,0]])
    z = ca.vertcat(*[p*i for i in range(8)])
    f = ca.Function("f", [x],[ztf(z)])
    L=[1,2,3]
    f_out = f(L)
    zt = f_out.full()
    zr = array([[L[0]*i,L[1]*i,L[2]*i] for i in range(8)])
    checkarray(self,zrf(zr),zt,name)
    return (zt,zrf(zr))

class MXtests(casadiTestCase):

  def setUp(self):
    self.pool=FunctionPool()
    self.pool.append(lambda x: ca.sqrt(x[0]),ca.sqrt,"sqrt")
    self.pool.append(lambda x: ca.sin(x[0]),ca.sin,"sin")
    self.pool.append(lambda x: ca.cos(x[0]),ca.cos,"cos")
    self.pool.append(lambda x: ca.tan(x[0]),ca.tan,"tan")
    self.pool.append(lambda x: ca.arctan(x[0]),ca.arctan,"arctan")
    self.pool.append(lambda x: ca.arcsin(x[0]),ca.arcsin,"arcsin")
    self.pool.append(lambda x: ca.arccos(x[0]),ca.arccos,"arccos")
    self.pool.append(lambda x: ca.exp(x[0]),ca.exp,"exp")
    self.pool.append(lambda x: ca.log(x[0]),ca.log,"log",flags={'nozero'})
    self.pool.append(lambda x: x[0]**0,lambda x : x**0,"x^0",flags={'nozero'})
    self.pool.append(lambda x: x[0]**1,lambda x : x**1,"^1")
    self.pool.append(lambda x: x[0]**(-2),lambda x : x**(-2),"^-2",flags={'nozero'})
    self.pool.append(lambda x: x[0]**(0.3),lambda x : x**(0.3),"^0.3")
    self.pool.append(lambda x: ca.floor(x[0]),ca.floor,"floor")
    self.pool.append(lambda x: ca.ceil(x[0]),ca.ceil,"ceil")
    self.Jpool=FunctionPool()
    self.Jpool.append(lambda x: ca.sqrt(x[0]),lambda x:ca.diag(1/(2.0*ca.sqrt(x))),"sqrt")
    self.Jpool.append(lambda x: ca.sin(x[0]),lambda x:ca.diag(ca.cos(x)),"sin")
    self.Jpool.append(lambda x: ca.cos(x[0]),lambda x:ca.diag(-ca.sin(x)),"cos")
    self.Jpool.append(lambda x: ca.tan(x[0]),lambda x:ca.diag(1.0/ca.cos(x)**2),"tan")
    self.Jpool.append(lambda x: ca.arctan(x[0]),lambda x:ca.diag( 1.0/(x**2+1)),"arctan")
    self.Jpool.append(lambda x: ca.arcsin(x[0]),lambda x:ca.diag( 1.0/ca.sqrt(1-x**2)),"arcsin")
    self.Jpool.append(lambda x: ca.arccos(x[0]),lambda x: ca.diag(-1.0/ca.sqrt(1-x**2)),"arccos")
    self.Jpool.append(lambda x: ca.exp(x[0]),lambda x: ca.diag(ca.exp(x)),"exp")
    self.Jpool.append(lambda x: ca.log(x[0]),lambda x: ca.diag(1.0/x),"log")
    self.Jpool.append(lambda x: x[0]**0,lambda x :ca.diag(zeros(x.shape)),"x^0")
    self.Jpool.append(lambda x: x[0]**1,lambda x : ca.diag(ones(x.shape)),"^1")
    self.Jpool.append(lambda x: x[0]**(-2),lambda x : ca.diag(-2.0/x**3),"^-2")
    self.Jpool.append(lambda x: x[0]**(0.3),lambda x :ca.diag( 0.3/x**0.7),"^0.3")
    self.matrixpool=FunctionPool()
    #self.matrixpool.append(lambda x: norm_2(x[0]),linalg.norm,"norm_2")
    #self.matrixpool.append(lambda x: norm_1(x[0]),lambda x: sum(sum(abs(x))),"norm_1")
    #self.matrixpool.append(lambda x: norm_inf(x[0]),lambda x: abs(matrix(x)).max(),"norm_inf")
    self.matrixbinarypool=FunctionPool()
    self.matrixbinarypool.append(lambda a: a[0]+a[1],lambda a: a[0]+a[1],"Matrix+Matrix")
    self.matrixbinarypool.append(lambda a: a[0]-a[1],lambda a: a[0]-a[1],"Matrix-Matrix")
    self.matrixbinarypool.append(lambda a: a[0]*a[1],lambda a: a[0]*a[1],"Matrix*Matrix")
    self.matrixbinarypool.append(lambda a: ca.fmax(a[0],a[1]),lambda a: ca.fmax(a[0],a[1]),"fmin")

    self.matrixbinarypool.append(lambda a: ca.fmin(a[0],a[1]),lambda a: ca.fmin(a[0],a[1]),"fmax")
    self.matrixbinarypool.append(lambda a: a[0] @ a[1].T,lambda a: numpy.dot(a[0],a[1].T),"mtimes(Matrix,Matrix.T)")
    self.matrixbinarypool.append(lambda a: ca.arctan2(a[0],a[1]),lambda a: ca.arctan2(a[0],a[1]),"arctan2")
    #self.matrixbinarypool.append(lambda a: inner_mul(a[0],trans(a[1])),lambda a: c.dot(a[0].T,a[1]),name="inner_mul(Matrix,Matrix)")
    self.matrixbinarypool.append(lambda a: a[0] @ a[1].T,lambda a: numpy.dot(a[0],a[1].T),"mtimes(Matrix,Matrix.T)")

  def test_MX1(self):
    self.message("MX constructor")
    x = ca.MX.sym("x",2,3)
    self.assertEqual(x.size1(),2,"MX fails to indicate its size1")
    self.assertEqual(x.size2(),3,"MX fails to indicate its size2")

  def test_MXvertcat(self):
    self.message("MX vertcat")
    x = ca.MX.sym("x",1,3)
    y = ca.MX.sym("y",1,3)
    z=ca.vertcat(*(x,y))
    self.assertEqual(z.size1(),2,"MX fails to indicate its size1")
    self.assertEqual(z.size2(),3,"MX fails to indicate its size2")

  def test_MX_fun1(self):
    self.message("MXFunction single input, single output")
    # check if x->2*x
    # evaluates correctly for x=3
    x = ca.MX.sym("x")
    y = 2*x
    f = ca.Function("f", [x],[y])
    self.assertEqual(f.n_in(),1,"Function fails to indicate correct number of inputs")
    self.assertEqual(f.n_out(),1,"Function fails to indicate correct number of outputs")
    f_out = f(3)
    yt = tuple(f_out.nonzeros())
    self.assertEqual(type(yt),TupleType,"Output of Function is expected to be tuple of floats")
    self.assertEqual(len(yt),1,"Output of Function was tuple of floats, as expected, but length is incorrect.")
    y=yt[0]
    self.assertEqual(type(y),float,"Output of Function is expected to be tuple of floats")
    self.assertAlmostEqual(y, 2*3,10)

  def test_MXfunction2(self):
    self.message("Function multi input, multi output")
      # check if [x,y]->[y+x,y*x]
    # evaluates correctly for x=3,y=7
    x = ca.MX.sym("x")
    y = ca.MX.sym("y")
    f = ca.Function("f", [x,y],[x+y,y*x])
    self.assertEqual(f.n_in(),2,"Function fails to indicate correct number of inputs")
    self.assertEqual(f.n_out(),2,"Function fails to indicate correct number of outputs")
    f_out = f(3, 7)
    zt1 = tuple(f_out[0].nonzeros())
    zt2 = tuple(f_out[1].nonzeros())
    self.assertEqual(type(zt1),TupleType,"Output of Function is expected to be tuple of floats")
    self.assertEqual(type(zt2),TupleType,"Output of Function is expected to be tuple of floats")
    self.assertEqual(len(zt1),1,"Output of Function was tuple of floats, as expected, but length is incorrect.")
    self.assertEqual(len(zt2),1,"Output of Function was tuple of floats, as expected, but length is incorrect.")
    z1=zt1[0]
    z2=zt2[0]
    self.assertEqual(type(z1),float,"Output of Function is expected to be tuple of floats")
    self.assertEqual(type(z2),float,"Output of Function is expected to be tuple of floats")
    self.assertAlmostEqual(z2, 21,10)
    self.assertAlmostEqual(z1, 10,10)



  def test_MXfunction3(self):
    self.message("Function single input, multi output (1)")
    # check if [x,y]->[y+x,y*x]
    # evaluates correctly for x=3,y=7
    # now with single input, multi output
    xy = ca.MX.sym("xy",2)
    f = ca.Function("f", [xy],[xy[0]+xy[1],xy[0]*xy[1]])
    self.assertEqual(f.n_in(),1,"Function fails to indicate correct number of inputs")
    self.assertEqual(f.n_out(),2,"Function fails to indicate correct number of outputs")
    f_out = f([3,7])
    zt1 = tuple(f_out[0].nonzeros())
    zt2 = tuple(f_out[1].nonzeros())
    self.assertEqual(type(zt1),TupleType,"Output of Function is expected to be tuple of floats")
    self.assertEqual(type(zt2),TupleType,"Output of Function is expected to be tuple of floats")
    self.assertEqual(len(zt1),1,"Output of Function was tuple of floats, as expected, but length is incorrect.")
    self.assertEqual(len(zt2),1,"Output of Function was tuple of floats, as expected, but length is incorrect.")
    z1=zt1[0]
    z2=zt2[0]
    self.assertEqual(type(z1),float,"Output of Function is expected to be tuple of floats")
    self.assertEqual(type(z2),float,"Output of Function is expected to be tuple of floats")
    self.assertAlmostEqual(z2, 21,10)
    self.assertAlmostEqual(z1, 10,10)

  def test_MXfunction3b(self):
    self.message("Function single input, multi output (2)")
    # check if [x,y]->[y+x,y*x]
    # evaluates correctly for x=3,y=7
    # now with single input, multi output
    xy = ca.MX.sym("xy",1,2)
    f = ca.Function("f", [xy],[xy[0,0]+xy[0,1],xy[0,0]*xy[0,1]])
    self.assertEqual(f.n_in(),1,"Function fails to indicate correct number of inputs")
    self.assertEqual(f.n_out(),2,"Function fails to indicate correct number of outputs")
    f_out = f([3,7])
    zt1 = f_out[0].full()
    zt2 = f_out[1].full()

    self.assertEqual(type(zt1),ndarray,"Output of Function is expected to be numpy.ndarray")
    self.assertEqual(zt1.shape[0],1,"Output of Function is of wrong shape.")
    self.assertEqual(zt1.shape[1],1,"Output of Function is of wrong shape.")

    self.assertEqual(type(zt2),ndarray,"Output of Function is expected to be numpy.ndarray")
    self.assertEqual(zt2.shape[0],1,"Output of Function is of wrong shape.")
    self.assertEqual(zt2.shape[1],1,"Output of Function is of wrong shape.")

    z1=zt1[0,0]
    z2=zt2[0,0]
    self.assertEqual(type(z1),numpy.float64,"Output of Function is expected to be numpy.ndarray of floats")
    self.assertEqual(type(z2),numpy.float64,"Output of Function is expected to be numpy.ndarray of floats")
    self.assertAlmostEqual(z2, 21,10)
    self.assertAlmostEqual(z1, 10,10)

  def test_MXfunction4(self):
    self.message("Function single input, single output , using vertcat")
    # check if [x,y]->[y+x,y*x]
    # evaluates correctly for x=3,y=7
    # now with single input, single output
    xy = ca.MX.sym("xy",2)
    z=ca.vertcat(*[xy[0]+xy[1],xy[0]*xy[1]])
    f = ca.Function("f", [xy],[z])
    self.assertEqual(f.n_in(),1,"Function fails to indicate correct number of inputs")
    self.assertEqual(f.n_out(),1,"Function fails to indicate correct number of outputs")
    f_out = f([3,7])
    zt=f_out.full()
    self.assertEqual(type(zt),ndarray,"Output of Function is expected to be numpy.ndarray")
    self.assertEqual(zt.shape[0],2,"Output of Function is of wrong shape.")
    self.assertEqual(zt.shape[1],1,"Output of Function is of wrong shape.")
    z1=zt[0,0]
    z2=zt[1,0]
    self.assertEqual(type(z1),numpy.float64,"Output of Function is expected to be numpy.ndarray of floats")
    self.assertEqual(type(z2),numpy.float64,"Output of Function is expected to be numpy.ndarray of floats")
    self.assertAlmostEqual(z2, 21,10)
    self.assertAlmostEqual(z1, 10,10)

  def test_MXfunction5(self):
    self.message("Function single input, single output , using horzcat")
    # check if [x,y]->[y+x,y*x]
    # evaluates correctly for x=3,y=7
    # now with single input, single output
    xy = ca.MX.sym("xy",2)
    z=ca.horzcat(*[xy[0]+xy[1],xy[0]*xy[1]])
    f = ca.Function("f", [xy],[z])
    self.assertEqual(f.n_in(),1,"Function fails to indicate correct number of inputs")
    self.assertEqual(f.n_out(),1,"Function fails to indicate correct number of outputs")
    f_out = f([3,7])
    zt = f_out.full()
    self.assertEqual(type(zt),ndarray,"Output of Function is expected to be numpy.ndarray")
    self.assertEqual(zt.shape[0],1,"Output of Function is of wrong shape.")
    self.assertEqual(zt.shape[1],2,"Output of Function is of wrong shape.")
    z1=zt[0,0]
    z2=zt[0,1]
    self.assertEqual(type(z1),numpy.float64,"Output of Function is expected to be numpy.ndarray of floats")
    self.assertEqual(type(z2),numpy.float64,"Output of Function is expected to be numpy.ndarray of floats")
    self.assertAlmostEqual(z2, 21,10)
    self.assertAlmostEqual(z1, 10,10)

  def test_which_depends_empty(self):
    for X in [ca.SX,ca.MX]:
      x=X.sym("x")

      for tr in [True,False]:
        for i in [0,1,2]:
          self.assertEqual(ca.which_depends(x,X(0,1),i,tr),[False]*(1 if tr else 0))

          self.assertEqual(ca.which_depends(X(0,1),x,i,tr),[False]*(0 if tr else 1))
          self.assertTrue(len(ca.which_depends(X(0,1),X(0,1),i,tr))==0)

  def test_issue83(self):
    x=ca.MX.sym("x")
    y=ca.MX.sym("y")

    z = x + y

    f = ca.Function("f", [x,y],[z])

    fc = f(ca.MX(3),y)

    g = ca.Function("g", [y],[fc])
    g_in = [7]
    g_out = g(g_in)

    self.assertAlmostEqual(g_out[0],10,10,"issue #83")

    fc = f(x,ca.MX(7))

    g = ca.Function("g", [x],[fc])
    g_in = [3]
    g_out = g(g_in)

    self.assertAlmostEqual(g_out[0],10,10,"issue #83")

  def test_identitySX(self):
    self.message("identity SXFunction")
    x = ca.SX.sym("x")
    f = ca.Function("f", [x],[x])
    f_in = [3]
    f_out = f(f_in)
    self.assertAlmostEqual(f_out[0,0], 3,10)

  def test_identityMX(self):
    self.message("identity Function")
    x = ca.MX.sym("x")
    f = ca.Function("f", [x],[x])
    f_in = [3]
    f_out = f(f_in)
    self.assertAlmostEqual(f_out[0,0], 3,10)

  def test_MXorder(self):
    self.message("Function order of non-zero elements")
    x = ca.MX.sym("x",2,3)
    f = ca.Function("f", [x],[x+x])

    self.assertEqual(f.n_in(),1,"Function fails to indicate correct number of inputs")
    self.assertEqual(f.n_out(),1,"Function fails to indicate correct number of outputs")
    L=[1,2,3,4,5,6]
    f_in = ca.DM(f.sparsity_in(0),L)
    f_out = f(f_in)
    zt = f_out.full()
    self.assertEqual(zt.shape[0],2,"Output of Function is of wrong shape.")
    self.assertEqual(zt.shape[1],3,"Output of Function is of wrong shape.")

    Lr=numpy.reshape(L,(2,3),'F')
    for i in range(2):
      for j in range(3):
        self.assertAlmostEqual(Lr[i,j]*2, zt[i,j],10)

  def test_trans(self):
    self.message("trans")
    a = ca.MX(0,1)
    b = a.T
    self.assertEqual(b.size1(),1)
    self.assertEqual(b.size2(),0)

  def test_MXtrans(self):
    self.message("trans(MX)")
    x = ca.MX.sym("x",2,3)
    z=x.T
    self.assertEqual(z.size1(),3,"Vec returns MX of wrong dimension")
    self.assertEqual(z.size2(),2,"Vec returns MX of wrong dimension")
    f = ca.Function("f", [x],[z])
    self.assertEqual(f.n_in(),1,"Function fails to indicate correct number of inputs")
    self.assertEqual(f.n_out(),1,"Function fails to indicate correct number of outputs")
    L=[1,2,3,4,5,6]
    f_in = ca.DM(f.sparsity_in(0),L)
    f_out = f(f_in)
    zt = f_out.full()

    ztr=numpy.reshape(zt,(3,2))
    Lr=numpy.reshape(L,(2,3),'F')
    for i in range(2):
      for j in range(3):
        self.assertAlmostEqual(Lr[i,j], ztr[j,i],10)

  def test_MXvec(self):

    u = ca.DM([[10*j+i for i in range(3)] for j in range(4) ])

    U = ca.MX.sym("u",u.shape)

    f = ca.Function("f", [U],[ca.vec(U)])
    f_out = f(u)

    self.checkarray(ca.vec(u),f_out,"vec")

  def test_MXreshape(self):
    self.message("reshape(MX)")
    x = ca.MX.sym("x",2,3)
    z=c.reshape(x,(1,6))
    self.assertEqual(z.size1(),1,"Vec returns MX of wrong dimension")
    self.assertEqual(z.size2(),6,"Vec returns MX of wrong dimension")
    f = ca.Function("f", [x],[z])
    self.assertEqual(f.n_in(),1,"Function fails to indicate correct number of inputs")
    self.assertEqual(f.n_out(),1,"Function fails to indicate correct number of outputs")
    L=[1,2,3,4,5,6]
    f_in = ca.DM(f.sparsity_in(0),L)
    f_out = f(f_in)
    zt = f_out.full()
    for i in range(len(L)):
      self.assertAlmostEqual(L[i], zt[0,i],10)

  def test_sparsity_cast(self):
    sp_source = ca.sparsify(ca.DM([[1, 0, 1],[0, 0, 1],[1, 0,0]])).sparsity()
    x   = ca.MX.sym("x",sp_source)
    xsx = ca.SX.sym("x",sp_source)



    sp = ca.sparsify(ca.DM([[1, 0, 1],[1, 0,1]])).sparsity()

    with self.assertInException("mismatch"):
      ca.reshape(x,sp)

    with self.assertInException("Mismatching"):
      ca.sparsity_cast(x,ca.horzcat(sp,sp))

    y = ca.sparsity_cast(x,sp)
    ysx = ca.sparsity_cast(xsx,sp)

    fsx = ca.Function("fsx",[xsx],[ysx])
    f = ca.Function("f",[x],[y])

    inp = ca.DM([[1, 0, 2],[0, 0, 3],[4, 0, 0]])
    self.checkfunction(f,fsx,inputs=[inp])
    self.check_codegen(f,inputs=[inp])
    y = fsx(inp)
    self.checkarray(y,ca.DM([[1, 0, 2],[4, 0, 3]]))
    
    x = ca.MX.sym("x",2)
    sp = ca.sparsify(ca.blockcat([[1,0],[0,1]])).sparsity()
    y = ca.sparsity_cast(ca.MX.sym("y",2),sp)

    ca.vec(y)


  def test_MXcompose(self):
    self.message("compositions of vec, trans, reshape with vertcat")
    checkMXoperations(self,lambda x: x,lambda x: x,'vertcat')
    checkMXoperations(self,lambda x: x.T,lambda x: x.T,'trans(vertcat)')
    checkMXoperations(self,lambda x: x.T.T,lambda x: x,'trans(trans(vertcat))')
    checkMXoperations(self,lambda x: ca.vec(x.T),lambda x: numpy.reshape(x,(numpy.prod(x.shape),1)),'vec(trans(vertcat))')
    checkMXoperations(self,lambda x: ca.vec(x).T,lambda x: numpy.reshape(x.T,(numpy.prod(x.shape),1)).T,'vec(trans(vertcat))')
    checkMXoperations(self,lambda x: c.reshape(x.T,(6,4)).T,lambda x: numpy.reshape(x,(4,6)),'reshape(vertcat)')
    checkMXoperations(self,lambda x: c.reshape(x,(6,4)).T,lambda x: numpy.reshape(x.T,(4,6)),'reshape(trans(vertcat))')
    checkMXoperations(self,lambda x: c.reshape(x.T,(6,4)),lambda x: numpy.reshape(x,(4,6)).T,'trans(reshape(vertcat))')

  def test_MXcompose2(self):
    self.message("compositions of vec, trans, reshape with horzcat")
    checkMXoperations2(self,lambda x: x,lambda x: x,'horzcat')
    checkMXoperations2(self,lambda x: x.T,lambda x: x.T,'trans(horzcat)')
    checkMXoperations2(self,lambda x: x.T.T,lambda x: x,'trans(trans(horzcat))')
    checkMXoperations2(self,lambda x: ca.vec(x.T),lambda x: numpy.reshape(x,(numpy.prod(x.shape),1)),'vec(trans(horzcat))')
    checkMXoperations2(self,lambda x: ca.vec(x).T,lambda x: numpy.reshape(x.T,(numpy.prod(x.shape),1)).T,'vec(trans(horzcat))')
    checkMXoperations2(self,lambda x: c.reshape(x.T,(6,4)).T,lambda x: numpy.reshape(x,(4,6)),'reshape(horzcat)')
    checkMXoperations2(self,lambda x: c.reshape(x,(6,4)).T,lambda x: numpy.reshape(x.T,(4,6)),'reshape(trans(horzcat))')
    checkMXoperations2(self,lambda x: c.reshape(x.T,(6,4)),lambda x: numpy.reshape(x,(4,6)).T,'trans(reshape(horzcat))')

  def test_MXcompose3(self):
    self.message("compositions of vec, trans, reshape with vertcat")
    checkMXoperations3(self,lambda x: x,lambda x: x,'snippet')
    checkMXoperations3(self,lambda x: x.T,lambda x: x.T,'trans(snippet)')
    checkMXoperations3(self,lambda x: x.T.T,lambda x: x,'trans(trans(snippet))')
    checkMXoperations3(self,lambda x: ca.vec(x.T),lambda x: numpy.reshape(x,(numpy.prod(x.shape),1)),'vec(trans(snippet))')
    checkMXoperations3(self,lambda x: ca.vec(x).T,lambda x: numpy.reshape(x.T,(numpy.prod(x.shape),1)).T,'vec(trans(snippet))')
    checkMXoperations3(self,lambda x: c.reshape(x.T,(6,4)).T,lambda x: numpy.reshape(x,(4,6)),'reshape(snippet)')
    checkMXoperations3(self,lambda x: c.reshape(x,(6,4)).T,lambda x: numpy.reshape(x.T,(4,6)),'reshape(trans(snippet))')
    checkMXoperations3(self,lambda x: c.reshape(x.T,(6,4)),lambda x: numpy.reshape(x,(4,6)).T,'trans(reshape(snippet))')

  def test_MXcompose4(self):
    self.message("compositions of horzcat + vertcat")
    checkMXoperations(self,lambda x: ca.vertcat(*[x]),lambda x: x,'vertcat(*vertcat)')
    checkMXoperations(self,lambda x: ca.vertcat(*[x,x*2]),lambda x: numpy.vstack((x,x*2)),'vertcat(*vertcat,vertcat)')
    checkMXoperations(self,lambda x: ca.horzcat(*[x]),lambda x: x,'horzcat(*vertcat)')
    checkMXoperations(self,lambda x: ca.horzcat(*[x,x*2]),lambda x: numpy.hstack((x,x*2)),'horzcat(*vertcat,vertcat)')

    checkMXoperations2(self,lambda x: ca.vertcat(*[x]),lambda x: x,'vertcat(*horzcat)')
    checkMXoperations2(self,lambda x: ca.vertcat(*[x,x*2]),lambda x: numpy.vstack((x,x*2)),'vertcat(*horzcat,horzcat)')
    checkMXoperations2(self,lambda x: ca.horzcat(*[x]),lambda x: x,'horzcat(*horzcat)')
    checkMXoperations2(self,lambda x: ca.horzcat(*[x,x*2]),lambda x: numpy.hstack((x,x*2)),'horzcat(*horzcat,horzcat)')

    checkMXoperations3(self,lambda x: ca.vertcat(*[x]),lambda x: x,'vertcat(*snippet)')
    checkMXoperations3(self,lambda x: ca.vertcat(*[x,x*2]),lambda x: numpy.vstack((x,x*2)),'vertcat(*snippet,snippet)')
    checkMXoperations3(self,lambda x: ca.horzcat(*[x]),lambda x: x,'horzcat(*snippet)')
    checkMXoperations3(self,lambda x: ca.horzcat(*[x,x*2]),lambda x: numpy.hstack((x,x*2)),'horzcat(*snippet,snippet)')

  @known_bug()  # Test refactoring, cf. #1436
  def test_MXslicingnew(self):
    self.message("MX slicing new")

    self.message(":dense")
    x = ca.MX.sym("x",3,2)
    x0=array([[0.738,0.2],[ 0.1,0.39 ],[0.99,0.999999]])
    self.numpyEvaluationCheck(lambda x:x[0][1,0],lambda x: x[1,0],[x],x0,'x[1,0]')
    self.numpyEvaluationCheck(lambda x:x[0][0,1],lambda x: x[0,1],[x],x0,'x[0,1]')
    self.numpyEvaluationCheck(lambda x:x[0][0,-1],lambda x: x[0,-1],[x],x0,'x[0,-1]')
    self.numpyEvaluationCheck(lambda x: x[0][:,0], lambda x: matrix(x)[:,0],[x],x0,name="x[:,0]")
    self.numpyEvaluationCheck(lambda x: x[0][:,1], lambda x: matrix(x)[:,1],[x],x0,name="x[:,1]")
    self.numpyEvaluationCheck(lambda x: x[0][1,:], lambda x: matrix(x)[1,:],[x],x0,name="x[1,:]")
    self.numpyEvaluationCheck(lambda x: x[0][0,:], lambda x: matrix(x)[0,:],[x],x0,name="x[0,:]")
    self.numpyEvaluationCheck(lambda x: x[0][-1,:], lambda x: matrix(x)[-1,:],[x],x0,name="x[-1,:]")
    self.numpyEvaluationCheck(lambda x: x[0][:,-2], lambda x: matrix(x)[:,-2],[x],x0,name="x[:,-2]")
    self.numpyEvaluationCheck(lambda x: x[0][0:-2,0:-1], lambda x: matrix(x)[0:-2,0:-1],[x],x0,name="x[0:-2,0:-1]")
    self.numpyEvaluationCheck(lambda x: x[0][0:2,0:2], lambda x: matrix(x)[0:2,0:2],[x],x0,name="x[0:2,0:2]")
    self.numpyEvaluationCheck(lambda x: x[0][[0,1],0:2], lambda x: matrix(x)[[0,1],0:2],[x],x0,name="x[[0,1],0:2]")
    self.numpyEvaluationCheck(lambda x: x[0].nz[[0,2,3]], lambda x: matrix([x[0,0],x[2,0],x[0,1]]).T,[x],x0,name="x[[0,2,3]]")

    self.numpyEvaluationCheck(lambda x: x[0].nz[1], lambda x: matrix(x.T.ravel()[1]).T,[x],x0,name="x[1] on dense matrix")
    self.numpyEvaluationCheck(lambda x: x[0].nz[-1], lambda x: matrix(x.ravel()[-1]).T,[x],x0,name="x[-1] on dense matrix")
    self.numpyEvaluationCheck(lambda x: x[0][[0,1],0:1],lambda x: x[[0,1],0:1],[x],x0,name='x[:,0:1]')
    self.numpyEvaluationCheck(lambda x: x[0][0:1,[0,1]],lambda x: x[0:1,[0,1]],[x],x0,name='x[0:1,:]')

    self.message(":sparse")

    sp=ca.Sparsity(4,3,[0,2,2,3],[1,2,1])
    x=ca.MX.sym("X",sp)
    sx0=[0.738,0.39,0.99]
    x0=ca.DM(ca.Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.39,0.99]).full()
    self.numpyEvaluationCheck(lambda x: x[0][0,0], lambda x: matrix(x)[0,0],[x],x0,name="x[0,0]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][1,0], lambda x: matrix(x)[1,0],[x],x0,name="x[1,0]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][0,1], lambda x: matrix(x)[0,1],[x],x0,name="x[1,0]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][0,-1], lambda x: matrix(x)[0,-1],[x],x0,name="x[0,-1]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][:,0], lambda x: matrix(x)[:,0],[x],x0,name="x[:,0]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][:,1], lambda x: matrix(x)[:,1],[x],x0,name="x[:,1]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][1,:], lambda x: matrix(x)[1,:],[x],x0,name="x[1,:]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][0,:], lambda x: matrix(x)[0,:],[x],x0,name="x[0,:]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][-1,:], lambda x: matrix(x)[-1,:],[x],x0,name="x[-1,:]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][:,-2], lambda x: matrix(x)[:,-2],[x],x0,name="x[:,-2]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][0:-2,0:-1], lambda x: matrix(x)[0:-2,0:-1],[x],x0,name="x[0:-2,0:-1]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][0:2,0:2], lambda x: matrix(x)[0:2,0:2],[x],x0,name="x[0:2,0:2]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0][[0,1],0:2], lambda x: matrix(x)[[0,1],0:2],[x],x0,name="x[[0,1],0:2]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0].nz[[2,1]], lambda x: matrix([x[1,2],x[2,0]]).T,[x],x0,name="x[[2,1]]")
    self.numpyEvaluationCheck(lambda x: x[0].nz[0:2], lambda x: matrix(sx0[0:2]).T,[x],x0,name="x[0:2] on dense matrix")
    self.numpyEvaluationCheck(lambda x: x[0].nz[1], lambda x: matrix(sx0[1]).T,[x],x0,name="x[1]",setx0=[sx0])
    self.numpyEvaluationCheck(lambda x: x[0].nz[-1], lambda x: matrix(sx0[-1]).T,[x],x0,name="x[-1]",setx0=[sx0])

  def test_mx_in(self):
    self.message("mx_out/mx_in")
    x=ca.MX.sym("x",2,3)
    f = ca.Function("f", [x],[3*x])
    x_in = f.mx_in()
    x_out = f.call(x_in)
    g = ca.Function("g", [x_in[0]],[6*x_out[0]])
    n=[1,2,3,4,5,6]
    f_in=ca.DM(f.sparsity_in(0),n)
    f_out = f(f_in)
    g_in=ca.DM(g.sparsity_in(0),n)
    g_out = g(g_in)
    checkarray(self,6*f_out.full(),g_out.full(),"slicing(trans)")

  def test_scalarMX(self):
      x=ca.MX.sym("x")
      x0=0.738
      self.numpyEvaluationCheckPool(self.pool,[x],x0,name="scalarMX")

      self.numpyEvaluationCheckPool(self.matrixpool,[x],x0,name="scalarMX")

  def test_MXJacobian(self):
    self.message("MX(1,1) unary operation, jacobian")
    self.Jpool=FunctionPool()
    self.message("SX(1,1) unary operation, jacobian")
    x=ca.MX.sym("x")
    x0=array([[0.738]])

    def fmod(f,x):
      return jacobian_old(f, 0, 0)

    self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="MX unary operations, jacobian",fmod=fmod)

  def test_MXJacobians(self):
      self.message("MX(3,1) unary operation, jacobian")
      x=ca.MX.sym("x",3,1)

      x0=array([0.738,0.9,0.3])

      def fmod(f,x):
        return jacobian_old(f, 0, 0)

      self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="MX unary operations, jacobian",fmod=fmod)

  def test_MXbinary(self):
      self.message("MX binary operations")
      x=ca.MX.sym("x",3,2)
      y=ca.MX.sym("x",3,2)
      x0=array([[0.738,0.2],[ 0.1,0.39 ],[0.99,0.999999]])
      y0=array([[1.738,0.6],[ 0.7,12 ],[0,-6]])
      self.numpyEvaluationCheckPool(self.matrixbinarypool,[x,y],[x0,y0],name="MX")

  def test_MXSparse(self):
      self.message("MX unary operations, sparse")
      sp=ca.Sparsity(4,3,[0,2,2,3],[1,2,1])

      x=ca.MX.sym("x",sp)
      if scipy_available:
        x0=ca.DM(ca.Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).sparse()

        self.numpyEvaluationCheckPool(self.pool,[x],array(x0.todense()),name="MX",setx0=x0,excludeflags={'nozero'})
        self.numpyEvaluationCheckPool(self.matrixpool,[x],array(x0.todense()),name="MX",setx0=x0)
      else:
        x0=ca.DM(ca.Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).full()

        self.numpyEvaluationCheckPool(self.pool,[x],x0,name="MX",setx0=x0,excludeflags={'nozero'})
        self.numpyEvaluationCheckPool(self.matrixpool,[x],x0,name="MX",setx0=x0)

  def test_MXbinarySparse(self):
      self.message("SX binary operations")
      spx=ca.Sparsity(4,3,[0,2,2,3],[1,2,1])
      spy=ca.Sparsity(4,3,[0,2,2,3],[0,2,3])
      xx=ca.MX.sym("x",spx)
      yy=ca.MX.sym("y",spy)
      if scipy_available:
        x0=ca.DM(ca.Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).sparse()
        y0=ca.DM(ca.Sparsity(4,3,[0,2,2,3],[0,2,3]),[1.738,0.7,-6]).sparse()

        self.numpyEvaluationCheckPool(self.matrixbinarypool,[xx,yy],[array(x0.todense()),array(y0.todense())],name="MX",setx0=[x0,y0])
      else:
        x0=ca.DM(ca.Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).full()
        y0=ca.DM(ca.Sparsity(4,3,[0,2,2,3],[0,2,3]),[1.738,0.7,-6]).full()

        self.numpyEvaluationCheckPool(self.matrixbinarypool,[xx,yy],[x0,y0],name="MX",setx0=[x0,y0])

  def test_symbolcheck(self):
    self.message("Check if non-symbolic inputs are caught")
    self.assertRaises(RuntimeError, lambda : ca.Function("f", [ca.MX(0)],[ca.MX.sym("x")]))

  def test_unite(self):
    self.message("unite operation")
    import numpy
    numpy.random.seed(42)
    xn = numpy.random.random((3,4))
    x=ca.MX(3,4)
    y=ca.MX.sym("x",3,4)
    z=ca.unite(x,y)
    f = ca.Function("f", [y],[z])
    f_out = f(xn)
    self.checkarray(f_out,xn,"unite dense")

    spx=ca.Sparsity(4,3,[0,2,2,3],[1,2,1])
    spy=ca.Sparsity(4,3,[0,1,2,3],[0,2,2])

    nx=ca.DM.zeros(spx)
    for k in range(nx.nnz()):
      nx.nz[k]= numpy.random.rand()
    ny=ca.DM.zeros(spy)
    for k in range(nx.nnz()):
      ny.nz[k]= numpy.random.rand()

    nxn = nx.full()
    nyn = ny.full()
    x=ca.MX.sym("x",spx)
    y=ca.MX.sym("y",spy)
    z=ca.unite(x,y)

    f = ca.Function("f", [x,y],[z])
    f_out = f(nx, ny)
    self.checkarray(f_out,nxn+nyn,"unite sparse")

  def test_imatrix_index(self):
    self.message("IM indexing")
    X = ca.MX.sym("x",2,2)
    Y = X.nz[np.array([[0,2],[1,1],[3,3]])]

    f = ca.Function("f", [X],[Y])
    f_out = f(ca.DM(f.sparsity_in(0),[1,2,3,4]))

    self.checkarray(f_out,array([[1,3],[2,2],[4,4]]),"IM indexing")

    Y = X[:,:]
    Y.nz[np.array([[0,2]])] = ca.DM([[9,8]])  # pyright: ignore[reportCallIssue,reportArgumentType]

    f = ca.Function("f", [X],[Y])
    f_in = ca.DM(f.sparsity_in(0),[1,2,3,4])
    f_out = f(f_in)

    self.checkarray(f_out,array([[9,8],[2,4]]),"IM indexing assignment")

  def test_subsass(self):
     self.message("Check subscripted assignment")

     X = ca.MX.sym("x",2,2)
     X[0,0]=ca.MX(5)
     X[0,0]=5
     X[:,0]=8

     x=ca.MX.sym("X",3,4)
     import numpy
     numpy.random.seed(42)
     xn = numpy.random.random((3,4))
     r = numpy.zeros((7,8))
     y=ca.MX.zeros(7,8)
     y[1:4,[2,4,6,7]]=x
     r[1:4,[2,4,6,7]]=xn
     fy = ca.Function("fy", [x],[y])
     fy_out = fy(xn)

     self.checkarray(fy_out,r,"subscripted assigment")

     y=ca.MX(7,8)
     y[1:4,[2,4,6,7]]=x
     r[1:4,[2,4,6,7]]=xn
     fy = ca.Function("fy", [x],[y])
     fy_out = fy(xn)
     self.checkarray(fy_out,r,"subscripted assigment")

     kl=[2,4,5,8]

     s=y.sparsity()
     for k in kl:
       r[s.row()[k],s.get_col()[k]]=1.0

     y.nz[kl]=ca.MX(1)
     fy = ca.Function("fy", [x],[y])
     fy_out = fy(xn)
     self.checkarray(fy_out,r,"subscripted assigment")

     y.nz[kl]=x.nz[[0,1,2,3]]
     s=y.sparsity()
     sx=x.sparsity()
     cnt=0
     for k in kl:
       r[s.row()[k],s.get_col()[k]]=xn[sx.row()[cnt],sx.get_col()[cnt]]
       cnt+=1
     fy = ca.Function("fy", [x],[y])
     fy_out = fy(xn)
     self.checkarray(fy_out,r,"subscripted assigment")

  def test_erase(self):
    self.message("Erase function")
    self.message(":dense")
    y=ca.MX.sym("Y",7,8)
    import numpy
    r=2*numpy.ones((7,8))
    r[1:4,[2,4,6,7]]=numpy.zeros((3,4))
    z = y *2
    z.erase([1,2,3],[2,4,6,7])
    f = ca.Function("f", [y],[z])
    f_in = [0]*f.n_in()  # type: list

    f_in[0]=ca.DM(f.sparsity_in(0),[1]*56)
    f_out = f(0)
    e = f_out
    self.checkarray(f_out,e,"erase") # fishy
    self.message(":sparse")

  def test_MXalgebraDense(self):
    self.message("Test some dense algebraic properties of matrices")
    # issue 96
    n = 3
    m = 4
    import numpy
    numpy.random.seed(42)
    A_ = numpy.random.random((m,n))
    A = ca.MX.sym("A",m,n)
    b_ = numpy.random.random((m,1))
    b = ca.MX.sym("b",m,1)
    C_ = numpy.random.random((m,m))
    C = ca.MX.sym("C",m,m)
    D_ = numpy.random.random((m,n))
    D = ca.MX.sym("D",m,n)
    e_ = numpy.random.random((m,1))
    e = ca.MX.sym("e",m,1)
    x_ = numpy.random.random((n,1))
    x = ca.MX.sym("x",n,1)

    Axb = A @ x+b
    Dxe = D @ x+e
    a = Axb.T @ C @ Dxe

    f = ca.Function("f", [x,A,b,C,D,e],[a])
    f_out = f(x_, A_, b_, C_, D_, e_)

    f_ = numpy.dot(numpy.dot((numpy.dot(A_,x_)+b_).T,C_),(numpy.dot(D_,x_)+e_))

    self.checkarray(f_out,f_,"evaluation")


    J_ = numpy.dot(numpy.dot((numpy.dot(D_,x_)+e_).T,C_.T),A_) + numpy.dot(numpy.dot((numpy.dot(A_,x_)+b_).T,C_),D_)

    for w in [0, 1]:
      f = ca.Function("f", [x,A,b,C,D,e], [a], {"ad_weight":w, "ad_weight_sp":w})
      J = jacobian_old(f, 0, 0)
      J_in = [0]*J.n_in()  # type: list

      J_in[0]=x_
      J_in[1]=A_
      J_in[2]=b_
      J_in[3]=C_
      J_in[4]=D_
      J_in[5]=e_
      J_out = J.call(J_in)

      self.checkarray(J_out[0],J_,"evaluation")

  def test_MXalgebraSparse(self):
    self.message("Test some sparse algebraic properties of matrices")
    if not(scipy_available):
      return
    # issue 97
    n = 3
    m = 4
    import numpy
    numpy.random.seed(42)

    def randsparsity(m,n):
      sp = ca.Sparsity(m,n)
      for i in range(int((n*m)/2)):
        sp.add_nz(numpy.random.randint(m),numpy.random.randint(n))
      return sp

    def gentest(m,n):
      As = randsparsity(m,n)
      A_ = ca.DM.zeros(As)
      for k in range(As.nnz()):
        A_.nz[k]= numpy.random.rand()
      A = ca.MX.sym("A",As)
      return (A_.sparse(),A)

    (A_,A)=gentest(m,n)
    (b_,b)=gentest(m,1)
    (C_,C)=gentest(m,m)
    (D_,D)=gentest(m,n)
    (e_,e)=gentest(m,1)
    x_ = numpy.random.random((n,1))
    x = ca.MX.sym("x",n,1)

    Axb = A @ x+b
    Dxe = D @ x+e
    a = Axb.T @ C @ Dxe

    f = ca.Function("f", [x,A,b,C,D,e],[a])
    f_out = f(x_, A_, b_, C_, D_, e_)


    Axb_ = A_*x_+b_
    Dxe_ = D_*x_+e_

    f_ = Axb_.T*C_*Dxe_

    self.checkarray(f_out,f_,"evaluation")


    J_ = (D_*x_+e_).T*C_.T*A_ + (A_*x_+b_).T*C_*D_

    for w in [0, 1]:
      f = ca.Function("f", [x,A,b,C,D,e], [a], {"ad_weight":w, "ad_weight_sp":w})

      J = jacobian_old(f, 0, 0)
      J_in = [0]*J.n_in()  # type: list

      J_in[0]=x_
      J_in[1]=A_
      J_in[2]=b_
      J_in[3]=C_
      J_in[4]=D_
      J_in[5]=e_
      J_out = J.call(J_in)

      self.checkarray(J_out[0],J_,"evaluation")

  #@unittest.skipIf(not(scipy_available))
  def test_MXalgebraSparseSparse(self):
    if not(scipy_available):
        return
    self.message("Test some sparse algebraic properties of matrices")
    n = 8
    m = 10
    import numpy
    numpy.random.seed(42)

    def randsparsity(m,n):
      sp = ca.Sparsity(m,n)
      for k in range(int((n*m)/2)):
        i = numpy.random.randint(m)
        j = numpy.random.randint(n)
        if not(i == int(m/2)):
          if n==1 or not(j == int(n/2)):
            sp.add_nz(i,j)
      return sp

    def gentest(m,n):
      As = randsparsity(m,n)
      A_ = ca.DM.zeros(As)
      for k in range(As.nnz()):
        A_.nz[k]= numpy.random.rand()
      A = ca.MX.sym("A",As)
      return (A_.sparse(),A)

    (A_,A)=gentest(m,n)
    (b_,b)=gentest(m,1)
    (C_,C)=gentest(m,m)
    (D_,D)=gentest(m,n)
    (e_,e)=gentest(m,1)
    x_ = numpy.random.random((n,1))
    x = ca.MX.sym("x",n,1)

    Axb = A @ x+b
    Dxe = D @ x+e
    a = Axb.T @ C @ Dxe

    f = ca.Function("f", [x,A,b,C,D,e],[a])
    f_out = f(x_, A_, b_, C_, D_, e_)


    Axb_ = A_*x_+b_
    Dxe_ = D_*x_+e_

    f_ = Axb_.T*C_*Dxe_

    self.checkarray(f_out,f_,"evaluation")


    J_ = (D_*x_+e_).T*C_.T*A_ + (A_*x_+b_).T*C_*D_

    for w in [0, 1]:
      f = ca.Function("f", [x,A,b,C,D,e], [a], {"ad_weight":w, "ad_weight_sp":w})
      J = jacobian_old(f, 0, 0)
      J_in = [0]*J.n_in()  # type: list

      J_in[0]=x_
      J_in[1]=A_
      J_in[2]=b_
      J_in[3]=C_
      J_in[4]=D_
      J_in[5]=e_
      J_out = J.call(J_in)

      self.checkarray(J_out[0],J_,"evaluation")


  def test_chaining(self):
    self.message("Chaining SX and MX together")
    x=ca.SX.sym("x")
    y=x**3
    f=ca.Function("f", [x],[y])
    J = jacobian_old(f, 0, 0)

    X=ca.MX.sym("X")
    F=ca.Function("F", [X], list(J(X)))

    x_=1.7
    F_out = F(x_)
    self.checkarray(F_out[0],3*x_**2,"Chaining eval")

  def test_issue107(self):
    self.message("Regression test for issue 107: +=")
    x=ca.MX.sym("x")
    y=ca.MX.sym("y")

    z=x
    z+=y

    self.assertTrue(x.is_symbolic())
    self.assertFalse(z.is_symbolic())

  def test_MXd_trivial(self):
    self.message("symbolic variables and constants jac")
    X =  ca.MX.sym("X",10)
    V =  ca.MX.sym("V")
    J = ca.jacobian(X,X)
    self.assertTrue(isinstance(J,ca.MX))
    self.assertEqual(J.nnz(),10)
    self.assertEqual(J.size1(),10)
    self.assertEqual(J.size2(),10)

    g = ca.Function("g", [],[J])
    [g_out] = g.call([])
    self.checkarray(g_out,eye(10),"unit matrix")
    g = ca.Function("g", [],[ca.jacobian(ca.MX.eye(3),X)])
    [g_out] = g.call([])
    self.checkarray(g_out,zeros((9,10)),"zero matrix")
    g = ca.Function("g", [],[ca.jacobian(X,V)])
    [g_out] = g.call([])
    self.checkarray(g_out,zeros((10,1)),"zero matrix")

    g = ca.Function("g", [],[ca.jacobian(ca.MX.eye(3),V)])
    [g_out] = g.call([])
    self.checkarray(g_out,zeros((9,1)),"zero matrix")

  def test_MXd_substractionl(self):
    self.message("substraction jac")
    V =  ca.MX.sym("V")
    X =  ca.MX.sym("X")
    g = ca.Function("g", [],[ca.jacobian(X-V, X)])
    [g_out] = g.call([])
    self.checkarray(g_out,ones((1,1)), "one")

    g = ca.Function("g", [],[ca.jacobian(X-V, V)])
    [g_out] = g.call([])
    self.checkarray(g_out,-ones((1,1)), "one")

    g = ca.Function("g", [],[ca.jacobian(V-X, X)])
    [g_out] = g.call([])
    self.checkarray(g_out,-ones((1,1)), "one")

    g = ca.Function("g", [],[ca.jacobian(V-X, V)])
    [g_out] = g.call([])
    self.checkarray(g_out,ones((1,1)),"one")

  def test_MXd_mapping(self):
    self.message("mapping jac")
    X = ca.MX.sym("X",3)
    Y = ca.MX.sym("Y",2)
    J = ca.jacobian(ca.vertcat(X,Y),X)
    JJ = ca.DM.ones(J.sparsity())
    self.checkarray(JJ,numpy.vstack((eye(3),zeros((2,3)))),"diag")
    J = ca.jacobian(ca.vertcat(X,Y),Y)
    JJ = ca.DM.ones(J.sparsity())
    self.checkarray(JJ,numpy.vstack((zeros((3,2)),eye(2))),"diag")

  def test_null(self):
    self.message("Function null")
    x = ca.MX.sym("x")

    f = ca.Function("f", [x],[x**2,ca.MX()])

    self.assertEqual(f.size1_out(1),0)
    self.assertEqual(f.size2_out(1),0)
    f_out = f(0)

    f = ca.Function("f", [x,ca.MX()],[x**2,ca.MX()])

    self.assertEqual(f.size1_out(1),0)
    self.assertEqual(f.size2_out(1),0)
    f_out = f(0,0)

    r = f(x,ca.MX())
    self.assertTrue(r[1].is_empty(True))

    r = f(ca.MX(),ca.MX())
    self.assertTrue(r[1].is_empty(True))

    #self.assertRaises(Exception,lambda : f([x,x],True))
    #self.assertRaises(Exception,lambda : f([[],[]],True))

  def test_issue184(self):
    self.message("Regression test issue #184")
    x = ca.MX.sym("x", 3)
    y = x[0:0]
    self.assertEqual(y.nnz(),0)

  def test_indexingOutOfBounds(self):
    self.message("Indexing out of bounds")
    y = ca.DM.zeros(4, 5)
    self.assertRaises(RuntimeError,lambda : y[12,0] )
    self.assertRaises(RuntimeError,lambda : y[12,12] )
    self.assertRaises(RuntimeError,lambda : y[0,12] )
    self.assertRaises(RuntimeError,lambda : y[12,:] )
    self.assertRaises(RuntimeError,lambda : y[12:15,0] )
    self.assertRaises(RuntimeError,lambda : y[:,12] )
    self.assertRaises(RuntimeError,lambda : y[0,12:15] )
    y[-1,2]
    self.assertRaises(RuntimeError,lambda : y[-12,2] )
    y[-3:-1,2]
    self.assertRaises(RuntimeError,lambda : y[-12:-9,2] )

    def test():
      y[12,0] = 0
    self.assertRaises(RuntimeError,test)
    def test():
      y[12,12] = 0
    self.assertRaises(RuntimeError,test)
    def test():
      y[0,12] = 0
    self.assertRaises(RuntimeError,test)
    def test():
      y[12,:] = 0
    self.assertRaises(RuntimeError,test)
    def test():
      y[12:15,0] = 0
    self.assertRaises(RuntimeError,test)
    def test():
      y[:,12] = 0
    self.assertRaises(RuntimeError,test)
    def test():
      y[0,12:15] = 0
    self.assertRaises(RuntimeError,test)
    y[-1,2] = 0
    def test():
      y[-12,2] = 0
    self.assertRaises(RuntimeError,test)
    y[-3:-1,2] = 0
    def test():
      y[-12:-9,2] = 0
    self.assertRaises(RuntimeError,test)

  def test_indexinglimits(self):
    self.message("Limits of indexing")
    y = ca.MX.sym("y", 3)
    self.assertRaises(RuntimeError,lambda : y[[0, 5]] )
    try:
      y[[0, 5]] = ca.MX.sym("a")
      self.assertTrue(False)
    except RuntimeError:
      pass
    y[[0, 2]]
    y[[0, 2]] = ca.MX.sym("a")

  def test_mtimes(self):
    A = ca.MX(ca.DM.ones((4,3)))
    B = ca.MX(ca.DM.ones((3,8)))
    C = ca.MX(ca.DM.ones((8,7)))

    self.assertRaises(RuntimeError,lambda : ca.mtimes([]))

    D = ca.mtimes([A])

    self.assertEqual(D.shape[0],4)
    self.assertEqual(D.shape[1],3)

    D = ca.mtimes([A,B])

    self.assertEqual(D.shape[0],4)
    self.assertEqual(D.shape[1],8)

    D = ca.mtimes([A,B,C])

    self.assertEqual(D.shape[0],4)
    self.assertEqual(D.shape[1],7)



  @memory_heavy()
  def test_mtimes_blas_plugin(self):
    import os, tempfile
    
    ca.DM.rng(1)
    
    m, k, n = 3, 4, 5

    def parse_sp(art):
      """ASCII art -> Sparsity. 'X' = nonzero, '.' (or anything else) = zero."""
      rows = [r.strip() for r in art.strip().split('\n') if r.strip()]
      return ca.sparsify(ca.DM([[1.0 if c == 'X' else 0.0 for c in r]
                          for r in rows])).sparsity()

    # A is m x k = 3 x 4. Compactible: rows {0,2} x cols {0,2}.
    A_pats = [
      ('dense', parse_sp("""
          XXXX
          XXXX
          XXXX""")),
      ('compactible', parse_sp("""
          X.X.
          ....
          X.X.""")),
      ('sparse', parse_sp("""
          X.XX
          .X.X
          XX..""")),
    ]
    # B is k x n = 4 x 5. Compactible: rows {0,2} x cols {1,3} -- the row set
    # matches A's compactible col set so the (compactible, compactible) combo
    # exercises the PseudoDenseMultiplication path.
    B_pats = [
      ('dense', parse_sp("""
          XXXXX
          XXXXX
          XXXXX
          XXXXX""")),
      ('compactible', parse_sp("""
          .X.X.
          .....
          .X.X.
          .....""")),
      ('sparse', parse_sp("""
          X.X.X
          .X..X
          XX.X.
          .X.XX""")),
    ]
    Adense_val = ca.DM.rand(m, k)
    Bdense_val = ca.DM.rand(k, n)

    plugins = ['reference']
      
    for opt in ['classic', 'blasfeo']:
      try:
        ca.load_blas(opt); plugins.append(opt)
      except Exception:
        pass

    # Per-plugin spelling of the dense matrix kernel. Reference emits the
    # canonical casadi runtime call; classic/blasfeo emit through the macro
    # / direct BLAS symbol.
    dense_kernel = {
        'reference': 'casadi_mtimes_dense(',
        'classic':   'CASADI_BLAS_DGEMM(',
        'blasfeo':   'blasfeo_blas_dgemm(',
    }

    # Link line for codegen+compile of each plugin. The "classic" line is
    # taken from CasadiMeta.lapack_libraries() so the test stays agnostic to
    # the build's BLAS choice (OpenBLAS / MKL / system BLAS / ...). It's
    # CMake's native LAPACK_LIBRARIES form: a semicolon-separated mix of
    # full paths and -l flags. helpers.py's check_codegen passes paths
    # through verbatim and prefixes bare names with -l.
    classic_libs = [t for t in ca.CasadiMeta.lapack_libraries().split(";") if t]
    plugin_extralibs = {
        'reference': [],
        'classic':   classic_libs,
        'blasfeo':   ['blasfeo'],
    }

    def predict_category(A_sp, B_sp, z_sp):
      if A_sp.is_dense() and B_sp.is_dense() and z_sp.is_dense():
        return 'dense'             # DenseMultiplication
      if A_sp.is_dense() and z_sp.is_dense():
        return 'dense_sparse'      # DenseSparseMultiplication
      ok_a, ar, ac = A_sp.is_compactible()
      ok_b, br, bc = B_sp.is_compactible()
      ok_z, zr, zc = z_sp.is_compactible()
      if (ok_a and ok_b and ok_z and z_sp.nnz() > 0
          and ac == br and ar == zr and bc == zc):
        return 'dense'             # PseudoDenseMultiplication
      return 'generic'             # base Multiplication

    Adense_sym = ca.MX.sym('A', m, k)
    Bdense_sym = ca.MX.sym('B', k, n)
    Zdense_sym = ca.MX.sym('Z', m, n)
    Zdense_val = ca.DM(np.random.randn(m, n))

    for name in plugins:
      for a_label, a_sp in A_pats:
        for b_label, b_sp in B_pats:
          for op_label in ['mtimes', 'mac']:
            tag = "%s/%s/%s/%s" % (op_label, name, a_label, b_label)
            A = ca.project(Adense_sym, a_sp)
            B = ca.project(Bdense_sym, b_sp)
            Amask = ca.densify(ca.DM(a_sp, 1.0))
            Bmask = ca.densify(ca.DM(b_sp, 1.0))

            if op_label == 'mtimes':
              # mtimes: Z is implicit MX::zeros with sparsity mul(A_sp, B_sp).
              # densify on output so f and f_ref share output sparsity;
              # otherwise checkfunction's random adjoint seeds differ in size
              # and forward/adjoint sensitivities aren't comparable.
              out_f   = ca.densify(ca.mtimes(A, B, name))
              out_ref = ca.mtimes(Adense_sym * Amask, Bdense_sym * Bmask, "reference")
              sym_in, val_in = [Adense_sym, Bdense_sym], [Adense_val, Bdense_val]
              z_sp = (ca.DM(a_sp, 1.0) @  ca.DM(b_sp, 1.0)).sparsity()
            else:
              # mac: Z is dense, supplied explicitly. Output is Z + A*B and
              # already dense, so no densify needed.
              out_f   = ca.mac(A, B, Zdense_sym, name)
              out_ref = ca.mac(Adense_sym * Amask, Bdense_sym * Bmask,
                            Zdense_sym, "reference")
              sym_in = [Adense_sym, Bdense_sym, Zdense_sym]
              val_in = [Adense_val, Bdense_val, Zdense_val]
              z_sp = ca.Sparsity.dense(m, n)

            f     = ca.Function('f',     sym_in, [out_f])
            f_ref = ca.Function('f_ref', sym_in, [out_ref])

            self.checkfunction(f, f_ref, inputs=val_in)

            # Also exercise the DM-level Matrix<double>::mtimes / ::mac path
            # directly. Project the dense values, call the operation on DMs,
            # and compare against the reference Function output (which uses
            # the all-dense kernel via "reference").
            A_dm = ca.project(Adense_val, a_sp)
            B_dm = ca.project(Bdense_val, b_sp)
            if op_label == 'mtimes':
              dm_out = ca.mtimes(A_dm, B_dm, name)
            else:
              dm_out = ca.mac(A_dm, B_dm, Zdense_val, name)
            self.checkarray(ca.densify(dm_out), f_ref(*val_in),
                            "%s (DM-level)" % tag)

            cat = predict_category(a_sp, b_sp, z_sp)
            if cat == 'dense':
              expected = dense_kernel[name]
            elif cat == 'dense_sparse':
              expected = 'casadi_mtimes_dense_sparse('  # plugin-agnostic
            else:
              expected = 'casadi_mtimes('               # plugin-agnostic

            # Is the expected kernel in the generated code?
            f.generate('f.c')
            with open('f.c', 'r') as codefile:
                code = codefile.read()
            self.assertIn(expected, code,
                "%s: expected kernel %s in generated code" % (tag, expected))

            # Does serialization preserve the expected kernel?
            f.save('f.casadi')
            fdeser = ca.Function.load('f.casadi')
            fdeser.generate('fdeser.c')
            with open('fdeser.c', 'r') as codefile:
                code = codefile.read()
            self.assertIn(expected, code,
                "%s: expected kernel %s in deserialized code" % (tag, expected))

            self.checkfunction_light(f, fdeser, inputs=val_in)

            extralibs = plugin_extralibs.get(name, [])
            if name == 'reference' or extralibs:
              if expected == dense_kernel[name] or name=="reference":
                self.check_codegen(f, inputs=val_in, extralibs=extralibs)


  def test_l1_blas_plugin(self):
    # Counterpart to test_mtimes_blas_plugin, but for the L1 ops dispatched
    # through GlobalOptions::setDefaultBlas: dot, norm_1, norm_2, plus axpy
    # (which appears in MapSum codegen, not directly in MX). Simpler input
    # space (vectors, no sparsity combinatorics) so we exercise more
    # operations instead.

    ca.DM.rng(2)

    plugins = ['reference']
    for opt in ['classic', 'blasfeo']:
      try:
        ca.load_blas(opt); plugins.append(opt)
      except Exception:
        pass

    classic_libs = [t for t in ca.CasadiMeta.lapack_libraries().split(";") if t]
    plugin_extralibs = {
        'reference': [],
        'classic':   classic_libs,
        'blasfeo':   ['blasfeo'],
    }

    # Per (plugin, op) -> substring that should appear in the generated C
    # for the L1 op's auxiliary block. None means "this plugin falls back
    # to reference for that op". The reference column documents the
    # baseline shape (inline static definition of the casadi_*<T> body).
    expected_kernel = {
      ('reference', 'dot'    ): 'casadi_real casadi_dot(',
      ('reference', 'norm_1' ): 'casadi_real casadi_norm_1(',
      ('reference', 'norm_2' ): 'casadi_real casadi_norm_2(',
      ('reference', 'axpy'   ): 'void casadi_axpy(',
      # copy is never plugin-dispatched (no BLAS analog for the
      # x==NULL -> zero-fill semantics); always emits the memcpy/memset
      # body from Blas::codegen_copy_aux regardless of default_blas.
      ('reference', 'copy'   ): 'memcpy(y, x',
      ('classic',   'dot'    ): 'CASADI_BLAS_DDOT',
      ('classic',   'norm_1' ): 'CASADI_BLAS_DASUM',
      ('classic',   'norm_2' ): 'CASADI_BLAS_DNRM2',
      ('classic',   'axpy'   ): 'CASADI_BLAS_DAXPY',
      ('classic',   'copy'   ): 'memcpy(y, x',
      ('blasfeo',   'dot'    ): 'blasfeo_blas_ddot',
      # blasfeo plugin only supplies daxpy/ddot; norm_1/norm_2 fall back to
      # the inline reference template (no namespaced extern emitted).
      ('blasfeo',   'norm_1' ): 'casadi_real casadi_norm_1(',
      ('blasfeo',   'norm_2' ): 'casadi_real casadi_norm_2(',
      ('blasfeo',   'axpy'   ): 'blasfeo_blas_daxpy',
      ('blasfeo',   'copy'   ): 'memcpy(y, x',
    }

    n = 7
    x_sym = ca.MX.sym('x', n)
    y_sym = ca.MX.sym('y', n)
    x_val = ca.DM.rand(n)
    y_val = ca.DM.rand(n)

    # Reference baseline (default_==0) for numeric comparison across plugins.
    # setDefaultBlas selects the plugin for both runtime dispatch and (when the
    # 'l1_blas' codegen option is on) the L1 codegen path; codegen ignores it
    # unless l1_blas=True.
    ca.GlobalOptions.setDefaultBlas("reference")

    # Op -> (output expression, expected reference numeric value)
    def build_ops():
      return {
        'dot':    ca.dot(x_sym, y_sym),
        'norm_1': ca.norm_1(x_sym),
        'norm_2': ca.norm_2(x_sym),
        # axpy: emitted inside MapSum codegen. Wrap a tiny inner Function
        # and reduce-sum its output across N reps.
      }

    ref_ops = build_ops()
    f_ref = ca.Function('f_ref', [x_sym, y_sym], list(ref_ops.values()),
                     ['x', 'y'], list(ref_ops.keys()))
    ref_out = f_ref(x_val, y_val)

    # Build the axpy carrier separately: a MapSum-shaped Function that the
    # MX codegen lowers using g.axpy(...) -> AUX_AXPY.
    inner_x = ca.MX.sym('x', n)
    inner = ca.Function('inner', [inner_x], [2.0 * inner_x])
    f_axpy_ref = inner.map(4, [False], [True])
    f_axpy_inputs = [ca.horzcat(*[x_val for _ in range(4)])]
    f_axpy_ref_out = f_axpy_ref(*f_axpy_inputs)

    try:
      for name in plugins:
        ca.GlobalOptions.setDefaultBlas(name)
        tag_base = "blas=%s" % name

        # Numerics: each L1 op should agree with the reference baseline.
        ops = build_ops()
        f = ca.Function('f', [x_sym, y_sym], list(ops.values()),
                     ['x', 'y'], list(ops.keys()))
        got = f(x_val, y_val)
        for i, op_name in enumerate(ops.keys()):
          self.checkarray(got[i], ref_out[i],
                          "%s op=%s numerical agreement" % (tag_base, op_name))

        # axpy via MapSum:
        f_axpy = inner.map(4, [False], [True])
        got_axpy = f_axpy(*f_axpy_inputs)
        self.checkarray(got_axpy, f_axpy_ref_out,
                        "%s op=axpy numerical agreement" % tag_base)

        # Codegen: per-op kernel string check. Generate a function whose
        # body uses *only* one op so the assertion isn't muddied by other
        # auxiliaries that pull in unrelated symbols (e.g. AUX_NORM_2 also
        # pulls AUX_DOT in the reference path).
        per_op_carriers = {
          'dot':    ca.Function('g_dot',    [x_sym, y_sym], [ca.dot(x_sym, y_sym)]),
          'norm_1': ca.Function('g_norm_1', [x_sym],        [ca.norm_1(x_sym)]),
          'norm_2': ca.Function('g_norm_2', [x_sym],        [ca.norm_2(x_sym)]),
          'axpy':   inner.map(4, [False], [True]),
          # Identity Function emits a casadi_copy from the input slot to the
          # output slot during codegen, pulling in AUX_COPY.
          'copy':   ca.Function('g_copy',   [x_sym],        [x_sym]),
        }
        # Codegen L1 routing is opt-in via the 'l1_blas' CodeGenerator option
        # (default off => reference); when on it follows the active default_blas.
        for op_name, g in per_op_carriers.items():
          g.generate('l1.c', {'l1_blas': True})
          with open('l1.c', 'r') as codefile:
            code = codefile.read()
          self.assertIn(expected_kernel[(name, op_name)], code,
              "%s op=%s: expected %r in generated code"
              % (tag_base, op_name, expected_kernel[(name, op_name)]))

        # Serialization round-trip: the BLAS choice is a session-wide
        # GlobalOption (not stored on the Function), so a deserialized
        # Function picks up whatever default_blas is active at codegen
        # time. Verify it round-trips to the SAME emission with the same
        # default_blas in effect.
        f.save('l1.casadi')
        fdeser = ca.Function.load('l1.casadi')
        self.checkfunction_light(f, fdeser, inputs=[x_val, y_val])

        # check_codegen end-to-end: compile + dlopen + evaluate. Skip for
        # plugins where the Python bindings or system libs don't resolve.
        extralibs = plugin_extralibs.get(name, [])
        if name == 'reference' or extralibs:
          self.check_codegen(f, inputs=[x_val, y_val], extralibs=extralibs,
                             opts={'l1_blas': True})
          # Also check the axpy carrier so the AUX_AXPY emission goes
          # through the same compile+run gauntlet.
          self.check_codegen(f_axpy, inputs=f_axpy_inputs, extralibs=extralibs,
                             opts={'l1_blas': True})
    finally:
      # Restore so other tests in the session aren't surprised.
      ca.GlobalOptions.setDefaultBlas("reference")


  def test_truth(self):
    self.message("Truth values")
    self.assertRaises(Exception, lambda : bool(ca.MX.sym("x")))
    #self.assertRaises(Exception, lambda : bool(MX.sym("x")>0))
    self.assertTrue(bool(ca.MX(1)))
    self.assertFalse(bool(ca.MX(0)))
    self.assertTrue(bool(ca.MX(0.2)))
    self.assertTrue(bool(ca.MX(-0.2)))
    self.assertRaises(Exception, lambda : bool(ca.MX(ca.DM([2.0,3]))))
    self.assertRaises(Exception, lambda : bool(ca.MX()))


  def test_MXbool(self):
    self.message("bool")

    xy = ca.MX.sym("x",2)
    x = xy[0]
    y = xy[1]

    f = ca.Function("f", [xy],[ca.vertcat(*[ca.logic_and(x,y),ca.logic_or(x,y),ca.logic_not(x)])])


    for t1 in [0,1]:
      for t2 in [0,1]:
        T1 = t1!=0
        T2 = t2!=0
        f_out = f([t1,t2])
        self.checkarray(f_out,ca.DM([T1 and T2,T1 or T2,not T1]),"bool(%d,%d): %s" % (t1,t2,str(f_out)))
  
  def test_short_circuiting_codegen(self):
    x = ca.MX.sym('x', 5)

    y = ca.if_else(x > 0, 0, x)

    f = ca.Function('f', [x], [y], ['x'], ['y'])
    
    self.check_codegen(f,inputs=[ca.vertcat(-200, 200, -100, 100, 0)])

    x = ca.MX.sym('x', 8)
    y = ca.MX.sym('y', 8)
    
    z = ca.logic_and(x,y)
    z2 = ca.logic_or(x,y)

    f = ca.Function('f', [x,y], [z,z2])
    
    self.check_codegen(f,inputs=[ca.vertcat(0, 1, 0, 1, 0, 1, 0, 1),ca.vertcat(0, 0, 1, 1, 0, 0, 1, 1)])



  def test_MXineq(self):
    self.message("SX ineq")

    xy = ca.MX.sym("x",2)
    x = xy[0]
    y = xy[1]


    f = ca.Function("f", [xy],[ca.vertcat(*[x<y,x<=y,x>=y,x==y,x!=y])])

    for t1 in [-10,0.1,0,1,10]:
      for t2 in [-10,0.1,0,1,10]:
        T1 = t1
        T2 = t2
        f_out = f([t1,t2])
        self.checkarray(f_out,ca.DM([T1 < T2,T1 <= T2, T1 >= T2, T1 == T2, T1 != T2]),"ineq(%d,%d)" % (t1,t2))

  def test_if_else_zero(self):
    x = ca.MX.sym("x")
    y = ca.if_else(x,5,0)
    f = ca.Function("f", [x],[y])
    f_in = 1
    f_out = f(f_in)
    self.assertTrue(f_out==5,"if_else_zero %s " % str(f_out))
    f_in = 0
    f_out = f(f_in)
    self.assertTrue(f_out==0,"if_else_zero")


  def test_if_else(self):
    x = ca.MX.sym("x")
    y = ca.if_else(x,1,2)
    f = ca.Function("f", [x],[y])
    f_in = 1
    f_out = f(f_in)
    self.assertTrue(f_out==1,"if_else")
    f_in = 0
    f_out = f(f_in)
    self.assertTrue(f_out==2,"if_else")

  def test_regression491(self):
    self.message("regression #491")
    u = ca.SX.sym("u")
    x = ca.SX.sym("x")

    F = ca.Function("F", [u,x],[u+1/x])

    U = ca.MX.sym("U")

    X = F(U,U)
    G = F(U,X)

    for kk in range(2):
      gfcn = 0
      if kk==0:
        gfcn = ca.Function("gfcn", [U], [G]).expand("e_gfcn", {"ad_weight":1})
      else:
        gfcn = ca.Function("gfcn", [U],[G], {"ad_weight":1})
      J = jacobian_old(gfcn, 0, 0)
      J_in = [0]*J.n_in()  # type: list

      J_in[0]=1
      J_out = J.call(J_in)
      self.assertAlmostEqual(J_out[0],1,9)

  def test_ticket(self):
    J = [] + ca.MX.sym("x")
    J = ca.MX.sym("x") + []

  def test_jacobian_tools(self):
    self.message("jacobian")

    X = ca.MX.sym("X")

    Y = ca.jacobian(X**2,X)

    f = ca.Function("f", [X], [Y])

    f_in=2.3
    f_out = f(f_in)

    self.assertAlmostEqual(f_out,4.6)

  def test_reshape(self):
    self.message("reshape")
    X = ca.MX.sym("X",10)

    i = ca.DM(ca.Sparsity.lower(3),list(range(6)))

    i.print_dense()
    print(i.T.nz[:])

    T = X.nz[i]

    q = T.T.nz[:]**2
    f = ca.Function("f", [X], [q])
    f_out = f(list(range(10)))

    self.checkarray(ca.DM([0,1,9,4,16,25]),f_out)

    Y = ca.MX.sym("Y",10)

    ff = ca.Function("ff", [Y],f.call([Y],True))
    ff_out = ff(list(range(10)))

    self.checkarray(ca.DM([0,1,9,4,16,25]),ff_out)

    J = ca.Function("J", [X],[ca.jacobian(q, X)])
    J_out = J(list(range(10)))

    i = ca.horzcat(*[ca.diag([0,2,4,6,8,10]),ca.DM.zeros(6,4)])
    i[[2,3],:] = i[[3,2],:]

    self.checkarray(i,J_out)

    q = (T.T).nz[:]**2
    J = ca.Function("J", [X],[ca.jacobian(q,X)])
    J_out = J(list(range(10)))

    i = ca.horzcat(*[ca.diag([0,2,4,6,8,10]),ca.DM.zeros(6,4)])
    i[[2,3],:] = i[[3,2],:]

    self.checkarray(i,J_out)

  def test_vertcat(self):
    self.message("vertcat")
    X = ca.MX.sym("X",10)

    T = ca.vertcat(*[X[4],X[2]])
    q = T**2
    f = ca.Function("f", [X],[q])
    f_in = [0]*f.n_in()  # type: list

    f_in[0]=list(range(10))
    f_out = f.call(f_in)

    self.checkarray(ca.DM([16,4]),f_out[0])

    Y = ca.MX.sym("Y",10)

    ff = ca.Function("ff", [Y],f.call([Y],True))
    ff_out = ff(list(range(10)))

    self.checkarray(ca.DM([16,4]),ff_out)
    J = ca.Function("J", [X],[ca.jacobian(q,X)])
    J_out = J(list(range(10)))

    i = ca.DM.zeros(2,10)
    i[0,4] = 8
    i[1,2] = 4

    self.checkarray(i,J_out)
    q = T**2
    J = ca.Function("J", [X],[ca.jacobian(q, X)])
    J_out = J(list(range(10)))

    self.checkarray(i,J_out)

  def test_blockcat(self):
    x = ca.MX.sym("x")

    y = ca.blockcat([[x,2*x],[3*x,4*x]])
    f = ca.Function("f", [x],[y])
    f_out = f(3)
    self.checkarray(f_out,ca.DM([[3,6],[9,12]]))


  def test_veccats(self):
    x= ca.MX.sym("x",2)
    self.assertTrue(hash(ca.vec(x))==hash(x))

  def test_constmxmtimes(self):
    0.1*ca.MX.ones(2)

  def test_is_regular(self):
    self.assertTrue(ca.MX(ca.DM([0,1])).is_regular())
    self.assertFalse(ca.MX(ca.DM([0,inf])).is_regular())
    with self.assertRaises(Exception):
      self.assertFalse(ca.MX.sym("x",2).is_regular())

  def test_diagcat(self):
    C = ca.diagcat(*[ca.MX(ca.DM(([[-1.4,-3.2],[-3.2,-28]]))),ca.DM([[15,-12,2.1],[-12,16,-3.8],[2.1,-3.8,15]]),1.8,-4.0])
    self.assertTrue(isinstance(C,ca.MX))
    r = ca.DM([[-1.4,-3.2,0,0,0,0,0],[-3.2,-28,0,0,0,0,0],[0,0,15,-12,2.1,0,0],[0,0,-12,16,-3.8,0,0],[0,0,2.1,-3.8,15,0,0],[0,0,0,0,0,1.8,0],[0,0,0,0,0,0,-4]])
    r = ca.sparsify(r)
    f = ca.Function("f", [],[C])
    f_out = f.call([])

    self.checkarray(f_out[0],r)
    
    cx = 3
    cy = 4
    for X in [ca.DM,ca.SX,ca.MX]:
        for nx in [2,1,0]:
            for ny in [2,1,0]:
                print(X,nx,ny)
                r = ca.diagcat(X(nx,ny),X.ones(cx,cy),X(nx,ny))
                self.assertEqual(r.shape[0],cx+2*nx)
                self.assertEqual(r.shape[1],cy+2*ny)

  def test_tril2symm(self):
    x = ca.MX.sym("x",ca.Sparsity.lower(3))
    f = ca.Function("f", [x],[ca.tril2symm(x)])
    f_out = f(ca.DM(f.sparsity_in(0),list(range(6))))
    self.checkarray(f_out,ca.DM([[0,1,2],[1,3,4],[2,4,5]]))

  def test_sparsity_indexing(self):
    self.message("sparsity")

    B_ = ca.DM([[1,2,3,4,5],[6,7,8,9,10]])
    B = ca.MX.sym("B",2,5)

    A = ca.DM([[1,1,0,0,0],[0,0,1,0,0]])
    A = ca.sparsify(A)
    sp = A.sparsity()
    import copy

    def meval(m):
      f = ca.Function("f", [B],[m])
      f_out = f(B_)
      return f_out

    self.checkarray(meval(B[sp]),ca.DM([[1,2,0,0,0],[0,0,8,0,0]]),"sparsity indexing")

    Bmod = copy.copy(B)
    Bmod[sp] = -4

    self.checkarray(meval(Bmod),ca.DM([[-4,-4,3,4,5],[6,7,-4,9,10]]),"sparsity indexing assignement")

    Bmod = copy.copy(B)
    Bmod[sp] = 2*B

    self.checkarray(meval(Bmod),ca.DM([[2,4,3,4,5],[6,7,16,9,10]]),"Imatrix indexing assignement")

    self.assertRaises(Exception, lambda : B[ca.Sparsity.dense(4,4)])

  def test_symvar(self):
    a = ca.MX.sym("a")
    b = ca.MX.sym("b")
    c = ca.MX.sym("c")
    e = ca.cos(a*b) + c
    w = ca.symvar(e)
    self.assertEqual(len(w),3)
    if ca.GlobalOptions.getSimplificationOnTheFly():
      self.assertTrue(ca.is_equal(w[0],a))
      self.assertTrue(ca.is_equal(w[1],b))
      self.assertTrue(ca.is_equal(w[2],c))

  def test_simplifications(self):
    for X in [ca.SX,ca.MX]:
        a = X.sym("a")
        b = X.sym("b")
        c = X.sym("c")
        
        A = X.sym("A",2,2)
        B = X.sym("B",2,2)
        
        args = [a,b,c,A,B]
        x = ca.sin(a)
        y = ca.sin(b)
        z = ca.sin(c)
        
        X = ca.sin(A)
        Y = ca.sin(B)
        
        ca.DM.rng(1)
        
        a0 = ca.DM.rand(a.sparsity())
        b0 = ca.DM.rand(b.sparsity())
        c0 = ca.DM.rand(c.sparsity())
        A0 = ca.DM.rand(A.sparsity())
        B0 = ca.DM.rand(B.sparsity())

        x0 = ca.sin(a0)
        y0 = ca.sin(b0)
        z0 = ca.sin(c0)
        X0 = ca.sin(A0)
        Y0 = ca.sin(B0)
        
        dx = x-y
          
        count = 1477
 
        for refcount,on_the_fly,genA,B in [
            (False,True,lambda x,y,z, X,Y :y*(x)+y*(1-x), y),
            (False,X is ca.MX,lambda x,y,z, X,Y :2*x-x, x),
            (False,True,lambda x,y,z, X,Y :4*x-3*x,x),
            (False,True,lambda x,y,z, X,Y :0.2*x+0.8*x,x),
            (False,True,lambda x,y,z, X,Y :(x + 0), x),
            (False,X is ca.MX,lambda x,y,z, X,Y :(0 + x), x),
            (False,True,lambda x,y,z, X,Y :(x + (-y)), x - y),
            (False,X is ca.MX,lambda x,y,z, X,Y :((-x) + y), y - x),
            (False,True,lambda x,y,z, X,Y :((0.5*x) + (0.5*x)), x),
            (False,True,lambda x,y,z, X,Y :((x/2) + (x/2)), x),
            (False,True,lambda x,y,z, X,Y :((x - y) + y), x),
            (False,True,lambda x,y,z, X,Y :(y + (x - y)), x),
            (False,True,lambda x,y,z, X,Y :(ca.sin(x)**2 + ca.cos(x)**2), 1),
            (False,True,lambda x,y,z, X,Y :(x - 0), x),
            (False,X is ca.MX,lambda x,y,z, X,Y :(0 - x), -x),
            (False,True,lambda x,y,z, X,Y :(x - x), 0),
            (False,True,lambda x,y,z, X,Y :(x - (-y)), x + y),
            (False,True,lambda x,y,z, X,Y :((x + y) - y), x),
            (False,True,lambda x,y,z, X,Y :((x + y) - x), y),
            (False,True,lambda x,y,z, X,Y :(x - (y + x)), -y),
            (False,True,lambda x,y,z, X,Y :(x - (x + y)), -y),
            (False,True,lambda x,y,z, X,Y :((-x) - y), -(x + y)),
            (False,True,lambda x,y,z, X,Y :(x * x), x**2),
            (False,True,lambda x,y,z, X,Y :(x * 3), 3 * x),
            (False,True,lambda x,y,z, X,Y :(x * 0), 0),
            (False,True,lambda x,y,z, X,Y :(0 * x), 0),
            (False,True,lambda x,y,z, X,Y :(x * 1), x),
            (False,X is ca.MX,lambda x,y,z, X,Y :(1 * x), x),
            (False,True,lambda x,y,z, X,Y :(x * -1), -x),
            (False,X is ca.MX,lambda x,y,z, X,Y :(-1 * x), -x),
            (False,True,lambda x,y,z, X,Y :(x * (1/y)), x/y),
            (False,X is ca.MX,lambda x,y,z, X,Y :((1/x) * y), y/x),
            (False,True,lambda x,y,z, X,Y :(5 * (0.2 * x)), x),
            (False,True,lambda x,y,z, X,Y :(0.5 * (2 * x)), x),
            (False,True,lambda x,y,z, X,Y :(5 * (x / 5)), x),
            (False,True,lambda x,y,z, X,Y :((2 / x) * x), 2),
            (False,True,lambda x,y,z, X,Y :(2*x)/x,2),
            (False,True,lambda x,y,z, X,Y :(5*x)/x,5),
            (False,True,lambda x,y,z, X,Y :((x) * (2 / x)), 2),
            (False,X is ca.MX,lambda x,y,z, X,Y :((-x) * y), -(x*y)),
            (False,True,lambda x,y,z, X,Y :(x * (-y)), -(x*y)),
            (False,X is ca.MX,lambda x,y,z, X,Y :(2*(0.5*x)),x),
            (False,X is ca.MX,lambda x,y,z, X,Y :2*(x*0.5),x),
            (False,X is ca.MX,lambda x,y,z, X,Y :(2*(0.5*X)),X),
            (False,X is ca.MX,lambda x,y,z, X,Y :2*(X*0.5),X),
            (False,True,lambda x,y,z, X,Y :(0 / x), 0),
            (False,True,lambda x,y,z, X,Y :(x / 1), x),
            (False,True,lambda x,y,z, X,Y :(x / -1), -x),
            (False,True,lambda x,y,z, X,Y :(x / x), 1),
            (False,True,lambda x,y,z, X,Y :(2*x / 2), x),
            (False,True,lambda x,y,z, X,Y :((x * y) / x), y),
            (False,X is ca.MX,lambda x,y,z, X,Y :((x * y) / y), x),
            (False,X is ca.MX,lambda x,y,z, X,Y :(1 / x), x**-1),
            (False,True,lambda x,y,z, X,Y :(x / (1/y)), x*y),
            (False,True,lambda x,y,z, X,Y :((2*x) / (2*y)), x/y),
            (False,True,lambda x,y,z, X,Y :((x/5)/0.2), x),
            (False,True,lambda x,y,z, X,Y :(x / (2*x)), 1.0/2),
            (False,True,lambda x,y,z, X,Y :((-x) / x), -1),
            (False,True,lambda x,y,z, X,Y :(x / (-x)), -1),
            (False,True,lambda x,y,z, X,Y :((-x)/(-x)), 1),
            (False,True,lambda x,y,z, X,Y :((x / y) / x), 1/y),
            (False,X is ca.MX,lambda x,y,z, X,Y :((-x) / y), -(x / y)),
            (False,True,lambda x,y,z, X,Y :(x / (-y)), -(x / y)),
            (False,True,lambda x,y,z, X,Y :(x ** 0), 1),
            (False,True,lambda x,y,z, X,Y :(x ** 2), x*x),
            (False,True,lambda x,y,z, X,Y :(x ** 3), x*(x*x)),
            (False,True,lambda x,y,z, X,Y :(x ** -3), 1/(x**3)),
            (False,True,lambda x,y,z, X,Y :(x ** 4), (x**2)*(x**2)),
            (False,True,lambda x,y,z, X,Y :(x ** 0.5), x**0.5),
            (False,False,lambda x,y,z, X,Y :(x ** y), x**y),
            (False,True,lambda x,y,z, X,Y :(x**2 >= 0), 1),
            (False,True,lambda x,y,z, X,Y :(2*x**2 >= x**2), 1),
            (False,True,lambda x,y,z, X,Y :(ca.fmin(x, float('inf'))), x),
            (False,True,lambda x,y,z, X,Y :(ca.fmin(float('inf'), x)), x),
            (False,True,lambda x,y,z, X,Y :(ca.fmin(-float('inf'), x)), -float('inf')),
            (False,True,lambda x,y,z, X,Y :(ca.fmin(x, x)), x),
            (False,True,lambda x,y,z, X,Y :(ca.fmax(-float('inf'), x)), x),
            (False,True,lambda x,y,z, X,Y :(ca.fmax(x, -float('inf'))), x),
            (False,True,lambda x,y,z, X,Y :(ca.fmax(x, float('inf'))), float('inf')),
            (False,True,lambda x,y,z, X,Y :(ca.fmax(x, x)), x),
            (False,True,lambda x,y,z, X,Y :x**2 < 0, 0),
            (False,True,lambda x,y,z, X,Y :(x == x), 1),
            (False,True,lambda x,y,z, X,Y :(x != x), 0),
            (False,True,lambda x,y,z, X,Y :ca.sqrt(x)**2,x),
            (False,True,lambda x,y,z, X,Y :ca.sqrt(X)**2,X),
            (False,True,lambda x,y,z, X,Y :(-x)**2,x**2),
            (False,X is ca.MX,lambda x,y,z, X,Y :ca.fabs(ca.fabs(x)),ca.fabs(x)),
            (False,X is ca.MX,lambda x,y,z, X,Y :ca.fabs(ca.sqrt(x)),ca.sqrt(x)),
            (False,X is ca.MX,lambda x,y,z, X,Y :ca.fabs(ca.exp(x)),ca.exp(x)),
            (False,X is ca.MX,lambda x,y,z, X,Y :ca.log(ca.exp(x)),x),
            (False,X is ca.MX,lambda x,y,z, X,Y :1/(1/x),x),
            (False,True,lambda x,y,z, X,Y :ca.sqrt(x**2),ca.fabs(x)),
            (False,X is ca.MX,lambda x,y,z, X,Y :ca.fabs(-x),ca.fabs(x)),
            (False,True,lambda x,y,z, X,Y :ca.fabs(x)**2,x**2),
            (False,X is ca.MX,lambda x,y,z, X,Y :ca.cos(-x),ca.cos(x)),
            (False,X is ca.MX,lambda x,y,z, X,Y :ca.cos(ca.fabs(x)),ca.cos(x)),
            (False,False,lambda x,y,z, X,Y :-(-x),x),
            (False,True,lambda x,y,z, X,Y :ca.cosh(x*0),x*0+1),
            (False,True,lambda x,y,z, X,Y :x/0.5,2*x),
            (False,True,lambda x,y,z, X,Y :x+x,2*x),
            (False,True,lambda x,y,z, X,Y :x<x,0),
            (False,True,lambda x,y,z, X,Y :2*x+x,3*x),
            (False,True,lambda x,y,z, X,Y :x+2*x,3*x),
            (False,X is ca.MX,lambda x,y,z, X,Y :2*x-x,x),
            (False,True,lambda x,y,z, X,Y :x-2*x,-x),
            (True,False,lambda x,y,z, X,Y :-(x-y),y-x),
               (True,False,lambda x,y,z, X,Y :-(x-y)+(x-y)**2, (-dx)+dx**2),
            (True,False,lambda x,y,z, X,Y :(x+y/z)*z,x*z + y),
            (True,False,lambda x,y,z, X,Y :(x/z+y)*z,x + y*z),
            (True,False,lambda x,y,z, X,Y :z*(x/z+y),x + z*y),
            (True,False,lambda x,y,z, X,Y :z*(x+y/z),z*x + y),
            (True,False,lambda x,y,z, X,Y :(x-y/z)*z,x*z - y),
            (True,False,lambda x,y,z, X,Y :(x/z-y)*z,x - y*z),
            (True,False,lambda x,y,z, X,Y :z*(x/z-y),x - z*y),
            (True,False,lambda x,y,z, X,Y :z*(x-y/z),z*x - y),
            ]:
          print(count)
          count =count+1
          A = genA(x,y,z,X,Y)
          f = ca.Function('f',args,[A])
          if refcount:
            self.assertNotEqual(str(A),str(B))
            print(A)
            f = f.transform()
            fref = ca.Function('f',args,[B])
            self.assertEqual(str(f.call(args,True,False)),str(fref.call(args,True,False)))
          else:
            self.assertEqual(str(A),str(B))
          if on_the_fly:
              ca.GlobalOptions.setSimplificationOnTheFly(False)
              A = genA(x,y,z,X,Y)
              ca.GlobalOptions.setSimplificationOnTheFly(True)
              self.assertNotEqual(str(A),str(B))

          Aval = f(a0,b0,c0,A0,B0)
          Aref = genA(x0,y0,z0,X0,Y0)
          print(Aval)
          print(Aref)
          self.checkarray(Aval,Aref)
          
        
        print(ca.DM(ca.Sparsity.upper(3),0).is_nonnegative())
        print(ca.DM(ca.Sparsity.upper(3),0.5).is_half())
        print(ca.MX(ca.DM(ca.Sparsity.upper(3),0.5)).is_half())
        print(ca.MX(ca.DM(ca.Sparsity.dense(3,3),0.5)).is_half())
        2*(ca.blockcat([[2,3],[6,8]])*A)


  
    
  @known_bug()
  def test_vertcat_empty(self):
    a = ca.MX(ca.DM(0,2))
    v = ca.vertcat(*[a,a])

    self.assertEqual(v.size1(),0)
    self.assertEqual(v.size2(),2)

    a = ca.MX(ca.DM(2,0))
    v = ca.vertcat(*[a,a])

    self.assertEqual(v.size1(),4)
    self.assertEqual(v.size2(),0)

  def test_jacobian_empty(self):
    x = ca.MX.sym("x",3)

    s = ca.jacobian(ca.DM(0,0),x).shape
    self.assertEqual(s[0],0)
    self.assertEqual(s[1],3)

    s = ca.jacobian(x,ca.MX.sym("x",0,4)).shape
    self.assertEqual(s[0],3)
    self.assertEqual(s[1],0)

  def test_mul_sparsity(self):

    N = 10
    x = ca.MX.sym("x",N,N)
    y = ca.MX.sym("y",N,N)

    x_ = self.randDM(N,N)
    y_ = self.randDM(N,N)

    filt = ca.Sparsity.diag(N)+ca.Sparsity.triplet(N,N,[1],[3])

    f = ca.Function("f", [x,y],[x @ y])
    f_in = (x_, y_)
    g = ca.Function("g", [x,y],[ca.mac(x,y,ca.MX.zeros(filt))])
    g_in = (x_, y_)

    f_out = f(*f_in)
    g_out = g(*g_in)

    self.checkarray(ca.DM.ones(filt),ca.DM.ones(g.sparsity_out(0)))

    self.checkarray(f_out[filt],g_out)

  def test_mul_zero_wrong(self):
    with self.assertRaises(RuntimeError):
      ca.MX.sym("X",4,5) @ ca.MX.zeros(3,2)

  def test_vertsplit(self):
    a = ca.MX.sym("X",ca.Sparsity.lower(5))
    v = ca.vertsplit(a,[0,2,4,5])

    Nr = int(5*6/2)

    f = ca.Function("f", [a],v)
    f_in = [0]*f.n_in()  # type: list

    f_in[0]=ca.DM(f.sparsity_in(0),list(range(Nr)))

    f_out = f.call(f_in)
    v = [f_out[i] for i in range(len(v))]

    self.assertEqual(len(v),3)
    self.checkarray(v[0],ca.DM([[0,0,0,0,0],[1,5,0,0,0]]))
    self.checkarray(v[1],ca.DM([[2,6,9,0,0],[3,7,10,12,0]]))
    self.checkarray(v[2],ca.DM([[4,8,11,13,14]]))

    v = ca.vertsplit(a)

    f = ca.Function("f", [a],v)
    f_in = [0]*f.n_in()  # type: list

    f_in[0]=ca.DM(f.sparsity_in(0),list(range(Nr)))

    f_out = f.call(f_in)
    v = [f_out[i] for i in range(len(v))]

    self.assertEqual(len(v),a.size1())
    self.checkarray(v[0],ca.DM([[0,0,0,0,0]]))
    self.checkarray(v[1],ca.DM([[1,5,0,0,0]]))
    self.checkarray(v[2],ca.DM([[2,6,9,0,0]]))
    self.checkarray(v[3],ca.DM([[3,7,10,12,0]]))
    self.checkarray(v[4],ca.DM([[4,8,11,13,14]]))

    v = ca.vertsplit(a,2)

    f = ca.Function("f", [a],v)
    f_in = [0]*f.n_in()  # type: list

    f_in[0]=ca.DM(f.sparsity_in(0),list(range(Nr)))

    f_out = f.call(f_in)
    v = [f_out[i] for i in range(len(v))]

    self.assertEqual(len(v),3)
    self.checkarray(v[0],ca.DM([[0,0,0,0,0],[1,5,0,0,0]]))
    self.checkarray(v[1],ca.DM([[2,6,9,0,0],[3,7,10,12,0]]))
    self.checkarray(v[2],ca.DM([[4,8,11,13,14]]))

    v = ca.vertsplit(a,[0,0,3,a.size1()])

    f = ca.Function("f", [a],v)
    f_in = [0]*f.n_in()  # type: list

    f_in[0]=ca.DM(f.sparsity_in(0),list(range(Nr)))

    f_out = f(*f_in)
    V = [f_out[i] for i in range(len(v))]

    self.assertEqual(len(v),3)
    self.assertEqual(v[0].size1(),0)
    self.assertEqual(v[0].size2(),5)  # why not 5?
    self.checkarray(V[1],ca.DM([[0,0,0,0,0],[1,5,0,0,0],[2,6,9,0,0]]))
    self.checkarray(V[2],ca.DM([[3,7,10,12,0],[4,8,11,13,14]]))

  def test_horzsplit(self):
    a = ca.MX.sym("X",ca.Sparsity.lower(5))
    v = ca.horzsplit(a,[0,2,4,5])

    Nr = int(5*6/2)
    f = ca.Function("f", [a],v)
    f_in = ca.DM(f.sparsity_in(0),list(range(Nr)))

    f_out = f(f_in)
    v = [f_out[i] for i in range(len(v))]
    self.assertEqual(len(v),3)
    self.checkarray(v[0],ca.DM([[0,0],[1,5],[2,6],[3,7],[4,8]]))
    self.checkarray(v[1],ca.DM([[0,0],[0,0],[9,0],[10,12],[11,13]]))
    self.checkarray(v[2],ca.DM([[0],[0],[0],[0],[14]]))

    v = ca.horzsplit(a)

    f = ca.Function("f", [a],v)
    f_in = ca.DM(f.sparsity_in(0),list(range(Nr)))
    f_out = f(f_in)
    v = [f_out[i] for i in range(len(v))]
    self.assertEqual(len(v),a.size1())
    self.checkarray(v[0],ca.DM([0,1,2,3,4]))
    self.checkarray(v[1],ca.DM([0,5,6,7,8]))
    self.checkarray(v[2],ca.DM([0,0,9,10,11]))
    self.checkarray(v[3],ca.DM([0,0,0,12,13]))
    self.checkarray(v[4],ca.DM([0,0,0,0,14]))

    v = ca.horzsplit(a,2)

    f = ca.Function("f", [a],v)
    f_in = ca.DM(f.sparsity_in(0),list(range(Nr)))
    f_out = f(f_in)
    v = [f_out[i] for i in range(len(v))]

    self.assertEqual(len(v),3)
    self.checkarray(v[0],ca.DM([[0,0],[1,5],[2,6],[3,7],[4,8]]))
    self.checkarray(v[1],ca.DM([[0,0],[0,0],[9,0],[10,12],[11,13]]))
    self.checkarray(v[2],ca.DM([[0],[0],[0],[0],[14]]))

    v = ca.horzsplit(a,[0,0,3,a.size2()])
    f = ca.Function("f", [a],v)
    f_in = ca.DM(f.sparsity_in(0),list(range(Nr)))
    f_out = f(f_in)
    V = [f_out[i] for i in range(len(v))]

    self.assertEqual(len(v),3)
    self.assertEqual(v[0].size1(),5)
    self.assertEqual(v[0].size2(),0)
    self.checkarray(V[1],ca.DM([[0,0,0],[1,5,0],[2,6,9],[3,7,10],[4,8,11]]))
    self.checkarray(V[2],ca.DM([[0,0],[0,0],[0,0],[12,0],[13,14]]))

  def test_blocksplit(self):
    a = ca.MX.sym("X",ca.Sparsity.lower(5))
    v = ca.blocksplit(a,[0,2,4,5],[0,1,3,5])

    Nr = int(5*6/2)
    fs = [ca.Function("fs", [a],vr) for vr in v]
    v = [fs[i](ca.DM(fs[i].sparsity_in(0),list(range(Nr)))) for i in range(3)]

    self.checkarray(v[0][0],ca.DM([0,1]))
    self.checkarray(v[0][1],ca.DM([[0,0],[5,0]]))
    self.checkarray(v[1][0],ca.DM([2,3]))
    self.checkarray(ca.blockcat(v),ca.DM(fs[0].sparsity_in(0),list(range(Nr))))  # pyright: ignore[reportCallIssue,reportArgumentType]
    
    
  def test_is_one_etc(self):
    for op in [lambda x: x, ca.MX, ca.SX]:
        self.assertTrue(op(ca.DM(ca.Sparsity.upper(3),0)).is_zero())
        self.assertTrue(op(ca.DM(ca.Sparsity.dense(3,4),0)).is_zero())
        self.assertTrue(op(ca.DM(3,4)).is_zero())
        self.assertFalse(op(ca.blockcat([[0,0],[0,7]])).is_zero())
        
        self.assertFalse(op(ca.DM(ca.Sparsity.upper(3),1)).is_one())
        self.assertTrue(op(ca.DM(ca.Sparsity.dense(3,4),1)).is_one())
        self.assertFalse(op(ca.blockcat([[1,1],[1,7]])).is_one())
       
        self.assertFalse(op(ca.DM(ca.Sparsity.upper(3),0.5)).is_half())
        self.assertTrue(op(ca.DM(ca.Sparsity.dense(3,4),0.5)).is_half())
        self.assertFalse(op(ca.blockcat([[0.5,0.5],[0.5,7]])).is_half())
         
        self.assertFalse(op(ca.DM(ca.Sparsity.upper(3),-1)).is_minus_one())
        self.assertTrue(op(ca.DM(ca.Sparsity.dense(3,4),-1)).is_minus_one())
        self.assertFalse(op(ca.blockcat([[-1,-1],[-1,7]])).is_minus_one())

        self.assertFalse(op(ca.DM(ca.Sparsity.upper(3),np.inf)).is_inf())
        self.assertTrue(op(ca.DM(ca.Sparsity.dense(3,4),np.inf)).is_inf())
        self.assertFalse(op(ca.blockcat([[np.inf,np.inf],[np.inf,7]])).is_inf())

        self.assertFalse(op(ca.DM(ca.Sparsity.upper(3),-np.inf)).is_minus_inf())
        self.assertTrue(op(ca.DM(ca.Sparsity.dense(3,4),-np.inf)).is_minus_inf())
        self.assertFalse(op(ca.blockcat([[-np.inf,-np.inf],[-np.inf,7]])).is_minus_inf())

        self.assertTrue(op(ca.DM(ca.Sparsity.upper(3),3)).is_integer())
        self.assertTrue(op(ca.DM(ca.Sparsity.dense(3,4),3)).is_integer())
        self.assertFalse(op(ca.blockcat([[3,3],[3,0.7]])).is_integer())
        self.assertTrue(op(ca.blockcat([[3,3],[3,7]])).is_integer())

        self.assertFalse(op(ca.DM(ca.Sparsity.upper(3),7)).is_value(7))
        self.assertTrue(op(ca.DM(ca.Sparsity.dense(3,4),7)).is_value(7))
        self.assertFalse(op(ca.blockcat([[7,7],[7,1]])).is_value(7))
        self.assertTrue(op(ca.DM(ca.Sparsity.upper(3),0)).is_value(0))
        self.assertTrue(op(ca.DM(ca.Sparsity.dense(3,4),0)).is_value(0))
        self.assertTrue(op(ca.DM(3,4)).is_value(0))
        self.assertFalse(op(ca.blockcat([[0,0],[0,7]])).is_value(0))

        self.assertFalse(op(ca.DM(ca.Sparsity.upper(3),1)).is_eye())
        self.assertTrue(op(ca.DM(ca.Sparsity.diag(3),1)).is_eye())
        self.assertFalse(op(ca.blockcat([[1,0],[0,1]])).is_eye()) # Why?
        
        
        self.assertFalse(op(ca.DM(ca.Sparsity.upper(3),-1)).is_nonnegative())
        self.assertTrue(op(ca.DM(ca.Sparsity.upper(3),1)).is_nonnegative())
        self.assertTrue(op(ca.DM(ca.Sparsity.dense(3,4),1)).is_nonnegative())
        self.assertTrue(op(ca.blockcat([[1,1],[1,7]])).is_nonnegative())
        self.assertFalse(op(ca.blockcat([[1,1],[1,-7]])).is_nonnegative())
        
  def test_mxnulloutput(self):
     a = ca.MX(5,0)
     b = ca.MX.sym("x",2)

     f = ca.Function("f", [b],[a])
     c = f(b)

     self.assertEqual(c.size1(),5)
     self.assertEqual(c.size2(),0)

     c = f.call([b],True)[0]

     self.assertEqual(c.size1(),5)
     self.assertEqual(c.size2(),0)

     a = ca.MX(0,0)
     b = ca.MX.sym("x",2)

     f = ca.Function("f", [b],[a])
     c = f(b)

     self.assertEqual(c.size1(),0)
     self.assertEqual(c.size2(),0)

     c = f.call([b],True)[0]

     self.assertEqual(c.size1(),0)
     self.assertEqual(c.size2(),0)

  def test_mxnull(self):
     a = ca.MX(5,0)
     b = ca.MX(0,3)

     c = a @ b

     self.assertEqual(c.nnz(),0)

     a = ca.MX(5,3)
     b = ca.MX(3,4)

     c = a @ b

     self.assertEqual(c.nnz(),0)

  def  test_mxnullop(self):
    c = ca.MX(0,0)
    x = ca.MX.sym("x",2,3)

    # https://github.com/casadi/casadi/issues/2628
    if swig4:
      with self.assertRaises(TypeError):
        d = x + c
    else:
      with self.assertRaises(RuntimeError):
        d = x + c

    if swig4:
      with self.assertRaises(TypeError):
        d = x / c
    else:
      with self.assertRaises(RuntimeError):
        d = x / c

  @slow()
  @memory_heavy()
  def test_MX_shapes(self):
      self.message("MX unary operations")

      #self.checkarray(DM(Sparsity.lower(4),1),DM(Sparsity.dense(4,4),1))

      for sp in [ca.Sparsity.dense(0,0),ca.Sparsity.dense(0,2),ca.Sparsity.dense(2,0),ca.Sparsity.dense(1,1),ca.Sparsity.dense(2,2), ca.Sparsity(4,3,[0,2,2,3],[1,2,1])]:
        for v in [0,1,0.2]:
          x_ = ca.DM(sp,v)

          xx = ca.MX.sym("x",sp.size1(),sp.size2())
          x=xx[sp]

          for (casadiop, numpyop,name, flags) in self.pool.zip():
            if 'nozero' in flags and v==0: continue
            r = casadiop([x])
            f = ca.Function("f", [xx],[r])
            rr = f(v)
            print(r,rr,numpyop(x_))
            self.checkarray(rr,numpyop(x_))

            a = ca.DM(f.sparsity_out(0),1)
            b = ca.DM.ones(ca.DM(numpyop(x_)).sparsity())

            c = b-a
            if c.nnz()>0:
              # At least as sparse as DM calculus
              self.assertTrue(min(c.nonzeros())>=0,str([sp,v,name]))

      for sp in [ca.Sparsity(1,1),ca.Sparsity.dense(1,1),ca.Sparsity(3,4),ca.Sparsity.dense(3,4), ca.Sparsity(4,3,[0,2,2,3],[1,2,1]).T]:
        for v1 in [0,1,0.2,-0.2]:
          x1_ = ca.DM(sp,v1)
          xx1 = ca.MX.sym("x",sp.size1(),sp.size2())
          x1=xx1[sp]
          xx1s = ca.SX.sym("x",sp.size1(),sp.size2())
          x1s=xx1s[sp]
          for sp2 in [ca.Sparsity(1,1),ca.Sparsity.dense(1,1),ca.Sparsity(3,4),ca.Sparsity.dense(3,4), ca.Sparsity(4,3,[0,2,2,3],[1,2,1]).T]:
            for v2 in [0,1,0.2,-0.2]:
              x2_ = ca.DM(sp2,v2)
              xx2 = ca.MX.sym("x",sp2.size1(),sp2.size2())
              x2=xx2[sp2]
              xx2s = ca.SX.sym("x",sp2.size1(),sp2.size2())
              x2s=xx2s[sp2]
              for (casadiop, numpyop,name, flags) in self.matrixbinarypool.zip():
                if ("mul" in name or "mtimes" in name) and (sp.numel()==1 or sp2.numel()==1): continue
                r = casadiop([x1,x2])
                f = ca.Function("f", [xx1,xx2],[r])
                f_in = [v1, v2]
                f_out = f(*f_in)
                g = ca.Function("g", [xx1,xx2],[r])
                g_out = g(v1, v2)

                self.checkarray(f_out,numpyop([x1_,x2_]),str([sp,sp2,v1,v2,x1_,x2_,name]))


                if "mul" not in name:
                  a = ca.DM.ones(f.sparsity_out(0))
                  b = ca.DM.ones(g.sparsity_out(0))

                  c = b-a
                  if c.nnz()>0:
                    # At least as sparse as DM calculus
                    self.assertTrue(min(c.nonzeros())>=0,str([sp,sp2,v1,v2,name,a,b]))

                if sp.nnz()>0 and sp2.nnz()>0 and v1!=0 and v2!=0:
                  self.checkfunction(f,g,inputs=f_in,hessian=False,failmessage=str([sp,sp2,v1,v2,x1_,x2_,name]))

  @memory_heavy()
  def test_MXConstant(self):
      self.message("MX unary operations, constant")

      #self.checkarray(DM(Sparsity.lower(4),1),DM(Sparsity.dense(4,4),1))

      for sp in [ca.Sparsity.dense(0,0),ca.Sparsity.dense(0,2),ca.Sparsity.dense(2,0),ca.Sparsity.dense(1,1),ca.Sparsity.dense(2,2), ca.Sparsity(4,3,[0,2,2,3],[1,2,1])]:
        for v in [0,1,0.2]:
          x_ = ca.DM(sp,v)

          x=ca.MX(sp,v)

          for (casadiop, numpyop,name, flags) in self.pool.zip():
            if 'nozero' in flags and (v==0 or not sp.is_dense()): continue
            r = casadiop([x])
            print(r)
            self.assertTrue(r.is_constant())

            self.checkarray(r.to_DM(),numpyop(x_),str([x_,name]))

            a = ca.DM.ones(r.to_DM().sparsity())
            b = ca.DM.ones(ca.DM(numpyop(x_)).sparsity())

            c = b-a
            if c.nnz()>0:
              # At least as sparse as DM calculus
              self.assertTrue(min(c.nonzeros())>=0,str([sp,v,name]))

      for sp in [ca.Sparsity.dense(1,1),ca.Sparsity(1,1),ca.Sparsity(3,4),ca.Sparsity.dense(3,4), ca.Sparsity(4,3,[0,2,2,3],[1,2,1]).T]:
        for v1 in [0,1,0.2,-0.2]:
          x1_ = ca.DM(sp,v1)
          x1=ca.MX(sp,v1)
          for sp2 in [ca.Sparsity.dense(1,1),ca.Sparsity(1,1),ca.Sparsity(3,4),ca.Sparsity.dense(3,4), ca.Sparsity(4,3,[0,2,2,3],[1,2,1]).T]:
            for v2 in [0,1,0.2,-0.2]:
              x2_ = ca.DM(sp2,v2)
              x2=ca.MX(sp2,v2)
              for (casadiop, numpyop,name, flags) in self.matrixbinarypool.zip():
                if ("mul" in name or "mtimes" in name) and (sp.numel()==1 or sp2.numel()==1): continue
                r = casadiop([x1,x2])
                f = ca.Function("f", [],[r]) # Should not be needed -> constant folding
                f_out = f.call([])


                self.checkarray(f_out[0],numpyop([x1_,x2_]),str([sp,sp2,v1,v2,name]))
                if "mul" not in name and "mtimes" not in name:
                  a = ca.DM.ones(f.sparsity_out(0))
                  b = ca.DM.ones(ca.DM(numpyop([x1_,x2_])).sparsity())

                  c = b-a
                  if c.nnz()>0:
                    # At least as sparse as DM calculus
                    self.assertTrue(min(c.nonzeros())>=0,str([sp,sp2,v1,v2,name]))

  def test_graph_substitute(self):
    x=ca.MX.sym("X",4,4)
    y=ca.MX.sym("Y",4,4)
    b=ca.MX.sym("B",4,4)

    c = x*y
    d = b*y
    f = c+d


    C = ca.MX.sym("C",4,4)
    f = ca.graph_substitute(f,[c],[C])

    F = ca.Function("F", [y,b,C],[f])
    F_out = F(1, 2, 3)

    self.checkarray(F_out,5*ca.DM.ones(4,4))

    D = ca.MX.sym("D",4,4)
    f = ca.graph_substitute(f,[d],[D])

    F = ca.Function("F", [D,C],[f])
    F_out = F(4, 5)

    self.checkarray(F_out,9*ca.DM.ones(4,4))
    
    z = ca.vertcat(c,d,c+d)
    
    f = ca.graph_substitute(z,[d],[D])
    
    F = ca.Function("F", [x,y,b,D],[f])
    F_ref = ca.Function("F", [x,y,b,D],[ca.vertcat(c,D,c+D)])
    self.checkfunction_light(F,F_ref,inputs=[3,5,7,12])

  def test_graph_substitute_function(self):
    x=ca.MX.sym("x")
    y=ca.MX.sym("y")


        
    f = ca.Function("forig",[x,y],[x*y],{"never_inline":True})
    g = ca.Function("fsubs",[x,y],[x*y],{"never_inline":True})
    
    sx = ca.sin(x)
    y2 = 2*y
    
    a = f(sx,y2)
    
    x2 = x**2
    
    z = ca.sqrt(x2*a)+a**2
    print(z)
    z2 = ca.graph_substitute(z,[x2],[x])
    print(z2)
        
    f = ca.Function("forig",[x,y],[x*y,x-y,x+y],{"never_inline":True})
    g = ca.Function("fsubs",[x,y],[x*y,x-y,x+y],{"never_inline":True})
    
    sx = ca.sin(x)
    y2 = 2*y
    
    [a,b,c] = f(sx,y2)
    
    z = ca.sqrt(a)+c**2
    print(z)
    z2 = ca.graph_substitute(z,[],[])
    print(z2)
    
    # Substitute a single output of a call
    z2 = ca.graph_substitute(z,[a],[g(sx,y2)[0]])
    
    print(z2)
    self.assertTrue("fsubs" in str(z2))
    self.assertTrue("forig" in str(z2))
    
    F_ref = ca.Function("F_ref",[x,y],[z])
    F = ca.Function("F",[x,y],[z2])
    
    self.checkfunction_light(F,F_ref,inputs=[3.1,7.1])

    # Substitute a single output of a call
    z2 = ca.graph_substitute(z2,[b],[np.nan])
    self.assertTrue("fsubs" in str(z2))
    self.assertTrue("forig" in str(z2))
    
    F_ref = ca.Function("F_ref",[x,y],[z])
    F = ca.Function("F",[x,y],[z2])
    
    self.checkfunction_light(F,F_ref,inputs=[3.1,7.1])

    # Substitute a single output of a call
    z2 = ca.graph_substitute(z2,[c],[g(sx,y2)[2]])
    self.assertTrue("fsubs" in str(z2))
    self.assertFalse("forig" in str(z2))
    
    F_ref = ca.Function("F_ref",[x,y],[z])
    F = ca.Function("F",[x,y],[z2])
    
    self.checkfunction_light(F,F_ref,inputs=[3.1,7.1])
    
    # Substitute the entire Multi-output call
    
    z2 = ca.graph_substitute(z,[a.dep(0)],[g(sx,y2)[0].dep(0)])
    
    self.assertTrue("fsubs" in str(z2))
    self.assertFalse("forig" in str(z2))

    F_ref = ca.Function("F_ref",[x,y],[z])
    F = ca.Function("F",[x,y],[z2])
    
    self.checkfunction_light(F,F_ref,inputs=[3.1,7.1])
    
    
  def test_graph_substitute_call(self):
    x = ca.MX.sym("x")
    y = ca.MX.sym("y")
    z = ca.MX.sym("z")

    f = ca.Function("f",[x,y,z],[x**2])

    ca.vertsplit

    sx = ca.sin(x)
    res = ca.graph_substitute([f(x+1,y,z),ca.sin(f(x+1,y,z)),f(sx,y,z),ca.sin(f(sx,y,z))],[sx],[ca.cos(x)])
    res = str(res)
    self.assertEqual(res.count("{0}"),4)

  def test_matrix_expand(self):
    n = 2
    a = ca.MX.sym("a",n,n)
    b = ca.MX.sym("b",n,n)
    c = ca.MX.sym("c",n,n)

    d = a+b
    e = d*c

    self.assertEqual(ca.n_nodes(e),6)

    t0 = ca.matrix_expand(e)

    self.assertEqual(ca.n_nodes(t0),5)

    t1 = ca.matrix_expand(e,[d])
    self.assertEqual(ca.n_nodes(t1),6)

    print(e,t0,t1)


    outs = []
    for x in [e,t0,t1]:
      f = ca.Function("f", [a,b,c],[x])

      f_in = [1.1, 2.2, 3.3]
      f_out = f(*f_in)

      outs.append(f_out)
      if len(outs)>1:
        self.checkarray(outs[0],outs[-1])

    print(outs)

  def test_kron(self):
    a = ca.sparsify(ca.DM([[1,0,6],[2,7,0]]))
    b = ca.sparsify(ca.DM([[1,0,0],[2,3,7],[0,0,9],[1,12,13]]))

    A = ca.MX.sym("A",a.sparsity())
    B = ca.MX.sym("B",b.sparsity())
    C = c.kron(A,B)

    f = ca.Function("f", [A,B],[C])
    f_in = [a, b]
    f_out = f(*f_in)

    c_ = f_out

    self.assertEqual(c_.size1(),a.size1()*b.size1())
    self.assertEqual(c_.size2(),a.size2()*b.size2())
    self.assertEqual(c_.nnz(),a.nnz()*b.nnz())

    self.checkarray(c_,numpy.kron(a,b))

    # Kron is its own MXNode: the graph is O(1) regardless of a's size.
    self.assertEqual(f.n_nodes(), 4)

    # Forward-AD rule kron(dA,B) + kron(A,dB) must hold exactly.
    dA = ca.MX.sym("dA", A.sparsity())
    dB = ca.MX.sym("dB", B.sparsity())
    fwd = ca.jtimes(C, ca.vertcat(ca.vec(A), ca.vec(B)), ca.vertcat(ca.vec(dA), ca.vec(dB)))
    expected = c.kron(dA, B) + c.kron(A, dB)
    g = ca.Function("g", [A,B,dA,dB], [fwd - expected])
    ca.DM.rng(0)
    diff = g(a, b, ca.DM.rand(A.sparsity()), ca.DM.rand(B.sparsity()))
    self.assertTrue(float(ca.norm_inf(diff)) < 1e-14)

    # Codegen + serialization round-trips
    self.check_codegen(f, inputs=f_in)
    self.check_serialize(f, inputs=f_in)
    self.checkfunction(f, f.expand(), inputs=f_in)
    
    A = ca.SX.sym("A",a.sparsity())
    B = ca.SX.sym("B",b.sparsity())
    C = c.kron(A,B)

    fSX = ca.Function("f", [A,B],[C])
    self.checkfunction(f, fSX, inputs=f_in)
  def test_kron_flavors(self):
    # Kron has 4 dispatch paths based on operand sparsity (DenseKron,
    # DenseSparseKron, SparseDenseKron, and the general sparse-sparse base).
    # Verify each path: eval against numpy.kron, codegen reload, serialize
    # round-trip, n_nodes stays O(1).
    import numpy
    ca.DM.rng(11)
    sp_d2 = ca.Sparsity.dense(2, 3)
    sp_d3 = ca.Sparsity.dense(3, 2)
    sp_s2 = ca.sparsify(ca.DM([[1, 0, 1], [0, 1, 0]])).sparsity()
    sp_s3 = ca.sparsify(ca.DM([[1, 0], [1, 1], [0, 1]])).sparsity()
    cases = [
        ("dense+dense",  sp_d2, sp_d3),
        ("dense+sparse", sp_d2, sp_s3),
        ("sparse+dense", sp_s2, sp_d3),
        ("sparse+sparse",sp_s2, sp_s3),
    ]
    for label, sp_a, sp_b in cases:
      A = ca.MX.sym("A", sp_a)
      B = ca.MX.sym("B", sp_b)
      K = c.kron(A, B)
      f = ca.Function("f", [A, B], [K])
      self.assertEqual(f.n_nodes(), 4, label + ": graph not collapsed")
      a_val = ca.DM.rand(sp_a)
      b_val = ca.DM.rand(sp_b)
      ref = numpy.kron(numpy.array(a_val), numpy.array(b_val))
      self.checkarray(f(a_val, b_val), ref, label)
      self.check_codegen(f, inputs=[a_val, b_val])
      self.check_serialize(f, inputs=[a_val, b_val])
      self.checkfunction(f, f.expand(), inputs=[a_val, b_val])

  @memory_heavy()
  def test_kron_contract_flavors(self):
    # KronContract has 4 dispatch paths (DenseKronContract,
    # DenseSparseKronContract, SparseDenseKronContract, and base). For each
    # combination of M-sparsity x X-sparsity x {inner, outer}: verify eval
    # against a dense numpy reference, codegen round-trip, serialize.
    import numpy
    ca.DM.rng(13)
    def ref_kron_contract(m_np, x_np, inner):
      mA_mB, nA_nB = m_np.shape
      xrow, xcol = x_np.shape
      if inner:
        mB, nB = xrow, xcol
        mA, nA = mA_mB // mB, nA_nB // nB
        y = numpy.zeros((mA, nA))
        for i in range(mA):
          for j in range(nA):
            for r in range(mB):
              for s in range(nB):
                y[i, j] += m_np[i*mB+r, j*nB+s] * x_np[r, s]
      else:
        mA, nA = xrow, xcol
        mB, nB = mA_mB // mA, nA_nB // nA
        y = numpy.zeros((mB, nB))
        for r in range(mB):
          for s in range(nB):
            for i in range(mA):
              for j in range(nA):
                y[r, s] += x_np[i, j] * m_np[i*mB+r, j*nB+s]
      return y

    sp_d2 = ca.Sparsity.dense(2, 3)
    sp_d3 = ca.Sparsity.dense(3, 2)
    sp_s2 = ca.sparsify(ca.DM([[1, 0, 1], [0, 1, 0]])).sparsity()
    sp_s3 = ca.sparsify(ca.DM([[1, 0], [1, 1], [0, 1]])).sparsity()
    combos = [
      ("dense_M+dense_X",   sp_d2, sp_d3),
      ("dense_M+sparse_X",  sp_d2, sp_s3),
      ("sparse_M+dense_X",  sp_s2, sp_d3),
      ("sparse_M+sparse_X", sp_s2, sp_s3),
    ]
    for label, sp_a, sp_b in combos:
      sp_m = c.kron(sp_a, sp_b)
      # inner mode: X has sp_b, output has sp_a (approximately)
      M_inner = ca.MX.sym("M", sp_m)
      X_inner = ca.MX.sym("X", sp_b)
      Yi = c.kron_contract(M_inner, X_inner, True)
      fi = ca.Function("fi", [M_inner, X_inner], [Yi])
      mv = ca.DM.rand(sp_m); xv = ca.DM.rand(sp_b)
      ref = ref_kron_contract(numpy.array(mv), numpy.array(xv), True)
      # Y may be sparse if both M and X are sparse; compare projected dense.
      self.checkarray(ca.DM(fi(mv, xv)), ca.DM(ref) * ca.DM.ones(Yi.sparsity()),
                      label + " inner")
      self.check_codegen(fi, inputs=[mv, xv])
      self.check_serialize(fi, inputs=[mv, xv])
      self.checkfunction(fi, fi.expand(), inputs=[mv, xv])
      # outer mode: X has sp_a, output has sp_b (approximately)
      M_outer = ca.MX.sym("M", sp_m)
      X_outer = ca.MX.sym("X", sp_a)
      Yo = c.kron_contract(M_outer, X_outer, False)
      fo = ca.Function("fo", [M_outer, X_outer], [Yo])
      mv2 = ca.DM.rand(sp_m); xv2 = ca.DM.rand(sp_a)
      ref2 = ref_kron_contract(numpy.array(mv2), numpy.array(xv2), False)
      self.checkarray(ca.DM(fo(mv2, xv2)), ca.DM(ref2) * ca.DM.ones(Yo.sparsity()),
                      label + " outer")
      self.check_codegen(fo, inputs=[mv2, xv2])
      self.check_serialize(fo, inputs=[mv2, xv2])
      self.checkfunction(fo, fo.expand(), inputs=[mv2, xv2])
     
  def test_kron_contract_empty_x(self):
    # Regression: DenseSparseKronContract used to write mA*nA zeros into a
    # 0-size output buffer when X was a sparse-with-no-nz pattern, corrupting
    # the heap. Sparsity::kron_contract correctly returns a (mA, nA)-shaped
    # but nnz=0 output for that case; the dense kernel must short-circuit.
    sp_x_empty = ca.Sparsity(2, 2)                  # 2x2 with 0 nz
    sp_m = ca.Sparsity.dense(4, 4)                  # mA=mB=2, nA=nB=2
    for inner in [True, False]:
      M = ca.MX.sym("M", sp_m)
      X = ca.MX.sym("X", sp_x_empty)
      Y = c.kron_contract(M, X, inner)
      self.assertEqual(Y.shape, (2, 2))
      self.assertEqual(Y.nnz(), 0)
      f = ca.Function("f", [M, X], [Y])
      ca.DM.rng(0)
      mv = ca.DM.rand(sp_m); xv = ca.DM(sp_x_empty)
      out = f(mv, xv)
      self.assertEqual(out.shape, (2, 2))
      self.assertEqual(out.nnz(), 0)
      self.check_codegen(f, inputs=[mv, xv])
      self.check_serialize(f, inputs=[mv, xv])

  def test_kron_contract(self):
    ca.DM.rng(7)
    a = ca.sparsify(ca.DM([[1,0,6],[2,7,0]]))                    # 2x3
    b = ca.sparsify(ca.DM([[1,0,0],[2,3,7],[0,0,9],[1,12,13]]))  # 4x3

    A = ca.MX.sym("A", a.sparsity())
    B = ca.MX.sym("B", b.sparsity())
    K = c.kron(A, B)                                       # 8x9
    M = ca.MX.sym("M", K.sparsity())

    # ---- forward eval: kron_contract reproduces the math definition ----
    # inner: Y[i,j] = sum_{r,s} M[i*mB+r, j*nB+s] * B[r,s]
    Yi = c.kron_contract(M, B, True)
    self.assertEqual(Yi.size1(), A.size1())
    self.assertEqual(Yi.size2(), A.size2())
    # outer: Y[r,s] = sum_{i,j} A[i,j] * M[i*mB+r, j*nB+s]
    Yo = c.kron_contract(M, A, False)
    self.assertEqual(Yo.size1(), B.size1())
    self.assertEqual(Yo.size2(), B.size2())

    f_inner = ca.Function("fi", [M, B], [Yi])
    f_outer = ca.Function("fo", [M, A], [Yo])

    # Numeric check vs a dense reference computation
    m_val = ca.DM.rand(K.sparsity())
    m_np = numpy.array(m_val)
    a_np = numpy.array(a); b_np = numpy.array(b)
    mB1, nB1 = b.size1(), b.size2()
    mA1, nA1 = a.size1(), a.size2()
    ref_inner = numpy.zeros((mA1, nA1))
    ref_outer = numpy.zeros((mB1, nB1))
    for i in range(mA1):
      for j in range(nA1):
        ref_inner[i,j] = (m_np[i*mB1:(i+1)*mB1, j*nB1:(j+1)*nB1] * b_np).sum()
    for r in range(mB1):
      for s in range(nB1):
        for i in range(mA1):
          for j in range(nA1):
            ref_outer[r,s] += a_np[i,j] * m_np[i*mB1+r, j*nB1+s]
    self.checkarray(f_inner(m_val, b), ref_inner)
    self.checkarray(f_outer(m_val, a), ref_outer)

    # ---- AD closure: graph stays O(1) regardless of A's size ----
    sizes = [3, 5, 10, 20]
    node_counts = []
    for n in sizes:
      Ax = ca.MX.sym("Ax", n, n)
      Bx = ca.MX.sym("Bx", 3, 3)
      Kx = c.kron(Ax, Bx)
      gA = ca.jacobian(ca.vec(Kx), ca.vec(Ax))
      node_counts.append(ca.Function("f", [Ax, Bx], [gA]).n_nodes())
    self.assertEqual(min(node_counts), max(node_counts),
                     "kron jacobian graph must be scale-invariant: %s" % node_counts)

    # ---- AD closure: hessian (second-order) also O(1) ----
    hess_counts = []
    for n in sizes:
      Ax = ca.MX.sym("Ax", n, n)
      Bx = ca.MX.sym("Bx", 2, 2)
      Kx = c.kron(Ax, Bx)
      obj = 0.5 * ca.sumsqr(Kx)
      x = ca.vertcat(ca.vec(Ax), ca.vec(Bx))
      H, _ = ca.hessian(obj, x)
      hess_counts.append(ca.Function("H", [Ax, Bx], [H]).n_nodes())
    self.assertEqual(min(hess_counts), max(hess_counts),
                     "kron hessian graph must be scale-invariant: %s" % hess_counts)

    # ---- AD correctness: SX-expand path matches MX-graph path ----
    # checkfunction_light compares f with f.expand() over the same inputs --
    # exercises every AD direction implicitly.
    A2 = ca.MX.sym("A2", 3, 4); B2 = ca.MX.sym("B2", 2, 5)
    M2 = ca.MX.sym("M2", c.kron(A2, B2).sparsity())
    for expr, args in [(c.kron_contract(M2, B2, True),  [M2, B2]),
                       (c.kron_contract(M2, A2, False), [M2, A2])]:
      f = ca.Function("f", args, [expr])
      ca.DM.rng(11)
      inputs = [ca.DM.rand(arg.sparsity()) for arg in args]
      self.checkfunction_light(f, f.expand(), inputs=inputs)

    # ---- forward + reverse AD via jtimes self-consistency for inner mode ----
    A3 = ca.MX.sym("A3", 3, 3); B3 = ca.MX.sym("B3", 2, 2)
    M3 = ca.MX.sym("M3", c.kron(A3, B3).sparsity())
    Y = c.kron_contract(M3, B3, True)
    # Closed-form chain rule: d(kron_contract(M, B, inner)) =
    #   kron_contract(dM, B, inner) + kron_contract(M, dB, inner)
    dM = ca.MX.sym("dM", M3.sparsity())
    dB = ca.MX.sym("dB", B3.sparsity())
    fwd = ca.jtimes(Y, ca.vertcat(ca.vec(M3), ca.vec(B3)), ca.vertcat(ca.vec(dM), ca.vec(dB)))
    expected = c.kron_contract(dM, B3, True) + c.kron_contract(M3, dB, True)
    g = ca.Function("g", [M3, B3, dM, dB], [fwd - expected])
    ca.DM.rng(3)
    diff = g(ca.DM.rand(M3.sparsity()), ca.DM.rand(B3.sparsity()),
            ca.DM.rand(M3.sparsity()), ca.DM.rand(B3.sparsity()))
    self.assertTrue(float(ca.norm_inf(diff)) < 1e-14)

    # ---- codegen + serialization round-trips for both modes ----
    Av = ca.DM.rand(A2.sparsity()); Bv = ca.DM.rand(B2.sparsity())
    Mv = ca.DM.rand(M2.sparsity())
    fi = ca.Function("fi", [M2, B2], [c.kron_contract(M2, B2, True)])
    fo = ca.Function("fo", [M2, A2], [c.kron_contract(M2, A2, False)])
    self.check_codegen(fi, inputs=[Mv, Bv])
    self.check_codegen(fo, inputs=[Mv, Av])
    self.check_serialize(fi, inputs=[Mv, Bv])
    self.check_serialize(fo, inputs=[Mv, Av])

  def test_project(self):
    x = ca.MX.sym("x",ca.Sparsity.lower(3))
    y = ca.project(x, ca.Sparsity.lower(3).T)

    f = ca.Function("f", [x],[y])
    f_in=ca.DM(f.sparsity_in(0),list(range(1,int(4*3/2+1))))
    f_out = f(f_in)

    self.checkarray(f_out,ca.DM([[1,0,0],[0,4,0],[0,0,6]]))
    self.checkarray(ca.DM.ones(f.sparsity_out(0)),ca.DM.ones(ca.Sparsity.lower(3).T))

    self.check_codegen(f,inputs=[ca.DM([[1,0,0],[0,4,0],[0,0,6]])])
    self.check_serialize(f,inputs=[ca.DM([[1,0,0],[0,4,0],[0,0,6]])])

  def test_project_simplify(self):
    x = ca.MX.sym("x",ca.Sparsity.lower(3))
    ca.DM.rng(1)
    X0 = ca.DM.rand(x.sparsity())
    
    y = ca.project(x, ca.Sparsity.upper(3))
    z = ca.project(y, ca.sparsify(ca.blockcat([[0,0,0],[1,1,1],[1,1,1]])).sparsity())
    
    f = ca.Function('f',[x],[z])
    fs = f.transform()
    self.assertTrue(fs.n_nodes()>3)

    self.checkfunction(f,fs,inputs=[X0])

    z = ca.project(y, ca.sparsify(ca.blockcat([[1,0,0],[0,0,0],[0,0,1]])).sparsity())
    
    f = ca.Function('f',[x],[z])
    fs = f.transform()

    self.assertTrue(fs.n_nodes(),3)
    
    fs.disp(True)

    self.checkfunction(f,fs,inputs=[X0])
    
  def test_repmat(self):
    a = ca.DM([[1,2],[3,4],[5,6]])
    self.checkarray(ca.repmat(a,2,3),ca.kron(ca.DM.ones(2,3),a))

  def test_transpose_bug(self):

    A = [[-26.9091,00,00,1,00,00,00,00,00,00,00,00,00,00,00],
     [00,-26.9091,00,00,1,00,00,00,00,00,00,00,00,00,00],
     [00,00,-26.9091,00,00,1,00,00,00,00,00,00,00,00,00],
     [00,00,00,-14.1393,-0.348109,-0.382033,3.6491,9.40422,-5.05449,3.03347,-0.00126949,-0.000414104,00,0.00498232,0.00495617],
     [00,00,00,-0.348109,-13.6315,-0.194202,-6.41868,-0.0189287,3.34725,1.86987,-0.000645326,-0.000210504,00,0.0025327,0.0025194],
     [00,00,00,-0.382033,-0.194202,-13.6677,10.1342,-13.0101,-0.890078,-8.99621,-0.000708215,-0.000231018,00,0.00277951,0.00276493],
     [00,00,00,00,00,00,-27.0487,0.36775,-0.277186,0.564995,0.424389,0.0235157,0.227774,00,00],
     [00,00,00,00,00,00,0.113776,-27.3241,0.159254,-0.520988,-0.0235157,0.424389,0.132128,00,00],
     [00,00,00,00,00,00,0.227472,-0.0735535,-26.9135,0.206825,-0.227774,-0.132128,0.424389,00,00],
     [00,00,00,00,00,00,0.332188,-1.02565,-0.0471482,-28.3499,0.132128,-0.227774,0.0235157,00,00],
     [00,00,00,-0.00126949,-0.000645326,-0.000708215,0.00763923,0.0016513,-0.0044564,-0.00284703,-0.119342,-0.000941301,-0.00160961,-1.62792e-05,0.00156863],
     [00,00,00,-0.000414104,-0.000210504,-0.000231018,0.0024919,0.000538651,-0.00145367,-0.000928698,0.000941301,-0.119493,0.000742542,-0.00155303,-7.50989e-06],
     [00,00,00,00,00,00,00,00,00,00,2.1684e-19,-2.1684e-19,-0.171654,-0.000868032,0.000868032],
     [00,00,00,0.00181259,0.000921407,0.0010112,-0.0109074,-0.00235775,0.00636292,0.00406504,0,-0.000567731,-0.000868032,-0.00087529,00],
     [00,00,00,0.00166876,0.000848291,0.000930959,-0.0100419,-0.00217066,0.005858,0.00374247,0.00052268,0,-0.000868032,00,-0.000874062]
     ]

    A = ca.sparsify(ca.DM(A))

    As = ca.MX.sym("As",A.sparsity())

    f = ca.Function("f", [As],[ca.densify(As.T),ca.densify(As).T,As.T,As,ca.densify(As)])

    f_out = f(A)

    self.checkarray(f_out[0],A.T)
    self.checkarray(f_out[1],A.T)
    self.checkarray(f_out[2],A.T)
    self.checkarray(f_out[3],A)
    self.checkarray(f_out[4],A)

  @requiresPlugin(ca.Linsol,"csparse")
  def test_bizarre_bug(self):

    A = [[-26.9091,00,00,1,00,00,00,00,00,00,00,00,00,00,00],
     [00,-26.9091,00,00,1,00,00,00,00,00,00,00,00,00,00],
     [00,00,-26.9091,00,00,1,00,00,00,00,00,00,00,00,00],
     [00,00,00,-14.1393,-0.348109,-0.382033,3.6491,9.40422,-5.05449,3.03347,-0.00126949,-0.000414104,00,0.00498232,0.00495617],
     [00,00,00,-0.348109,-13.6315,-0.194202,-6.41868,-0.0189287,3.34725,1.86987,-0.000645326,-0.000210504,00,0.0025327,0.0025194],
     [00,00,00,-0.382033,-0.194202,-13.6677,10.1342,-13.0101,-0.890078,-8.99621,-0.000708215,-0.000231018,00,0.00277951,0.00276493],
     [00,00,00,00,00,00,-27.0487,0.36775,-0.277186,0.564995,0.424389,0.0235157,0.227774,00,00],
     [00,00,00,00,00,00,0.113776,-27.3241,0.159254,-0.520988,-0.0235157,0.424389,0.132128,00,00],
     [00,00,00,00,00,00,0.227472,-0.0735535,-26.9135,0.206825,-0.227774,-0.132128,0.424389,00,00],
     [00,00,00,00,00,00,0.332188,-1.02565,-0.0471482,-28.3499,0.132128,-0.227774,0.0235157,00,00],
     [00,00,00,-0.00126949,-0.000645326,-0.000708215,0.00763923,0.0016513,-0.0044564,-0.00284703,-0.119342,-0.000941301,-0.00160961,-1.62792e-05,0.00156863],
     [00,00,00,-0.000414104,-0.000210504,-0.000231018,0.0024919,0.000538651,-0.00145367,-0.000928698,0.000941301,-0.119493,0.000742542,-0.00155303,-7.50989e-06],
     [00,00,00,00,00,00,00,00,00,00,2.1684e-19,-2.1684e-19,-0.171654,-0.000868032,0.000868032],
     [00,00,00,0.00181259,0.000921407,0.0010112,-0.0109074,-0.00235775,0.00636292,0.00406504,0,-0.000567731,-0.000868032,-0.00087529,00],
     [00,00,00,0.00166876,0.000848291,0.000930959,-0.0100419,-0.00217066,0.005858,0.00374247,0.00052268,0,-0.000868032,00,-0.000874062]
     ]

    A = ca.sparsify(ca.DM(A))

    b = ca.DM(list(range(15)))
    H = 5

    Bs = ca.MX.sym("Bs",b.sparsity())
    As = ca.MX.sym("As",A.sparsity())

    Ast = As.T

    r= ca.Function("r", [As,Bs],[ca.solve(Ast,Bs,"csparse")])
    R= ca.Function("R", [As,Bs],[ca.solve(ca.densify(Ast),Bs,"csparse")])

    r_in = (A, b)
    r_out = r(*r_in)
    R_out = R(*r_in)

    r.sparsity_out(0).spy()
    R.sparsity_out(0).spy()

    self.checkarray(R_out,numpy.linalg.solve(A.T,b))
    self.checkarray(r_out,R_out)

  def test_depends_on(self):
    a = ca.MX.sym("a")
    b = ca.MX.sym("b")

    self.assertTrue(ca.depends_on(a**2,a))
    self.assertTrue(ca.depends_on(a,a))
    self.assertFalse(ca.depends_on(0,a))

    ab = ca.vertcat(*(a,b))
    self.assertTrue(ca.depends_on(a**2,ab))
    self.assertTrue(ca.depends_on(a,ab))
    self.assertFalse(ca.depends_on(0,ab))
    self.assertTrue(ca.depends_on(b**2,ab))
    self.assertTrue(ca.depends_on(b,ab))
    self.assertTrue(ca.depends_on(a**2+b**2,ab))
    self.assertTrue(ca.depends_on(a+b,ab))
    self.assertTrue(ca.depends_on(ca.vertcat(*[0,a]),a))
    self.assertTrue(ca.depends_on(ca.vertcat(*[a,0]),a))
    self.assertFalse(ca.depends_on(ca.vertcat(*[0,b]),a))
    self.assertFalse(ca.depends_on(ca.vertcat(*[b,0]),a))
    self.assertTrue(ca.depends_on(ca.vertcat(*[a**2,b**2]),ab))
    self.assertTrue(ca.depends_on(ca.vertcat(*[a,0]),ab))
    self.assertTrue(ca.depends_on(ca.vertcat(*[0,b]),ab))
    self.assertTrue(ca.depends_on(ca.vertcat(*[b,0]),ab))
    self.assertFalse(ca.depends_on(ca.vertcat(*[0,0]),ab))

  def test_vertcat_simp(self):
    x = ca.MX.sym("x",10)
    y = ca.MX.sym("y")
    z = ca.MX.sym("z")
    x_ = ca.DM(list(range(10)))
    y_ = ca.DM([20])
    z_ = ca.DM([30])

    def evalvertcat(*a):
      f = ca.Function("f", [x,y,z],[ca.vertcat(*a)])
      f_in = [0]*f.n_in()  # type: list

      f_in[0]=x_
      f_in[1]=y_
      f_in[2]=z_
      return f(*f_in)

    self.checkarray(evalvertcat(*ca.vertsplit(x)),x_)
    self.checkarray(evalvertcat(*ca.vertsplit(x)+[y]),ca.vertcat(*[x_,y_]))
    self.checkarray(evalvertcat(*[z]+ca.vertsplit(x)+[y] + ca.vertsplit(x)+[z]),ca.vertcat(*[z_,x_,y_,x_,z_]))
    self.checkarray(evalvertcat(*ca.vertsplit(x)[:-1]),x_[:-1])
    self.checkarray(evalvertcat(*ca.vertsplit(x)[:-1]+[y]),ca.vertcat(*[x_[:-1],y_]))
    self.checkarray(evalvertcat(*[z]+ca.vertsplit(x)[:-1]+[y] + ca.vertsplit(x)[:-1]+[z]),ca.vertcat(*[z_,x_[:-1],y_,x_[:-1],z_]))
    self.checkarray(evalvertcat(*ca.vertsplit(x)[1:]),x_[1:])
    self.checkarray(evalvertcat(*ca.vertsplit(x)[1:]+[y]),ca.vertcat(*[x_[1:],y_]))
    self.checkarray(evalvertcat(*[z]+ca.vertsplit(x)[1:]+[y] + ca.vertsplit(x)[1:]+[z]),ca.vertcat(*[z_,x_[1:],y_,x_[1:],z_]))
    g = ca.vertsplit(x)[5:]+ca.vertsplit(x)[:5]
    self.checkarray(evalvertcat(*g),ca.vertcat(*[x_[5:],x_[:5]]))
    self.checkarray(evalvertcat(*g+[y]),ca.vertcat(*[x_[5:],x_[:5],y_]))
    self.checkarray(evalvertcat(*[z]+g+[y] + g+[z]),ca.vertcat(*[z_,x_[5:],x_[:5],y_,x_[5:],x_[:5],z_]))



    w = ca.vertsplit(x,2)
    r = builtins.sum([ca.vertsplit(i) for i in w],[])

    self.checkarray(evalvertcat(*r),x_)

    w = ca.vertsplit(x,2)
    r = builtins.sum([ca.vertsplit(i)+[y] for i in w],[])
    print("vertcat:", r)
    print("result:", ca.vertcat(*r))

    w = ca.vertsplit(x,2)
    r = builtins.sum([ca.vertsplit(i) for i in w],[])
    print("vertcat:", r)
    print("result:", ca.vertcat(*r+[y]))

    self.assertTrue(ca.is_equal(ca.vertcat(*ca.vertsplit(x)),x))

  def test_horzcat_simp(self):
    x = ca.MX.sym("x",1,10)
    y = ca.MX.sym("y")
    z = ca.MX.sym("z")
    x_ = ca.DM(list(range(10))).T
    y_ = ca.DM([20])
    z_ = ca.DM([30])

    def evalhorzcat(*a):
      f = ca.Function("f", [x,y,z],[ca.horzcat(*a)])
      return f(x_, y_, z_)

    self.checkarray(evalhorzcat(*ca.horzsplit(x)),x_)
    self.checkarray(evalhorzcat(*ca.horzsplit(x)+[y]),ca.horzcat(*[x_,y_]))
    self.checkarray(evalhorzcat(*[z]+ca.horzsplit(x)+[y] + ca.horzsplit(x)+[z]),ca.horzcat(*[z_,x_,y_,x_,z_]))
    self.checkarray(evalhorzcat(*ca.horzsplit(x)[:-1]),x_[0,:-1])
    self.checkarray(evalhorzcat(*ca.horzsplit(x)[:-1]+[y]),ca.horzcat(*[x_[0,:-1],y_]))
    self.checkarray(evalhorzcat(*[z]+ca.horzsplit(x)[:-1]+[y] + ca.horzsplit(x)[:-1]+[z]),ca.horzcat(*[z_,x_[0,:-1],y_,x_[0,:-1],z_]))
    self.checkarray(evalhorzcat(*ca.horzsplit(x)[1:]),x_[0,1:])
    self.checkarray(evalhorzcat(*ca.horzsplit(x)[1:]+[y]),ca.horzcat(*[x_[0,1:],y_]))
    self.checkarray(evalhorzcat(*[z]+ca.horzsplit(x)[1:]+[y] + ca.horzsplit(x)[1:]+[z]),ca.horzcat(*[z_,x_[0,1:],y_,x_[0,1:],z_]))
    g = ca.horzsplit(x)[5:]+ca.horzsplit(x)[:5]
    self.checkarray(evalhorzcat(*g),ca.horzcat(*[x_[0,5:],x_[0,:5]]))
    self.checkarray(evalhorzcat(*g+[y]),ca.horzcat(*[x_[0,5:],x_[0,:5],y_]))
    self.checkarray(evalhorzcat(*[z]+g+[y] + g+[z]),ca.horzcat(*[z_,x_[0,5:],x_[0,:5],y_,x_[0,5:],x_[0,:5],z_]))


    w = ca.horzsplit(x,2)
    r = builtins.sum([ca.horzsplit(i) for i in w],[])

    self.checkarray(evalhorzcat(*r),x_)

    w = ca.horzsplit(x,2)
    r = builtins.sum([ca.horzsplit(i)+[y] for i in w],[])
    print("vertcat:", r)
    print("result:", ca.horzcat(*r))

    w = ca.horzsplit(x,2)
    r = builtins.sum([ca.horzsplit(i) for i in w],[])
    print("vertcat:", r)
    print("result:", ca.horzcat(*r+[y]))

    self.assertTrue(ca.is_equal(ca.horzcat(*ca.horzsplit(x)),x))

  def test_vertsplit_simp(self):

    dvars = [ca.MX.sym("abcdefghijklm"[i]) for i in range(5) ]
    dvars_ = list(range(5))

    zz = ca.MX.sym("zz",2)
    zz_ = ca.DM([11,12])
    y = ca.MX.sym("y")
    z = ca.MX.sym("z")
    y_ = ca.DM([20])
    z_ = ca.DM([30])

    aa = ca.MX.sym("aa",5)
    aa_ = list(range(100,105))

    def evalvertsplit(a,*args):
      print(ca.vertsplit(a,*args))
      f = ca.Function("f", dvars+[y,z,zz,aa],ca.vertsplit(a,*args))
      f_in = [0]*f.n_in()  # type: list
      for i in range(5):
        f_in[i]=dvars_[i]
      f_in[5+0]=y_
      f_in[5+1]=z_
      f_in[5+2]=zz_
      f_in[5+3]=aa_
      f_out = f.call(f_in)
      return [f_out[i] for i in range(f.n_out())]

    s= evalvertsplit(ca.vertcat(*[y]+dvars+[z]))
    self.checkarray(s[0],y_)
    for i in range(5):
      self.checkarray(s[1+i],dvars_[i])
    self.checkarray(s[6],z_)

    s= evalvertsplit(ca.vertcat(*[y]+dvars+[z]),2)

    self.checkarray(s[0],ca.vertcat(*[y_,dvars_[0]]))
    self.checkarray(s[1],ca.vertcat(*[dvars_[1],dvars_[2]]))
    self.checkarray(s[2],ca.vertcat(*[dvars_[3],dvars_[4]]))
    self.checkarray(s[3],ca.vertcat(*[z_]))

    s= evalvertsplit(ca.vertcat(*[y,zz,z,zz]),2)

    self.checkarray(s[0],ca.vertcat(*[y_,zz_[0]]))
    self.checkarray(s[1],ca.vertcat(*[zz_[1],z_]))
    self.checkarray(s[2],zz_)

    s= evalvertsplit(ca.vertcat(*[y,zz,z,zz]),3)

    self.checkarray(s[0],ca.vertcat(*[y_,zz_[0],zz_[1]]))
    self.checkarray(s[1],ca.vertcat(*[z_,zz_[0],zz_[1]]))

    s= evalvertsplit(ca.vertcat(*[zz,zz]),2)
    self.checkarray(s[0],zz_)
    self.checkarray(s[1],zz_)

    s= evalvertsplit(ca.vertcat(*[zz]+dvars))
    self.checkarray(s[0],zz_[0])
    self.checkarray(s[1],zz_[1])

    for i in range(5):
      self.checkarray(s[2+i],dvars_[i])

    s= evalvertsplit(ca.vertcat(*dvars+[aa]),5)
    self.checkarray(s[0],ca.DM(dvars_))
    self.checkarray(s[1],ca.DM(aa_))

    s= evalvertsplit(ca.vertcat(*dvars+[aa]),4)
    self.checkarray(s[0],ca.DM(dvars_[:4]))
    self.checkarray(s[1],ca.DM([dvars_[-1]]+aa_[:3]))
    self.checkarray(s[2],ca.DM(aa_[3:]))

    s= evalvertsplit(ca.vertcat(*dvars+[aa]),6)
    self.checkarray(s[0],ca.DM(dvars_+[aa_[0]]))
    self.checkarray(s[1],ca.DM(aa_[1:]))

    for i in range(5):
      self.assertTrue(ca.is_equal(ca.vertsplit(ca.vertcat(*dvars))[i],dvars[i]))

  def test_horzsplit_simp(self):

    dvars = [ca.MX.sym("abcdefghijklm"[i]) for i in range(5) ]
    dvars_ = list(range(5))

    zz = ca.MX.sym("zz",1,2)
    zz_ = ca.DM([11,12]).T
    y = ca.MX.sym("y")
    z = ca.MX.sym("z")
    y_ = ca.DM([20])
    z_ = ca.DM([30])

    aa = ca.MX.sym("aa",1,5)
    aa_ = list(range(100,105))

    def evalhorzsplit(a,*args):
      print(ca.horzsplit(a,*args))
      f = ca.Function("f", dvars+[y,z,zz,aa],ca.horzsplit(a,*args))
      f_in = [0]*f.n_in()  # type: list
      for i in range(5):
        f_in[i]=dvars_[i]
      f_in[5+0]=y_
      f_in[5+1]=z_
      f_in[5+2]=zz_
      f_in[5+3]=aa_
      f_out = f(*f_in)
      return [f_out[i] for i in range(f.n_out())]

    s= evalhorzsplit(ca.horzcat(*[y]+dvars+[z]))
    self.checkarray(s[0],y_)
    for i in range(5):
      self.checkarray(s[1+i],dvars_[i])
    self.checkarray(s[6],z_)

    s= evalhorzsplit(ca.horzcat(*[y]+dvars+[z]),2)

    self.checkarray(s[0],ca.vertcat(*[y_,dvars_[0]]).T)
    self.checkarray(s[1],ca.vertcat(*[dvars_[1],dvars_[2]]).T)
    self.checkarray(s[2],ca.vertcat(*[dvars_[3],dvars_[4]]).T)
    self.checkarray(s[3],ca.vertcat(*[z_]).T)

    s= evalhorzsplit(ca.horzcat(*[y,zz,z,zz]),2)

    self.checkarray(s[0],ca.vertcat(*[y_,zz_[0,0]]).T)
    self.checkarray(s[1],ca.vertcat(*[zz_[0,1],z_]).T)
    self.checkarray(s[2],zz_)

    s= evalhorzsplit(ca.horzcat(*[y,zz,z,zz]),3)

    self.checkarray(s[0],ca.vertcat(*[y_,zz_[0,0],zz_[0,1]]).T)
    self.checkarray(s[1],ca.vertcat(*[z_,zz_[0,0],zz_[0,1]]).T)

    s= evalhorzsplit(ca.horzcat(*[zz,zz]),2)
    self.checkarray(s[0],zz_)
    self.checkarray(s[1],zz_)

    s= evalhorzsplit(ca.horzcat(*[zz]+dvars))
    self.checkarray(s[0],zz_[0,0])
    self.checkarray(s[1],zz_[0,1])

    for i in range(5):
      self.checkarray(s[2+i],dvars_[i])

    s= evalhorzsplit(ca.horzcat(*dvars+[aa]),5)
    self.checkarray(s[0],ca.DM(dvars_).T)
    self.checkarray(s[1],ca.DM(aa_).T)

    s= evalhorzsplit(ca.horzcat(*dvars+[aa]),4)
    self.checkarray(s[0],ca.DM(dvars_[:4]).T)
    self.checkarray(s[1],ca.DM([dvars_[-1]]+aa_[:3]).T)
    self.checkarray(s[2],ca.DM(aa_[3:]).T)

    s= evalhorzsplit(ca.horzcat(*dvars+[aa]),6)
    self.checkarray(s[0],ca.DM(dvars_+[aa_[0]]).T)
    self.checkarray(s[1],ca.DM(aa_[1:]).T)

    for i in range(5):
      self.assertTrue(ca.is_equal(ca.horzsplit(ca.horzcat(*dvars))[i],dvars[i]))

  def test_vertsplit_derivative(self):
    m = ca.MX.sym("X",10)

    f = ca.Function("f", [m],[ca.vertsplit(m)[0]])

    f.reverse(1)

  def test_MX_const_sp(self):
    x = ca.MX.sym("x",4,1)

    sp = ca.Sparsity.triplet(3,3,[0,1,2,2],[0,0,1,2])

    f = ca.Function("f", [x],[x.nz[ca.DM(sp,list(range(sp.nnz())))]])

    g = ca.Function("g", [x],[ca.MX(sp,x)])

    f_in = [0]*f.n_in()  # type: list

    f_in[0]=ca.DM(list(range(1,5)))

    self.checkfunction(f,g,inputs=f_in)
    self.check_codegen(f,inputs=f_in)
    self.check_codegen(g,inputs=f_in)

  def test_reshape_sp(self):
    x = ca.MX.sym("x",4,1)

    f = ca.Function("f", [x],[x.reshape((2,2))])

    sx = ca.SX.sym("x",4,1)

    g = ca.Function("g", [sx],[sx.reshape((2,2))])

    f_in = [0]*f.n_in()  # type: list

    f_in[0]=ca.DM(list(range(1,5)))

    self.checkfunction(f,g,inputs=f_in)
    self.check_codegen(f,inputs=f_in)
    self.check_codegen(g,inputs=f_in)
    
  def test_issue1041(self):
    x = ca.MX.sym("x",2)

    y = ca.vertsplit(x,[0,1,2])[1]

    f = ca.Function("f", [x],[y])

    H = hessian_old(f, 0, 0)

  def test_bug_1042(self):

    x = ca.MX.sym('x',2,1)

    mf = ca.Function("mf", [x],[x*x[0,0]])

    mfunction = mf.expand('expand_'+mf.name())

    mfg = mf.reverse(1)

    mfunctiong = mfunction.reverse(1)

    f_in = [0, 5, ca.DM([1,2])]

    self.checkfunction(mfg,mfunctiong,inputs=f_in)

  def test_bug_1042bis(self):
    x = ca.MX.sym('x',2,1)
    a = ca.MX.sym("ax",2,1)
    i1 = x[0,0]
    z = i1*x
    i3 = i1*a
    i3= c.dot(x,a)
    d = ca.Function("d", [x,a],[z,i3])

    dx = d.expand('expand_'+d.name())
    dx_in = [0]*dx.n_in()  # type: list

    dx_in[0]=ca.DM([1,2])
    dx_in[1]=ca.DM([3,4])

    self.checkfunction(d,dx,inputs=dx_in)

  def test_bug_1042tris(self):
    x = ca.MX.sym('x',2,1)
    a = ca.MX.sym("ax",2,1)
    d = ca.Function("d", [x,a],[c.dot(x,a)])

    dx = d.expand('expand_'+d.name())
    dx_in = [0]*dx.n_in()  # type: list

    dx_in[0]=ca.DM([1,2])
    dx_in[1]=ca.DM([3,4])

    self.checkfunction(d,dx,inputs=dx_in)

  def test_bug_1046(self):
    x = ca.MX.sym('x',1,1)
    y = ca.MX.sym('y',1,1)
    z = ca.jacobian(x,y)

    self.assertTrue(z.nnz()==0)

  def test_expand_free(self):
    x = ca.MX.sym("x")
    y = ca.MX.sym("y")
    with self.assertInException("free"):
      f = ca.Function('f',[x],[x+y])

    f = ca.Function('f',[x],[x+y],{"allow_free": True})
    g = ca.Function('g',[x,y],[f(x)])
    ge = g.expand()
    self.checkfunction(g,g,inputs=[1,2])

  def test_singularcat(self):

    for c in [ca.MX,ca.SX,ca.DM]:
      x0 = c.zeros(10,0)

      x1s = ca.vertsplit(x0, [0,5,10])

      for x in x1s:
        self.checkarray(x.shape,(5,0))


      x2 = ca.vertcat(*x1s)
      self.checkarray(x2.shape,(10,0))

      x2 = ca.vertcat(*[c.zeros(0,0)] + x1s + [c.zeros(0,0)])
      self.checkarray(x2.shape,(10,0))

      x0 = c.zeros(0,1)
      x2 = ca.vertcat(x0,c.zeros(0,0),x0)
      self.checkarray(x2.shape,(0,1))

    for c in [ca.MX,ca.SX,ca.DM]:
      x0 = c.zeros(0,10)

      x1s = ca.horzsplit(x0, [0,5,10])

      for x in x1s:
        self.checkarray(x.shape,(0,5))

      x2 = ca.horzcat(*x1s)
      self.checkarray(x2.shape,(0,10))

      x2 = ca.horzcat(*[c.zeros(0,0)] + x1s + [c.zeros(0,0)])
      self.checkarray(x2.shape,(0,10))

      x0 = c.zeros(1,0)
      x2 = ca.horzcat(x0,c.zeros(0,0),x0)
      self.checkarray(x2.shape,(1,0))

    for c in [ca.MX,ca.SX,ca.DM]:
      x0 = c.zeros(10,0)

      x1s = ca.vertsplit(x0, [0,5,10])

      x0 = c.zeros(0,10)
      x1st = ca.horzsplit(x0, [0,5,10])

      x2 = ca.diagcat(*x1s)
      self.checkarray(x2.shape,(10,0))

      x2 = ca.diagcat(*[c.zeros(0,0)] + x1s + [c.zeros(0,0)])
      self.checkarray(x2.shape,(10,0))

      x2 = ca.diagcat(*x1st)
      self.checkarray(x2.shape,(0,10))

      x2 = ca.diagcat(*[c.zeros(0,0)] + x1st + [c.zeros(0,0)])
      self.checkarray(x2.shape,(0,10))

      x2 = ca.diagcat(*x1s+x1st)
      self.checkarray(x2.shape,(10,10))

      x2 = ca.diagcat(*[c.zeros(0,0)] + x1s+x1st + [c.zeros(0,0)])
      self.checkarray(x2.shape,(10,10))

  def test_exprjac(self):


    for fun in [lambda x: x[0]*x[1],lambda x: x[0]*ca.sin(x[1])]:
      def hessian1(f, x): return ca.hessian(f, x)[0]
      for op in [c.gradient, ca.jacobian, hessian1]:
        print(fun, op)
        x = ca.MX.sym("x",2)
        print(fun(x))
        print(op(fun(x),x))
        f = ca.Function("f", [x],[op(fun(x),x)])

        x = ca.SX.sym("x",2)
        fr = ca.Function("fr", [x],[op(fun(x),x)])

        self.checkfunction(f,fr,inputs=[0])
        self.check_codegen(f,inputs=[0])

  @memory_heavy()
  def test_einstein(self):
        dim_b = [3, 3]
        dim_a = [3, 3]




        def free(a):
            return set([i for i in a if i<0])


        def dim_map(dim_a, dim_b, dim_c, ind_a, ind_b, ind_c):
          dim_map = {}
          for dim, ind in [(dim_a,ind_a), (dim_b,ind_b)]:
            for i,j in enumerate(ind):
              if j<0:
                if j in dim_map:
                  assert dim_map[j]==dim[i]
                else:
                  dim_map[j]=dim[i]
          return dim_map

        def sub2ind(dims, sub):
          ret = []
          for i in range(len(dims)):
            ret.append(sub % dims[i])
            sub/= dims[i]
          return ret

        def ind2sub(dims, ind):
          ret=0
          cumprod = 1
          for i in range(len(dims)):
            ret+= ind[i]*cumprod
            cumprod*= dims[i]
          return ret

        def combine(dim_a, dim_b, ind_a, ind_b, ind_c):
          dim_map = {}
          for dim, ind in [(dim_a,ind_a), (dim_b,ind_b)]:
            for i,j in enumerate(ind):
              if j<0:
                if j in dim_map:
                  assert dim_map[j]==dim[i]
                else:
                  dim_map[j]=dim[i]
          return [dim_map[i] for i in ind_c]


        def subs(e,old,n):
          return [ n[old.index(i)] if i in old else i for i in e ]


        def einstein_tests(dim_a, dim_b, dim_c, ind_a, ind_b, ind_c):


          A = ca.MX.sym("A", np.prod(dim_a))
          B = ca.MX.sym("B", np.prod(dim_b))
          C = ca.MX.sym("C", int(np.prod(dim_c)),1)

          A_sx = ca.SX.sym("A", np.prod(dim_a))
          B_sx = ca.SX.sym("B", np.prod(dim_b))
          C_sx = ca.SX.sym("C", int(np.prod(dim_c)),1)

          np.random.seed(0)

          A_ = np.random.random(A.shape)
          B_ = np.random.random(B.shape)
          C_ = np.random.random(C.shape)



          out = ca.einstein(A, B, C, dim_a, dim_b, dim_c, ind_a, ind_b, ind_c)
          out_sx = ca.einstein(A_sx, B_sx, C_sx, dim_a, dim_b, dim_c, ind_a, ind_b, ind_c)
          f = ca.Function('f',[A,B,C],[out],{"ad_weight_sp": 0})
          frev = ca.Function('f',[A,B,C],[out],{"ad_weight_sp": 1})
          fr = ca.Function('fr',[A,B,C],[my_einstein(A, B, C, dim_a, dim_b, dim_c, ind_a, ind_b, ind_c)])
          f_sx = ca.Function('f',[A_sx,B_sx,C_sx],[out_sx])

          fsx = f.expand()
          self.checkfunction(f, fr, inputs=[A_,B_,C_])
          self.checkfunction(fsx, fr, inputs=[A_,B_,C_])
          self.check_codegen(f, inputs=[A_,B_,C_])
          self.check_serialize(f, inputs=[A_,B_,C_])
          self.checkfunction(fsx, f_sx, inputs=[A_,B_,C_])

          for i in range(3):
            self.check_sparsity(f.jac_sparsity(0, i), fsx.jac_sparsity(0, i))
            self.check_sparsity(frev.jac_sparsity(0, i), fsx.jac_sparsity(0, i))


        def my_einstein(A, B, C, dim_a, dim_b, dim_c, ind_a, ind_b, ind_c):
          try:
            R = ca.DM(C)
          except:
            R = ca.MX(C)

          d = dim_map(dim_a, dim_b, dim_c, ind_a, ind_b, ind_c)

          dk = list(d.keys())

          for val in itertools.product(*[range(d[k]) for k in dk]):
            ind_a_ = subs(ind_a, dk, val)
            ind_b_ = subs(ind_b, dk, val)
            ind_c_ = subs(ind_c, dk, val)


            ai = ind2sub(dim_a, ind_a_)
            bi = ind2sub(dim_b, ind_b_)
            ci = ind2sub(dim_c, ind_c_)

            R[ci]+= A[ai]*B[bi]

          return R

        ind = [ [0, -1], [1, -1], [0, 1], [-1, -2], [-2, -1] ]

        for ind_a in ind:
            for ind_b in ind:
                for ind_c in itertools.permutations(list(free(ind_a) ^ free(ind_b))):
                  dim_c = combine(dim_a, dim_b, ind_a, ind_b, ind_c)
                  einstein_tests(dim_a, dim_b, dim_c, ind_a, ind_b, ind_c)

        einstein_tests([2,4,3], [2,5,3], [5, 4], [-1, -2, -3], [-1, -4, -3], [-4, -2])

  def test_sparsity_operation(self):
    L = [ca.MX(ca.Sparsity(1,1)),ca.MX(ca.Sparsity(2,1)), ca.MX.sym("x",1,1), ca.MX.sym("x", ca.Sparsity(1,1)), ca.DM(1), ca.DM(ca.Sparsity(1,1),1), ca.DM(ca.Sparsity(2,1),1), ca.DM(ca.Sparsity.dense(2,1),1)]

    for a in L:
      for b in L:
        c = a*b

        if a.nnz()==0 or b.nnz()==0:
          self.assertTrue(c.nnz()==0)
        else:
          self.assertTrue(c.nnz()>0)

  @requires_linsol("lapackqr")
  def test_solve(self):
    print(123)
    N = 3
    nrhs = 50
    np.random.seed(0)
    A = np.random.random((3,3))
    B = np.random.random((3,50))


    C = ca.solve(A, B, "lapackqr", {"max_nrhs": 50})
    C1 = ca.solve(A, B, "lapackqr", {"max_nrhs": 20})
    C2 = ca.solve(A, B, "lapackqr")

    self.checkarray(C, C1)
    self.checkarray(C1, C2)
  def test_sparse_lt(self):
    x = ca.MX.sym("x",ca.Sparsity.lower(5))
    self.assertEqual((x>0).nnz(),5*6/2)
    self.assertEqual((x>=0).nnz(),5*5)
  def test_inv(self):
   np.random.seed(0)

   for X in [ca.SX, ca.MX]:
     A  = X.sym("x",3,3)
     Av = np.random.random((3,3))
     f = ca.Function('f',[A],[ca.inv(A),ca.inv(ca.DM(Av)),A.__mpower__(-1), ca.DM(Av).__mpower__(-1)])
     out = f(Av)
     for o in out:
       self.checkarray(o, np.linalg.inv(np.array(Av)))


  def test_interp1d(self):
    v = [7,3,4,-3]
    vs = ca.MX.sym("x",4,1)

    xq = [10,0,1,2,4,8,7,5,1.5]
    x = [1,2,4,8]
    F = ca.Function("f",[vs],[ca.interp1d(x,vs,xq)])

    self.checkarray(F(v),np.interp(xq,x,v))

    F = ca.Function("f",[vs],[ca.interp1d(x,vs,xq,"floor")])

    self.checkarray(F(v),ca.DM([-3,7,7,3,4,-3,4,4,7]))

    F = ca.Function("f",[vs],[ca.interp1d(x,vs,xq,"ceil")])

    self.checkarray(F(v),ca.DM([-3,7,7,3,4,-3,-3,-3,3]))

    v = [7,3,4,-3]
    vs = ca.MX.sym("x",4,1)

    xq = [10,0,1,2,3,4,3.5,3.25,1.5]
    x = [1,2,3,4]
    F = ca.Function("f",[vs],[ca.interp1d(x,vs,xq,"linear",True)])

    self.checkarray(F(v),np.interp(xq,x,v))

  @memory_heavy()
  def test_bilin_etc(self):
    x = ca.MX.sym("x",3,3)
    y = ca.MX.sym("y",3,1)
    z = ca.MX.sym("z",3,1)

    import numpy
    numpy.random.seed(42)
    x0 = numpy.random.random((3,3))
    y0 = numpy.random.random((3,1))
    z0 = numpy.random.random((3,1))

    for e in [(ca.bilin(x,y,z),y.T @ x @ z),(ca.rank1(x,0.3,y,z),x+0.3*(y @ z.T))]:
      f = ca.Function('f',[x,y,z],[e[0]])
      fref = ca.Function('fref',[x,y,z],[e[1]])
      f_sx = f.expand()
      self.checkfunction(f,fref,inputs=[x0,y0,z0])
      self.checkfunction(f_sx,fref,inputs=[x0,y0,z0])
      self.check_codegen(f,inputs=[x0,y0,z0])
      
  def test_bilin_short(self):
    for X in [ca.SX,ca.MX]:
        x = X.sym("x",3,3)
        y = X.sym("y",3,1)

        import numpy
        numpy.random.seed(42)
        x0 = numpy.random.random((3,3))
        y0 = numpy.random.random((3,1))
        
        f = ca.Function("f",[x,y],[ca.bilin(x,y)])
        fref = ca.Function("f",[x,y],[ca.bilin(x,y,y)])
        
        self.checkfunction(f,fref,inputs=[x0,y0])
        self.check_codegen(f,inputs=[x0,y0])
        self.check_codegen(fref,inputs=[x0,y0])

  def test_det_shape(self):
    X = ca.MX.sym("x",2,3)
    with self.assertRaises(RuntimeError):
      ca.det(X)
    X = ca.MX.sym("x",3,3)
    ca.det(X)


  def test_det(self):
    X = ca.MX.sym("x",3,3)
    x = ca.SX.sym("x",3,3)

    import numpy
    numpy.random.seed(42)
    x0 = numpy.random.random((3,3))

    f = ca.Function('f',[x],[ca.det(x)])
    F = ca.Function('F',[X],[ca.det(X)])
    self.checkfunction(f,F,inputs=[x0])
    self.check_codegen(F,inputs=[x0])

  @memory_heavy()
  def test_det_solvers(self):
    # Determinant through the linear-solver plugins, across sizes and a sparse pattern
    numpy.random.seed(1)

    solvers = ["qr"]
    if ca.has_linsol("csparse"): solvers.append("csparse")

    # Dense cases of several sizes
    dense_cases = [numpy.random.random((n,n)) for n in [1,2,3,4,5]]

    # A genuinely sparse (but nonsingular) pattern
    S = numpy.random.random((5,5))
    S[0,2]=S[0,4]=S[2,0]=S[3,1]=S[4,2]=0
    S += 3*numpy.eye(5)

    for A in dense_cases:
      n = A.shape[0]
      ref = numpy.linalg.det(A)
      X = ca.MX.sym("x", n, n)
      x = ca.SX.sym("x", n, n)

      # DM (default cofactor) and SX (cofactor) match numpy
      self.checkarray(ca.det(ca.DM(A)), ref, digits=9)
      self.checkarray(ca.evalf(ca.det(ca.SX(ca.DM(A)))), ref, digits=9)

      # SX symbolic determinant via the symbolicqr plugin's static method
      xs = ca.SX.sym("xs", n, n)
      self.checkarray(ca.evalf(ca.substitute(ca.det(xs, "symbolicqr"), xs, ca.DM(A))), ref, digits=9)

      f = ca.Function('f',[x],[ca.det(x)])
      for ls in solvers:
        # MX node with this solver: symbolic value + AD consistency vs cofactor
        F = ca.Function('F',[X],[ca.det(X, ls)])
        self.checkfunction(f, F, inputs=[A], digits=8)
        # DM numeric through the plugin
        self.checkarray(ca.det(ca.DM(A), ls), ref, digits=9)
        # Codegen for solvers that support it (qr)
        if ls=="qr":
          self.check_codegen(F, inputs=[A])

      # Options dict is forwarded to the solver
      self.checkarray(ca.det(ca.DM(A), "qr", {"eps":1e-13}), ref, digits=9)

    # Sparse determinant (numeric DM path) for every available plugin
    refS = numpy.linalg.det(S)
    self.checkarray(ca.det(ca.DM(S)), refS, digits=9)
    for ls in solvers:
      self.checkarray(ca.det(ca.sparsify(ca.DM(S)), ls), refS, digits=9)

    # Symbolic plugin determinant is only provided by symbolicqr; a plugin
    # without the exposed static method raises a clear error.
    with self.assertInException("does not provide a symbolic determinant"):
      ca.det(ca.SX.sym("x",2,2), "qr")

    # csparse exposes a numeric det but not the symbolic (static) one. This must
    # be a clean error, NOT a segfault from an uninitialised Exposed.det pointer.
    if ca.has_linsol("csparse"):
      with self.assertInException("does not provide a symbolic determinant"):
        ca.det(ca.SX.sym("x",2,2), "csparse")

  def test_det_symbolicqr_btf(self):
    # The symbolicqr determinant exploits block-triangular structure (BTF):
    # det = sign(perms) * prod(det(diagonal blocks)).
    numpy.random.seed(2)

    def sxdet(Anp):
      Asp = ca.sparsify(ca.DM(Anp))
      x = ca.SX.sym("x", Asp.sparsity())
      f = ca.Function("f", [x], [ca.det(x, "symbolicqr")])
      return f(Asp), f.n_nodes()

    # Block diagonal (two 4x4 blocks) must equal numpy det and be cheaper than dense
    n = 8
    B = numpy.zeros((n, n))
    B[:4, :4] = numpy.random.random((4, 4))
    B[4:, 4:] = numpy.random.random((4, 4))
    vb, nodes_b = sxdet(B)
    self.checkarray(vb, numpy.linalg.det(B), digits=9)

    Adense = numpy.random.random((n, n))
    xd = ca.SX.sym("x", n, n)
    nodes_d = ca.Function("f", [xd], [ca.det(xd, "symbolicqr")]).n_nodes()
    self.assertTrue(nodes_b < nodes_d)  # sparsity actually exploited

    # Block lower-triangular [[A11,0],[A21,A22]]: det = det(A11)*det(A22), so the
    # result must be structurally INDEPENDENT of the coupling block A21. BTF drops
    # A21 entirely; a single global QR would factorize it into the expression
    # (giving a false structural dependence). This is the unambiguous BTF win.
    L = numpy.random.random((6, 6)) + 3*numpy.eye(6)
    L[:3, 3:] = 0                                   # zero upper-right -> block lower-tri
    Lsp = ca.sparsify(ca.DM(L))
    xL = ca.SX.sym("x", Lsp.sparsity())
    dL = ca.det(xL, "symbolicqr")
    self.checkarray(ca.evalf(ca.substitute(dL, xL, Lsp)), numpy.linalg.det(L), digits=9)
    rows = Lsp.sparsity().row(); cols = Lsp.sparsity().get_col()
    coupling = [k for k in range(Lsp.nnz()) if rows[k] >= 3 and cols[k] < 3]
    self.assertTrue(len(coupling) > 0)              # there really is a coupling block
    self.assertEqual(ca.jacobian(dL, xL)[:, coupling].nnz(), 0)  # det ignores it

    # Fully triangular -> all 1x1 blocks -> product of the diagonal entries
    U = numpy.triu(numpy.random.random((5, 5)) + numpy.eye(5))
    vU, nodes_U = sxdet(U)
    self.checkarray(vU, numpy.linalg.det(U), digits=9)
    self.assertTrue(nodes_U < 30)  # collapses to a tiny product

    # Structurally singular (a blank row) -> determinant is structurally zero
    Sg = numpy.random.random((4, 4)); Sg[2, :] = 0
    self.checkarray(sxdet(Sg)[0], 0, digits=12)

  def test_mtimes_mismatch_segfault(self):
    with self.assertInException("incompatible dimensions"):
      ca.DM(ca.Sparsity.lower(5)) @ ca.MX.sym('x',100)

  def test_monitor(self):
    x = ca.MX.sym("x")
    y = ca.sqrt(x.monitor("hey"))

    f = ca.Function('f',[x],[y])
    with capture_stdout() as out:
      f(3)
    print(out)
    
    
    ref = "hey:\n3.0000000000000000e+00\n"
    self.assertTrue(out[0]==ref)

    if args.run_slow:
        with self.assertInException("Parsing error"):
            self.check_codegen(f,inputs=[3],main=True)
                    
        with open("f_out.txt","r") as f_out:
            out2 = f_out.read()
            
        print(out2)
            
        self.assertTrue(ref in out2)

    x = ca.MX.sym("x",2)
    y = ca.sqrt(x.monitor("hey"))

    f = ca.Function('f',[x],[y])
    with capture_stdout() as out:
      f(3)
    print(out)
    ref = "hey:\n2x1: [3.0000000000000000e+00, 3.0000000000000000e+00]\n"
    self.assertTrue(out[0]==ref)

    if args.run_slow:
        with self.assertInException("Parsing error"):
            self.check_codegen(f,inputs=[3],main=True)
                    
        with open("f_out.txt","r") as f_out:
            out2 = f_out.read()
            
        print(out2)
            
        self.assertTrue(ref in out2)

  def test_dump(self):
    import os
    import glob as globmod

    # Cleanup any leftover dump files
    for f in globmod.glob("myvar.*.mtx"):
      os.remove(f)

    # Scalar test
    x = ca.MX.sym("x")
    y = ca.sqrt(x.dump("myvar"))

    f = ca.Function('f', [x], [y])

    f(9)

    self.assertTrue(os.path.exists("myvar.000000.mtx"))
    val = ca.DM.from_file("myvar.000000.mtx")
    self.checkarray(val, 9)

    # Second call increments counter
    f(16)
    self.assertTrue(os.path.exists("myvar.000001.mtx"))
    val = ca.DM.from_file("myvar.000001.mtx")
    self.checkarray(val, 16)

    # Cleanup
    for fname in globmod.glob("myvar.*.mtx"):
      os.remove(fname)

    # Sparse matrix test
    x = ca.MX.sym("x", ca.Sparsity.lower(2))
    y = 2*x.dump("mysp")
    f = ca.Function('f', [x], [y])

    inp = ca.sparsify(ca.DM([[1, 0], [3, 4]]))
    f(inp)

    self.assertTrue(os.path.exists("mysp.000000.mtx"))
    val = ca.DM.from_file("mysp.000000.mtx")
    self.checkarray(val, inp)

    for fname in globmod.glob("mysp.*.mtx"):
      os.remove(fname)

    # Test with dir option
    dump_dir = "dump_test_dir"
    if not os.path.exists(dump_dir):
      os.makedirs(dump_dir)
    x = ca.MX.sym("x")
    y = x.dump("myvar", {"dir": dump_dir})
    f = ca.Function('f', [x], [y])

    f(42)

    dump_path = os.path.join(dump_dir, "myvar.000000.mtx")
    self.assertTrue(os.path.exists(dump_path))
    val = ca.DM.from_file(dump_path)
    self.checkarray(val, 42)

    for fname in globmod.glob(os.path.join(dump_dir, "myvar.*.mtx")):
      os.remove(fname)
    os.rmdir(dump_dir)

    if args.run_slow and "ghc-filesystem" in ca.CasadiMeta.feature_list():
      import shutil
      for d in ["mydump", "mydump_codegen"]:
        if os.path.exists(d):
          shutil.rmtree(d)

      # Codegen test
      x = ca.MX.sym("x")
      y = ca.sqrt(x.dump("x",{"dir": "mydump"}))
      f = ca.Function('f', [x], [y])
      self.check_codegen(f, inputs=[9],std="c99",main=True,opts={"dump_dir_suffix": "_codegen"})

      self.file_equal("mydump/x.000000.mtx", "mydump_codegen/x.000000.mtx")

      for d in ["mydump", "mydump_codegen"]:
        if os.path.exists(d):
          shutil.rmtree(d)

  def test_codegen_specials(self):
    x = ca.MX.sym("x")
    y = ca.MX.sym("y")

    for z in [ x**2, ca.if_else(y>0,2*x,x*y), ca.fmin(x,y), ca.fmax(x,y), ca.sign(x*y)]:
      f = ca.Function('f',[x,y],[z])
      self.check_codegen(f,inputs=[1,2])
      self.check_codegen(f,inputs=[1,-2])


  @memory_heavy()
  def test_getsetnonzeros(self):
    import numpy
    numpy.random.seed(42)
    for S in [ca.Sparsity.lower(5),ca.Sparsity.dense(5,5)]:
      M = ca.MX.sym("X",S.nnz())

      m = ca.sin(ca.MX(S,M))
      for i in [0,2]:
        for ind in [(i,i),(1,i),(i,1),(slice(None),i),(i,slice(None)),(slice(i,i+2),slice(i,i+2)),(slice(i,i+2),slice(None)),([i,i],[0,0]),([],[])]:
          E = m.__getitem__(ind)
          e = ca.cos(E)
          e = ca.dot(e,e)
          f = ca.Function('f',[M],[e])
          self.checkfunction(f,f.expand(),inputs=[ numpy.random.random((S.nnz(),1))])
          self.check_codegen(f,inputs=[ numpy.random.random((S.nnz(),1))],digits=13)

          mc = m+0

          Y = ca.MX.sym("y",E.nnz())
          y = ca.MX(E.sparsity(),Y)
          mc.__setitem__(ind,y)
          e = ca.cos(mc)
          e = ca.dot(e,e)

          f = ca.Function('f',[M,Y],[e])
          self.checkfunction(f,f.expand(),inputs=[ numpy.random.random((S.nnz(),1)), numpy.random.random((E.nnz(),1))])
          self.check_serialize(f,inputs=[ numpy.random.random((S.nnz(),1)), numpy.random.random((E.nnz(),1))])
          self.check_codegen(f,inputs=[ numpy.random.random((S.nnz(),1)), numpy.random.random((E.nnz(),1))],digits=13)

  def test_evalf(self):
    x = ca.MX.sym("x")

    p = ca.MX.sym("p")
    f = ca.Function('f',[x],[ca.sin(x)])
    y = f.call([p],False,True)[0]
    y = ca.substitute(y,p,3)

    with self.assertInException("not defined"):
      y.to_DM()
    self.checkarray(ca.evalf(y),ca.sin(3))
    with self.assertInException("since variables [x] are free"):
      ca.evalf(x)


  def test_mmin(self):

    def mx_eval(f,X,x):
      y = X.sym("y",x.sparsity())
      return ca.Function('f',[y],f.call([y],X is ca.MX,False))

    for X in [ca.SX,ca.MX]:
      for mod in [lambda f:f, lambda f:f.expand(), lambda f: mx_eval(f,X,x)]:
        x = X.sym("x",2,2)

        f = mod(ca.Function('f',[x],[ca.mmin(x),ca.mmax(x)]))


        [mmin_res,mmax_res] = f(ca.DM([[-2,3],[-3,5]]))

        self.checkarray(mmin_res, -3)
        self.checkarray(mmax_res, 5)

        x = X.sym("x",ca.Sparsity.lower(2))

        f = mod(ca.Function('f',[x],[ca.mmin(x),ca.mmax(x)]))

        [mmin_res,mmax_res] = f(ca.sparsify(ca.DM([[-2,0],[-3,-5]])))

        self.checkarray(mmin_res, -5)
        self.checkarray(mmax_res, 0)

        [mmin_res,mmax_res] = f(ca.sparsify(ca.DM([[2,0],[3,5]])))

        self.checkarray(mmin_res, 0)
        self.checkarray(mmax_res, 5)

        x = X.sym("x",ca.Sparsity(2,2))
        f = mod(ca.Function('f',[x],[ca.mmin(x),ca.mmax(x)]))
        [mmin_res,mmax_res] = f(ca.DM(2,2))

        self.checkarray(mmin_res, 0)
        self.checkarray(mmax_res, 0)

        x = X.sym("x",ca.Sparsity(2,0))
        f = mod(ca.Function('f',[x],[ca.mmin(x),ca.mmax(x)]))
        [mmin_res,mmax_res] = f(ca.DM(2,0))

        self.assertTrue(mmin_res.is_empty())
        self.assertTrue(mmax_res.is_empty())

        x = X.sym("x",ca.Sparsity(0,0))
        f = mod(ca.Function('f',[x],[ca.mmin(x),ca.mmax(x)]))
        [mmin_res,mmax_res] = f(ca.DM(0,0))

        self.assertTrue(mmin_res.is_empty())
        self.assertTrue(mmax_res.is_empty())

    for m in [ca.mmin,ca.mmax]:
      x = ca.MX.sym("X",2)
      f = ca.Function("f",[x],[m(x)])
      #self.checkfunction(f,f.expand(),inputs=[[0.2,0.3]])

      #J = f.jacobian()
      #print(J([2,2],0))
      f.expand().disp(True)
      J = f.expand().jacobian()
      f.expand().jacobian().disp(True)
      print(J([2,2],0))
      self.checkfunction(f,f.expand(),inputs=[[2,2]])
      self.check_codegen(f,inputs=[[2,2]],std="c99")

  def test_doc_expression_tools(self):
    self.assertTrue("Given a repeated matrix, computes the sum of repeated parts." in (ca.repsum.__doc__ or ""))

  def test_densify(self):
    I = ca.sparsify(ca.DM([[0,1,0,1],[1,1,0,1],[0,1,1,0]])).sparsity()

    x = ca.MX.sym("x",I)

    f = ca.Function("f",[x],[ca.densify(x)])

    a = ca.DM([[0,1,0,2],[3,4,0,5],[0,6,7,0]])

    y = f(ca.sparsify(a))
    self.assertTrue(y.sparsity()==a.sparsity())
    self.checkarray(y,a)

    self.check_codegen(f,inputs=[a])


    x = ca.MX.sym("x",3,4)

    f = ca.Function("f",[x],[x[I]])

    a = ca.DM([[0,1,0,2],[3,4,0,5],[0,6,7,0]])
    b = ca.sparsify(a)
    a = ca.DM([[9,1,9,2],[3,4,9,5],[9,6,7,9]])

    y = f(a)
    self.assertTrue(y.sparsity()==b.sparsity())
    self.checkarray(y,b)

    self.check_codegen(f,inputs=[a])
    self.check_serialize(f,inputs=[a])

  def test_constant_from_file(self):

    open("test.txt", "w").write("5.6 1e-5")
    with self.assertInException("Failed to read a double from 'test.txt'. Expected 3 doubles."):
      x = ca.MX(ca.Sparsity.lower(2), "test.txt")

    open("test.txt", "w").write("5.6 abc 1e-5")
    with self.assertInException("Failed to read a double from 'test.txt'. Expected 3 doubles."):
      x = ca.MX(ca.Sparsity.lower(2), "test.txt")

    with self.assertInException("Cannot open file 'testQWHAL567p.txt'."):
      x = ca.MX(ca.Sparsity.lower(2), "testQWHAL567p.txt")

    open("test.txt", "w").write("5.6 1e-5 -12")
    x = ca.MX(ca.Sparsity.lower(2), "test.txt")

    a = ca.MX.sym("a")

    f = ca.Function("f",[a],[a*x])
    self.checkarray(f(2),ca.DM(ca.Sparsity.lower(2),[2*5.6,2*1e-5,-2*12]))

    self.check_codegen(f,inputs=[2])

    with self.assertInException("Not defined for ConstantFile"):
      x.to_DM()

  def test_nonzeros_param(self):
    x = ca.MX.sym("x",2)
    y = ca.MX.sym("y")

    for ad_weight_sp in [0,1]:
      opts = {"helper_options":{"ad_weight_sp":ad_weight_sp}}
      z = x.nz[y]
      J = ca.jacobian(z,x,opts)
      self.assertTrue(J.size()==(1,2))
      self.assertTrue(J.nnz()==2)
      J = ca.jacobian(z,y,opts)
      self.assertTrue(J.size()==(1,1))
      self.assertTrue(J.nnz()==0)

      q = ca.MX.sym("q")
      z = ca.MX(x)
      z.nz[y] = q

      J = ca.jacobian(z,x,opts)
      self.assertTrue(J.size()==(2,2))
      self.assertTrue(J.nnz()==4)
      J = ca.jacobian(z,y,opts)
      self.assertTrue(J.size()==(2,1))
      self.assertTrue(J.nnz()==0)
      J = ca.jacobian(z,q,opts)
      self.assertTrue(J.size()==(2,1))
      self.assertTrue(J.nnz()==2)

  def test_print_instructions(self):
    for X in [ca.MX,ca.SX]:
        x = X.sym("x",ca.Sparsity.diag(5))
        y = X.sym("y",5,5)
        z = X.sym("z",5)
        w = X.sym("w")
        g = ca.Function('g',[z],[z[0]+z[2],z[1]+z[4]],{"never_inline":True})

        [a,b] = g((y+x) @ z)
        f = ca.Function('f',[x,y,z,w],[ca.sin(a*(7*w)+6)*b],{"print_instructions": True, "print_canonical":True})
        ca.DM.rng(1)
        inputs = [ca.DM.rand(f.sparsity_in(i)) for i in range(f.n_in())]
        with capture_stdout() as out:
            f(*inputs)
        out = out[0].split("\n")[:-1]
        
        self.assertTrue("inputs" in out[0])
        
        
        if args.run_slow:
            with self.assertInException("Parsing error"):
                self.check_codegen(f,inputs,main=True)
                        
            with open("f_out.txt","r") as f_out:
                out2 = f_out.read().split("\n")[:-2]
            
            self.assertEqual(len(out),len(out2))
            for a,b in zip(out,out2):
                self.assertEqual(a,b)


    results = {}
    results_slow = {}
    for X in [ca.MX,ca.SX]:
        x = X.sym("x")
        y = X.sym("y")
        f = ca.Function('f',[x,y],[x+y],{"print_instructions": True, "print_canonical":True})
        f.disp(True)
        ca.DM.rng(1)
        inputs = [ca.DM.rand(f.sparsity_in(i)) for i in range(f.n_in())]
        with capture_stdout() as out:
            f(*inputs)
        out = out[0].split("\n")[:-1]
        results[X] = out
        
        self.assertTrue("inputs" in out[0])
        
        
        if args.run_slow:
            with self.assertInException("Parsing error"):
                self.check_codegen(f,inputs,main=True)
                        
            with open("f_out.txt","r") as f_out:
                out2 = f_out.read().split("\n")[:-2]
                results_slow[X] = out2
            
            
            self.assertEqual(len(out),len(out2))
            for a,b in zip(out,out2):
                self.assertEqual(a,b)

    self.assertEqual(len(results[ca.SX]),len(results[ca.MX]))
    for a,b in zip(results[ca.SX],results[ca.MX]):
        self.assertEqual(a,b)
    if args.run_slow:
        self.assertEqual(len(results_slow[ca.SX]),len(results_slow[ca.MX]))
        for a,b in zip(results_slow[ca.SX],results_slow[ca.MX]):
            self.assertEqual(a,b)

  def test_low(self):
    v = ca.MX.sym("v",5)
    p0 = [-1,0,0.2,1,1.3,2,2.5,3,3.5,4,4.5]
    p = ca.MX.sym("p",len(p0))


    f = ca.Function("f",[v,p],[ca.low(v, p)])

    inputs = [list(range(5)),p0]
    self.check_codegen(f,inputs=inputs)
    self.check_serialize(f,inputs=inputs)

    self.checkarray(f(*inputs),ca.DM([0, 0, 0, 1, 1, 2, 2, 3, 3, 3, 3]))

  def test_linear_interpn(self):
    N = 4
    n_dim=1
    nq = 3
    n_out = 2
    grid = []
    x0 = []
    for i in range(n_dim):
      g = sorted(list(np.random.random(N)))
      grid.append(g)
      x0.append(np.linspace(g[0],g[-1],nq+2)[1:-1])

    x0r = list(zip(*x0))

    D = np.random.random([n_out]+[N]*n_dim)

    d = D.ravel(order='F')

    strides = [n_out]
    for i in range(n_dim):
      strides.append(strides[-1]*N)

    F = ca.interpolant('F', 'linear', grid, d)

    ref = ca.vcat([F(x).T for x in x0r])

    x = [ca.MX.sym("x%d" % i,1,N) for i in range(n_dim)]
    v = ca.MX.sym("v",N**n_dim*n_out)
    xq = [ca.MX.sym("xq%d" % i,nq) for i in range(n_dim)]

    e = ca.MX.interpn_linear(x,v,xq)

    f = ca.Function('f',[ca.hcat(x),v,ca.vcat(xq)],[e])
    inputs = [ca.hcat([ca.hcat(e) for e in grid]),d,ca.vcat(x0)]
    r = f(*inputs)
    self.checkarray(r, ref)
    self.check_codegen(f,inputs=inputs,std="c99")
    self.check_serialize(f,inputs=inputs)

  def test_paramslice(self):
    N = 4
    M = 6
    A = ca.DM.rand(N,M)
    As = ca.MX.sym("A", A.shape)
    ASs = ca.MX.sym("AS",ca.tril(As.sparsity()))
    print(ASs.sparsity().spy())
    vr=list(range(N))
    vc=list(range(M))
    vk=list(range(M*N))
    for r in [1, slice(1,3), slice(0,N,2), slice(None)]:
      rs = ca.MX.sym("r",A[r,:].shape[0])
      for c in [1, slice(1,3), slice(0,M,2), slice(None)]:
        cs = ca.MX.sym("c",A[:,c].shape[1])
        for rr in [r, rs]:
          for cc in [c, cs]:
            f = ca.Function('f',[As,rs,cs],[As[rr,cc]])
            self.checkarray(f(A,vr[r],vc[c]),A[r,c])
            if isinstance(rr,ca.MX) or isinstance(cc,ca.MX):
              with self.assertInException("Parametric slicing only supported for dense matrices"):
                ASs[rr,cc]
    for k in [1, slice(1,3), slice(0,N,2), slice(0,2*N,2), slice(None)]:
      ks = ca.MX.sym("r",A[k].shape[0])
      for kk in [k, ks]:
        f = ca.Function('f',[As,ks],[As[kk]])
        self.checkarray(f(A,vk[k]),A[k])
        with self.assertInException("Parametric slicing only supported for dense matrices"):
          ASs[ks]

  def test_mapping(self):
    A = ca.MX.sym("A",4,4)
    i = ca.DM([[0,3],[1,2]])
    self.checkarray(i,A[i].mapping())


  @memory_heavy()
  def test_convexify(self):
    A = ca.diagcat(1,2,-1,ca.blockcat([[1.2,1.3],[1.3,4]]),ca.sparsify(ca.blockcat([[0,1,0],[1,4,7],[0,7,9]])),ca.DM(2,2))

    margin = 1e-7

    np.random.seed(0)
    p = np.random.permutation(A.shape[0])


    [D,V] = np.linalg.eigh(np.array(A))
    Dr= ca.fmax(abs(D),1e-7)
    Dc= ca.fmax(D,1e-7)  # pyright: ignore[reportCallIssue,reportArgumentType]

    Ar_ref = V @ ca.diag(Dr) @ V.T  # pyright: ignore[reportCallIssue,reportArgumentType]
    Ac_ref = V @ ca.diag(Dc) @ V.T  # pyright: ignore[reportCallIssue,reportArgumentType]
    As = ca.MX.sym("As",A.sparsity())

    for opts,ref in [({"strategy":"eigen-reflect"},Ar_ref), ({"strategy":"eigen-clip"},Ac_ref), ({"strategy":"regularize"},A+(4+margin)*ca.DM.eye(A.shape[0]))]:


      ops = [lambda e: e, lambda e: ca.triu(e), lambda e: ca.tril(e),
                 lambda e: e[p,p], lambda e: ca.triu(e[p,p]), lambda e: ca.tril(e[p,p])]
      if "regularize" in str(opts): ops = [lambda e: e, lambda e: e[p,p]]
      for op in ops:

        f= ca.Function("f",[As],[ca.convexify(op(A),opts)])

        self.checkarray(f(A),op(ref),digits=8)

        self.check_serialize(f,inputs=[A])

        self.check_codegen(f,inputs=[A])

  def test_convexify_bugs(self):
    for eps in [0,1e-100,1]:
      q = 4
      Q = 10
      for s in [-4,0.1,0,0.1,4]:
        for trans in [lambda e: e, lambda e: ca.sparsify(e)]:
          A = trans(ca.blockcat([[q,s,0+eps],[s,Q,0],[0+eps,0,2]]))
          print(A)
          self.assertTrue(np.all(np.linalg.eigh(A)[0]>0))
          Ac = ca.evalf(ca.convexify(A,{"strategy":"eigen-reflect"}))
          self.checkarray(A,Ac,digits=8)


  def test_logsumexp(self):
    x = ca.MX.sym("x",3)

    f_ref = ca.Function("f_ref",[x],[ca.log(ca.exp(x[0])+ca.exp(x[1])+ca.exp(x[2]))])
    f = ca.Function("f",[x],[ca.logsumexp(x)])

    self.checkfunction(f,f_ref,inputs=[ca.vertcat(1.1,1.3,1.7)])
    self.check_codegen(f,inputs=[ca.vertcat(1.1,1.3,1.7)])
    self.checkfunction(f,f_ref,inputs=[ca.vertcat(1.1,1.3,1.3)])
    self.checkfunction(f,f_ref,inputs=[ca.vertcat(1.3,1.3,1.3)])
    

    # Avoid overflow
    res = f(ca.vertcat(100,1000,10000))
    self.checkarray(res,10000)

  def test_empty_broadcast(self):
    for nc in [0,2]:
      res = ca.atan2(ca.MX.sym("c",nc,1),ca.MX.sym("t",nc,3))

      self.assertEqual(res.shape[0],nc)
      self.assertEqual(res.shape[1],3)
  
      with self.assertInException("Dimension mismatch"):
        res = ca.atan2(ca.MX.sym("c",nc,2),ca.MX.sym("t",nc,3))

      res = ca.atan2(ca.MX.sym("c",nc,2),ca.MX.sym("t",nc,4))
      self.assertEqual(res.shape[0],nc)
      self.assertEqual(res.shape[1],4)

      res = ca.atan2(ca.MX.sym("c",nc,4),ca.MX.sym("t",nc,2))
      self.assertEqual(res.shape[0],nc)
      self.assertEqual(res.shape[1],4)

  def test_horzcat_sparsity(self):
    x = ca.MX.sym("x",2,2)

    for y in [ca.MX.sym("y",ca.Sparsity(2,2)), ca.MX.sym("y",ca.Sparsity.lower(2))]:

      z = ca.sin(ca.horzcat(x,ca.sin(y),x))
      z_alt = ca.sin(ca.sparsity_cast(ca.vertcat(ca.vec(x),ca.vec(ca.sin(y)),ca.vec(x)),z.sparsity()))

      f = ca.Function('f',[x,y],[z])
      f_alt = ca.Function('f_alt',[x,y],[z])

      for F in [f,f_alt]:
        self.assertEqual(F.call([x,y],True,False)[0].nnz(), z.nnz())
        self.assertEqual(F.call([x,y],False,True)[0].nnz(), z.nnz())
        self.assertEqual(F.call([x,x**2],True,False)[0].nnz(), z.nnz())
        self.assertEqual(F.call([x,x**2],False,True)[0].nnz(), z.nnz())

    self.checkfunction(f,f.expand(),inputs=[ca.DM([[1,2],[3,4]]),ca.DM([[5,6],[7,8]])])
    self.checkfunction(f,f_alt,inputs=[ca.DM([[1,2],[3,4]]),ca.DM([[5,6],[7,8]])])
    self.check_codegen(f,inputs=[ca.DM([[1,2],[3,4]]),ca.DM([[5,6],[7,8]])])

  def test_sparsity_cast_ad(self):
    #issue 3164

    y = ca.MX.sym("y",3)

    xy = ca.vertcat(y[0],y[1])
    w = ca.MX(ca.sparsify(ca.DM([1,1,0])).sparsity().T,xy) @ y
    
    f = ca.Function('f',[y],[w])
    self.checkfunction(f,f.expand(), inputs=[ca.vertcat(1.1,1.3,1.7)])
    self.check_codegen(f,inputs=[ca.vertcat(1.1,1.3,1.7)])


  def test_fractional_slicing(self):

    t = ca.MX.sym("x")

    f = ca.Function('f',[t],[ca.MX(ca.DM([1,2]))[t]])

    self.assertEqual(f(0),1)
    self.assertEqual(f(1),2)
    self.assertEqual(f(0.5),1)
    self.assertEqual(f(0.9),1)
    
    f = ca.Function('f',[t],[ca.MX(ca.DM([[1,2],[3,4]]))[:,t]])

    self.checkarray(f(0),ca.DM([1,3]))
    self.checkarray(f(1),ca.DM([2,4]))
    self.checkarray(f(0.5),ca.DM([1,3]))
    self.checkarray(f(0.9),ca.DM([1,3]))
    
  def test_primitives(self):
  
    def experiments():
        x = ca.MX.sym("x",3)
        y = ca.MX.sym("y",3)
        yield ca.vertcat(x,y)
        yield ca.horzcat(x,y)
        yield ca.diagcat(x,y)
        yield ca.vertcat(x,y,ca.MX(0,1))
        yield ca.vertcat(x,y).reshape((2,3))
        yield ca.sparsity_cast(ca.vertcat(x,y), ca.Sparsity.diag(6))
        
    for z in experiments():
        print(z)
        self.assertTrue(z.is_valid_input())
        
        print(z.primitives())
        f = ca.Function('f',[z],z.primitives())
        ca.DM.rng(1)
        a = ca.DM.rand(z.sparsity())
        
        ref = f(a)
        res = z.split_primitives(a)
        
        self.assertEqual(len(ref),len(res))  # pyright: ignore[reportArgumentType]
        for ea,eb in zip(ref,res):
            self.checkarray(ea, eb)
        
        a_back = z.join_primitives(res)
        self.checkarray(a,a_back)
        
        zsx = ca.SX.sym("x",z.sparsity())
        f2 = ca.Function('f',[zsx],z.split_primitives(zsx))
        
        self.checkfunction_light(f,f2,inputs=[a])
        self.check_codegen(f,inputs=[a])
        self.check_codegen(f2,inputs=[a])
        
        f3 = ca.Function('f',[zsx],[z.join_primitives(z.split_primitives(zsx))])
        
        res = f3(a)
        self.checkarray(a,res)
        
        
  def test_cache_output_node(self):
    x = ca.MX.sym("x")
    y = ca.MX.sym("y")
    f = ca.Function("f",[x,y],[x+y,x/y])
    
    [a,b] = f(ca.sin(x),ca.cos(y))
    base = a.dep(0)
    self.assertTrue(hash(a),hash(base.get_output(0)))
    self.assertTrue(hash(b),hash(base.get_output(1)))
    
    
  def test_naked_call_function(self):
    x = ca.MX.sym("x")
    y = ca.MX.sym("y")
    f = ca.Function("f",[x,y],[x+y,x/y])
    
    [a,b] = f(ca.sin(x),ca.cos(y))
    
    

    base = a.dep(0)

    with self.assertInException("get_output"):
        g = ca.Function("g",[x,y],[base])

    with self.assertInException("get_output"):
        3*base
        
  def test_substitute(self):
    x = ca.MX.sym("x")
    y = ca.MX.sym("y")
    w = ca.MX.sym("w")
    

    xy = x*y
    
    z = xy+w
    z2 = ca.substitute([z],[w],[3])[0]
    print(xy)
    print(z2.dep(1))
    self.assertNotEqual(hash(xy),hash(z2.dep(1)))

  def test_extract_parametric(self):
    def test_equal(a,b):
        if not isinstance(a,list) and not isinstance(a,tuple):
          a = [a]
        if not isinstance(b,list) and not isinstance(b,tuple):
          b = [b]
        f1 = ca.Function('f',[x,y,p],a)
        f2 = ca.Function('f',[x,y,p],b)
        ca.DM.rng(1)
        args = [ca.DM.rand(f1.sparsity_in(i)) for i in range(f1.n_in())]
        
        for a,b in zip(f1.call(args,False,False),f2.call(args,False,False)):
           if np.linalg.norm(a-b,1)>1e-12:
              return False
        return True

        #return cse(a-b).is_zero()  
    
    for X in [ca.SX,ca.MX]:


      x = X.sym("x")
      y = X.sym("y")
      p = X.sym("p")
      
      expr = 2*x*p+2*y*p**2
      expr_ret,symbols,parametric = ca.extract_parametric(expr,p, {"extract_trivial": True})
      print(expr_ret,parametric)
      self.assertTrue(test_equal(parametric,[p,p**2]))
      
      expr = ca.vertcat(2*x+3*y+8*p,3*x+2*y+4,7)
      
      expr_ret,symbols,parametric = ca.extract_parametric(expr,p, {"extract_trivial": True})
      
      print(expr_ret,parametric)
      
      self.assertTrue(test_equal(parametric[0],8*p))
      
      

        
    for X in [ca.SX,ca.MX]:


      x = X.sym("x")
      y = X.sym("y")
      p = X.sym("p")

      def tests():
        cp = ca.cos(p)
        expr = (x-ca.sqrt(p))*(y+cp**2)+cp/(x+cp)/(y-p)
        
        yield expr
        
        yield ca.vertcat(expr,p**2)
        
      for expr in tests():

        expr_ret,symbols,parametric = ca.extract_parametric(expr,p, {"extract_trivial": True})
        
        print(expr_ret,symbols,parametric)
        
        self.assertFalse(ca.depends_on(expr_ret,p))
        
        expr_recreated = ca.substitute([expr_ret],symbols,parametric)[0]
        self.assertTrue(test_equal(expr,expr_recreated))

  def test_separate_linear(self):
    for X in [ca.SX,ca.MX]:
        x = X.sym("x")
        y = X.sym("y")
        p = X.sym("p")

        z = 3*x+5*y+x*y+p*x+ca.cos(p)


        for expr, ref_const, ref_lin, ref_nonlin in [
                    [x,0 , x, 0],
                    [x+p, p, x, 0],
                    [x+y, 0, x+y, 0],
                    [x-p, -p, x, 0],
                    [p-x, p, -x, 0],
                    [2*x, 0, 2*x, 0],
                    [-x, 0, -x, 0],  
                    [3*x+y, 0, 3*x+y, 0],
                    [ca.cos(p),ca.cos(p),0,0],
                    [3*x+y+ca.cos(p), ca.cos(p), 3*x+y, 0],
                    [3*x+7, 7, 3*x, 0],
                    [p*(3*x+7), p*7, p*(3*x),0],
                    [y*(3*x+7), 0, 7*y, y*(3*x)],
                    [ca.cos(7+3*x), 0, 0, ca.cos(7+3*x)],
                    [ca.cos(3*x+7), 0, 0, ca.cos(3*x+7)],
                    [ca.cos(7*p),ca.cos(7*p),0,0],
                    [ca.cos(p+x),0,0,ca.cos(p+x)],
                    [x*ca.cos(p)+7,7,x*ca.cos(p),0],
                    [y*(x*ca.cos(p)+7),0,7*y,y*x*ca.cos(p)],
                    [p*(x+y+x*y), 0, p*(x+y), p*(x*y)],
                    [ca.cos(x+y+x*y), 0, 0, ca.cos(x+y+x*y)],
                    [p*(p+ca.cos(p)*x), p**2, p*ca.cos(p)*x, 0],
                    [(p+x+x*y)*(ca.cos(p)+7*x+3*x*y), p*ca.cos(p), p*(7*x)+x*ca.cos(p),p*(3*x*y)+x*y*(ca.cos(p)+7*x+3*x*y)+7*x**2++3*x**2*y ],
                    [(p+x+x*y)/ca.cos(p), p/ca.cos(p), x/ca.cos(p), x*y/ca.cos(p) ],
                    [(p+x+x*y)/(ca.cos(p)+7*x), 0,0,(p+x+x*y)/(ca.cos(p)+7*x) ],
                    [(p+x+x*y)/(ca.cos(p)+7*x+3*x*y), 0,0,(p+x+x*y)/(ca.cos(p)+7*x+3*x*y) ]  
                    ]:

            print("ref",expr,ref_const,ref_lin,ref_nonlin)
            
            [expr_const,expr_lin,expr_nonlin] = ca.separate_linear(expr,ca.vertcat(x,y),p)
            print("actual",expr_const,expr_lin,expr_nonlin)
            
            
            def test_equal(a,b):
                f1 = ca.Function('f',[x,y,p],[a])
                f2 = ca.Function('f',[x,y,p],[b])
                ca.DM.rng(1)
                args = [ca.DM.rand(f1.sparsity_in(i)) for i in range(f1.n_in())]
                return abs(f1(*args)-f2(*args))<1e-12
                
                #return cse(a-b).is_zero() 

            assert test_equal(expr,expr_const+expr_lin+expr_nonlin)

            assert test_equal(expr_const,ref_const)
            assert test_equal(expr_lin,ref_lin)
            assert test_equal(expr_nonlin,ref_nonlin)
            
            
            assert ca.is_linear(expr_lin,ca.vertcat(x,y))
            assert not ca.depends_on(expr_const,ca.vertcat(x,y))
            
        x = X.sym("x",2)
        y = X.sym("y",2)
        p = X.sym("p",2)

        for expr, ref_const, ref_lin, ref_nonlin in [
                    [ca.dot(p,x),0 , ca.dot(p,x), 0],
                    [ca.dot(p+x+x*y,ca.cos(p)+7*x+3*x*y), ca.dot(p,ca.cos(p)), ca.dot(p,7*x)+ca.dot(x,ca.cos(p)), ca.dot(p,3*x*y)+ca.dot(x,7*x+3*x*y)+ca.dot(x*y,ca.cos(p)+7*x+3*x*y)]
                    ]:

            print("ref",expr,ref_const,ref_lin,ref_nonlin)
            
            [expr_const,expr_lin,expr_nonlin] = ca.separate_linear(expr,ca.vertcat(x,y),p)
            print("actual",expr_const,expr_lin,expr_nonlin)
            
            
            def test_equal(a,b):
                f1 = ca.Function('f',[x,y,p],[a])
                f2 = ca.Function('f',[x,y,p],[b])
                ca.DM.rng(1)
                args = [ca.DM.rand(f1.sparsity_in(i)) for i in range(f1.n_in())]
                return abs(f1(*args)-f2(*args))<1e-12
                
                #return cse(a-b).is_zero() 

            assert test_equal(expr,expr_const+expr_lin+expr_nonlin)

            assert test_equal(expr_const,ref_const)
            assert test_equal(expr_lin,ref_lin)
            assert test_equal(expr_nonlin,ref_nonlin)


            assert ca.is_linear(expr_lin,ca.vertcat(x,y))
            assert not ca.depends_on(expr_const,ca.vertcat(x,y))

        x = X.sym("x",2,2)
        y = X.sym("y",2,2)
        p = X.sym("p",2,2)

        for expr, ref_const, ref_lin, ref_nonlin in [
                    [(p+x+x*y).T,p.T,x.T,(x*y).T],
                    [p @ x,0 , p @ x, 0],
                    [(p+x+x*y) @ (ca.cos(p)+7*x+3*x*y), p @ ca.cos(p), p @ (7*x)+x @ ca.cos(p), p @ (3*x*y)+x @ (7*x+3*x*y)+x*y @ (ca.cos(p)+7*x+3*x*y)]
                    ]:

            print("ref",expr,ref_const,ref_lin,ref_nonlin)
            
            [expr_const,expr_lin,expr_nonlin] = ca.separate_linear(expr,ca.veccat(x,y),p)
            print("actual",expr_const,expr_lin,expr_nonlin)
            
            
            def test_equal(a,b):
                f1 = ca.Function('f',[x,y,p],[a])
                f2 = ca.Function('f',[x,y,p],[b])
                ca.DM.rng(1)
                args = [ca.DM.rand(f1.sparsity_in(i)) for i in range(f1.n_in())]
                return np.linalg.norm(f1(*args)-f2(*args),1)<1e-12
                
                #return cse(a-b).is_zero() 

            assert test_equal(expr,expr_const+expr_lin+expr_nonlin)

            assert test_equal(expr_const,ref_const)
            assert test_equal(expr_lin,ref_lin)
            assert test_equal(expr_nonlin,ref_nonlin)


            assert ca.is_linear(expr_lin,ca.veccat(x,y))
            assert not ca.depends_on(expr_const,ca.veccat(x,y))

  def test_separate_linear_mat(self):

        x = ca.MX.sym("x",3)
        y = ca.MX.sym("y",3,3)
        z = ca.MX.sym("z",2)
        p = ca.MX.sym("p")
        
        A = ca.DM([[0,1,3],[2,0,0],[0,9,0]])
        Asp = ca.sparsify(A)
        

        B = ca.DM([[1.7,1.1],[6,0],[0,1.2]])
        Bsp = ca.sparsify(B)
        [p1,p2] = ca.vertsplit(9*ca.vertcat(3,3*x,y @ x),[0,2,7])
        [p1c,p2c] = ca.vertsplit(9*ca.vertcat(3,ca.DM.zeros(3,1),ca.DM.zeros(3,1)),[0,2,7])
        [p1l,p2l] = ca.vertsplit(9*ca.vertcat(0,3*x,ca.DM.zeros(3,1)),[0,2,7])
        [p1n,p2n] = ca.vertsplit(9*ca.vertcat(0,ca.DM.zeros(3,1),y @ x),[0,2,7])
        E = (9*ca.vertcat(3,3*x,y @ x)).nz[[2,3,0,4,6,6]]
        Ec = (9*ca.vertcat(3,ca.DM.zeros(3,1),ca.DM.zeros(3,1))).nz[[2,3,0,4,6,6]]
        El = (9*ca.vertcat(0,3*x,ca.DM.zeros(3,1))).nz[[2,3,0,4,6,6]]
        En = (9*ca.vertcat(0,ca.DM.zeros(3,1),y @ x)).nz[[2,3,0,4,6,6]]

        
        for expr, ref_const, ref_lin, ref_nonlin in [
                    [A @ x+B @ z,ca.DM.zeros(3,1) , A @ x, B @ z],
                    [Asp @ x+Bsp @ z,ca.DM.zeros(3,1) , Asp @ x, Bsp @ z],
                    [ca.vertcat(3*x+y @ x,9*z),ca.DM.zeros(5,1),ca.vertcat(3*x,ca.DM.zeros(2,1)), ca.vertcat(y @ x,9*z) ],
                    [7*ca.vertcat(3*x+y @ x,9*z),ca.DM.zeros(5,1),7*ca.vertcat(3*x,ca.DM.zeros(2,1)), 7*ca.vertcat(y @ x,9*z) ],
                    [9*p2, 9*p2c, 9*p2l, 9*p2n],
                    [E,Ec,El,En]
                    ]:

            print("separate_linear(", expr, ")")
            print("ref")
            print("  const",ref_const)
            print("  lin",ref_lin)
            print("  nonlin",ref_nonlin)
            
            [expr_const,expr_lin,expr_nonlin] = ca.separate_linear(expr,[x,y],[p])
            print("actual")
            print("  const",expr_const)
            print("  lin",expr_lin)
            print("  nonlin",expr_nonlin)
            
            def test_equal(a,b):
                f1 = ca.Function('f',[x,y,z,p],[a])
                f2 = ca.Function('f',[x,y,z,p],[b])
                ca.DM.rng(1)
                args = [ca.DM.rand(f1.sparsity_in(i)) for i in range(f1.n_in())]
                self.checkarray(f1(*args),f2(*args))

            test_equal(expr,expr_const+expr_lin+expr_nonlin)

            test_equal(expr_const,ref_const)
            test_equal(expr_lin,ref_lin)
            test_equal(expr_nonlin,ref_nonlin)
            
            
            assert ca.is_linear(expr_lin,ca.veccat(x,y))
            assert not ca.depends_on(expr_const,ca.veccat(x,y))
            

  def test_extract_parametric_opts(self):
      for X in [ca.SX,ca.MX]:
          x = X.sym("x")
          y = X.sym("y")
          p = X.sym("p")
          
          expr = 2*x*p+2*y*p**2
          expr_ret,symbols,parametric = ca.extract_parametric(expr,p, {"extract_trivial": True})
          
          self.assertEqual(len(symbols),2)
          
          self.assertEqual(symbols[0].name(),"e_0")
          self.assertEqual(symbols[1].name(),"e_1")
          
          expr = 2*x*p+2*y*p**2
          expr_ret,symbols,parametric = ca.extract_parametric(expr,p, {"offset": 5, "prefix": "foo_", "suffix": "bar", "extract_trivial": True})
          
          self.assertEqual(symbols[0].name(),"foo_5bar")
          self.assertEqual(symbols[1].name(),"foo_6bar")
          
          print(symbols)
          print(parametric)
          
          
          expr = 2*x*p+2*y*p**2
          expr_ret,symbols,parametric = ca.extract_parametric(expr,p, {"extract_trivial": False})
          self.assertEqual(len(symbols),1)
          self.assertTrue("sq(p)" in str(parametric))
          print(parametric)
      
          expr = 2*x*p+2*y*p**2
          expr_ret,symbols,parametric = ca.extract_parametric(expr,p, {"extract_trivial": True})
          self.assertEqual(len(symbols),2)
          self.assertTrue("sq(p)" in str(parametric))
          print(parametric)

  def test_weakref(self):
    x = ca.MX.sym("x")
    
    y = ca.WeakRef(x)
    
    self.assertTrue(y.alive())
    
    x = 0
    
    import gc
    gc.collect()
    
    self.assertFalse(y.alive())
    
    
    x = ca.MX.sym("x")
    y = ca.WeakRef(x)
    xx = y.shared()
    xx = 0
    
    import gc
    gc.collect()
    
    self.assertTrue(y.alive())
    
  def test_nonzeros(self):
    
    A = ca.MX.sym("A",3,3)

    x = ca.MX.sym("x")
    y = ca.MX.sym("y")
    z = ca.MX.sym("z")
    T = ca.MX.sym("T",ca.Sparsity.upper(3))
    w = ca.MX.sym("w",3,3)
    
    ca.DM.rng(1)
    
    for expr in [
        ca.vertcat(x,y,z),
        ca.veccat(x,y,z,w),
        ca.veccat(x,T,z,w),
        2*x,
        4*T*x,
        ca.vertcat(x,ca.vertcat(y,z))
    ]:
        fref = ca.Function('fref',[x,y,z,w,A,T],[expr.nz[:]])
        f = ca.Function('f',[x,y,z,w,A,T],[ca.vcat(expr.nonzeros())])
        
        self.checkfunction_light(f,fref,inputs=[ca.DM.rand(f.sparsity_in(i)) for i in range(f.n_in())])

    for expr in [
        ca.vertcat(x,y,z),
        ca.vertcat(x,ca.vertcat(y,z))
    ]:
        s = str(expr.nonzeros())
        self.assertTrue("[MX(x), MX(y), MX(z)]" in s)
        
  def test_printme_codegen(self):
    for X in [ca.SX,ca.MX]:
        x = X.sym("x")
        
        f = ca.Function("f",[x],[(x**2).printme(2)+6])
        self.check_codegen(f,inputs=[3])
        
  @memory_heavy()
  def test_norms(self):
  
    for f in [lambda x,X: ca.vertcat(x[0],X(1,1),x[1:]),
              lambda x,X: ca.diag(x)
                ]:
        for norm in [ca.norm_fro, ca.norm_1, ca.norm_inf]:
            print("norm",norm)
            x = ca.SX.sym("x",5)
            e = ca.log(norm(f(ca.sin(x),ca.SX)))
            fsx = ca.Function('f',[x],[e])    
            
            x = ca.MX.sym("x",5)
            e = ca.log(norm(f(ca.sin(x),ca.MX)))
            fmx = ca.Function('f',[x],[e])
            for args in [ca.vertcat(-0.3,0.3,0.21,0.17,0),ca.vertcat(-1,-2,-7,-7,-7)]:
                self.checkfunction(fsx,fmx,inputs=[args])
            args = [ca.vertcat(-0.3,0.3,0.21,0.17,0)]
            self.check_codegen(fmx,inputs=args)
      
  def test_contains(self):
    x = ca.MX.sym("x",2,2)
    y = ca.MX.sym("y")
    z = ca.MX.sym("z",3,3)
    
    e = x*y
    self.assertTrue(ca.contains([x,y,z],x))
    self.assertFalse(ca.contains([x,y],z))
    self.assertTrue(ca.contains_any([x,y],[y,z]))
    self.assertFalse(ca.contains_all([x,y],[y,z]))
    self.assertTrue(ca.contains_any([x,y],[x,y]))
    self.assertTrue(ca.contains_all([x,y],[x,y]))
    self.assertTrue(ca.contains([e,x],e))
   
  def test_issue4041(self):
    
    x = ca.MX.sym("x")
    y = ca.MX.sym("y")

    f = ca.Function('f',[x,y],[x**2],{"never_inline":True})

    expr = f(x,y)
    
    self.assertFalse(ca.depends_on(expr,y))
    self.assertTrue(ca.contains(ca.symvar(expr),y))
 
    f = ca.Function('f',[x,y],[x**2],{"always_inline":True})

    expr = f(x,y)
    
    self.assertFalse(ca.depends_on(expr,y))
    self.assertFalse(ca.contains(ca.symvar(expr),y))
      
  def test_symmetric_jacobian(self):
    x = ca.MX.sym('x', 4)
    u = ca.MX.sym('u', 2)

    h = x[1]*5*x[2]

    nh = 1
    lam_h = ca.MX.sym('lam_h', nh, 1)

    adj_ux = ca.densify(ca.jtimes(h, ca.vertcat(u, x), lam_h, True))
    hess_ux = ca.jacobian(adj_ux, ca.vertcat(u, x), {'symmetric': True})
    
    adj_ux = ca.jtimes(h, ca.vertcat(u, x), lam_h, True)
    with self.assertInException("Symmetry exploitation"):
        hess_ux = ca.jacobian(adj_ux, ca.vertcat(u, x), {'symmetric': True})
        
        
  def test_pow_zero(self):
  
    for X in [ca.SX,ca.MX]:

        x = X.sym('x', ca.Sparsity.diag(4))
        self.assertTrue((x**0).is_dense())
        self.assertFalse((x**0.3).is_dense())
        self.assertFalse((x**1).is_dense())
        self.assertFalse((x**2).is_dense())
        self.assertFalse((x**2.3).is_dense())
        self.assertTrue((x**(-2)).is_dense())
        p = X.sym("p")
        
        self.assertTrue((x**p).is_dense())


        x=X.sym("x")

        self.checkarray(ca.evalf(ca.substitute(ca.vertcat(x**2,x**1,x**0),x,0)),ca.vertcat(0,0,1))
        self.checkarray(ca.evalf(ca.substitute(x**ca.DM([2,1,0]),x,0)),ca.vertcat(0,0,1))

  def test_pow(self):
    x = ca.MX.sym("x")
    r = ca.evalf(ca.substitute(x**ca.vertcat(2,1,0),x,0))
    self.checkarray(r,ca.vertcat(0,0,1))
    
    y = ca.MX(4,1)
    self.assertEqual((y**0).nnz(),4)
    y = ca.MX(4,1)
    y[1] = x
    self.assertEqual((y**0).nnz(),4)
    y = ca.MX(4,1)
    y[1] = 0
    self.assertEqual((y**0).nnz(),4)
    y = ca.MX(4,1)
    y[1] = 1
    self.assertEqual((y**0).nnz(),4)

  def test_sum(self):
    test_cases = [ca.DM(4), ca.horzcat(1,2),ca.vertcat(1,2),ca.blockcat([[1,2],[3,4]]),ca.blockcat([[1,2,3],[3,4,5]]),ca.blockcat([[1,2],[3,4],[6,9]])]

    # np.sum on casadi types now stays in the casadi type system: the
    # axis-aware result is 2-D (matching casadi's column convention)
    # whereas numpy collapses to 1-D.  Reshape numpy's reference to
    # match the casadi shape for comparison.
    for e in test_cases:
        ref_total = np.sum(np.array(e))
        self.checkarray(np.array(np.sum(e)).reshape(()), ref_total)
        # axis=0: casadi shape (1, n);  numpy shape (n,)
        res = np.array(np.sum(e, 0))
        ref = np.sum(np.array(e), 0).reshape(res.shape)
        self.checkarray(res, ref)
        # axis=1: casadi shape (m, 1);  numpy shape (m,)
        res = np.array(np.sum(e, 1))
        ref = np.sum(np.array(e), 1).reshape(res.shape)
        self.checkarray(res, ref)

        with self.assertInException("axis 2 is out of bound"):
            np.sum(e,2)

  def test_jtimes_empty(self):
  
    with self.assertInException("Ambiguous"):
        x = ca.SX.sym("x",2,0)
        y = ca.jtimes(5,x,ca.DM.ones(2,5))
    with self.assertInException("Ambiguous"):
        x = ca.SX.sym("x",2,3)
        y = ca.jtimes(ca.DM.ones(4,0),x,ca.DM.ones(4,6),True)
    
    # non-transpose
    for X in [ca.SX,ca.MX]:
        for nx in [2,1,0]:
            for mx in [3,1,0]:
                x = X.sym("x",nx,mx)
                for ny in [3,1,0]:
                    for my in [4,1,0]:
                        for bm in [1,3,0]:
                            if bm!=1 and my==0: continue # ambiguous construction
                            y = ca.DM.rand(ny,nx) @  x @ ca.DM.rand(mx,my)
                            assert y.shape==(ny,my)
                            yb = ca.DM.rand(ny,my*bm)
                            
                            r = ca.jtimes(y,x,yb,True)
                            r_ref = ca.hcat([(ca.jacobian(y,ca.vec(x)).T @ ca.vec(e)).reshape(x.shape) for e in ca.horzsplit(yb,max(1,my))])
                            self.assertTrue(r.shape==(nx,mx*bm))
                            f = ca.Function('f',[x],[r])
                            f_ref = ca.Function('f',[x],[r_ref])
                            self.checkfunction_light(f,f_ref,inputs=[ca.DM.rand(nx,mx)])
        for nx in [2,1,0]:
            for mx in [3,1,0]:
                x = X.sym("x",nx,mx)
                for ny in [3,1,0]:
                    for my in [4,1,0]:
                        for bm in [1,3,0]:
                            if bm!=1 and mx==0: continue # ambiguous construction
                            y = ca.DM.rand(ny,nx) @  x @ ca.DM.rand(mx,my)
                            assert y.shape==(ny,my)
                            dx = ca.DM.rand(nx,mx*bm)

                            r = ca.jtimes(y,x,dx)
                            r_ref = ca.hcat([(ca.jacobian(y,ca.vec(x)) @ ca.vec(e)).reshape(y.shape) for e in ca.horzsplit(dx,max(1,mx))])
                            
                            self.assertTrue(r.shape==(ny,my*bm))
                            f = ca.Function('f',[x],[r])
                            f_ref = ca.Function('f',[x],[r_ref])
                            self.checkfunction_light(f,f_ref,inputs=[ca.DM.rand(nx,mx)])

  def test_issue2934(self):
    M = ca.MX.sym("X",5,5)
    Y = ca.MX.sym("y",2,2)

    args = [M,Y]

    m = ca.MX(M)


    m[[0,0],[0,0]] = Y

    print(m)
    f = ca.Function('f',[M,Y],[m])
    print(f.call([M,Y],True,False)[0])
    self.assertTrue("((2.*X)[-1, -1, -1, 0] = y)" in str(f.call([2*M,Y],True,False)[0]))

  def test_issue2935(self):


    M = ca.MX(ca.DM.rand(5,5))
    S = M.sparsity()
    m = M

    Y = ca.MX.sym("y",2,2)
    E = Y.sparsity()
    y = Y

    m[[0,0],[0,0]] = y

    e = ca.cos(m)
    e = ca.dot(e,e)

    eb = ca.MX.sym("eb")
    ge = ca.jtimes(e,Y,eb,True)

    f = ca.Function('f',[eb,Y],[ge])

    e1 = f.call([eb,Y],True,False)[0] # without eval_mx
    e2 = f.call([eb,Y+1e-300],True,False)[0] # with eval_mx

    print(e1)
    print(e2)

    f1 = ca.Function('f1',[eb,Y],[ge])
    f2 = ca.Function('f2',[eb,Y],[e2])

    args = [1,numpy.random.random(E.shape)]

    print(f1(*args))
    print(f2(*args))

    assert f1(*args).sparsity()==f2(*args).sparsity()
          
                    
  def test_extract(self):
    x = ca.MX.sym('x',2)
    y = ca.MX.sym('y',2)
    
    u = ca.vertcat(x,y)
    ca.DM.rng(1)
    u0 = ca.DM.rand(4,1)
    f = ca.Function('f',[x,y],[ca.sin(x+y)],{"never_inline":True})
    
    g = ca.Function('g',[x],[ca.cos(3*ca.sum(x))],{"never_inline":True})
    
    e = g(f(x/y,2*x*y)*x)/(x+y)

    res = ca.extract([e],{"lift_shared":False,"lift_calls":True})
    
    [vexpr,v,vdef] = res
    
    print(vexpr)
    print(v)
    print(vdef)

    self.assertEqual(len(v), 5)
    self.assertEqual(len(vdef), len(v))
    self.assertEqual(len(vexpr), 1)
    #assert()
    # Check inverted operation
    f_ref = ca.Function('f',[u],[e])
    print("e",e)
    e_reconstructed = ca.substitute_inplace(v, vdef, [e])[1]
    print(e_reconstructed)
    f = ca.Function('f',[u],e_reconstructed)
    self.checkarray(f_ref(u0),f(u0))
    
    print(e)
    
    print(ca.substitute_inplace(v, vdef, vexpr))
    
    vexpr = ca.vcat(vexpr)
    v = ca.vcat(v)
    vdef = ca.vcat(vdef)
    
    
    Jf_ref = ca.Function('J_ref',[u],[ca.jacobian(e,u)])
    
    J_v_u = ca.solve(ca.DM.eye(vdef.size1())-ca.jacobian(vdef,v),ca.jacobian(vdef,u))
    J = ca.jacobian(vexpr,u) + ca.jacobian(vexpr,v) @  J_v_u

    [vexpr,v,vdef] = res
    
    J = ca.substitute_inplace(v, vdef, [J])[1][0]
    
    Jf = ca.Function('J2',[u],[J])
    

    self.checkarray(Jf_ref(u0),Jf(u0))

  def test_const_folding_on_the_fly(self):
    x = ca.MX.sym('x')

    for a in [2,7]:
        for b in [2,7]:
            ref = (a*b)*x
            if a==7 and b==7:
                # Should probably not happen, but it is present in 3.7.0
                self.assertEqual(str(ref),str(a*(b*x)))
                self.assertEqual(str(ref),str((b*x)*a))
            else:
                # Dangerous: no caching of constants
                self.assertNotEqual(str(ref),str(a*(b*x)))
                self.assertNotEqual(str(ref),str((b*x)*a))

  def test_cse_reuse(self):
    for X in [ca.SX,ca.MX]:
        x = ca.MX.sym("x")
        y = ca.sin(ca.MX.sym("y"))
        z = ca.MX.sym("z")
        
        e = (x+y/z)
        r = ca.vertcat(ca.cse(e),2*y)
        print(r)
        self.assertTrue("@1=sin(y)" in str(r))
      
  def test_is_slice2(self):

    m = 13
    n = 11

    x = ca.MX.sym("x",m,n)
    ca.DM.rng(1)
    A = ca.DM.rand(m,n)


    for step_col in [1,2,3]:
        for step_row in [1,2,3]:
            nz = []
            
            cslice = range(3,7,step_col)
            for col in cslice:
                rslice = range(1,9,step_row)
                for row in rslice:
                    nz.append(row+m*col)

            assert(ca.is_slice2(nz))
            
            # adversarial insert
            for j in range(len(nz)):
                nz_mod = list(nz)
                if j==0:
                    val = nz[0]-1
                else:
                    val = int(ca.floor((nz[j]+nz[j-1])/2))
                nz_mod.insert(j,val)
                assert(not ca.is_slice2(nz_mod))
            r = x[rslice,cslice]
            assert(";" in str(r))
            f = ca.Function('f',[x],[r])
            self.checkarray(f(A),A[rslice,cslice])
  
  def test_is_diff_sparsity_propagation(self):
    for X in [ca.MX,ca.SX]:
        x = X.sym("x")
        u = X.sym("u")
        q = X.sym("q")
        
        is_diff_in = [True, True, False]
        is_diff_out = [True, False]

        F = ca.Function('F', [x, u, q], [ca.sin(u+x+q), u*q],
                     {"never_inline": True,
                      "is_diff_in": is_diff_in,
                      "is_diff_out": is_diff_out})
        Fs = [F]
        if args.run_slow:
            res = self.check_codegen(F,inputs=[1,2,3])
            Fs.append(res["F"])
        
            F = ca.Function('F', [x, u, q], [ca.sin(u+x+q), u*q])
            
            res = self.check_codegen(F,inputs=[1,2,3],with_jac_sparsity=True,with_forward=True,with_reverse=True)
            Fs.append(ca.external("F",res["libname"],{"is_diff_in": is_diff_in, "is_diff_out": is_diff_out}))

            res = self.check_codegen(F,inputs=[1,2,3])
            Fs.append(ca.external("F",res["libname"],{"is_diff_in": is_diff_in, "is_diff_out": is_diff_out}))
                 
        for F in Fs:

            for X2 in [ca.MX,ca.SX]:
                X0 = X2.sym("X0")
                U0 = X2.sym("U0")
                Q0 = X2.sym("Q0")
                [xn, qn] = F(5*X0, 6*U0, 7*Q0)
                w = ca.vertcat(X0, U0, Q0)
                
                for ad_weight_sp in [0,1]:
                    
                    H = ca.Function('H',[w],[2*xn,3*qn],{"ad_weight_sp": ad_weight_sp})
                    # xn should depend on X0 and U0 only (not Q0)
                    self.assertEqual(H.sparsity_jac(0,0).nnz(), 2)
                    # qn is non-diff output: no dependencies
                    self.assertEqual(H.sparsity_jac(0,1).nnz(), 0)

  def test_set_precision(self):
    # Issue #4326: numeric MX constants honor DM::set_precision
    # (the value lives in a DM-like node, so MX printing piggybacks
    # on DM's stream precision settings).
    pi_val = 3.141592653589793239
    mx_scalar = ca.MX(pi_val)
    mx_dm = ca.MX(ca.DM([1.234567890123, 9.876543210987]))
    try:
      self.assertEqual(str(mx_scalar), "3.14159")
      ca.DM.set_precision(12)
      self.assertEqual(str(mx_scalar), "3.14159265359")
      self.assertEqual(str(mx_dm), "[1.23456789012, 9.87654321099]")
      ca.DM.set_precision(4)
      ca.DM.set_scientific(True)
      self.assertTrue("e" in str(mx_scalar).lower())
    finally:
      ca.DM.set_precision(6)
      ca.DM.set_scientific(False)

if __name__ == '__main__':
    unittest.main()
