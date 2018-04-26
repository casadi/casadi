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
from numpy import random, array, linalg, matrix, zeros, ones, ndarray, eye
import unittest
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
      zr=reshape(zr,(zt.shape));
    self.assertEqual(zt.shape[0],zr.shape[0],"%s dimension error. Got %s, expected %s" % (name,str(zt.shape),str(zr.shape)))
    self.assertEqual(len(zt.shape),len(zr.shape),"%s dimension error. Got %s, expected %s" % (name,str(zt.shape),str(zr.shape)))
    self.assertEqual(zt.shape[1],zr.shape[1],"%s dimension error. Got %s, expected %s" % (name,str(zt.shape),str(zr.shape)))
    for i in range(zr.shape[0]):
      for j in range(zr.shape[1]):
        self.assertAlmostEqual(zt[i,j],zr[i,j],10,"%s evaluation error. %s <-> %s" % (name, str(zt),str(zr)))

def checkMXoperations(self,ztf,zrf,name):
    x = MX.sym("x",1,3)
    z=vertcat(*[x*(i+1) for i in range(8)])
    f = Function("f", [x],[ztf(z)])
    L=[1,2,3]
    f_out = f(L)
    zt = f_out.full()
    zr = array([[L[0]*(i+1),L[1]*(i+1),L[2]*(i+1)] for i in range(8)])
    checkarray(self,zrf(zr),zt,name)
    return (zt,zrf(zr))

def checkMXoperations2(self,ztf,zrf,name):
    x = MX.sym("x",3,1)
    z = horzcat(*[x*i for i in range(8)])
    f = Function("f", [x],[ztf(z)])
    L=[1,2,3]
    f_out = f(L)
    zt = f_out.full()
    zr = array([[L[0]*i,L[1]*i,L[2]*i] for i in range(8)]).T
    checkarray(self,zrf(zr),zt,name)
    return zt

def checkMXoperations3(self,ztf,zrf,name):
    x = MX.sym("x",3,1)
    p = horzcat(*[x[0,0],x[1,0],x[2,0]])
    z = vertcat(*[p*i for i in range(8)])
    f = Function("f", [x],[ztf(z)])
    L=[1,2,3]
    f_out = f(L)
    zt = f_out.full()
    zr = array([[L[0]*i,L[1]*i,L[2]*i] for i in range(8)])
    checkarray(self,zrf(zr),zt,name)
    return (zt,zrf(zr))

class MXtests(casadiTestCase):

  def setUp(self):
    self.pool=FunctionPool()
    self.pool.append(lambda x: sqrt(x[0]),sqrt,"sqrt")
    self.pool.append(lambda x: sin(x[0]),sin,"sin")
    self.pool.append(lambda x: cos(x[0]),cos,"cos")
    self.pool.append(lambda x: tan(x[0]),tan,"tan")
    self.pool.append(lambda x: arctan(x[0]),arctan,"arctan")
    self.pool.append(lambda x: arcsin(x[0]),arcsin,"arcsin")
    self.pool.append(lambda x: arccos(x[0]),arccos,"arccos")
    self.pool.append(lambda x: exp(x[0]),exp,"exp")
    self.pool.append(lambda x: log(x[0]),log,"log",flags={'nozero'})
    self.pool.append(lambda x: x[0]**0,lambda x : x**0,"x^0",flags={'nozero'})
    self.pool.append(lambda x: x[0]**1,lambda x : x**1,"^1")
    self.pool.append(lambda x: x[0]**(-2),lambda x : x**(-2),"^-2",flags={'nozero'})
    self.pool.append(lambda x: x[0]**(0.3),lambda x : x**(0.3),"^0.3")
    self.pool.append(lambda x: floor(x[0]),floor,"floor")
    self.pool.append(lambda x: ceil(x[0]),ceil,"ceil")
    self.Jpool=FunctionPool()
    self.Jpool.append(lambda x: sqrt(x[0]),lambda x:diag(1/(2.0*sqrt(x))),"sqrt")
    self.Jpool.append(lambda x: sin(x[0]),lambda x:diag(cos(x)),"sin")
    self.Jpool.append(lambda x: cos(x[0]),lambda x:diag(-sin(x)),"cos")
    self.Jpool.append(lambda x: tan(x[0]),lambda x:diag(1.0/cos(x)**2),"tan")
    self.Jpool.append(lambda x: arctan(x[0]),lambda x:diag( 1.0/(x**2+1)),"arctan")
    self.Jpool.append(lambda x: arcsin(x[0]),lambda x:diag( 1.0/sqrt(1-x**2)),"arcsin")
    self.Jpool.append(lambda x: arccos(x[0]),lambda x: diag(-1.0/sqrt(1-x**2)),"arccos")
    self.Jpool.append(lambda x: exp(x[0]),lambda x: diag(exp(x)),"exp")
    self.Jpool.append(lambda x: log(x[0]),lambda x: diag(1.0/x),"log")
    self.Jpool.append(lambda x: x[0]**0,lambda x :diag(zeros(x.shape)),"x^0")
    self.Jpool.append(lambda x: x[0]**1,lambda x : diag(ones(x.shape)),"^1")
    self.Jpool.append(lambda x: x[0]**(-2),lambda x : diag(-2.0/x**3),"^-2")
    self.Jpool.append(lambda x: x[0]**(0.3),lambda x :diag( 0.3/x**0.7),"^0.3")
    self.matrixpool=FunctionPool()
    #self.matrixpool.append(lambda x: norm_2(x[0]),linalg.norm,"norm_2")
    #self.matrixpool.append(lambda x: norm_1(x[0]),lambda x: sum(sum(abs(x))),"norm_1")
    #self.matrixpool.append(lambda x: norm_inf(x[0]),lambda x: abs(matrix(x)).max(),"norm_inf")
    self.matrixbinarypool=FunctionPool()
    self.matrixbinarypool.append(lambda a: a[0]+a[1],lambda a: a[0]+a[1],"Matrix+Matrix")
    self.matrixbinarypool.append(lambda a: a[0]-a[1],lambda a: a[0]-a[1],"Matrix-Matrix")
    self.matrixbinarypool.append(lambda a: a[0]*a[1],lambda a: a[0]*a[1],"Matrix*Matrix")
    self.matrixbinarypool.append(lambda a: fmax(a[0],a[1]),lambda a: fmax(a[0],a[1]),"fmin")

    self.matrixbinarypool.append(lambda a: fmin(a[0],a[1]),lambda a: fmin(a[0],a[1]),"fmax")
    self.matrixbinarypool.append(lambda a: mtimes(a[0],a[1].T),lambda a: numpy.dot(a[0],a[1].T),"mtimes(Matrix,Matrix.T)")
    self.matrixbinarypool.append(lambda a: arctan2(a[0],a[1]),lambda a: arctan2(a[0],a[1]),"arctan2")
    #self.matrixbinarypool.append(lambda a: inner_mul(a[0],trans(a[1])),lambda a: c.dot(a[0].T,a[1]),name="inner_mul(Matrix,Matrix)")
    self.matrixbinarypool.append(lambda a: mtimes(a[0],a[1].T),lambda a: numpy.dot(a[0],a[1].T),"mtimes(Matrix,Matrix.T)")

  def test_MX1(self):
    self.message("MX constructor")
    x = MX.sym("x",2,3)
    self.assertEqual(x.size1(),2,"MX fails to indicate its size1")
    self.assertEqual(x.size2(),3,"MX fails to indicate its size2")

  def test_MXvertcat(self):
    self.message("MX vertcat")
    x = MX.sym("x",1,3)
    y = MX.sym("y",1,3)
    z=vertcat(*(x,y))
    self.assertEqual(z.size1(),2,"MX fails to indicate its size1")
    self.assertEqual(z.size2(),3,"MX fails to indicate its size2")

  def test_MX_fun1(self):
    self.message("MXFunction single input, single output")
    # check if x->2*x
    # evaluates correctly for x=3
    x = MX.sym("x")
    y = 2*x
    f = Function("f", [x],[y])
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
    x = MX.sym("x")
    y = MX.sym("y")
    f = Function("f", [x,y],[x+y,y*x])
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
    xy = MX.sym("xy",2)
    f = Function("f", [xy],[xy[0]+xy[1],xy[0]*xy[1]])
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
    xy = MX.sym("xy",1,2)
    f = Function("f", [xy],[xy[0,0]+xy[0,1],xy[0,0]*xy[0,1]])
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
    xy = MX.sym("xy",2)
    z=vertcat(*[xy[0]+xy[1],xy[0]*xy[1]])
    f = Function("f", [xy],[z])
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
    xy = MX.sym("xy",2)
    z=horzcat(*[xy[0]+xy[1],xy[0]*xy[1]])
    f = Function("f", [xy],[z])
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
    for X in [SX,MX]:
      x=X.sym("x")

      for tr in [True,False]:
        for i in [0,1,2]:
          self.assertEqual(which_depends(x,X(0,1),i,tr),[False]*(1 if tr else 0))

          self.assertEqual(which_depends(X(0,1),x,i,tr),[False]*(0 if tr else 1))
          self.assertTrue(len(which_depends(X(0,1),X(0,1),i,tr))==0)

  def test_issue83(self):
    x=MX.sym("x")
    y=MX.sym("y")

    z = x + y

    f = Function("f", [x,y],[z])

    fc = f(MX(3),y)

    g = Function("g", [y],[fc])
    g_in = [7]
    g_out = g(g_in)

    self.assertAlmostEqual(g_out[0],10,10,"issue #83")

    fc = f(x,MX(7))

    g = Function("g", [x],[fc])
    g_in = [3]
    g_out = g(g_in)

    self.assertAlmostEqual(g_out[0],10,10,"issue #83")

  def test_identitySX(self):
    self.message("identity SXFunction")
    x = SX.sym("x")
    f = Function("f", [x],[x])
    f_in = [3]
    f_out = f(f_in)
    self.assertAlmostEqual(f_out[0,0], 3,10)

  def test_identityMX(self):
    self.message("identity Function")
    x = MX.sym("x")
    f = Function("f", [x],[x])
    f_in = [3]
    f_out = f(f_in)
    self.assertAlmostEqual(f_out[0,0], 3,10)

  def test_MXorder(self):
    self.message("Function order of non-zero elements")
    x = MX.sym("x",2,3)
    f = Function("f", [x],[x+x])

    self.assertEqual(f.n_in(),1,"Function fails to indicate correct number of inputs")
    self.assertEqual(f.n_out(),1,"Function fails to indicate correct number of outputs")
    L=[1,2,3,4,5,6]
    f_in = DM(f.sparsity_in(0),L)
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
    a = MX(0,1)
    b = a.T
    self.assertEqual(b.size1(),1)
    self.assertEqual(b.size2(),0)

  def test_MXtrans(self):
    self.message("trans(MX)")
    x = MX.sym("x",2,3)
    z=x.T
    self.assertEqual(z.size1(),3,"Vec returns MX of wrong dimension")
    self.assertEqual(z.size2(),2,"Vec returns MX of wrong dimension")
    f = Function("f", [x],[z])
    self.assertEqual(f.n_in(),1,"Function fails to indicate correct number of inputs")
    self.assertEqual(f.n_out(),1,"Function fails to indicate correct number of outputs")
    L=[1,2,3,4,5,6]
    f_in = DM(f.sparsity_in(0),L)
    f_out = f(f_in)
    zt = f_out.full()

    ztr=numpy.reshape(zt,(3,2))
    Lr=numpy.reshape(L,(2,3),'F')
    for i in range(2):
      for j in range(3):
        self.assertAlmostEqual(Lr[i,j], ztr[j,i],10)

  def test_MXvec(self):

    u = DM([[10*j+i for i in range(3)] for j in range(4) ])

    U = MX.sym("u",u.shape)

    f = Function("f", [U],[vec(U)])
    f_out = f(u)

    self.checkarray(vec(u),f_out,"vec")

  def test_MXreshape(self):
    self.message("reshape(MX)")
    x = MX.sym("x",2,3)
    z=c.reshape(x,(1,6))
    self.assertEqual(z.size1(),1,"Vec returns MX of wrong dimension")
    self.assertEqual(z.size2(),6,"Vec returns MX of wrong dimension")
    f = Function("f", [x],[z])
    self.assertEqual(f.n_in(),1,"Function fails to indicate correct number of inputs")
    self.assertEqual(f.n_out(),1,"Function fails to indicate correct number of outputs")
    L=[1,2,3,4,5,6]
    f_in = DM(f.sparsity_in(0),L)
    f_out = f(f_in)
    zt = f_out.full()
    for i in range(len(L)):
      self.assertAlmostEqual(L[i], zt[0,i],10)

  def test_MXcompose(self):
    self.message("compositions of vec, trans, reshape with vertcat")
    checkMXoperations(self,lambda x: x,lambda x: x,'vertcat')
    checkMXoperations(self,lambda x: x.T,lambda x: x.T,'trans(vertcat)')
    checkMXoperations(self,lambda x: x.T.T,lambda x: x,'trans(trans(vertcat))')
    checkMXoperations(self,lambda x: vec(x.T),lambda x: numpy.reshape(x,(numpy.prod(x.shape),1)),'vec(trans(vertcat))')
    checkMXoperations(self,lambda x: vec(x).T,lambda x: numpy.reshape(x.T,(numpy.prod(x.shape),1)).T,'vec(trans(vertcat))')
    checkMXoperations(self,lambda x: c.reshape(x.T,(6,4)).T,lambda x: numpy.reshape(x,(4,6)),'reshape(vertcat)')
    checkMXoperations(self,lambda x: c.reshape(x,(6,4)).T,lambda x: numpy.reshape(x.T,(4,6)),'reshape(trans(vertcat))')
    checkMXoperations(self,lambda x: c.reshape(x.T,(6,4)),lambda x: numpy.reshape(x,(4,6)).T,'trans(reshape(vertcat))')

  def test_MXcompose2(self):
    self.message("compositions of vec, trans, reshape with horzcat")
    checkMXoperations2(self,lambda x: x,lambda x: x,'horzcat')
    checkMXoperations2(self,lambda x: x.T,lambda x: x.T,'trans(horzcat)')
    checkMXoperations2(self,lambda x: x.T.T,lambda x: x,'trans(trans(horzcat))')
    checkMXoperations2(self,lambda x: vec(x.T),lambda x: numpy.reshape(x,(numpy.prod(x.shape),1)),'vec(trans(horzcat))')
    checkMXoperations2(self,lambda x: vec(x).T,lambda x: numpy.reshape(x.T,(numpy.prod(x.shape),1)).T,'vec(trans(horzcat))')
    checkMXoperations2(self,lambda x: c.reshape(x.T,(6,4)).T,lambda x: numpy.reshape(x,(4,6)),'reshape(horzcat)')
    checkMXoperations2(self,lambda x: c.reshape(x,(6,4)).T,lambda x: numpy.reshape(x.T,(4,6)),'reshape(trans(horzcat))')
    checkMXoperations2(self,lambda x: c.reshape(x.T,(6,4)),lambda x: numpy.reshape(x,(4,6)).T,'trans(reshape(horzcat))')

  def test_MXcompose3(self):
    self.message("compositions of vec, trans, reshape with vertcat")
    checkMXoperations3(self,lambda x: x,lambda x: x,'snippet')
    checkMXoperations3(self,lambda x: x.T,lambda x: x.T,'trans(snippet)')
    checkMXoperations3(self,lambda x: x.T.T,lambda x: x,'trans(trans(snippet))')
    checkMXoperations3(self,lambda x: vec(x.T),lambda x: numpy.reshape(x,(numpy.prod(x.shape),1)),'vec(trans(snippet))')
    checkMXoperations3(self,lambda x: vec(x).T,lambda x: numpy.reshape(x.T,(numpy.prod(x.shape),1)).T,'vec(trans(snippet))')
    checkMXoperations3(self,lambda x: c.reshape(x.T,(6,4)).T,lambda x: numpy.reshape(x,(4,6)),'reshape(snippet)')
    checkMXoperations3(self,lambda x: c.reshape(x,(6,4)).T,lambda x: numpy.reshape(x.T,(4,6)),'reshape(trans(snippet))')
    checkMXoperations3(self,lambda x: c.reshape(x.T,(6,4)),lambda x: numpy.reshape(x,(4,6)).T,'trans(reshape(snippet))')

  def test_MXcompose4(self):
    self.message("compositions of horzcat + vertcat")
    checkMXoperations(self,lambda x: vertcat(*[x]),lambda x: x,'vertcat(*vertcat)')
    checkMXoperations(self,lambda x: vertcat(*[x,x*2]),lambda x: numpy.vstack((x,x*2)),'vertcat(*vertcat,vertcat)')
    checkMXoperations(self,lambda x: horzcat(*[x]),lambda x: x,'horzcat(*vertcat)')
    checkMXoperations(self,lambda x: horzcat(*[x,x*2]),lambda x: numpy.hstack((x,x*2)),'horzcat(*vertcat,vertcat)')

    checkMXoperations2(self,lambda x: vertcat(*[x]),lambda x: x,'vertcat(*horzcat)')
    checkMXoperations2(self,lambda x: vertcat(*[x,x*2]),lambda x: numpy.vstack((x,x*2)),'vertcat(*horzcat,horzcat)')
    checkMXoperations2(self,lambda x: horzcat(*[x]),lambda x: x,'horzcat(*horzcat)')
    checkMXoperations2(self,lambda x: horzcat(*[x,x*2]),lambda x: numpy.hstack((x,x*2)),'horzcat(*horzcat,horzcat)')

    checkMXoperations3(self,lambda x: vertcat(*[x]),lambda x: x,'vertcat(*snippet)')
    checkMXoperations3(self,lambda x: vertcat(*[x,x*2]),lambda x: numpy.vstack((x,x*2)),'vertcat(*snippet,snippet)')
    checkMXoperations3(self,lambda x: horzcat(*[x]),lambda x: x,'horzcat(*snippet)')
    checkMXoperations3(self,lambda x: horzcat(*[x,x*2]),lambda x: numpy.hstack((x,x*2)),'horzcat(*snippet,snippet)')

  @known_bug()  # Test refactoring, cf. #1436
  def test_MXslicingnew(self):
    self.message("MX slicing new")

    self.message(":dense")
    x = MX.sym("x",3,2)
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

    sp=Sparsity(4,3,[0,2,2,3],[1,2,1])
    x=MX.sym("X",sp)
    sx0=[0.738,0.39,0.99]
    x0=DM(Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.39,0.99]).full()
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
    x=MX.sym("x",2,3)
    f = Function("f", [x],[3*x])
    x_in = f.mx_in()
    x_out = f.call(x_in)
    g = Function("g", [x_in[0]],[6*x_out[0]])
    n=[1,2,3,4,5,6]
    f_in=DM(f.sparsity_in(0),n)
    f_out = f(f_in)
    g_in=DM(g.sparsity_in(0),n)
    g_out = g(g_in)
    checkarray(self,6*f_out.full(),g_out.full(),"slicing(trans)")

  def test_scalarMX(self):
      x=MX.sym("x")
      x0=0.738
      self.numpyEvaluationCheckPool(self.pool,[x],x0,name="scalarMX")

      self.numpyEvaluationCheckPool(self.matrixpool,[x],x0,name="scalarMX")

  def test_MXJacobian(self):
    self.message("MX(1,1) unary operation, jacobian")
    self.Jpool=FunctionPool()
    self.message("SX(1,1) unary operation, jacobian")
    x=MX.sym("x")
    x0=array([[0.738]])

    def fmod(f,x):
      J=f.jacobian_old(0, 0)
      return J

    self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="MX unary operations, jacobian",fmod=fmod)

  def test_MXJacobians(self):
      self.message("MX(3,1) unary operation, jacobian")
      x=MX.sym("x",3,1)

      x0=array([0.738,0.9,0.3])

      def fmod(f,x):
        J=f.jacobian_old(0, 0)
        return J

      self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="MX unary operations, jacobian",fmod=fmod)

  def test_MXbinary(self):
      self.message("MX binary operations")
      x=MX.sym("x",3,2)
      y=MX.sym("x",3,2)
      x0=array([[0.738,0.2],[ 0.1,0.39 ],[0.99,0.999999]])
      y0=array([[1.738,0.6],[ 0.7,12 ],[0,-6]])
      self.numpyEvaluationCheckPool(self.matrixbinarypool,[x,y],[x0,y0],name="MX")

  def test_MXSparse(self):
      self.message("MX unary operations, sparse")
      sp=Sparsity(4,3,[0,2,2,3],[1,2,1])

      x=MX.sym("x",sp)
      if scipy_available:
        x0=DM(Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).sparse()

        self.numpyEvaluationCheckPool(self.pool,[x],array(x0.todense()),name="MX",setx0=x0,excludeflags={'nozero'})
        self.numpyEvaluationCheckPool(self.matrixpool,[x],array(x0.todense()),name="MX",setx0=x0)
      else:
        x0=DM(Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).full()

        self.numpyEvaluationCheckPool(self.pool,[x],x0,name="MX",setx0=x0,excludeflags={'nozero'})
        self.numpyEvaluationCheckPool(self.matrixpool,[x],x0,name="MX",setx0=x0)

  def test_MXbinarySparse(self):
      self.message("SX binary operations")
      spx=Sparsity(4,3,[0,2,2,3],[1,2,1])
      spy=Sparsity(4,3,[0,2,2,3],[0,2,3])
      xx=MX.sym("x",spx)
      yy=MX.sym("y",spy)
      if scipy_available:
        x0=DM(Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).sparse()
        y0=DM(Sparsity(4,3,[0,2,2,3],[0,2,3]),[1.738,0.7,-6]).sparse()

        self.numpyEvaluationCheckPool(self.matrixbinarypool,[xx,yy],[array(x0.todense()),array(y0.todense())],name="MX",setx0=[x0,y0])
      else:
        x0=DM(Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).full()
        y0=DM(Sparsity(4,3,[0,2,2,3],[0,2,3]),[1.738,0.7,-6]).full()

        self.numpyEvaluationCheckPool(self.matrixbinarypool,[xx,yy],[x0,y0],name="MX",setx0=[x0,y0])

  def test_symbolcheck(self):
    self.message("Check if non-symbolic inputs are caught")
    self.assertRaises(RuntimeError, lambda : Function("f", [MX(0)],[MX.sym("x")]))

  def test_unite(self):
    self.message("unite operation")
    import numpy
    numpy.random.seed(42)
    xn = numpy.random.random((3,4))
    x=MX(3,4)
    y=MX.sym("x",3,4)
    z=unite(x,y)
    f = Function("f", [y],[z])
    f_in = [0]*f.n_in();f_in[0]=xn
    f_out = f(*f_in)
    self.checkarray(f_out,xn,"unite dense")

    spx=Sparsity(4,3,[0,2,2,3],[1,2,1])
    spy=Sparsity(4,3,[0,1,2,3],[0,2,2])

    nx=DM.zeros(spx)
    for k in range(nx.nnz()):
      nx.nz[k]= numpy.random.rand()
    ny=DM.zeros(spy)
    for k in range(nx.nnz()):
      ny.nz[k]= numpy.random.rand()

    nxn = nx.full()
    nyn = ny.full()
    x=MX.sym("x",spx)
    y=MX.sym("y",spy)
    z=unite(x,y)

    f = Function("f", [x,y],[z])
    f_in = [0]*f.n_in();f_in[0]=nx
    f_in[1]=ny
    f_out = f(*f_in)
    self.checkarray(f_out,nxn+nyn,"unite sparse")

  def test_imatrix_index(self):
    self.message("IM indexing")
    X = MX.sym("x",2,2)
    Y = X.nz[IM([[0,2],[1,1],[3,3]])]

    f = Function("f", [X],[Y])
    f_in = [0]*f.n_in();f_in[0]=DM(f.sparsity_in(0),[1,2,3,4])
    f_out = f(*f_in)

    self.checkarray(f_out,array([[1,3],[2,2],[4,4]]),"IM indexing")

    Y = X[:,:]
    Y.nz[IM([[0,2]])] = DM([[9,8]])

    f = Function("f", [X],[Y])
    f_in = DM(f.sparsity_in(0),[1,2,3,4])
    f_out = f(f_in)

    self.checkarray(f_out,array([[9,8],[2,4]]),"IM indexing assignment")

  def test_subsass(self):
     self.message("Check subscripted assignment")

     X = MX.sym("x",2,2)
     X[0,0]=MX(5)
     X[0,0]=5
     X[:,0]=8

     x=MX.sym("X",3,4)
     import numpy
     numpy.random.seed(42)
     xn = numpy.random.random((3,4))
     r = numpy.zeros((7,8))
     y=MX.zeros(7,8)
     y[1:4,[2,4,6,7]]=x
     r[1:4,[2,4,6,7]]=xn
     fy = Function("fy", [x],[y])
     fy_in = [0]*fy.n_in();fy_in[0]=xn
     fy_out = fy(*fy_in)

     self.checkarray(fy_out,r,"subscripted assigment")

     y=MX(7,8)
     y[1:4,[2,4,6,7]]=x
     r[1:4,[2,4,6,7]]=xn
     fy = Function("fy", [x],[y])
     fy_out = fy(xn)
     self.checkarray(fy_out,r,"subscripted assigment")

     kl=[2,4,5,8]

     s=y.sparsity()
     for k in kl:
       r[s.row()[k],s.get_col()[k]]=1.0

     y.nz[kl]=MX(1)
     fy = Function("fy", [x],[y])
     fy_out = fy(xn)
     self.checkarray(fy_out,r,"subscripted assigment")

     y.nz[kl]=x.nz[[0,1,2,3]]
     s=y.sparsity()
     sx=x.sparsity()
     cnt=0
     for k in kl:
       r[s.row()[k],s.get_col()[k]]=xn[sx.row()[cnt],sx.get_col()[cnt]]
       cnt+=1
     fy = Function("fy", [x],[y])
     fy_out = fy(xn)
     self.checkarray(fy_out,r,"subscripted assigment")

  def test_erase(self):
    self.message("Erase function")
    self.message(":dense")
    y=MX.sym("Y",7,8)
    import numpy
    r=2*numpy.ones((7,8))
    r[1:4,[2,4,6,7]]=numpy.zeros((3,4))
    z = y *2
    z.erase([1,2,3],[2,4,6,7])
    f = Function("f", [y],[z])
    f_in = [0]*f.n_in();f_in[0]=DM(f.sparsity_in(0),[1]*56)
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
    A = MX.sym("A",m,n)
    b_ = numpy.random.random((m,1))
    b = MX.sym("b",m,1)
    C_ = numpy.random.random((m,m))
    C = MX.sym("C",m,m)
    D_ = numpy.random.random((m,n))
    D = MX.sym("D",m,n)
    e_ = numpy.random.random((m,1))
    e = MX.sym("e",m,1)
    x_ = numpy.random.random((n,1))
    x = MX.sym("x",n,1)

    Axb = mtimes(A,x)+b
    Dxe = mtimes(D,x)+e
    a = mtimes(mtimes(Axb.T,C),Dxe)

    f = Function("f", [x,A,b,C,D,e],[a])
    f_out = f(x_, A_, b_, C_, D_, e_)

    f_ = numpy.dot(numpy.dot((numpy.dot(A_,x_)+b_).T,C_),(numpy.dot(D_,x_)+e_))

    self.checkarray(f_out,f_,"evaluation")


    J_ = numpy.dot(numpy.dot((numpy.dot(D_,x_)+e_).T,C_.T),A_) + numpy.dot(numpy.dot((numpy.dot(A_,x_)+b_).T,C_),D_)

    for w in [0, 1]:
      f = Function("f", [x,A,b,C,D,e], [a], {"ad_weight":w, "ad_weight_sp":w})
      J = f.jacobian_old(0, 0)
      J_in = [0]*J.n_in();J_in[0]=x_
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
      sp = Sparsity(m,n)
      for i in range(int((n*m)/2)):
        sp.add_nz(numpy.random.randint(m),numpy.random.randint(n))
      return sp

    def gentest(m,n):
      As = randsparsity(m,n)
      A_ = DM.zeros(As)
      for k in range(As.nnz()):
        A_.nz[k]= numpy.random.rand()
      A = MX.sym("A",As)
      return (A_.sparse(),A)

    (A_,A)=gentest(m,n)
    (b_,b)=gentest(m,1)
    (C_,C)=gentest(m,m)
    (D_,D)=gentest(m,n)
    (e_,e)=gentest(m,1)
    x_ = numpy.random.random((n,1))
    x = MX.sym("x",n,1)

    Axb = mtimes(A,x)+b
    Dxe = mtimes(D,x)+e
    a = mtimes(mtimes(Axb.T,C),Dxe)

    f = Function("f", [x,A,b,C,D,e],[a])
    f_in = [0]*f.n_in();f_in[0]=x_
    f_in[1]=A_
    f_in[2]=b_
    f_in[3]=C_
    f_in[4]=D_
    f_in[5]=e_
    f_out = f(*f_in)


    Axb_ = A_*x_+b_
    Dxe_ = D_*x_+e_

    f_ = Axb_.T*C_*Dxe_

    self.checkarray(f_out,f_,"evaluation")


    J_ = (D_*x_+e_).T*C_.T*A_ + (A_*x_+b_).T*C_*D_

    for w in [0, 1]:
      f = Function("f", [x,A,b,C,D,e], [a], {"ad_weight":w, "ad_weight_sp":w})

      J = f.jacobian_old(0, 0)
      J_in = [0]*J.n_in();J_in[0]=x_
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
      sp = Sparsity(m,n)
      for k in range(int((n*m)/2)):
        i = numpy.random.randint(m)
        j = numpy.random.randint(n)
        if not(i == int(m/2)):
          if n==1 or not(j == int(n/2)):
            sp.add_nz(i,j)
      return sp

    def gentest(m,n):
      As = randsparsity(m,n)
      A_ = DM.zeros(As)
      for k in range(As.nnz()):
        A_.nz[k]= numpy.random.rand()
      A = MX.sym("A",As)
      return (A_.sparse(),A)

    (A_,A)=gentest(m,n)
    (b_,b)=gentest(m,1)
    (C_,C)=gentest(m,m)
    (D_,D)=gentest(m,n)
    (e_,e)=gentest(m,1)
    x_ = numpy.random.random((n,1))
    x = MX.sym("x",n,1)

    Axb = mtimes(A,x)+b
    Dxe = mtimes(D,x)+e
    a = mtimes(mtimes(Axb.T,C),Dxe)

    f = Function("f", [x,A,b,C,D,e],[a])
    f_in = [0]*f.n_in();f_in[0]=x_
    f_in[1]=A_
    f_in[2]=b_
    f_in[3]=C_
    f_in[4]=D_
    f_in[5]=e_
    f_out = f(*f_in)


    Axb_ = A_*x_+b_
    Dxe_ = D_*x_+e_

    f_ = Axb_.T*C_*Dxe_

    self.checkarray(f_out,f_,"evaluation")


    J_ = (D_*x_+e_).T*C_.T*A_ + (A_*x_+b_).T*C_*D_

    for w in [0, 1]:
      f = Function("f", [x,A,b,C,D,e], [a], {"ad_weight":w, "ad_weight_sp":w})
      J = f.jacobian_old(0, 0)
      J_in = [0]*J.n_in();J_in[0]=x_
      J_in[1]=A_
      J_in[2]=b_
      J_in[3]=C_
      J_in[4]=D_
      J_in[5]=e_
      J_out = J.call(J_in)

      self.checkarray(J_out[0],J_,"evaluation")


  def test_chaining(self):
    self.message("Chaining SX and MX together")
    x=SX.sym("x")
    y=x**3
    f=Function("f", [x],[y])
    J=f.jacobian_old(0, 0)

    X=MX.sym("X")
    F=Function("F", [X], J(X))

    x_=1.7
    F_out = F(x_)
    self.checkarray(F_out[0],3*x_**2,"Chaining eval")

  def test_issue107(self):
    self.message("Regression test for issue 107: +=")
    x=MX.sym("x")
    y=MX.sym("y")

    z=x
    z+=y

    self.assertTrue(x.is_symbolic())
    self.assertFalse(z.is_symbolic())

  def test_MXd_trivial(self):
    self.message("symbolic variables and constants jac")
    X =  MX.sym("X",10)
    V =  MX.sym("V")
    J = jacobian(X,X)
    self.assertTrue(isinstance(J,MX))
    self.assertEqual(J.nnz(),10)
    self.assertEqual(J.size1(),10)
    self.assertEqual(J.size2(),10)

    g = Function("g", [],[J])
    [g_out] = g.call([])
    self.checkarray(g_out,eye(10),"unit matrix")
    g = Function("g", [],[jacobian(MX.eye(3),X)])
    [g_out] = g.call([])
    self.checkarray(g_out,zeros((9,10)),"zero matrix")
    g = Function("g", [],[jacobian(X,V)])
    [g_out] = g.call([])
    self.checkarray(g_out,zeros((10,1)),"zero matrix")

    g = Function("g", [],[jacobian(MX.eye(3),V)])
    [g_out] = g.call([])
    self.checkarray(g_out,zeros((9,1)),"zero matrix")

  def test_MXd_substractionl(self):
    self.message("substraction jac")
    V =  MX.sym("V")
    X =  MX.sym("X")
    g = Function("g", [],[jacobian(X-V, X)])
    [g_out] = g.call([])
    self.checkarray(g_out,ones((1,1)), "one")

    g = Function("g", [],[jacobian(X-V, V)])
    [g_out] = g.call([])
    self.checkarray(g_out,-ones((1,1)), "one")

    g = Function("g", [],[jacobian(V-X, X)])
    [g_out] = g.call([])
    self.checkarray(g_out,-ones((1,1)), "one")

    g = Function("g", [],[jacobian(V-X, V)])
    [g_out] = g.call([])
    self.checkarray(g_out,ones((1,1)),"one")

  def test_MXd_mapping(self):
    self.message("mapping jac")
    X = MX.sym("X",3)
    Y = MX.sym("Y",2)
    J = jacobian(vertcat(X,Y),X)
    JJ = DM.ones(J.sparsity())
    self.checkarray(JJ,numpy.vstack((eye(3),zeros((2,3)))),"diag")
    J = jacobian(vertcat(X,Y),Y)
    JJ = DM.ones(J.sparsity())
    self.checkarray(JJ,numpy.vstack((zeros((3,2)),eye(2))),"diag")

  def test_null(self):
    self.message("Function null")
    x = MX.sym("x")

    f = Function("f", [x],[x**2,MX()])

    self.assertEqual(f.size1_out(1),0)
    self.assertEqual(f.size2_out(1),0)
    f_out = f(0)

    f = Function("f", [x,MX()],[x**2,MX()])

    self.assertEqual(f.size1_out(1),0)
    self.assertEqual(f.size2_out(1),0)
    f_out = f(0,0)

    r = f(x,MX())
    self.assertTrue(r[1].is_empty(True))

    r = f(MX(),MX())
    self.assertTrue(r[1].is_empty(True))

    #self.assertRaises(Exception,lambda : f([x,x],True))
    #self.assertRaises(Exception,lambda : f([[],[]],True))

  def test_issue184(self):
    self.message("Regression test issue #184")
    x = MX.sym("x", 3)
    y = x[0:0]
    self.assertEqual(y.nnz(),0)

  def test_indexingOutOfBounds(self):
    self.message("Indexing out of bounds")
    y = DM.zeros(4, 5)
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
    y = casadi.MX.sym("y", 3)
    self.assertRaises(RuntimeError,lambda : y[[0, 5]] )
    try:
      y[[0, 5]] = MX.sym("a")
      self.assertTrue(False)
    except RuntimeError:
      pass
    y[[0, 2]]
    y[[0, 2]] = MX.sym("a")

  def test_mtimes(self):
    A = MX(DM.ones((4,3)))
    B = MX(DM.ones((3,8)))
    C = MX(DM.ones((8,7)))

    self.assertRaises(RuntimeError,lambda : mtimes([]))

    D = mtimes([A])

    self.assertEqual(D.shape[0],4)
    self.assertEqual(D.shape[1],3)

    D = mtimes([A,B])

    self.assertEqual(D.shape[0],4)
    self.assertEqual(D.shape[1],8)

    D = mtimes([A,B,C])

    self.assertEqual(D.shape[0],4)
    self.assertEqual(D.shape[1],7)

  def test_truth(self):
    self.message("Truth values")
    self.assertRaises(Exception, lambda : bool(MX.sym("x")))
    #self.assertRaises(Exception, lambda : bool(MX.sym("x")>0))
    self.assertTrue(bool(MX(1)))
    self.assertFalse(bool(MX(0)))
    self.assertTrue(bool(MX(0.2)))
    self.assertTrue(bool(MX(-0.2)))
    self.assertRaises(Exception, lambda : bool(MX(DM([2.0,3]))))
    self.assertRaises(Exception, lambda : bool(MX()))


  def test_MXbool(self):
    self.message("bool")

    xy = MX.sym("x",2)
    x = xy[0]
    y = xy[1]

    f = Function("f", [xy],[vertcat(*[logic_and(x,y),logic_or(x,y),logic_not(x)])])


    for t1 in [0,1]:
      for t2 in [0,1]:
        T1 = t1!=0
        T2 = t2!=0
        f_out = f([t1,t2])
        self.checkarray(f_out,DM([T1 and T2,T1 or T2,not T1]),"bool(%d,%d): %s" % (t1,t2,str(f_out)))

  def test_MXineq(self):
    self.message("SX ineq")

    xy = MX.sym("x",2)
    x = xy[0]
    y = xy[1]


    f = Function("f", [xy],[vertcat(*[x<y,x<=y,x>=y,x==y,x!=y])])

    for t1 in [-10,0.1,0,1,10]:
      for t2 in [-10,0.1,0,1,10]:
        T1 = t1
        T2 = t2
        f_out = f([t1,t2])
        self.checkarray(f_out,DM([T1 < T2,T1 <= T2, T1 >= T2, T1 == T2, T1 != T2]),"ineq(%d,%d)" % (t1,t2))

  def test_if_else_zero(self):
    x = MX.sym("x")
    y = if_else(x,5,0)
    f = Function("f", [x],[y])
    f_in = 1
    f_out = f(f_in)
    self.assertTrue(f_out==5,"if_else_zero %s " % str(f_out))
    f_in = 0
    f_out = f(f_in)
    self.assertTrue(f_out==0,"if_else_zero")


  def test_if_else(self):
    x = MX.sym("x")
    y = if_else(x,1,2)
    f = Function("f", [x],[y])
    f_in = 1
    f_out = f(f_in)
    self.assertTrue(f_out==1,"if_else")
    f_in = 0
    f_out = f(f_in)
    self.assertTrue(f_out==2,"if_else")

  def test_regression491(self):
    self.message("regression #491")
    u = SX.sym("u")
    x = SX.sym("x")

    F = Function("F", [u,x],[u+1/x])

    U = MX.sym("U")

    X = F(U,U)
    G = F(U,X)

    for kk in range(2):
      gfcn = 0
      if kk==0:
        gfcn = Function("gfcn", [U], [G]).expand("e_gfcn", {"ad_weight":1})
      else:
        gfcn = Function("gfcn", [U],[G], {"ad_weight":1})
      J = gfcn.jacobian_old(0, 0)
      J_in = [0]*J.n_in();J_in[0]=1
      J_out = J.call(J_in)
      self.assertAlmostEqual(J_out[0],1,9)

  def test_ticket(self):
    J = [] + MX.sym("x")
    J = MX.sym("x") + []

  def test_jacobian_tools(self):
    self.message("jacobian")

    X = MX.sym("X")

    Y = jacobian(X**2,X)

    f = Function("f", [X], [Y])

    f_in=2.3
    f_out = f(f_in)

    self.assertAlmostEqual(f_out,4.6)

  def test_reshape(self):
    self.message("reshape")
    X = MX.sym("X",10)

    i = IM(Sparsity.lower(3),list(range(6)))

    i.print_dense()
    print(i.T.nz[:])

    T = X.nz[i]

    q = T.T.nz[:]**2
    f = Function("f", [X], [q])
    f_out = f(list(range(10)))

    self.checkarray(IM([0,1,9,4,16,25]),f_out)

    Y = MX.sym("Y",10)

    ff = Function("ff", [Y],f.call([Y],True))
    ff_out = ff(list(range(10)))

    self.checkarray(IM([0,1,9,4,16,25]),ff_out)

    J = Function("J", [X],[jacobian(q, X)])
    J_out = J(list(range(10)))

    i = horzcat(*[diag([0,2,4,6,8,10]),IM.zeros(6,4)])
    i[[2,3],:] = i[[3,2],:]

    self.checkarray(i,J_out)

    q = (T.T).nz[:]**2
    J = Function("J", [X],[jacobian(q,X)])
    J_in = [0]*J.n_in();J_in[0]=list(range(10))
    J_out = J(*J_in)

    i = horzcat(*[diag([0,2,4,6,8,10]),IM.zeros(6,4)])
    i[[2,3],:] = i[[3,2],:]

    self.checkarray(i,J_out)

  def test_vertcat(self):
    self.message("vertcat")
    X = MX.sym("X",10)

    T = vertcat(*[X[4],X[2]])
    q = T**2
    f = Function("f", [X],[q])
    f_in = [0]*f.n_in();f_in[0]=list(range(10))
    f_out = f.call(f_in)

    self.checkarray(IM([16,4]),f_out[0])

    Y = MX.sym("Y",10)

    ff = Function("ff", [Y],f.call([Y],True))
    ff_in = [0]*ff.n_in();ff_in[0]=list(range(10))
    ff_out = ff(*ff_in)

    self.checkarray(IM([16,4]),ff_out)
    J = Function("J", [X],[jacobian(q,X)])
    J_in = [0]*J.n_in();J_in[0]=list(range(10))
    J_out = J(*J_in)

    i = IM.zeros(2,10)
    i[0,4] = 8
    i[1,2] = 4

    self.checkarray(i,J_out)
    q = T**2
    J = Function("J", [X],[jacobian(q, X)])
    J_in = [0]*J.n_in();J_in[0]=list(range(10))
    J_out = J(*J_in)

    self.checkarray(i,J_out)

  def test_blockcat(self):
    x = MX.sym("x")

    y = blockcat([[x,2*x],[3*x,4*x]])
    f = Function("f", [x],[y])
    f_out = f(3)
    self.checkarray(f_out,DM([[3,6],[9,12]]))


  def test_veccats(self):
    x= MX.sym("x",2)
    self.assertTrue(hash(vec(x))==hash(x))

  def test_constmxmtimes(self):
    0.1*MX.ones(2)

  def test_is_regular(self):
    self.assertTrue(MX(DM([0,1])).is_regular())
    self.assertFalse(MX(DM([0,inf])).is_regular())
    with self.assertRaises(Exception):
      self.assertFalse(MX.sym("x",2).is_regular())

  def test_diagcat(self):
    C = diagcat(*[MX(DM(([[-1.4,-3.2],[-3.2,-28]]))),DM([[15,-12,2.1],[-12,16,-3.8],[2.1,-3.8,15]]),1.8,-4.0])
    self.assertTrue(isinstance(C,MX))
    r = DM([[-1.4,-3.2,0,0,0,0,0],[-3.2,-28,0,0,0,0,0],[0,0,15,-12,2.1,0,0],[0,0,-12,16,-3.8,0,0],[0,0,2.1,-3.8,15,0,0],[0,0,0,0,0,1.8,0],[0,0,0,0,0,0,-4]])
    r = sparsify(r)
    f = Function("f", [],[C])
    f_out = f.call([])

    self.checkarray(f_out[0],r)

  def test_tril2symm(self):
    x = MX.sym("x",Sparsity.lower(3))
    f = Function("f", [x],[tril2symm(x)])
    f_in = [0]*f.n_in();f_in[0]=DM(f.sparsity_in(0),list(range(6)))
    f_out = f(*f_in)
    self.checkarray(f_out,DM([[0,1,2],[1,3,4],[2,4,5]]))

  def test_sparsity_indexing(self):
    self.message("sparsity")

    B_ = DM([[1,2,3,4,5],[6,7,8,9,10]])
    B = MX.sym("B",2,5)

    A = IM([[1,1,0,0,0],[0,0,1,0,0]])
    A = sparsify(A)
    sp = A.sparsity()
    import copy

    def meval(m):
      f = Function("f", [B],[m])
      f_in = [0]*f.n_in();f_in[0]=B_
      f_out = f(*f_in)
      return f_out

    self.checkarray(meval(B[sp]),DM([[1,2,0,0,0],[0,0,8,0,0]]),"sparsity indexing")

    Bmod = copy.copy(B)
    Bmod[sp] = -4

    self.checkarray(meval(Bmod),DM([[-4,-4,3,4,5],[6,7,-4,9,10]]),"sparsity indexing assignement")

    Bmod = copy.copy(B)
    Bmod[sp] = 2*B

    self.checkarray(meval(Bmod),DM([[2,4,3,4,5],[6,7,16,9,10]]),"Imatrix indexing assignement")

    self.assertRaises(Exception, lambda : B[Sparsity.dense(4,4)])

  def test_symvar(self):
    a = MX.sym("a")
    b = MX.sym("b")
    c = MX.sym("c")
    e = cos(a*b) + c
    w = symvar(e)
    self.assertEqual(len(w),3)
    if GlobalOptions.getSimplificationOnTheFly():
      self.assertTrue(is_equal(w[0],a))
      self.assertTrue(is_equal(w[1],b))
      self.assertTrue(is_equal(w[2],c))

  @known_bug()
  def test_vertcat_empty(self):
    a = MX(DM(0,2))
    v = vertcat(*[a,a])

    self.assertEqual(v.size1(),0)
    self.assertEqual(v.size2(),2)

    a = MX(DM(2,0))
    v = vertcat(*[a,a])

    self.assertEqual(v.size1(),4)
    self.assertEqual(v.size2(),0)

  def test_jacobian_empty(self):
    x = MX.sym("x",3)

    s = jacobian(DM(0,0),x).shape
    self.assertEqual(s[0],0)
    self.assertEqual(s[1],3)

    s = jacobian(x,MX.sym("x",0,4)).shape
    self.assertEqual(s[0],3)
    self.assertEqual(s[1],0)

  def test_mul_sparsity(self):

    N = 10
    x = MX.sym("x",N,N)
    y = MX.sym("y",N,N)

    x_ = self.randDM(N,N)
    y_ = self.randDM(N,N)

    filt = Sparsity.diag(N)+Sparsity.triplet(N,N,[1],[3])

    f = Function("f", [x,y],[mtimes(x,y)])
    f_in = (x_, y_)
    g = Function("g", [x,y],[mac(x,y,MX.zeros(filt))])
    g_in = (x_, y_)

    f_out = f(*f_in)
    g_out = g(*g_in)

    self.checkarray(IM.ones(filt),IM.ones(g.sparsity_out(0)))

    self.checkarray(f_out[filt],g_out)

  def test_mul_zero_wrong(self):
    with self.assertRaises(RuntimeError):
      mtimes(MX.sym("X",4,5),MX.zeros(3,2))

  def test_vertsplit(self):
    a = MX.sym("X",Sparsity.lower(5))
    v = vertsplit(a,[0,2,4,5])

    Nr = int(5*6/2)

    f = Function("f", [a],v)
    f_in = [0]*f.n_in();f_in[0]=DM(f.sparsity_in(0),list(range(Nr)))

    f_out = f.call(f_in)
    v = [f_out[i] for i in range(len(v))]

    self.assertEqual(len(v),3)
    self.checkarray(v[0],DM([[0,0,0,0,0],[1,5,0,0,0]]))
    self.checkarray(v[1],DM([[2,6,9,0,0],[3,7,10,12,0]]))
    self.checkarray(v[2],DM([[4,8,11,13,14]]))

    v = vertsplit(a)

    f = Function("f", [a],v)
    f_in = [0]*f.n_in();f_in[0]=DM(f.sparsity_in(0),list(range(Nr)))

    f_out = f.call(f_in)
    v = [f_out[i] for i in range(len(v))]

    self.assertEqual(len(v),a.size1())
    self.checkarray(v[0],DM([[0,0,0,0,0]]))
    self.checkarray(v[1],DM([[1,5,0,0,0]]))
    self.checkarray(v[2],DM([[2,6,9,0,0]]))
    self.checkarray(v[3],DM([[3,7,10,12,0]]))
    self.checkarray(v[4],DM([[4,8,11,13,14]]))

    v = vertsplit(a,2)

    f = Function("f", [a],v)
    f_in = [0]*f.n_in();f_in[0]=DM(f.sparsity_in(0),list(range(Nr)))

    f_out = f.call(f_in)
    v = [f_out[i] for i in range(len(v))]

    self.assertEqual(len(v),3)
    self.checkarray(v[0],DM([[0,0,0,0,0],[1,5,0,0,0]]))
    self.checkarray(v[1],DM([[2,6,9,0,0],[3,7,10,12,0]]))
    self.checkarray(v[2],DM([[4,8,11,13,14]]))

    v = vertsplit(a,[0,0,3,a.size1()])

    f = Function("f", [a],v)
    f_in = [0]*f.n_in();f_in[0]=DM(f.sparsity_in(0),list(range(Nr)))

    f_out = f(*f_in)
    V = [f_out[i] for i in range(len(v))]

    self.assertEqual(len(v),3)
    self.assertEqual(v[0].size1(),0)
    self.assertEqual(v[0].size2(),5)  # why not 5?
    self.checkarray(V[1],DM([[0,0,0,0,0],[1,5,0,0,0],[2,6,9,0,0]]))
    self.checkarray(V[2],DM([[3,7,10,12,0],[4,8,11,13,14]]))

  def test_horzsplit(self):
    a = MX.sym("X",Sparsity.lower(5))
    v = horzsplit(a,[0,2,4,5])

    Nr = int(5*6/2)
    f = Function("f", [a],v)
    f_in = DM(f.sparsity_in(0),list(range(Nr)))

    f_out = f(f_in)
    v = [f_out[i] for i in range(len(v))]
    self.assertEqual(len(v),3)
    self.checkarray(v[0],DM([[0,0],[1,5],[2,6],[3,7],[4,8]]))
    self.checkarray(v[1],DM([[0,0],[0,0],[9,0],[10,12],[11,13]]))
    self.checkarray(v[2],DM([[0],[0],[0],[0],[14]]))

    v = horzsplit(a)

    f = Function("f", [a],v)
    f_in = DM(f.sparsity_in(0),list(range(Nr)))
    f_out = f(f_in)
    v = [f_out[i] for i in range(len(v))]
    self.assertEqual(len(v),a.size1())
    self.checkarray(v[0],DM([0,1,2,3,4]))
    self.checkarray(v[1],DM([0,5,6,7,8]))
    self.checkarray(v[2],DM([0,0,9,10,11]))
    self.checkarray(v[3],DM([0,0,0,12,13]))
    self.checkarray(v[4],DM([0,0,0,0,14]))

    v = horzsplit(a,2)

    f = Function("f", [a],v)
    f_in = DM(f.sparsity_in(0),list(range(Nr)))
    f_out = f(f_in)
    v = [f_out[i] for i in range(len(v))]

    self.assertEqual(len(v),3)
    self.checkarray(v[0],DM([[0,0],[1,5],[2,6],[3,7],[4,8]]))
    self.checkarray(v[1],DM([[0,0],[0,0],[9,0],[10,12],[11,13]]))
    self.checkarray(v[2],DM([[0],[0],[0],[0],[14]]))

    v = horzsplit(a,[0,0,3,a.size2()])
    f = Function("f", [a],v)
    f_in = DM(f.sparsity_in(0),list(range(Nr)))
    f_out = f(f_in)
    V = [f_out[i] for i in range(len(v))]

    self.assertEqual(len(v),3)
    self.assertEqual(v[0].size1(),5)
    self.assertEqual(v[0].size2(),0)
    self.checkarray(V[1],DM([[0,0,0],[1,5,0],[2,6,9],[3,7,10],[4,8,11]]))
    self.checkarray(V[2],DM([[0,0],[0,0],[0,0],[12,0],[13,14]]))

  def test_blocksplit(self):
    a = MX.sym("X",Sparsity.lower(5))
    v = blocksplit(a,[0,2,4,5],[0,1,3,5])

    Nr = int(5*6/2)
    fs = [Function("fs", [a],vr) for vr in v]
    v = [fs[i](DM(fs[i].sparsity_in(0),list(range(Nr)))) for i in range(3)]

    self.checkarray(v[0][0],DM([0,1]))
    self.checkarray(v[0][1],DM([[0,0],[5,0]]))
    self.checkarray(v[1][0],DM([2,3]))
    self.checkarray(blockcat(v),DM(fs[0].sparsity_in(0),list(range(Nr))))

  def test_mxnulloutput(self):
     a = MX(5,0)
     b = MX.sym("x",2)

     f = Function("f", [b],[a])
     c = f(b)

     self.assertEqual(c.size1(),5)
     self.assertEqual(c.size2(),0)

     c = f.call([b],True)[0]

     self.assertEqual(c.size1(),5)
     self.assertEqual(c.size2(),0)

     a = MX(0,0)
     b = MX.sym("x",2)

     f = Function("f", [b],[a])
     c = f(b)

     self.assertEqual(c.size1(),0)
     self.assertEqual(c.size2(),0)

     c = f.call([b],True)[0]

     self.assertEqual(c.size1(),0)
     self.assertEqual(c.size2(),0)

  def test_mxnull(self):
     a = MX(5,0)
     b = MX(0,3)

     c = mtimes(a,b)

     self.assertEqual(c.nnz(),0)

     a = MX(5,3)
     b = MX(3,4)

     c = mtimes(a,b)

     self.assertEqual(c.nnz(),0)

  def  test_mxnullop(self):
    c = MX(0,0)
    x = MX.sym("x",2,3)

    with self.assertRaises(RuntimeError):
      d = x + c

    with self.assertRaises(RuntimeError):
      d = x / c

  @slow()
  @memory_heavy()
  def test_MX_shapes(self):
      self.message("MX unary operations")

      #self.checkarray(DM(Sparsity.lower(4),1),DM(Sparsity.dense(4,4),1))

      for sp in [Sparsity.dense(0,0),Sparsity.dense(0,2),Sparsity.dense(2,0),Sparsity.dense(1,1),Sparsity.dense(2,2), Sparsity(4,3,[0,2,2,3],[1,2,1])]:
        for v in [0,1,0.2]:
          x_ = DM(sp,v)

          xx = MX.sym("x",sp.size1(),sp.size2())
          x=xx[sp]

          for (casadiop, numpyop,name, flags) in self.pool.zip():
            if 'nozero' in flags and v==0: continue
            r = casadiop([x])
            f = Function("f", [xx],[r])
            rr = f(v)
            self.checkarray(rr,numpyop(x_))

            a = IM(f.sparsity_out(0),1)
            b = IM.ones(DM(numpyop(x_)).sparsity())

            c = b-a
            if c.nnz()>0:
              # At least as sparse as DM calculus
              self.assertTrue(min(c.nonzeros())>=0,str([sp,v,name]))

      for sp in [Sparsity(1,1),Sparsity.dense(1,1),Sparsity(3,4),Sparsity.dense(3,4), Sparsity(4,3,[0,2,2,3],[1,2,1]).T]:
        for v1 in [0,1,0.2,-0.2]:
          x1_ = DM(sp,v1)
          xx1 = MX.sym("x",sp.size1(),sp.size2())
          x1=xx1[sp]
          xx1s = SX.sym("x",sp.size1(),sp.size2())
          x1s=xx1s[sp]
          for sp2 in [Sparsity(1,1),Sparsity.dense(1,1),Sparsity(3,4),Sparsity.dense(3,4), Sparsity(4,3,[0,2,2,3],[1,2,1]).T]:
            for v2 in [0,1,0.2,-0.2]:
              x2_ = DM(sp2,v2)
              xx2 = MX.sym("x",sp2.size1(),sp2.size2())
              x2=xx2[sp2]
              xx2s = SX.sym("x",sp2.size1(),sp2.size2())
              x2s=xx2s[sp2]
              for (casadiop, numpyop,name, flags) in self.matrixbinarypool.zip():
                if ("mul" in name or "mtimes" in name) and (sp.numel()==1 or sp2.numel()==1): continue
                r = casadiop([x1,x2])
                f = Function("f", [xx1,xx2],[r])
                f_in = [v1, v2]
                f_out = f(*f_in)
                g = Function("g", [xx1,xx2],[r])
                g_out = g(v1, v2)

                self.checkarray(f_out,numpyop([x1_,x2_]),str([sp,sp2,v1,v2,x1_,x2_,name]))


                if "mul" not in name:
                  a = IM.ones(f.sparsity_out(0))
                  b = IM.ones(g.sparsity_out(0))

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

      for sp in [Sparsity.dense(0,0),Sparsity.dense(0,2),Sparsity.dense(2,0),Sparsity.dense(1,1),Sparsity.dense(2,2), Sparsity(4,3,[0,2,2,3],[1,2,1])]:
        for v in [0,1,0.2]:
          x_ = DM(sp,v)

          x=MX(sp,v)

          for (casadiop, numpyop,name, flags) in self.pool.zip():
            if 'nozero' in flags and (v==0 or not sp.is_dense()): continue
            r = casadiop([x])
            print(r)
            self.assertTrue(r.is_constant())

            self.checkarray(r.to_DM(),numpyop(x_),str([x_,name]))

            a = IM.ones(r.to_DM().sparsity())
            b = IM.ones(DM(numpyop(x_)).sparsity())

            c = b-a
            if c.nnz()>0:
              # At least as sparse as DM calculus
              self.assertTrue(min(c.nonzeros())>=0,str([sp,v,name]))

      for sp in [Sparsity.dense(1,1),Sparsity(1,1),Sparsity(3,4),Sparsity.dense(3,4), Sparsity(4,3,[0,2,2,3],[1,2,1]).T]:
        for v1 in [0,1,0.2,-0.2]:
          x1_ = DM(sp,v1)
          x1=MX(sp,v1)
          for sp2 in [Sparsity.dense(1,1),Sparsity(1,1),Sparsity(3,4),Sparsity.dense(3,4), Sparsity(4,3,[0,2,2,3],[1,2,1]).T]:
            for v2 in [0,1,0.2,-0.2]:
              x2_ = DM(sp2,v2)
              x2=MX(sp2,v2)
              for (casadiop, numpyop,name, flags) in self.matrixbinarypool.zip():
                if ("mul" in name or "mtimes" in name) and (sp.numel()==1 or sp2.numel()==1): continue
                r = casadiop([x1,x2])
                f = Function("f", [],[r]) # Should not be needed -> constant folding
                f_out = f.call([])


                self.checkarray(f_out[0],numpyop([x1_,x2_]),str([sp,sp2,v1,v2,name]))
                if "mul" not in name and "mtimes" not in name:
                  a = IM.ones(f.sparsity_out(0))
                  b = IM.ones(DM(numpyop([x1_,x2_])).sparsity())

                  c = b-a
                  if c.nnz()>0:
                    # At least as sparse as DM calculus
                    self.assertTrue(min(c.nonzeros())>=0,str([sp,sp2,v1,v2,name]))

  def test_graph_substitute(self):
    x=MX.sym("X",4,4)
    y=MX.sym("Y",4,4)
    b=MX.sym("B",4,4)

    c = x*y
    d = b*y
    f = c+d


    C = MX.sym("C",4,4)
    f = graph_substitute(f,[c],[C])

    F = Function("F", [y,b,C],[f])
    F_out = F(1, 2, 3)

    self.checkarray(F_out,5*DM.ones(4,4))

    D = MX.sym("D",4,4)
    f = graph_substitute(f,[d],[D])

    F = Function("F", [D,C],[f])
    F_out = F(4, 5)

    self.checkarray(F_out,9*DM.ones(4,4))


  def test_matrix_expand(self):
    n = 2
    a = MX.sym("a",n,n)
    b = MX.sym("b",n,n)
    c = MX.sym("c",n,n)

    d = a+b
    e = d*c

    self.assertEqual(n_nodes(e),6)

    t0 = matrix_expand(e)

    self.assertEqual(n_nodes(t0),5)

    t1 = matrix_expand(e,[d])
    self.assertEqual(n_nodes(t1),6)

    print(e,t0,t1)


    outs = []
    for x in [e,t0,t1]:
      f = Function("f", [a,b,c],[x])

      f_in = [1.1, 2.2, 3.3]
      f_out = f(*f_in)

      outs.append(f_out)
      if len(outs)>1:
        self.checkarray(outs[0],outs[-1])

    print(outs)

  def test_kron(self):
    a = sparsify(DM([[1,0,6],[2,7,0]]))
    b = sparsify(DM([[1,0,0],[2,3,7],[0,0,9],[1,12,13]]))

    A = MX.sym("A",a.sparsity())
    B = MX.sym("B",b.sparsity())
    C = c.kron(A,B)

    f = Function("f", [A,B],[C])
    f_in = [a, b]
    f_out = f(*f_in)

    c_ = f_out

    self.assertEqual(c_.size1(),a.size1()*b.size1())
    self.assertEqual(c_.size2(),a.size2()*b.size2())
    self.assertEqual(c_.nnz(),a.nnz()*b.nnz())

    self.checkarray(c_,numpy.kron(a,b))

  def test_project(self):
    x = MX.sym("x",Sparsity.lower(3))
    y = project(x, Sparsity.lower(3).T)

    f = Function("f", [x],[y])
    f_in=DM(f.sparsity_in(0),list(range(1,int(4*3/2+1))))
    f_out = f(f_in)

    self.checkarray(f_out,DM([[1,0,0],[0,4,0],[0,0,6]]))
    self.checkarray(IM.ones(f.sparsity_out(0)),IM.ones(Sparsity.lower(3).T))

  def test_repmat(self):
    a = DM([[1,2],[3,4],[5,6]])
    self.checkarray(repmat(a,2,3),kron(DM.ones(2,3),a))

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

    A = sparsify(DM(A))

    As = MX.sym("As",A.sparsity())

    f = Function("f", [As],[densify(As.T),densify(As).T,As.T,As,densify(As)])

    f_in = [0]*f.n_in()
    f_in[0]=A
    f_out = f(*f_in)

    self.checkarray(f_out[0],A.T)
    self.checkarray(f_out[1],A.T)
    self.checkarray(f_out[2],A.T)
    self.checkarray(f_out[3],A)
    self.checkarray(f_out[4],A)

  @requiresPlugin(Linsol,"csparse")
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

    A = sparsify(DM(A))

    b = DM(list(range(15)))
    H = 5

    Bs = MX.sym("Bs",b.sparsity())
    As = MX.sym("As",A.sparsity())

    Ast = As.T

    r= Function("r", [As,Bs],[solve(Ast,Bs,"csparse")])
    R= Function("R", [As,Bs],[solve(densify(Ast),Bs,"csparse")])

    r_in = (A, b)
    r_out = r(*r_in)
    R_out = R(*r_in)

    r.sparsity_out(0).spy()
    R.sparsity_out(0).spy()

    self.checkarray(R_out,numpy.linalg.solve(A.T,b))
    self.checkarray(r_out,R_out)

  def test_depends_on(self):
    a = MX.sym("a")
    b = MX.sym("b")

    self.assertTrue(depends_on(a**2,a))
    self.assertTrue(depends_on(a,a))
    self.assertFalse(depends_on(0,a))

    ab = vertcat(*(a,b))
    self.assertTrue(depends_on(a**2,ab))
    self.assertTrue(depends_on(a,ab))
    self.assertFalse(depends_on(0,ab))
    self.assertTrue(depends_on(b**2,ab))
    self.assertTrue(depends_on(b,ab))
    self.assertTrue(depends_on(a**2+b**2,ab))
    self.assertTrue(depends_on(a+b,ab))
    self.assertTrue(depends_on(vertcat(*[0,a]),a))
    self.assertTrue(depends_on(vertcat(*[a,0]),a))
    self.assertFalse(depends_on(vertcat(*[0,b]),a))
    self.assertFalse(depends_on(vertcat(*[b,0]),a))
    self.assertTrue(depends_on(vertcat(*[a**2,b**2]),ab))
    self.assertTrue(depends_on(vertcat(*[a,0]),ab))
    self.assertTrue(depends_on(vertcat(*[0,b]),ab))
    self.assertTrue(depends_on(vertcat(*[b,0]),ab))
    self.assertFalse(depends_on(vertcat(*[0,0]),ab))

  def test_vertcat_simp(self):
    x = MX.sym("x",10)
    y = MX.sym("y")
    z = MX.sym("z")
    x_ = DM(list(range(10)))
    y_ = DM([20])
    z_ = DM([30])

    def evalvertcat(*a):
      f = Function("f", [x,y,z],[vertcat(*a)])
      f_in = [0]*f.n_in();f_in[0]=x_
      f_in[1]=y_
      f_in[2]=z_
      return f(*f_in)

    self.checkarray(evalvertcat(*vertsplit(x)),x_)
    self.checkarray(evalvertcat(*vertsplit(x)+[y]),vertcat(*[x_,y_]))
    self.checkarray(evalvertcat(*[z]+vertsplit(x)+[y] + vertsplit(x)+[z]),vertcat(*[z_,x_,y_,x_,z_]))
    self.checkarray(evalvertcat(*vertsplit(x)[:-1]),x_[:-1])
    self.checkarray(evalvertcat(*vertsplit(x)[:-1]+[y]),vertcat(*[x_[:-1],y_]))
    self.checkarray(evalvertcat(*[z]+vertsplit(x)[:-1]+[y] + vertsplit(x)[:-1]+[z]),vertcat(*[z_,x_[:-1],y_,x_[:-1],z_]))
    self.checkarray(evalvertcat(*vertsplit(x)[1:]),x_[1:])
    self.checkarray(evalvertcat(*vertsplit(x)[1:]+[y]),vertcat(*[x_[1:],y_]))
    self.checkarray(evalvertcat(*[z]+vertsplit(x)[1:]+[y] + vertsplit(x)[1:]+[z]),vertcat(*[z_,x_[1:],y_,x_[1:],z_]))
    g = vertsplit(x)[5:]+vertsplit(x)[:5]
    self.checkarray(evalvertcat(*g),vertcat(*[x_[5:],x_[:5]]))
    self.checkarray(evalvertcat(*g+[y]),vertcat(*[x_[5:],x_[:5],y_]))
    self.checkarray(evalvertcat(*[z]+g+[y] + g+[z]),vertcat(*[z_,x_[5:],x_[:5],y_,x_[5:],x_[:5],z_]))



    w = vertsplit(x,2)
    r = builtins.sum([vertsplit(i) for i in w],[])

    self.checkarray(evalvertcat(*r),x_)

    w = vertsplit(x,2)
    r = builtins.sum([vertsplit(i)+[y] for i in w],[])
    print("vertcat:", r)
    print("result:", vertcat(*r))

    w = vertsplit(x,2)
    r = builtins.sum([vertsplit(i) for i in w],[])
    print("vertcat:", r)
    print("result:", vertcat(*r+[y]))

    self.assertTrue(is_equal(vertcat(*vertsplit(x)),x))

  def test_horzcat_simp(self):
    x = MX.sym("x",1,10)
    y = MX.sym("y")
    z = MX.sym("z")
    x_ = DM(list(range(10))).T
    y_ = DM([20])
    z_ = DM([30])

    def evalhorzcat(*a):
      f = Function("f", [x,y,z],[horzcat(*a)])
      return f(x_, y_, z_)

    self.checkarray(evalhorzcat(*horzsplit(x)),x_)
    self.checkarray(evalhorzcat(*horzsplit(x)+[y]),horzcat(*[x_,y_]))
    self.checkarray(evalhorzcat(*[z]+horzsplit(x)+[y] + horzsplit(x)+[z]),horzcat(*[z_,x_,y_,x_,z_]))
    self.checkarray(evalhorzcat(*horzsplit(x)[:-1]),x_[0,:-1])
    self.checkarray(evalhorzcat(*horzsplit(x)[:-1]+[y]),horzcat(*[x_[0,:-1],y_]))
    self.checkarray(evalhorzcat(*[z]+horzsplit(x)[:-1]+[y] + horzsplit(x)[:-1]+[z]),horzcat(*[z_,x_[0,:-1],y_,x_[0,:-1],z_]))
    self.checkarray(evalhorzcat(*horzsplit(x)[1:]),x_[0,1:])
    self.checkarray(evalhorzcat(*horzsplit(x)[1:]+[y]),horzcat(*[x_[0,1:],y_]))
    self.checkarray(evalhorzcat(*[z]+horzsplit(x)[1:]+[y] + horzsplit(x)[1:]+[z]),horzcat(*[z_,x_[0,1:],y_,x_[0,1:],z_]))
    g = horzsplit(x)[5:]+horzsplit(x)[:5]
    self.checkarray(evalhorzcat(*g),horzcat(*[x_[0,5:],x_[0,:5]]))
    self.checkarray(evalhorzcat(*g+[y]),horzcat(*[x_[0,5:],x_[0,:5],y_]))
    self.checkarray(evalhorzcat(*[z]+g+[y] + g+[z]),horzcat(*[z_,x_[0,5:],x_[0,:5],y_,x_[0,5:],x_[0,:5],z_]))


    w = horzsplit(x,2)
    r = builtins.sum([horzsplit(i) for i in w],[])

    self.checkarray(evalhorzcat(*r),x_)

    w = horzsplit(x,2)
    r = builtins.sum([horzsplit(i)+[y] for i in w],[])
    print("vertcat:", r)
    print("result:", horzcat(*r))

    w = horzsplit(x,2)
    r = builtins.sum([horzsplit(i) for i in w],[])
    print("vertcat:", r)
    print("result:", horzcat(*r+[y]))

    self.assertTrue(is_equal(horzcat(*horzsplit(x)),x))

  def test_vertsplit_simp(self):

    dvars = [MX.sym("abcdefghijklm"[i]) for i in range(5) ]
    dvars_ = list(range(5))

    zz = MX.sym("zz",2)
    zz_ = DM([11,12])
    y = MX.sym("y")
    z = MX.sym("z")
    y_ = DM([20])
    z_ = DM([30])

    aa = MX.sym("aa",5)
    aa_ = list(range(100,105))

    def evalvertsplit(a,*args):
      print(vertsplit(a,*args))
      f = Function("f", dvars+[y,z,zz,aa],vertsplit(a,*args))
      f_in = [0]*f.n_in()
      for i in range(5):
        f_in[i]=dvars_[i]
      f_in[5+0]=y_
      f_in[5+1]=z_
      f_in[5+2]=zz_
      f_in[5+3]=aa_
      f_out = f.call(f_in)
      return [f_out[i] for i in range(f.n_out())]

    s= evalvertsplit(vertcat(*[y]+dvars+[z]))
    self.checkarray(s[0],y_)
    for i in range(5):
      self.checkarray(s[1+i],dvars_[i])
    self.checkarray(s[6],z_)

    s= evalvertsplit(vertcat(*[y]+dvars+[z]),2)

    self.checkarray(s[0],vertcat(*[y_,dvars_[0]]))
    self.checkarray(s[1],vertcat(*[dvars_[1],dvars_[2]]))
    self.checkarray(s[2],vertcat(*[dvars_[3],dvars_[4]]))
    self.checkarray(s[3],vertcat(*[z_]))

    s= evalvertsplit(vertcat(*[y,zz,z,zz]),2)

    self.checkarray(s[0],vertcat(*[y_,zz_[0]]))
    self.checkarray(s[1],vertcat(*[zz_[1],z_]))
    self.checkarray(s[2],zz_)

    s= evalvertsplit(vertcat(*[y,zz,z,zz]),3)

    self.checkarray(s[0],vertcat(*[y_,zz_[0],zz_[1]]))
    self.checkarray(s[1],vertcat(*[z_,zz_[0],zz_[1]]))

    s= evalvertsplit(vertcat(*[zz,zz]),2)
    self.checkarray(s[0],zz_)
    self.checkarray(s[1],zz_)

    s= evalvertsplit(vertcat(*[zz]+dvars))
    self.checkarray(s[0],zz_[0])
    self.checkarray(s[1],zz_[1])

    for i in range(5):
      self.checkarray(s[2+i],dvars_[i])

    s= evalvertsplit(vertcat(*dvars+[aa]),5)
    self.checkarray(s[0],DM(dvars_))
    self.checkarray(s[1],DM(aa_))

    s= evalvertsplit(vertcat(*dvars+[aa]),4)
    self.checkarray(s[0],DM(dvars_[:4]))
    self.checkarray(s[1],DM([dvars_[-1]]+aa_[:3]))
    self.checkarray(s[2],DM(aa_[3:]))

    s= evalvertsplit(vertcat(*dvars+[aa]),6)
    self.checkarray(s[0],DM(dvars_+[aa_[0]]))
    self.checkarray(s[1],DM(aa_[1:]))

    for i in range(5):
      self.assertTrue(is_equal(vertsplit(vertcat(*dvars))[i],dvars[i]))

  def test_horzsplit_simp(self):

    dvars = [MX.sym("abcdefghijklm"[i]) for i in range(5) ]
    dvars_ = list(range(5))

    zz = MX.sym("zz",1,2)
    zz_ = DM([11,12]).T
    y = MX.sym("y")
    z = MX.sym("z")
    y_ = DM([20])
    z_ = DM([30])

    aa = MX.sym("aa",1,5)
    aa_ = list(range(100,105))

    def evalhorzsplit(a,*args):
      print(horzsplit(a,*args))
      f = Function("f", dvars+[y,z,zz,aa],horzsplit(a,*args))
      f_in = [0]*f.n_in()
      for i in range(5):
        f_in[i]=dvars_[i]
      f_in[5+0]=y_
      f_in[5+1]=z_
      f_in[5+2]=zz_
      f_in[5+3]=aa_
      f_out = f(*f_in)
      return [f_out[i] for i in range(f.n_out())]

    s= evalhorzsplit(horzcat(*[y]+dvars+[z]))
    self.checkarray(s[0],y_)
    for i in range(5):
      self.checkarray(s[1+i],dvars_[i])
    self.checkarray(s[6],z_)

    s= evalhorzsplit(horzcat(*[y]+dvars+[z]),2)

    self.checkarray(s[0],vertcat(*[y_,dvars_[0]]).T)
    self.checkarray(s[1],vertcat(*[dvars_[1],dvars_[2]]).T)
    self.checkarray(s[2],vertcat(*[dvars_[3],dvars_[4]]).T)
    self.checkarray(s[3],vertcat(*[z_]).T)

    s= evalhorzsplit(horzcat(*[y,zz,z,zz]),2)

    self.checkarray(s[0],vertcat(*[y_,zz_[0,0]]).T)
    self.checkarray(s[1],vertcat(*[zz_[0,1],z_]).T)
    self.checkarray(s[2],zz_)

    s= evalhorzsplit(horzcat(*[y,zz,z,zz]),3)

    self.checkarray(s[0],vertcat(*[y_,zz_[0,0],zz_[0,1]]).T)
    self.checkarray(s[1],vertcat(*[z_,zz_[0,0],zz_[0,1]]).T)

    s= evalhorzsplit(horzcat(*[zz,zz]),2)
    self.checkarray(s[0],zz_)
    self.checkarray(s[1],zz_)

    s= evalhorzsplit(horzcat(*[zz]+dvars))
    self.checkarray(s[0],zz_[0,0])
    self.checkarray(s[1],zz_[0,1])

    for i in range(5):
      self.checkarray(s[2+i],dvars_[i])

    s= evalhorzsplit(horzcat(*dvars+[aa]),5)
    self.checkarray(s[0],DM(dvars_).T)
    self.checkarray(s[1],DM(aa_).T)

    s= evalhorzsplit(horzcat(*dvars+[aa]),4)
    self.checkarray(s[0],DM(dvars_[:4]).T)
    self.checkarray(s[1],DM([dvars_[-1]]+aa_[:3]).T)
    self.checkarray(s[2],DM(aa_[3:]).T)

    s= evalhorzsplit(horzcat(*dvars+[aa]),6)
    self.checkarray(s[0],DM(dvars_+[aa_[0]]).T)
    self.checkarray(s[1],DM(aa_[1:]).T)

    for i in range(5):
      self.assertTrue(is_equal(horzsplit(horzcat(*dvars))[i],dvars[i]))

  def test_vertsplit_derivative(self):
    m = MX.sym("X",10)

    f = Function("f", [m],[vertsplit(m)[0]])

    f.reverse(1)

  def test_MX_const_sp(self):
    x = MX.sym("x",4,1)

    sp = Sparsity.triplet(3,3,[0,1,2,2],[0,0,1,2])

    f = Function("f", [x],[x.nz[IM(sp,list(range(sp.nnz())))]])

    g = Function("g", [x],[MX(sp,x)])

    f_in = [0]*f.n_in();f_in[0]=DM(list(range(1,5)))

    self.checkfunction(f,g,inputs=f_in)

  def test_reshape_sp(self):
    x = MX.sym("x",4,1)

    f = Function("f", [x],[x.reshape((2,2))])

    sx = SX.sym("x",4,1)

    g = Function("g", [sx],[sx.reshape((2,2))])

    f_in = [0]*f.n_in();f_in[0]=DM(list(range(1,5)))

    self.checkfunction(f,g,inputs=f_in)

  def test_issue1041(self):
    x = MX.sym("x",2)

    y = vertsplit(x,[0,1,2])[1]

    f = Function("f", [x],[y])

    H = f.hessian_old(0, 0)

  def test_bug_1042(self):

    x = MX.sym('x',2,1)

    mf = Function("mf", [x],[x*x[0,0]])

    mfunction = mf.expand('expand_'+mf.name())

    mfg = mf.reverse(1)

    mfunctiong = mfunction.reverse(1)

    f_in = [0, 5, DM([1,2])]

    self.checkfunction(mfg,mfunctiong,inputs=f_in)

  def test_bug_1042bis(self):
    x = MX.sym('x',2,1)
    a = MX.sym("ax",2,1)
    i1 = x[0,0]
    z = i1*x
    i3 = i1*a
    i3= c.dot(x,a)
    d = Function("d", [x,a],[z,i3])

    dx = d.expand('expand_'+d.name())
    dx_in = [0]*dx.n_in();dx_in[0]=DM([1,2])
    dx_in[1]=DM([3,4])

    self.checkfunction(d,dx,inputs=dx_in)

  def test_bug_1042tris(self):
    x = MX.sym('x',2,1)
    a = MX.sym("ax",2,1)
    d = Function("d", [x,a],[c.dot(x,a)])

    dx = d.expand('expand_'+d.name())
    dx_in = [0]*dx.n_in();dx_in[0]=DM([1,2])
    dx_in[1]=DM([3,4])

    self.checkfunction(d,dx,inputs=dx_in)

  def test_bug_1046(self):
    x = MX.sym('x',1,1)
    y = MX.sym('y',1,1)
    z = jacobian(x,y)

    self.assertTrue(z.nnz()==0)

  def test_singularcat(self):

    for c in [MX,SX,DM]:
      x0 = c.zeros(10,0)

      x1s = vertsplit(x0, [0,5,10])

      for x in x1s:
        self.checkarray(x.shape,(5,0))


      x2 = vertcat(*x1s)
      self.checkarray(x2.shape,(10,0))

      x2 = vertcat(*[c.zeros(0,0)] + x1s + [c.zeros(0,0)])
      self.checkarray(x2.shape,(10,0))
      
      x0 = c.zeros(0,1)
      x2 = vertcat(x0,c.zeros(0,0),x0)
      self.checkarray(x2.shape,(0,1))

    for c in [MX,SX,DM]:
      x0 = c.zeros(0,10)

      x1s = horzsplit(x0, [0,5,10])

      for x in x1s:
        self.checkarray(x.shape,(0,5))

      x2 = horzcat(*x1s)
      self.checkarray(x2.shape,(0,10))

      x2 = horzcat(*[c.zeros(0,0)] + x1s + [c.zeros(0,0)])
      self.checkarray(x2.shape,(0,10))

      x0 = c.zeros(1,0)
      x2 = horzcat(x0,c.zeros(0,0),x0)
      self.checkarray(x2.shape,(1,0))

    for c in [MX,SX,DM]:
      x0 = c.zeros(10,0)

      x1s = vertsplit(x0, [0,5,10])

      x0 = c.zeros(0,10)
      x1st = horzsplit(x0, [0,5,10])

      x2 = diagcat(*x1s)
      self.checkarray(x2.shape,(10,0))

      x2 = diagcat(*[c.zeros(0,0)] + x1s + [c.zeros(0,0)])
      self.checkarray(x2.shape,(10,0))

      x2 = diagcat(*x1st)
      self.checkarray(x2.shape,(0,10))

      x2 = diagcat(*[c.zeros(0,0)] + x1st + [c.zeros(0,0)])
      self.checkarray(x2.shape,(0,10))

      x2 = diagcat(*x1s+x1st)
      self.checkarray(x2.shape,(10,10))

      x2 = diagcat(*[c.zeros(0,0)] + x1s+x1st + [c.zeros(0,0)])
      self.checkarray(x2.shape,(10,10))

  def test_exprjac(self):


    for fun in [lambda x: x[0]*x[1],lambda x: x[0]*sin(x[1])]:
      def hessian1(f, x): return hessian(f, x)[0]
      for op in [c.gradient, jacobian, hessian1]:
        print(fun, op)
        x = MX.sym("x",2)
        print(fun(x))
        print(op(fun(x),x))
        f = Function("f", [x],[op(fun(x),x)])

        x = SX.sym("x",2)
        fr = Function("fr", [x],[op(fun(x),x)])

        self.checkfunction(f,fr,inputs=[0])

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


          A = MX.sym("A", np.product(dim_a))
          B = MX.sym("B", np.product(dim_b))
          C = MX.sym("C", int(np.product(dim_c)),1)

          A_sx = SX.sym("A", np.product(dim_a))
          B_sx = SX.sym("B", np.product(dim_b))
          C_sx = SX.sym("C", int(np.product(dim_c)),1)

          np.random.seed(0)

          A_ = np.random.random(A.shape)
          B_ = np.random.random(B.shape)
          C_ = np.random.random(C.shape)



          out = casadi.einstein(A, B, C, dim_a, dim_b, dim_c, ind_a, ind_b, ind_c)
          out_sx = casadi.einstein(A_sx, B_sx, C_sx, dim_a, dim_b, dim_c, ind_a, ind_b, ind_c)
          f = Function('f',[A,B,C],[out],{"ad_weight_sp": 0})
          frev = Function('f',[A,B,C],[out],{"ad_weight_sp": 1})
          fr = Function('fr',[A,B,C],[my_einstein(A, B, C, dim_a, dim_b, dim_c, ind_a, ind_b, ind_c)])
          f_sx = Function('f',[A_sx,B_sx,C_sx],[out_sx])

          fsx = f.expand()
          self.checkfunction(f, fr, inputs=[A_,B_,C_])
          self.checkfunction(fsx, fr, inputs=[A_,B_,C_])
          self.check_codegen(f, inputs=[A_,B_,C_])
          self.checkfunction(fsx, f_sx, inputs=[A_,B_,C_])

          for i in range(3):
            self.check_sparsity(f.sparsity_jac(i, 0), fsx.sparsity_jac(i, 0))
            self.check_sparsity(frev.sparsity_jac(i, 0), fsx.sparsity_jac(i, 0))


        def my_einstein(A, B, C, dim_a, dim_b, dim_c, ind_a, ind_b, ind_c):
          try:
            R = DM(C)
          except:
            R = MX(C)

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
    L = [MX(Sparsity(1,1)),MX(Sparsity(2,1)), MX.sym("x",1,1), MX.sym("x", Sparsity(1,1)), DM(1), DM(Sparsity(1,1),1), DM(Sparsity(2,1),1), DM(Sparsity.dense(2,1),1)]

    for a in L:
      for b in L:
        c = a*b

        if a.nnz()==0 or b.nnz()==0:
          self.assertTrue(c.nnz()==0)
        else:
          self.assertTrue(c.nnz()>0)

  @requiresPlugin(Linsol,"lapackqr")
  def test_solve(self):
    N = 3
    nrhs = 50
    np.random.seed(0)
    A = np.random.random((3,3))
    B = np.random.random((3,50))


    C = solve(A, B, "lapackqr", {"max_nrhs": 50})
    C1 = solve(A, B, "lapackqr", {"max_nrhs": 20})
    C2 = solve(A, B, "lapackqr")

    self.checkarray(C, C1)
    self.checkarray(C1, C2)
  def test_sparse_lt(self):
    x = MX.sym("x",Sparsity.lower(5))
    self.assertEqual((x>0).nnz(),5*6/2)
    self.assertEqual((x>=0).nnz(),5*5)
  def test_inv(self):
   np.random.seed(0)

   for X in [SX, MX]:
     A  = X.sym("x",3,3)
     Av = np.random.random((3,3))
     f = Function('f',[A],[inv(A),inv(DM(Av)),A.__mpower__(-1), DM(Av).__mpower__(-1)])
     out = f(Av)
     for o in out:
       self.checkarray(o, np.linalg.inv(np.array(Av)))
    

  def test_interp1d(self):
    v = [7,3,4,-3]
    vs = MX.sym("x",4,1)

    xq = [10,0,1,2,4,8,7,5,1.5]
    x = [1,2,4,8]
    F = Function("f",[vs],[interp1d(x,vs,xq)])

    self.checkarray(F(v),np.interp(xq,x,v))

    F = Function("f",[vs],[interp1d(x,vs,xq,"floor")])

    self.checkarray(F(v),DM([-3,7,7,3,4,-3,4,4,7]))

    F = Function("f",[vs],[interp1d(x,vs,xq,"ceil")])

    self.checkarray(F(v),DM([-3,7,7,3,4,-3,-3,-3,3]))

    v = [7,3,4,-3]
    vs = MX.sym("x",4,1)

    xq = [10,0,1,2,3,4,3.5,3.25,1.5]
    x = [1,2,3,4]
    F = Function("f",[vs],[interp1d(x,vs,xq,"linear",True)])

    self.checkarray(F(v),np.interp(xq,x,v))
    
  def test_bilin_etc(self):
    x = MX.sym("x",3,3)
    y = MX.sym("y",3,1)
    z = MX.sym("z",3,1)
    
    import numpy
    numpy.random.seed(42)
    x0 = numpy.random.random((3,3))
    y0 = numpy.random.random((3,1))
    z0 = numpy.random.random((3,1))
    
    for e in [(bilin(x,y,z),mtimes(mtimes(y.T,x),z)),(rank1(x,0.3,y,z),x+0.3*mtimes(y,z.T))]:
      f = Function('f',[x,y,z],[e[0]])
      fref = Function('fref',[x,y,z],[e[1]])
      f_sx = f.expand()    
      self.checkfunction(f,fref,inputs=[x0,y0,z0])
      self.checkfunction(f_sx,fref,inputs=[x0,y0,z0])
      self.check_codegen(f,inputs=[x0,y0,z0])

  def test_det_shape(self):
    X = MX.sym("x",2,3)
    with self.assertRaises(RuntimeError):
      det(X)
    X = MX.sym("x",3,3)
    det(X)
    
    
  @known_bug()
  def test_det(self):
    X = MX.sym("x",3,3)
    x = SX.sym("x",3,3)

    import numpy
    numpy.random.seed(42)
    x0 = numpy.random.random((3,3))
    
    f = Function('f',[x],[det(x)])
    F = Function('F',[X],[det(X)])
    self.checkfunction(f,F,inputs=[x0])

  def test_mtimes_mismatch_segfault(self):
    with self.assertInException("incompatible dimensions"):
      mtimes(DM(Sparsity.lower(5)),MX.sym('x',100))
    
  def test_monitor(self):
    x = MX.sym("x")
    y = sqrt(x.monitor("hey"))

    f = Function('f',[x],[y])
    with capture_stdout() as out:
      f(3)
    self.assertTrue(out[0]=="hey:\n[3]\n")

    with capture_stdout() as out2:
      self.check_codegen(f,inputs=[3])
    if args.run_slow:
      self.assertTrue(out2[0]=="hey:\n[3]\n")

  def test_codegen_specials(self):
    x = MX.sym("x")
    y = MX.sym("y")

    for z in [ x**2, if_else(y>0,2*x,x*y), fmin(x,y), fmax(x,y), sign(x*y)]:
      f = Function('f',[x,y],[z])
      self.check_codegen(f,inputs=[1,2])
      self.check_codegen(f,inputs=[1,-2])


  @memory_heavy()    
  def test_getsetnonzeros(self):
    import numpy
    numpy.random.seed(42)
    for S in [Sparsity.lower(5),Sparsity.dense(5,5)]:
      M = MX.sym("X",S.nnz())
      
      m = sin(MX(S,M))
      for i in [0,2]:
        for ind in [(i,i),(1,i),(i,1),(slice(None),i),(i,slice(None)),(slice(i,i+2),slice(i,i+2)),(slice(i,i+2),slice(None)),([i,i],[0,0]),([],[])]:
          E = m.__getitem__(ind)
          e = cos(E)
          e = dot(e,e)
          f = Function('f',[M],[e])
          self.checkfunction(f,f.expand(),inputs=[ numpy.random.random((S.nnz(),1))])
        
          mc = m+0
          
          Y = MX.sym("y",E.nnz())
          y = MX(E.sparsity(),Y)
          mc.__setitem__(ind,y)
          e = cos(mc)
          e = dot(e,e)
          
          f = Function('f',[M,Y],[e])
          self.checkfunction(f,f.expand(),inputs=[ numpy.random.random((S.nnz(),1)), numpy.random.random((E.nnz(),1))])

  def test_evalf(self):
    x = MX.sym("x")

    p = MX.sym("p")
    f = Function('f',[x],[sin(x)])
    y = f.call([p],False,True)[0]
    y = substitute(y,p,3)

    with self.assertInException("not defined"):
      y.to_DM()
    self.checkarray(evalf(y),sin(3))
    with self.assertInException("since variables [x] are free"):
      evalf(x)


  def test_mmin(self):

    def mx_eval(f,X,x):
      y = X.sym("y",x.sparsity())
      return Function('f',[y],f.call([y],X is MX,False))

    for X in [SX,MX]:
      for mod in [lambda f:f, lambda f:f.expand(), lambda f: mx_eval(f,X,x)]:
        x = X.sym("x",2,2)

        f = mod(Function('f',[x],[mmin(x),mmax(x)]))


        [mmin_res,mmax_res] = f(DM([[-2,3],[-3,5]]))

        self.checkarray(mmin_res, -3)
        self.checkarray(mmax_res, 5)

        x = X.sym("x",Sparsity.lower(2))

        f = mod(Function('f',[x],[mmin(x),mmax(x)]))

        [mmin_res,mmax_res] = f(sparsify(DM([[-2,0],[-3,-5]])))

        self.checkarray(mmin_res, -5)
        self.checkarray(mmax_res, 0)

        [mmin_res,mmax_res] = f(sparsify(DM([[2,0],[3,5]])))

        self.checkarray(mmin_res, 0)
        self.checkarray(mmax_res, 5)

        x = X.sym("x",Sparsity(2,2))
        f = mod(Function('f',[x],[mmin(x),mmax(x)]))
        [mmin_res,mmax_res] = f(DM(2,2))

        self.checkarray(mmin_res, 0)
        self.checkarray(mmax_res, 0)

        x = X.sym("x",Sparsity(2,0))
        f = mod(Function('f',[x],[mmin(x),mmax(x)]))
        [mmin_res,mmax_res] = f(DM(2,0))

        self.assertTrue(mmin_res.is_empty())
        self.assertTrue(mmax_res.is_empty())

        x = X.sym("x",Sparsity(0,0))
        f = mod(Function('f',[x],[mmin(x),mmax(x)]))
        [mmin_res,mmax_res] = f(DM(0,0))

        self.assertTrue(mmin_res.is_empty())
        self.assertTrue(mmax_res.is_empty())

  def test_doc_expression_tools(self):
    self.assertTrue("Given a repeated matrix, computes the sum of repeated parts." in repsum.__doc__)

if __name__ == '__main__':
    unittest.main()
