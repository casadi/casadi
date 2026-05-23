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

import casadi as c
import numpy
from numpy import random, array, linalg, matrix, zeros, ones
import unittest
from types import *
from helpers import *
import casadi as ca
import numpy as np
from numpy import inf, pi

scipy_available = True
try:
	from scipy.sparse import csr_matrix
except:
	scipy_available = False

class SXtests(casadiTestCase):

  def setUp(self):
    self.pool=FunctionPool()
    self.pool.append(lambda x: ca.sqrt(x[0]),ca.sqrt,"sqrt")
    self.pool.append(lambda x: ca.sin(x[0]),ca.sin,"sin")
    self.pool.append(lambda x: ca.cos(x[0]),ca.cos,"cos")
    self.pool.append(lambda x: ca.tan(x[0]),ca.tan,"tan")
    self.pool.append(lambda x: ca.fabs(x[0]),ca.fabs,"fabs")
    self.pool.append(lambda x: ca.sign(x[0]),ca.sign,"sign")
    self.pool.append(lambda x: ca.arctan(x[0]),ca.arctan,"arctan")
    self.pool.append(lambda x: ca.arcsin(x[0]),ca.arcsin,"arcsin")
    self.pool.append(lambda x: ca.arccos(x[0]),ca.arccos,"arccos")
    self.pool.append(lambda x: ca.exp(x[0]),ca.exp,"exp")
    self.pool.append(lambda x: ca.log(x[0]),ca.log,"log")
    self.pool.append(lambda x: x[0]**0,lambda x : x**0,"x^0",flags={'nozero'})
    self.pool.append(lambda x: x[0]**1,lambda x : x**1,"^1")
    self.pool.append(lambda x: x[0]**(-2),lambda x : x**(-2),"^-2",flags={'nozero'})
    self.pool.append(lambda x: x[0]**(0.3),lambda x : x**(0.3),"^0.3")
    self.pool.append(lambda x: ca.floor(x[0]),ca.floor,"floor")
    self.pool.append(lambda x: ca.ceil(x[0]),ca.ceil,"ceil")
    self.pool.append(lambda x: ca.log1p(x[0]),np.log1p,"log1p")
    self.pool.append(lambda x: ca.expm1(x[0]),np.expm1,"log1p")
    self.Jpool=FunctionPool()
    self.Jpool.append(lambda x: ca.sqrt(x[0]),lambda x:ca.diag(1/(2.0*ca.sqrt(x))),"sqrt")
    self.Jpool.append(lambda x: ca.sin(x[0]),lambda x:ca.diag(ca.cos(x)),"sin")
    self.Jpool.append(lambda x: ca.fabs(x[0]),lambda x:ca.diag(ca.sign(x)),"fabs")
    self.Jpool.append(lambda x: ca.sign(x[0]),lambda x:ca.diag(x*0),"fabs")
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
    self.Jpool.append(lambda x: ca.log1p(x[0]),lambda x :ca.diag( 1.0/(1+x)),"log1p")
    self.Jpool.append(lambda x: ca.expm1(x[0]),lambda x :ca.diag( ca.exp(x)),"log1p")
    self.matrixpool=FunctionPool()
    self.matrixpool.append(lambda x: ca.norm_2(x[0]),linalg.norm,"norm_2")
    self.matrixbinarypool=FunctionPool()
    self.matrixbinarypool.append(lambda a: a[0]+a[1],lambda a: a[0]+a[1],"Matrix+Matrix")
    self.matrixbinarypool.append(lambda a: a[0]-a[1],lambda a: a[0]-a[1],"Matrix-Matrix")
    self.matrixbinarypool.append(lambda a: a[0]*a[1],lambda a: a[0]*a[1],"Matrix*Matrix")
    self.matrixbinarypool.append(lambda a: ca.fmax(a[0],a[1]),lambda a: ca.fmax(a[0],a[1]),"fmin")
    self.matrixbinarypool.append(lambda a: ca.fmin(a[0],a[1]),lambda a: ca.fmin(a[0],a[1]),"fmax")
    self.matrixbinarypool.append(lambda a: ca.hypot(a[0],a[1]),lambda a: np.hypot(a[0],a[1]),"hypot")
    #self.matrixbinarypool.append(lambda a: dot(a[0],trans(a[1])),lambda a: dot(a[0].T,a[1]),name="dot(Matrix,Matrix)")
    self.matrixbinarypool.append(lambda a: a[0] @ a[1].T,lambda a: np.dot(a[0],a[1].T),"dot(Matrix,Matrix.T)")

    #self.pool.append(lambda x: erf(x[0]),erf,"erf") # numpy has no erf

  def test_equivalence(self):
    x = ca.SX.sym("x")
    y = ca.SX.sym("y")    
    for expr,equiv in [(ca.log1p(x),ca.log(1+x)),(ca.expm1(x),ca.exp(x)-1),(ca.hypot(x,y),ca.sqrt(x**2+y**2))]:
      f=ca.Function("f",[x,y],[expr])
      equiv_f = ca.Function("equiv_f",[x,y],[equiv])
      self.checkfunction(f,equiv_f,inputs=[1.3,1.7])

  def test_scalarSX(self):
      x=ca.SX.sym("x")
      x0=0.738

      self.numpyEvaluationCheckPool(self.pool,[x],x0,name="scalarSX")

  def test_gradient(self):
      self.message("jacobian of SX**number")
      x=ca.SX.sym("x");
      x0=1;
      p=3 # increase to 20 to showcase ticket #56
      y=x**p;
      dx=ca.jacobian(y,x);
      dxr=p;
      self.evaluationCheck([dx],dxr,[x],x0,name="jacobian");
      dxr=1
      for i in list(range(p)):
        y=ca.jacobian(y,x)
        dxr=dxr*(p-i)


      self.evaluationCheck([y],dxr,[x],x0,name="recursive jacobian");

  def test_gradient2(self):
      self.message("jacobian of SX**SX")
      x=ca.SX.sym("x");
      p=ca.SX.sym("p");
      x0=1;
      p0=3 # increase to 20 to showcase ticket #56
      y=x**p;
      dx=ca.jacobian(y,x);
      #print dx
      dxr=p0;
      self.evaluationCheck([dx],dxr,[x,p],[x0,p0],name="jacobian");

      dxr=1
      for i in list(range(p0)):
        y=ca.jacobian(y,x)
        dxr=dxr*(p0-i)

      self.evaluationCheck([y],dxr,[x,p],[x0,p0],name="jacobian");

  def test_SXJacobian(self):
      self.message("SX(1,1) unary operation, jacobian")
      x=ca.SX.sym("x")
      x0=array([[0.738]])

      def fmod(f,x):
        return jacobian_old(f, 0, 0)

      self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="SX unary operations, jacobian",fmod=fmod)

  def test_SXJac(self):
      self.message("SX(1,1) unary operation, jac")
      x=ca.SX.sym("x")
      x0=array([[0.738]])

      def fmod(f,x):
          y = f.call(x)
          J = ca.Function('J', x, [ca.jacobian(y[0],x[0])])
          return J

      self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="SX unary operations, jac",fmod=fmod)

  def test_SXJacobians(self):
      self.message("SX(3,1) unary operation, jacobian")
      x=ca.SX.sym("x",3)
      x0=array([0.738,0.9,0.3])

      def fmod(f,x):
        return jacobian_old(f, 0, 0)

      self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="SX unary operations, jacobian",fmod=fmod)

  def test_SXJacobians2(self):
      self.message("SX(1,3) unary operation, jacobian")
      x=ca.SX.sym("x",1,3)

      x0=array([0.738,0.9,0.3])

      def fmod(f,x):
        return jacobian_old(f, 0, 0)

      self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="SX unary operations, jacobian",fmod=fmod)

  def test_SX(self):
      self.message("SX unary operations")
      x=ca.SX.sym("x",3,2)
      x0=array([[0.738,0.2],[ 0.1,0.39 ],[0.99,0.999999]])

      self.numpyEvaluationCheckPool(self.pool,[x],x0,name="SX")

      x=ca.SX.sym("x",3,3)
      x0=array([[0.738,0.2,0.3],[ 0.1,0.39,-6 ],[0.99,0.999999,-12]])
      #self.numpyEvaluationCheck(lambda x: c.det(x[0]), lambda   x: linalg.det(x),[x],x0,name="det(SX)")
      self.numpyEvaluationCheck(lambda x: ca.SX(c.det(x[0])), lambda   x: linalg.det(x),[x],x0,name="det(SX)")
      self.numpyEvaluationCheck(lambda x: c.inv(x[0]), lambda   x: linalg.inv(x),[x],x0,name="inv(SX)")

  def test_SXSparse(self):
      self.message("SX unary operations, sparse")
      x=ca.SX.sym("x")
      y=ca.SX.sym("y")
      z=ca.SX.sym("z")
      x=ca.SX(ca.Sparsity(4,3,[0,2,2,3],[1,2,1]),ca.vertcat(*[x,y,z]))
      if scipy_available:
        x0=ca.DM(ca.Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).sparse()

        self.numpyEvaluationCheckPool(self.pool,[x],array(x0.todense()),name="SX",setx0=x0,excludeflags={'nozero'})
      else:
        x0=ca.DM(ca.Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).full()

        self.numpyEvaluationCheckPool(self.pool,[x],x0,name="SX",setx0=x0,excludeflags={'nozero'})

  def test_SXbinary(self):
      self.message("SX binary operations")
      x=ca.SX.sym("x",3,2)
      y=ca.SX.sym("x",3,2)
      x0=array([[0.738,0.2],[ 0.1,0.39 ],[0.99,0.999999]])
      y0=array([[1.738,0.6],[ 0.7,12 ],[0,-6]])
      self.numpyEvaluationCheckPool(self.matrixbinarypool,[x,y],[x0,y0],name="SX")
      self.assertRaises(RuntimeError, lambda : x @ y)

  def test_SXbinary_codegen(self):
      self.message("SX binary operations")
      x=ca.SX.sym("x",4,2)
      y=ca.SX.sym("x",4,2)
      x0=array([[0.738,0.2],[ 0.1,0.39 ],[0.99,0.999999],[1,2]])
      y0=array([[1.738,0.6],[ 0.7,12 ],[0,-6],[1,2]])
      for f in self.matrixbinarypool.casadioperators:
        f = ca.Function('f',[x,y],[f([x,y])])
        self.check_codegen(f,inputs=[x0,y0],std="c99")
        self.check_codegen(f,inputs=[x0,y0],std="c89")

  def test_SXbinary_diff(self):
      self.message("SX binary operations")
      x=ca.SX.sym("x",4,2)
      y=ca.SX.sym("x",4,2)
      dx=ca.SX.sym("x",4,2)
      dy=ca.SX.sym("x",4,2)
      x0=array([[0.738,0.2],[ 0.1,0.39 ],[0.99,0.999999],[1,2]])
      y0=array([[1.738,0.6],[ 0.7,12 ],[0,-6],[1,2]])
      for f,name in zip(self.matrixbinarypool.casadioperators,self.matrixbinarypool.names):
        f = ca.Function('f',[x,y],[f([x,y])])
        f.disp(True)
        df = ca.Function('f',[x,y,dx,dy],[ca.jtimes(f(x,y),x,dx)+ca.jtimes(f(x,y),y,dy)])
        eps = 1e-7
        ca.DM.rng(0)
        dx0 = ca.DM.rand(x.sparsity())
        dy0 = dx0 if name in ["fmin","fmax"] else ca.DM.rand(y.sparsity())
        dy0 = dx0 if name in ["fmin","fmax"] else ca.DM.rand(y.sparsity())
        df_experiment = (f(x0+eps*dx0,y0+eps*dy0)-f(x0,y0))/eps
        print(df_experiment,df(x0,y0,dx0,dy0))
        self.checkarray(df(x0,y0,dx0,dy0),df_experiment,digits=4)

  def test_DMbinary(self):
      self.message("SX binary operations")
      x=ca.SX.sym("x",3,2)
      y=ca.SX.sym("x",3,2)
      x0=array([[0.738,0.2],[ 0.1,0.39 ],[0.99,0.999999]])
      y0=array([[1.738,0.6],[ 0.7,12 ],[0,-6]])
      for f,fr,label,flags in self.matrixbinarypool.zip():
        self.checkarray(f(ca.vertcat(*[x0,y0])),fr(ca.vertcat(*[x0,y0])),label)

  def test_SXbinarySparse(self):
      self.message("SX binary operations")
      x=ca.SX.sym("x")
      y=ca.SX.sym("y")
      z=ca.SX.sym("z")
      x2=ca.SX.sym("x2")
      y2=ca.SX.sym("y2")
      z2=ca.SX.sym("z2")
      xx=ca.SX(ca.Sparsity(4,3,[0,2,2,3],[1,2,1]),ca.vertcat(*[x,y,z]))
      yy=ca.SX(ca.Sparsity(4,3,[0,2,2,3],[0,2,3]),ca.vertcat(*[x2,z2,y2]))

      if scipy_available:
        x0=ca.DM(ca.Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).sparse()
        y0=ca.DM(ca.Sparsity(4,3,[0,2,2,3],[0,2,3]),[1.738,0.7,-6]).sparse()

        self.numpyEvaluationCheckPool(self.matrixbinarypool,[xx,yy],[array(x0.todense()),array(y0.todense())],name="SX",setx0=[x0,y0])
      else:
        x0=ca.DM(ca.Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).full()
        y0=ca.DM(ca.Sparsity(4,3,[0,2,2,3],[0,2,3]),[1.738,0.7,-6]).full()

        self.numpyEvaluationCheckPool(self.matrixbinarypool,[xx,yy],[x0,y0],name="SX",setx0=[x0,y0])
      self.assertRaises(RuntimeError, lambda : xx @ yy)


  @known_bug()  # Test refactoring, cf. #1436
  def test_SXslicing(self):
      self.message("SX slicing/indexing")
      x=ca.SX.sym("x",3,2)
      x0=array([[0.738,0.2],[ 0.1,0.39 ],[0.99,0.999999]])

      self.message(":dense")
      self.numpyEvaluationCheck(lambda x: ca.SX(x[0][0,0]), lambda x: matrix(x)[0,0],[x],x0,name="x[0,0]")
      self.numpyEvaluationCheck(lambda x: ca.SX(x[0][1,0]), lambda x: matrix(x)[1,0],[x],x0,name="x[1,0]")
      self.numpyEvaluationCheck(lambda x: ca.SX(x[0][0,1]), lambda x: matrix(x)[0,1],[x],x0,name="x[1,0]")
      self.numpyEvaluationCheck(lambda x: ca.SX(x[0][0,-1]), lambda x: matrix(x)[0,-1],[x],x0,name="x[0,-1]")
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

      myarray=array([0,2,3])
      mylist=list(myarray)
      #self.numpyEvaluationCheck(lambda x: x[0][mylist], lambda x: matrix([x[0,0],x[1,0],x[1,1]]).T,[x],x0,name="x[[0,2,3]]")
      self.numpyEvaluationCheck(lambda x: x[0].nz[0:2], lambda x: matrix(x.T.ravel()[0:2]).T,[x],x0,name="x[0:2] on dense matrix")
      self.numpyEvaluationCheck(lambda x: x[0].nz[1], lambda x: matrix(x.T.ravel()[1]).T,[x],x0,name="x[1]")
      self.numpyEvaluationCheck(lambda x: x[0].nz[-1], lambda x: matrix(x.ravel()[-1]).T,[x],x0,name="x[-1]")

      self.message(":sparse")

      x=ca.SX(ca.Sparsity(4,3,[0,2,2,3],[1,2,1]),ca.vertcat(*[ca.SX.sym("x"),ca.SX.sym("y"),ca.SX.sym("z")]))
      sx0=[0.738,0.39,0.99]
      x0=ca.DM(ca.Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.39,0.99]).full()
      self.numpyEvaluationCheck(lambda x: ca.SX(x[0][0,0]), lambda x: matrix(x)[0,0],[x],x0,name="x[0,0]",setx0=[sx0])
      self.numpyEvaluationCheck(lambda x: ca.SX(x[0][0,0]), lambda x: matrix(x)[0,0],[x],x0,name="x[0,0]",setx0=[sx0])
      self.numpyEvaluationCheck(lambda x: ca.SX(x[0][1,0]), lambda x: matrix(x)[1,0],[x],x0,name="x[1,0]",setx0=[sx0])
      self.numpyEvaluationCheck(lambda x: ca.SX(x[0][0,1]), lambda x: matrix(x)[0,1],[x],x0,name="x[1,0]",setx0=[sx0])
      self.numpyEvaluationCheck(lambda x: ca.SX(x[0][0,-1]), lambda x: matrix(x)[0,-1],[x],x0,name="x[0,-1]",setx0=[sx0])
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


  def test_SX1(self):
    self.message("SXFunction evaluation")
    fun=lambda x,y: [x+y,x*y,x**2+y**3]
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")
    f=ca.Function("f", [ca.vertcat(*[x,y])],[ca.vertcat(*fun(x,y))])
    L=[2,3]
    f_in = [0]*f.n_in()  # type: list

    f_in[0]=L
    f_out = f.call(f_in)
    z=f_out[0].full().squeeze()
    zr=fun(*L)
    for i in range(3):
      self.assertAlmostEqual(z[i], zr[i],10,'SXfunction output in correct')
    self.message("SXFunction jacobian evaluation")
    J = jacobian_old(f, 0, 0)
    J_in = [0]*J.n_in()  # type: list

    J_in[0]=L
    J_out = J.call(J_in)
    Jr=np.array([[1,1],[3,2],[4,27]])
    self.checkarray(J_out[0],Jr,"SXfunction jacobian evaluates incorrectly")

  def test_SX2(self):
    self.message("SXFunction evalution 2")
    fun = lambda x,y: [3-ca.sin(x*x)-y, ca.sqrt(y)*x]
    # variables
    x = ca.SX.sym("x")
    y = ca.SX.sym("y")

    # Create function
    f = fun(x,y)
    if ca.GlobalOptions.getSimplificationOnTheFly():
      self.assertEqual(str(f),'[SX(((3-sin(sq(x)))-y)), SX((sqrt(y)*x))]','SX representation is wrong '+str(f))
    else:
      self.assertEqual(str(f),'[SX(((3-sin((x*x)))-y)), SX((sqrt(y)*x))]','SX representation is wrong '+str(f))
    fcn = ca.Function("fcn", [ca.vertcat(*[x,y])],[ca.vertcat(*f)])

    self.assertEqual(repr(fcn),'Function(fcn:(i0[2])->(o0[2]) SXFunction)','SX representation is wrong')

    # Pass inputs
    L=[2,3]
    fcn_in = [0]*fcn.n_in()  # type: list

    fcn_in[0]=L

    # Evaluate numerically
    fcn_out = fcn.call(fcn_in)

    # Get the results
    res = tuple(fcn_out[0].nonzeros())
    self.assertAlmostEqual(res[0], fun(*L)[0],10,'SXfunction evaluation wrong')
    self.assertAlmostEqual(res[1], fun(*L)[1],10,'SXfunction evaluation wrong')

  def test_const_folding_on_the_fly(self):
    x = ca.SX.sym('x')

    for a in [2,7]:
        for b in [2,7]:
            ref = (a*b)*x
            self.assertEqual(str(ref),str(a*(b*x)))
            self.assertEqual(str(ref),str((b*x)*a))

  def test_SX_func(self):
    self.message("Function constructors")
    x0=ca.SX.sym("x")
    x1=ca.SX.sym("x")
    x2=ca.SX.sym("x")
    x3=ca.SX.sym("x")
    x4=ca.SX.sym("x")
    x5=ca.SX.sym("x")
    x6=ca.SX.sym("x")
    y=ca.SX.sym("y",2,3)

    f=ca.Function("f", [y],[y])
    self.checkarray(f.size_in(0),(2,3),"Function constructors")
    self.checkarray(f.size_out(0),(2,3),"Function constructors")

    self.assertRaises(TypeError if systemswig else NotImplementedError,lambda: ca.Function("f", y,[y,y]))  # pyright: ignore[reportCallIssue,reportArgumentType]
    self.assertRaises(TypeError if systemswig else NotImplementedError,lambda: ca.Function("f", x0,[x0,x1]))  # pyright: ignore[reportCallIssue,reportArgumentType]

  def test_evalfail(self):
    self.message("eval fail test")
    x = ca.SX.sym("x",2,2)
    f = ca.Function("f", [x], [x])
    self.assertRaises(TypeError if systemswig else NotImplementedError,lambda: f.call(x))  # pyright: ignore[reportCallIssue,reportArgumentType]

  def test_SXconversion(self):
    self.message("Conversions from and to SX")
    y=ca.SX.sym("y")
    x=ca.SX.sym("x",3,3)
    ca.SX(y)
    ca.SX(x)
    c.det(x)
    y=array(ca.DM(x))
    c.det(y)

  def test_SXbool(self):
    self.message("bool")

    x = ca.SX.sym("x")
    y = ca.SX.sym("y")

    f = ca.Function("f", [ca.vertcat(*[x,y])],[ca.vertcat(*[ca.logic_and(x,y),ca.logic_or(x,y),ca.logic_not(x)])])


    for t1 in [0,1]:
      for t2 in [0,1]:
        T1 = t1!=0
        T2 = t2!=0
        f_in = [0]*f.n_in()  # type: list

        f_in[0]=[t1,t2]
        f_out = f.call(f_in)
        self.checkarray(f_out[0],ca.DM([T1 and T2,T1 or T2,not T1]),"bool(%d,%d): %s" % (t1,t2,str(f_out[0])))

  def test_SXineq(self):
    self.message("SX ineq")

    x = ca.SX.sym("x")
    y = ca.SX.sym("y")

    f = ca.Function("f", [ca.vertcat(*[x,y])],[ca.vertcat(*[x<y,x<=y,x>=y,x==y,x!=y])])


    for t1 in [-10,0.1,0,1,10]:
      for t2 in [-10,0.1,0,1,10]:
        T1 = t1
        T2 = t2
        f_in = [0]*f.n_in()  # type: list

        f_in[0]=[t1,t2]
        f_out = f.call(f_in)
        self.checkarray(f_out[0],ca.DM([T1 < T2,T1 <= T2, T1 >= T2, T1 == T2, T1 != T2]),"ineq(%d,%d)" % (t1,t2))



  def test_SX_func2(self):
    self.message("SXmatrix typemaps constructors")
    ca.simplify(ca.SX.sym("x"))
    list = [    ("number",2.3, (1,1)),
                ("SX", ca.SX.sym("x"), (1,1))
    ];
    for name, arg,shape in list:
      self.message(":" + name)
      i=c.transpose(c.transpose(arg))
      self.assertEqual(i.shape[0],shape[0],"shape mismatch")
      self.assertEqual(i.shape[1],shape[1],"shape mismatch")
      ca.SX(arg).is_empty()

  def test_SX_func3(self):
    self.message("vector(SXmatrix) typemaps constructors")
    y=ca.SX.sym("y")
    x=ca.SX.sym("x",3,1)
    ca.vertcat(*[x,x])
    ca.vertcat(*[y,y])
    ca.vertcat(*[x,[]])

  def test_eval(self):
    self.message("Function eval")
    x=ca.SX.sym("x",2,2)
    y=ca.SX.sym("y",2,2)
    f  = ca.Function("f", [x,y], [x*y])
    f(x,y)

  def test_symbolcheck(self):
    self.message("Check if non-symbolic inputs are caught")
    self.assertRaises(RuntimeError, lambda : ca.Function("f", [ca.SX(0)],[ca.SX.sym("x")]))

  def test_sparseconstr(self):
    self.message("Check sparsity constructors")
    self.checkarray(ca.DM.ones(ca.Sparsity.lower(3)).full(),np.array([[1,0,0],[1,1,0],[1,1,1]]),"tril")
    self.checkarray(ca.DM.ones(ca.Sparsity.diag(3)).full(),np.array([[1,0,0],[0,1,0],[0,0,1]]),"diag")

  def test_subsassignment(self):
    self.message("Check subscripted assignment")

    import numpy
    numpy.random.seed(42)
    xn = numpy.random.random((3,4))

    x=ca.DM(xn)

    y=ca.DM(7,8)
    z = numpy.zeros((7,8))
    y[0,0]=12; z[0,0] = 12
    self.checkarray(y,z,"scalar assignment")
    z[1:4,[2,4,5,6]]=xn
    y[1:4,[2,4,5,6]]=x
    self.checkarray(y,z,"range assignment")

    kl=[2,4,5,8]
    y.nz[kl]=1.0
    s=y.sparsity()
    for k in kl:
      z[s.row()[k],s.get_col()[k]]=1.0
    self.checkarray(y,z,"nonzero scalar assignment")
    y.nz[kl]=ca.DM(kl)

    cnt=0
    for k in kl:
      z[s.row()[k],s.get_col()[k]]=kl[cnt]
      cnt+=1
    self.checkarray(y,z,"nonzero range assignment")

  @skip(not ca.GlobalOptions.getSimplificationOnTheFly())
  def test_substitute(self):
    self.message("Basic symbolic algebra: substitute")
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")
    z = ca.cos(x)*y
    self.assertTrue(ca.depends_on(z,y))
    self.assertTrue(ca.depends_on(z,x))
    w = ca.substitute(z,x,0)
    self.assertTrue(w.is_symbolic())
    self.assertTrue(ca.depends_on(w,y))
    self.assertFalse(ca.depends_on(w,x))
    self.assertTrue(ca.is_equal(w,y))
    r=w-y
    self.assertFalse(r.is_symbolic())
    self.assertTrue(r.is_zero())
    self.assertEqual(float(r),0)
    self.assertEqual(float(r),0)
    y = ca.SX.sym("y",2)
    y = ca.substitute(y+6,y,0)
    self.assertEqual(int(y[0]),6)
    self.assertEqual(int(y[1]),6)

  def test_primitivefunctions(self):
    self.message("Primitive functions")
    x=ca.SX.sym("x")

    nums = [-2,-1.5,-1,-0.5,-0.25,0,0.25,0.5,1,1.5,2]

    def test(fun,comment,nums,reference):
      self.message(":"+comment)
      f = ca.Function("f", [x],[fun(x)])
      for n,r in zip(nums,reference):
        f_in = [0]*f.n_in()  # type: list

        f_in[0]=n
        f_out = f.call(f_in)
        self.assertEqual(f_out[0][0],r)

    test(ca.sign,"sign",nums,[-1,-1,-1,-1,-1,0,1,1,1,1,1])
    test(ca.heaviside,"heaviside",nums,[0,0,0,0,0,0.5,1,1,1,1,1])
    test(ca.ramp,"ramp",nums,[0,0,0,0,0,0,0.25,0.50,1,1.5,2])
    test(ca.rectangle,"rectangle",nums,[0,0,0,0.5,1,1,1,0.5,0,0,0])
    test(ca.triangle,"triangle",nums,[0,0,0,0.5,0.75,1,0.75,0.5,0,0,0])


  def test_taylor(self):
    self.message("univariate taylor expansion")
    x=ca.SX.sym("x")

    if ca.GlobalOptions.getSimplificationOnTheFly():
      self.assertTrue(ca.is_equal(ca.taylor(ca.sin(x),x),x))

    a_=0.13
    x_=0.15

    a = ca.SX.sym("a")

    def test(e,r):
      f = ca.Function("f", [x,a],[e])
      f_in = [0]*f.n_in()  # type: list

      f_in[0]=x_
      f_in[1]=a_
      f_out = f.call(f_in)
      self.assertAlmostEqual(f_out[0][0],r,10)

    test(ca.taylor(ca.sin(x),x,a,0),ca.sin(a_))
    test(ca.taylor(ca.sin(x),x,a,1),ca.sin(a_)+ca.cos(a_)*(x_-a_))
    test(ca.taylor(ca.sin(x),x,a,2),ca.sin(a_)+ca.cos(a_)*(x_-a_)-(ca.sin(a_)*(x_-a_)**2)/2.0)
    test(ca.taylor(ca.sin(x),x,a,3),ca.sin(a_)+ca.cos(a_)*(x_-a_)-(ca.sin(a_)*(x_-a_)**2)/2.0-(ca.cos(a_)*(x_-a_)**3)/6.0)

    M=ca.blockcat([[a*ca.sin(x),a*ca.cos(x)],[ca.exp(a*x),a*x**2],[ca.cos(x),0]])
    f = ca.Function("f", [x,a],[ca.taylor(M,x)])
    f_in = [0]*f.n_in()  # type: list

    f_in[0]=x_
    f_in[1]=a_
    f_out = f.call(f_in)
    self.checkarray(f_out[0],np.array([[x_*a_,a_],[1+a_*x_,0],[1,0]]),"taylor on dense matrices")

  def test_null(self):
    self.message("Function null")
    x = ca.SX.sym("x")

    f = ca.Function("f", [x],[x**2,[]])
    f_out = f.call([0])
    self.assertTrue(f_out[1].is_empty())

    f = ca.Function("f", [x,[]],[x**2,[]])
    f_out = f.call([0,0])
    self.assertTrue(f_out[1].is_empty())
    f_out = f.call([0,0])

    r = f.call([x,[]])
    self.assertTrue(r[1].is_empty())

    r = f.call([x,[]])
    self.assertTrue(r[1].is_empty())

    r = f.call([x,ca.SX(0,1)])
    self.assertTrue(r[1].is_empty())

    r = f.call([x,ca.SX(1,0)])
    self.assertTrue(r[1].is_empty())

    #self.assertRaises(Exception,lambda : f([x,x]))
    #self.assertRaises(Exception,lambda : f([[],[]]))

  def test_mtaylor(self):
    self.message("multivariate taylor expansions")
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")
    a=ca.SX.sym("a")
    b=ca.SX.sym("b")

    a_=0.13
    x_=0.15

    b_=0.73
    y_=0.75

    def test(e,r):
      f = ca.Function("f", [ca.vertcat(*[x,y]),ca.vertcat(*[a,b])],[e])
      f_in = [0]*f.n_in()  # type: list

      f_in[0]=[x_,y_]
      f_in[1]=[a_,b_]
      f_out = f.call(f_in)
      self.assertAlmostEqual(f_out[0][0],r,10)

    test(ca.mtaylor(ca.sin(x+y),ca.vertcat(*[x,y]),ca.vertcat(*[a,b]),0),ca.sin(a_+b_))
    test(ca.mtaylor(ca.sin(x+y),ca.vertcat(*[x,y]),ca.vertcat(*[a,b]),1),ca.sin(a_+b_)+(ca.cos(b_+a_)*(x_-a_)+ca.cos(b_+a_)*(y_-b_)))

    def sol(x,y,a,b):
      return ca.sin(b+a)+(ca.cos(b+a)*(x-a)+ca.cos(b+a)*(y-b))-(ca.sin(b+a)*(x-a)**2+2*ca.sin(b+a)*(y-b)*(x-a)+ca.sin(b+a)*(y-b)**2)/2

    test(ca.mtaylor(ca.sin(x+y),ca.vertcat(*[x,y]),ca.vertcat(*[a,b]),2),sol(x_,y_,a_,b_))

    def sol(x,y,a,b):
      return ca.sin(b+a)+(ca.cos(b+a)*(x-a)+ca.cos(b+a)*(y-b))-(ca.sin(b+a)*(x-a)**2+2*ca.sin(b+a)*(y-b)*(x-a)+ca.sin(b+a)*(y-b)**2)/2-(ca.cos(b+a)*(x-a)**3+3*ca.cos(b+a)*(y-b)*(x-a)**2+3*ca.cos(b+a)*(y-b)**2*(x-a)+ca.cos(b+a)*(y-b)**3)/6

    test(ca.mtaylor(ca.sin(x+y),ca.vertcat(*[x,y]),ca.vertcat(*[a,b]),3),sol(x_,y_,a_,b_))

    def sol(x,y,a,b):
      return (-2*ca.sin(b+a)*(x-a)*(y-b)-ca.sin(b+a)*(x-a)**2)/2+ca.cos(b+a)*(y-b)-(ca.cos(b+a)*(x-a)**3)/6+ca.cos(b+a)*(x-a)+ca.sin(b+a)

    test(ca.mtaylor(ca.sin(x+y),ca.vertcat(*[x,y]),ca.vertcat(*[a,b]),3,[1,2]),sol(x_,y_,a_,b_))

    test(ca.mtaylor(ca.sin(x+y),ca.vertcat(*[x,y]),ca.vertcat(*[0,0]),4,[1,2]),(-3*x_**2*y_-x_**3)/6+y_+x_)

  def test_issue107(self):
    self.message("Regression test for issue 107: +=")
    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    z=x
    z+=y

    self.assertTrue(x.is_symbolic())
    self.assertFalse(z.is_symbolic())

    x=ca.SX.sym("x")
    y=ca.SX.sym("y")

    z=x
    z+=y

    self.assertTrue(x.is_symbolic())
    self.assertFalse(z.is_symbolic())

  def test_evalchecking(self):
    x = ca.SX.sym("x",1,5)

    y = ca.SX.sym("y",1,3)
    z = ca.SX.sym("z",5,1)
    q = ca.SX.sym("z",1,6)

    f = ca.Function("f", [x],[x**2])

    self.assertRaises(RuntimeError, lambda : f(y))
    self.assertRaises(RuntimeError, lambda : f(q))
    f(z)

  def test_indexinglimits(self):
    self.message("Limits of indexing")
    y = ca.SX.sym("y", 3)
    self.assertRaises(RuntimeError,lambda : y[[0, 5]] )
    try:
      y[[0, 5]] = ca.SX.sym("a")
      self.assertTrue(False)
    except RuntimeError:
      pass
    y[[0, 2]]
    y[[0, 2]] = ca.SX.sym("a")

  def test_issue181(self):
    self.message("Regression test #181")
    x = ca.SX.sym("x")
    #self.assertRaises(TypeError,lambda : SX([x,None]))  # FIXME: this is leaking memory
    self.assertRaises(TypeError if systemswig else NotImplementedError,lambda: ca.Function("f", [[x], [None]], [[2 * x]]))  # pyright: ignore[reportCallIssue,reportArgumentType]

  @known_bug()  # Not implemented
  def test_is_equal(self):
    self.message("equivalent")
    x = ca.SX.sym("x")
    a = x*x
    b = x*x
    self.assertTrue(a.is_equal(b,1))  # pyright: ignore[reportAttributeAccessIssue]

  @skip(not ca.GlobalOptions.getSimplificationOnTheFly())
  def test_SXsimplifications(self):
    self.message("simplifications")
    x = ca.SX.sym("x")

    ops = []
    def temp(x):
      y = 0.5*x
      return y+y

    ops.append(temp)

    def temp(x):
      y = x/2
      return y+y

    ops.append(temp)

    def temp(x):
      y = x*0.5
      return y+y

    ops.append(temp)


    def temp(x):
      y = x*x
      return ((-y)/y)*(-x)

    ops.append(temp)

    ops.append(lambda x: ((-(x*x))/(x*x))*(-x))

    #ops.append(lambda x: ((-x*x)/(x*x))*(-x))

    def temp(x):
      y = x*x
      return (y/(-y))*(-x)

    ops.append(temp)

    def temp(x):
      y = x*x
      return ((-y)/(-y))*(x)

    ops.append(temp)
    ops.append(lambda x: (x-x) + x)
    ops.append(lambda x: ((x*x)-(x*x)) + x)
    ops.append(lambda x: 4*(0.25*x))
    ops.append(lambda x: 4*(x*0.25))
    ops.append(lambda x: 4*(0.25*x))
    ops.append(lambda x: 4*(x*0.25))
    ops.append(lambda x: (0.25*x)*4)
    ops.append(lambda x: (x*0.25)*4)
    ops.append(lambda x: (4*x)/4)
    ops.append(lambda x: 4*(x/4))
    ops.append(lambda x: (x/4)/0.25)
    ops.append(lambda x: x*(((4/x)*x)/4))
    ops.append(lambda x: x*((x*(2/x))/2))
    ops.append(lambda x: x*(((2*x)/x)/2))
    ops.append(lambda x: x*((x/(2*x))*2))
    ops.append(lambda x: x+0)
    ops.append(lambda x: 0+x)
    ops.append(lambda x: x-0)
    ops.append(lambda x: 0-(-x))
    ops.append(lambda x: x*1)
    ops.append(lambda x: 1*x)
    ops.append(lambda x: 1*(x*1))
    ops.append(lambda x: (1*x)*1)
    ops.append(lambda x: (0.5*x)+(0.5*x))
    ops.append(lambda x: (x/2)+(x/2))
    ops.append(lambda x: (x*0.5)+(0.5*x))
    ops.append(lambda x: (ca.SX(4)-ca.SX(4))+x)


    y = ca.SX.sym("x")

    ops.append(lambda x: ((x+y)-(y+x))+x)
    ops.append(lambda x: ((x*y)-(y*x))+x)
    ops.append(lambda x: ((-x)-(-x))+x)

    for op in ops:
      y = op(x)
      f = ca.Function("f", [x],[y])
      f_in = [0]*f.n_in()  # type: list

      f_in[0]=0.3
      f_out = f.call(f_in)
      self.checkarray(f_out[0],array(ca.DM(op(0.3))),"simplifications")
      self.assertEqual(str(y),"x")

      y = op(-x)
      f = ca.Function("f", [x],[y])
      f_in = [0]*f.n_in()  # type: list

      f_in[0]=0.3
      f_out = f.call(f_in)
      self.checkarray(f_out[0],array(ca.DM(op(-0.3))),"simplifications")
      self.assertEqual(str(y),"(-x)")

  def test_truth(self):
    self.message("Truth values")
    self.assertRaises(Exception, lambda : bool(ca.SX.sym("x")))
    self.assertRaises(Exception, lambda : bool(ca.SX.sym("x")>0))
    self.assertTrue(bool(ca.SX(1)))
    self.assertFalse(bool(ca.SX(0)))
    self.assertTrue(bool(ca.SX(0.2)))
    self.assertTrue(bool(ca.SX(-0.2)))
    self.assertRaises(Exception, lambda : bool(ca.SX.sym("x")))
    self.assertRaises(Exception, lambda : bool(ca.SX.sym("x")>0))
    self.assertTrue(bool(ca.SX(ca.SX(1))))
    self.assertFalse(bool(ca.SX(ca.SX(0))))
    self.assertTrue(bool(ca.SX(ca.SX(0.2))))
    self.assertTrue(bool(ca.SX(ca.SX(-0.2))))
    self.assertRaises(Exception, lambda : bool(ca.SX([2.0,3])))

  def test_if_else(self):
    x = ca.SX.sym("x")
    y = ca.if_else(x,1,2)
    f = ca.Function("f", [x],[y])
    f_in = [0]*f.n_in()  # type: list

    f_in[0]=1
    f_out = f.call(f_in)
    self.assertTrue(f_out[0]==1,"if_else")
    f_in = [0]*f.n_in()  # type: list

    f_in[0]=0
    f_out = f.call(f_in)
    self.assertTrue(f_out[0]==2,"if_else")

    x0 = 2.1
    y = ca.if_else(x>1,x**2,x**3)
    f = ca.Function("f", [x],[y])
    f_in = [0]*f.n_in()  # type: list

    f_in[0]=x0
    f_out = f.call(f_in)
    self.checkarray(f_out[0],x0**2,"if_else sens")

    x0 = -2.1
    f_in = [0]*f.n_in()  # type: list

    f_in[0]=x0
    f_out = f.call(f_in)
    self.checkarray(f_out[0],x0**3,"if_else sens")


  def test_is_regular(self):
    x = ca.SX.sym("x")

    self.assertTrue(ca.SX(0).is_regular())
    self.assertFalse(ca.SX(inf).is_regular())
    with self.assertRaises(Exception):
      self.assertTrue(x.nz[0])

    self.assertTrue(ca.SX(ca.DM([0,1])).is_regular())
    self.assertFalse(ca.SX(ca.DM([0,inf])).is_regular())
    self.assertFalse(ca.vertcat(*[x,inf]).is_regular())
    with self.assertRaises(Exception):
      self.assertFalse(ca.vertcat(*[x,x]).is_regular())


  def test_symvar(self):
    a = ca.SX.sym("a")
    b = ca.SX.sym("b")
    c = ca.SX.sym("c")
    e = ca.cos(a*b) + c
    w = ca.symvar(e)
    self.assertEqual(len(w),3)
    if ca.GlobalOptions.getSimplificationOnTheFly():
      self.assertTrue(ca.is_equal(w[0],a))
      self.assertTrue(ca.is_equal(w[1],b))
      self.assertTrue(ca.is_equal(w[2],c))

  def test_poly_coeff(self):
    x =ca.SX.sym("x")
    a= ca.SX.sym("a")
    c=ca.SX.sym("c")
    p=ca.poly_coeff(12*x**4+x**2+a*x+c,x)
    self.assertTrue(ca.is_equal(p[0],12))
    self.assertTrue(ca.is_equal(p[1],0))
    self.assertTrue(ca.is_equal(p[2],1))
    self.assertTrue(ca.is_equal(p[3],a))
    self.assertTrue(ca.is_equal(p[4],c))

    p=ca.poly_coeff((x-a)*(x+a),x)
    self.assertTrue(ca.is_equal(p[0],1))
    self.assertTrue(ca.is_equal(p[1],0))

  def test_poly_roots(self):

    p = ca.SX.sym("[a,b]")
    r = ca.poly_roots(p)

    f = ca.Function("f", [p],[r])
    f_in = [0]*f.n_in()  # type: list
    f_in[0]=ca.DM([2,7])
    a_ = f_in[0][0]
    b_ = f_in[0][1]
    f_out = f.call(f_in)
    f_out[0]
    self.checkarray(f_out[0],ca.vertcat(*[-b_/a_]))

    p = ca.SX.sym("[a,b]")
    r = ca.poly_roots(ca.vertcat(*[p,0]))

    f = ca.Function("f", [p],[r])
    f_in = [0]*f.n_in()  # type: list

    f_in[0]=ca.DM([2,7])
    a_ = f_in[0][0]
    b_ = f_in[0][1]
    f_out = f.call(f_in)
    f_out[0]
    self.checkarray(f_out[0],ca.vertcat(*[-b_/a_,0]))

    p = ca.SX.sym("[a,b,c]")
    r = ca.poly_roots(p)

    f = ca.Function("f", [p],[r])
    f_in = [0]*f.n_in()  # type: list

    f_in[0]=ca.DM([1.13,7,3])
    a_ = f_in[0][0]
    b_ = f_in[0][1]
    c_ = f_in[0][2]
    d = b_**2-4*a_*c_
    f_out = f.call(f_in)
    x0 = (-b_-ca.sqrt(d))/2/a_
    x1 = (-b_+ca.sqrt(d))/2/a_
    f_out[0]
    self.checkarray(f_out[0],ca.vertcat(*[x0,x1]))

    p = ca.SX.sym("[a,b,c,d]")
    r = ca.poly_roots(p)

    f = ca.Function("f", [p],[r])
    f_in = [0]*f.n_in()  # type: list

    f_in[0]=ca.DM([11,1.3,-1.7,0.1])
    f_out = f.call(f_in)
    f_out[0]
    self.checkarray(f_out[0],ca.DM([0.298028,-0.479787,0.0635774]),digits=5)

    p = ca.SX.sym("[a,b,c,d,e]")
    r = ca.poly_roots(p)

    f = ca.Function("f", [p],[r])
    f_in = [0]*f.n_in()  # type: list

    f_in[0]=ca.DM([3,6,-123,  -126,1080])
    f_out = f.call(f_in)
    f_out[0]
    self.checkarray(f_out[0],ca.DM([5,3,-4,-6]),digits=5)

  def test_eig_symbolic(self):
    x = ca.SX.sym("x",2,2)
    f = ca.Function("f", [x],[ca.eig_symbolic(x)])
    f_in = [0]*f.n_in()  # type: list

    f_in[0]=ca.DM([[2,0.1],[0.3,0.7]])
    f_out = f.call(f_in)
    self.checkarray(f_out[0],ca.DM([0.67732,2.02268]),digits=5)


    x = ca.SX.sym("x",2)
    f = ca.Function("f", [x],[ca.eig_symbolic(c.diag(x))])
    f_in = [0]*f.n_in()  # type: list

    f_in[0]=ca.DM([3,7])
    f_out = f.call(f_in)
    self.checkarray(f_out[0],f_in[0])


    x = ca.SX.sym("x",5)
    f = ca.Function("f", [x],[ca.eig_symbolic(c.diag(x))])
    f_in = [0]*f.n_in()  # type: list

    f_in[0]=ca.DM([3,7,2,1,6])
    f_out = f.call(f_in)
    self.checkarray(f_out[0],f_in[0])

    x = ca.SX.sym("x",2,2)
    y = ca.SX.sym("y",2)
    f = ca.Function("f", [x,y],[ca.eig_symbolic(ca.diagcat(*[x,c.diag(y)]))])
    f_in = [0]*f.n_in()  # type: list

    f_in[0]=ca.DM([[2,0.1],[0.3,0.7]])
    f_in[1]=[3,7]
    f_out = f.call(f_in)
    self.checkarray(f_out[0],ca.DM([0.67732,2.02268,3,7]),digits=5)

    x = ca.SX.sym("x",3,3)
    x[2,0] = 0
    x[1,0] = 0

    x = ca.sparsify(x)

    e = ca.eig_symbolic(x)

    f = ca.Function("f", [x],[e])
    f_in = [0]*f.n_in()  # type: list

    f_in[0]=ca.DM(f.sparsity_in(0),list(range(1,8)))
    f_in[0].print_dense()
    f_out = f.call(f_in)
    self.checkarray(f_out[0],ca.DM([1,-0.29150,10.29150]),digits=5)


    x = ca.SX.sym("x",3,3)
    x[2,0] = 0
    x[1,0] = 0
    x[2,1] = 0

    x = ca.sparsify(x)

    e = ca.eig_symbolic(x)

    f = ca.Function("f", [x],[e])
    f_in = [0]*f.n_in()  # type: list

    f_in[0]=ca.DM(f.sparsity_in(0),list(range(1,7)))
    f_in[0].print_dense()
    f_out = f.call(f_in)
    self.checkarray(f_out[0],ca.DM([1,3,6]),digits=5)

    x = ca.SX.sym("x",ca.Sparsity.upper(5))

    f = ca.Function("f", [x],[ca.eig_symbolic(x)])
    fin = ca.DM(x.sparsity(),0)
    fin[ca.Sparsity.diag(5)] = c.diag(list(range(5)))
    self.checkarray(f(fin), ca.DM(list(range(5))))

  def test_jacobian_empty(self):
    x = ca.SX.sym("x",3)

    s = ca.jacobian(ca.DM(0,0),x).shape
    self.assertEqual(s[0],0)
    self.assertEqual(s[1],3)

    s = ca.jacobian(x,ca.SX.sym("x",0,4)).shape
    self.assertEqual(s[0],3)
    self.assertEqual(s[1],0)

  def test_empty_SX(self):
    s = ca.SX([]).shape
    self.assertEqual(s[0],0)
    self.assertEqual(s[1],1)
    x = ca.vertcat(*(ca.SX.sym("x"),ca.SX([])))

    ca.SX(ca.Sparsity(3,3), [])

    ca.SX(ca.Sparsity(3,3), [5])
    with self.assertInException("fully sparse"):
      ca.SX(ca.Sparsity(3,3), [5,4])

  def test_mul_sparsity(self):

    N = 10
    x = ca.SX.sym("x",N,N)
    y = ca.SX.sym("y",N,N)

    x_ = self.randDM(N,N)
    y_ = self.randDM(N,N)

    filt = ca.Sparsity.diag(N)+ca.Sparsity.triplet(N,N,[1],[3])

    f = ca.Function("f", [x,y],[x @ y])
    f_in = [0]*f.n_in()  # type: list

    f_in[0]=x_
    f_in[1]=y_
    g = ca.Function("g", [x,y],[ca.mac(x,y,ca.SX.zeros(filt))])
    g_in = [0]*g.n_in()  # type: list

    g_in[0]=x_
    g_in[1]=y_

    f_out = f.call(f_in)
    g_out = g.call(g_in)

    self.checkarray(ca.DM.ones(filt),ca.DM.ones(g.sparsity_out(0)))

    self.checkarray(f_out[0][filt],g_out[0])

  @skip(platform_arch==32)
  @memory_heavy()
  def test_large_hessian(self):
    A = ca.Sparsity.from_file("../data/apoa1-2.mtx")


    H = ca.DM(A,list(range(A.nnz())))
    H = H + H.T

    H = H[:20000,:20000]

    x = ca.SX.sym("x",H.size1())

    f = ca.Function("f", [x],[ca.mtimes([x.T,H,x])], {'verbose':True})
    H *= 2

    h = hessian_old(f, 0, 0)
    h_out = h.call([0])

    self.assertTrue(h.sparsity_out(0)==H.sparsity())

    self.checkarray(h_out[0].nonzeros(),H.nonzeros())

  def test_mxnulloutput(self):
     a = ca.SX(5,0)
     b = ca.SX.sym("x",2)
     bm = ca.MX.sym("x",2)

     f = ca.Function("f", [b],[a])
     c = f(bm)

     self.assertEqual(c.size1(),5)
     self.assertEqual(c.size2(),0)

     c = f(b)

     self.assertEqual(c.size1(),5)
     self.assertEqual(c.size2(),0)

     a = ca.SX(0,0)

     f = ca.Function("f", [b],[a])

     c = f(bm)

     self.assertEqual(c.size1(),0)
     self.assertEqual(c.size2(),0)

     c = f(b)

     self.assertEqual(c.size1(),0)
     self.assertEqual(c.size2(),0)

  def test_mxnull(self):
     a = ca.SX(5,0)
     b = ca.SX(0,3)

     c = a @ b

     self.assertEqual(c.nnz(),0)

     a = ca.SX(5,3)
     b = ca.SX(3,4)

     c = a @ b

     self.assertEqual(c.nnz(),0)

  def  test_mxnullop(self):
    c = ca.SX(0,0)
    x = ca.SX.sym("x",2,3)

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

  def test_copysign(self):
    x = ca.SX.sym("x")
    y = ca.SX.sym("y")
    z = ca.copysign(x,y)

    f = ca.Function("f", [x,y],[z])

    f_in = [0]*f.n_in()  # type: list

    f_in[0]=2
    f_in[1]=0.5
    f_out = f.call(f_in)
    self.checkarray(f_out[0],ca.DM([2]))

    f_in = [0]*f.n_in()  # type: list

    f_in[0]=2
    f_in[1]=-0.5
    f_out = f.call(f_in)
    self.checkarray(f_out[0],ca.DM([-2]))

    f_in = [0]*f.n_in()  # type: list

    f_in[0]=-2
    f_in[1]=0.5
    f_out = f.call(f_in)
    self.checkarray(f_out[0],ca.DM([2]))

    f_in = [0]*f.n_in()  # type: list

    f_in[0]=-2
    f_in[1]=-0.5
    f_out = f.call(f_in)
    self.checkarray(f_out[0],ca.DM([-2]))

    f_in = [0]*f.n_in()  # type: list

    f_in[0]=2
    f_in[1]=0
    f_out = f.call(f_in)
    self.checkarray(f_out[0],ca.DM([2]))

    J = jacobian_old(f, 0, 0)

    J_in = [0]*J.n_in()  # type: list

    J_in[0]=2
    J_in[1]=0.5
    J_out = J.call(J_in)
    self.checkarray(J_out[0],ca.DM([1]))

    J_in = [0]*J.n_in()  # type: list

    J_in[0]=2
    J_in[1]=-0.5
    J_out = J.call(J_in)
    self.checkarray(J_out[0],ca.DM([-1]))

    J_in = [0]*J.n_in()  # type: list

    J_in[0]=-2
    J_in[1]=0.5
    J_out = J.call(J_in)
    self.checkarray(J_out[0],ca.DM([1]))

    J_in = [0]*J.n_in()  # type: list

    J_in[0]=-2
    J_in[1]=-0.5
    J_out = J.call(J_in)
    self.checkarray(J_out[0],ca.DM([-1]))

    J_in = [0]*J.n_in()  # type: list

    J_in[0]=2
    J_in[1]=0
    J_out = J.call(J_in)
    self.checkarray(J_out[0],ca.DM([1]))

    J = jacobian_old(f, 1, 0)

    J_in = [0]*J.n_in()  # type: list

    J_in[0]=2
    J_in[1]=0.5
    J_out = J.call(J_in)
    self.checkarray(J_out[0],ca.DM([0]))

    J_in = [0]*J.n_in()  # type: list

    J_in[0]=2
    J_in[1]=-0.5
    J_out = J.call(J_in)
    self.checkarray(J_out[0],ca.DM([0]))

    J_in = [0]*J.n_in()  # type: list

    J_in[0]=-2
    J_in[1]=0.5
    J_out = J.call(J_in)
    self.checkarray(J_out[0],ca.DM([0]))

    J_in = [0]*J.n_in()  # type: list

    J_in[0]=-2
    J_in[1]=-0.5
    J_out = J.call(J_in)
    self.checkarray(J_out[0],ca.DM([0]))

    J_in = [0]*J.n_in()  # type: list

    J_in[0]=2
    J_in[1]=0
    J_out = J.call(J_in)
    self.checkarray(J_out[0],ca.DM([0]))

  def test_depends_on(self):
    a = ca.SX.sym("a")
    b = ca.SX.sym("b")

    self.assertTrue(ca.depends_on(a**2,a))
    self.assertTrue(ca.depends_on(a,a))
    self.assertFalse(ca.depends_on(0,a))
    self.assertTrue(ca.depends_on(a**2,ca.vertcat(*[a,b])))
    self.assertTrue(ca.depends_on(a,ca.vertcat(*[a,b])))
    self.assertFalse(ca.depends_on(0,ca.vertcat(*[a,b])))
    self.assertTrue(ca.depends_on(b**2,ca.vertcat(*[a,b])))
    self.assertTrue(ca.depends_on(b,ca.vertcat(*[a,b])))
    self.assertTrue(ca.depends_on(a**2+b**2,ca.vertcat(*[a,b])))
    self.assertTrue(ca.depends_on(a+b,ca.vertcat(*[a,b])))
    self.assertTrue(ca.depends_on(ca.vertcat(*[0,a]),a))
    self.assertTrue(ca.depends_on(ca.vertcat(*[a,0]),a))
    self.assertTrue(ca.depends_on(ca.vertcat(*[a**2,b**2]),ca.vertcat(*[a,b])))
    self.assertTrue(ca.depends_on(ca.vertcat(*[a,0]),ca.vertcat(*[a,b])))
    self.assertTrue(ca.depends_on(ca.vertcat(*[0,b]),ca.vertcat(*[a,b])))
    self.assertTrue(ca.depends_on(ca.vertcat(*[b,0]),ca.vertcat(*[a,b])))
    self.assertFalse(ca.depends_on(ca.vertcat(*[0,0]),ca.vertcat(*[a,b])))

  @requires("is_smooth")
  def test_is_smooth(self):
    x = ca.SX.sym("a",2,2)
    import warnings
    with warnings.catch_warnings():
      warnings.simplefilter("error",DeprecationWarning)
      with self.assertRaises(Exception):
        is_smooth(x)  # pyright: ignore[reportUndefinedVariable]
      warnings.simplefilter("ignore")
      is_smooth(x)  # pyright: ignore[reportUndefinedVariable]

  def test_which_depends(self):
    for X in [ca.SX,ca.MX]:
      x = X.sym("x")
      y = X.sym("y")

      p = X.sym("p")

      e = ca.vertcat(0,x,y,p,2*p**3,x*y,x*p,ca.sin(x),ca.cos(y),ca.sqrt(x+y),p*p*x,x*y*p)

      self.checkarray(ca.which_depends(e, ca.vertcat(x,y),2,True),[0, 0, 0, 0,0, 1, 0, 1, 1, 1, 0, 1])
      self.checkarray(ca.which_depends(e, ca.vertcat(x,y),1,True),[0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1])

      z =X.sym("z")
      e = ca.vertcat(x*p,x+y)
      self.checkarray(ca.which_depends(e, ca.vertcat(x,y,p,z),2,False),[True, False, True, False])
      self.checkarray(ca.which_depends(e, ca.vertcat(x,y,p,z),1,False),[True, True, True, False])

      e = ca.vertcat(x*p,x+z*y)
      self.checkarray(ca.which_depends(e, ca.vertcat(x,y,p),2,False),[True, False, True])
      self.checkarray(ca.which_depends(e, ca.vertcat(x,y,p),1,False),[True, True, True])

      e = ca.vertcat(x*p,x+z*y)
      self.checkarray(ca.which_depends(e, ca.vertcat(x,y,p,z),2,False),[True, True, True, True])
      self.checkarray(ca.which_depends(e, ca.vertcat(x,y,p,z),1,False),[True, True, True, True])

      e = ca.vertcat(ca.sin(x+y)+p)
      self.checkarray(ca.which_depends(e, ca.vertcat(x,y,p,z),2,False),[True, True, False, False])
      self.checkarray(ca.which_depends(e, ca.vertcat(x,y,p,z),1,False),[True, True, True, False])

      e = ca.vertcat(ca.sin(x)*p**2,y**2)
      #self.checkarray(which_depends(e, vertcat(x,y,p),3,True),[True, False])
      #self.checkarray(which_depends(e, vertcat(x,y,p),3,False),[True, False, True])
      self.checkarray(ca.which_depends(e, ca.vertcat(x,y,p),2,True),[True, True])
      self.checkarray(ca.which_depends(e, ca.vertcat(x,y,p),2,False),[True, True, True])

      e = ca.vertcat(x**2*p,y)
      #self.checkarray(which_depends(e, vertcat(x,y,p),3,True),[True, False])
      #self.checkarray(which_depends(e, vertcat(x,y,p),3,False),[True, False, False])

      self.checkarray(ca.which_depends(e, ca.vertcat(x,y,p),2,True),[True, False])
      self.checkarray(ca.which_depends(e, ca.vertcat(x,y,p),2,False),[True, False, True])

  def test_if_else_zero_sens(self):

    for X in [ca.SX]:
      x=X.sym('x')


      a = 1+3*x+ca.sqrt(3*x)*x+7*x
      b = 1+2*x+ca.sin(2*x)*x +x
      z = ca.if_else(x>0,a,b)*x

      f = ca.Function("f",[x],[z,ca.jacobian(z,x)])
      fa = ca.Function("f",[x],[a*x,ca.jacobian(a*x,x)])
      fb = ca.Function("f",[x],[b*x,ca.jacobian(b*x,x)])

      for i,j in zip(f([3]),fa([3])):
        self.checkarray(i,j)

      for i,j in zip(f([-3]),fb([-3])):
        self.checkarray(i,j)


      f = ca.Function("f",[x],[z])
      fa = ca.Function("f",[x],[a*x])
      fb = ca.Function("f",[x],[b*x])

      self.checkfunction(f,fa,inputs=[3])
      self.checkfunction(f,fb,inputs=[-3],evals=1)

  def test_pw_const(self):
      t= ca.SX.sym("t")

      e = ca.pw_const(t, [0,2,3],[7,1,3,5])

      E = ca.Function("E",[t],[e])

      self.checkarray(E(-2),7)
      self.checkarray(E(-1),7)
      self.checkarray(E(0),1)
      self.checkarray(E(1),1)
      self.checkarray(E(1.9999),1)
      self.checkarray(E(2),3)
      self.checkarray(E(2.5),3)
      self.checkarray(E(3),5)
      self.checkarray(E(10),5)

  def test_pw_lin(self):
      t= ca.SX.sym("t")

      e = ca.pw_lin(t, [0,2,3,5], [7,1,3,2])

      E = ca.Function("E",[t],[e])

      self.checkarray(E(-2),13)
      self.checkarray(E(-1),10)
      self.checkarray(E(0),7)
      self.checkarray(E(1),4)
      self.checkarray(E(2),1)
      self.checkarray(E(2.5),2)
      self.checkarray(E(3),3)
      self.checkarray(E(4),2.5)
      self.checkarray(E(5),2)
      self.checkarray(E(7),1)

  def test_numpy_error(self):
      x = ca.SX.sym("x",3)
      # numpy.foo on symbolic values is gated behind the casadi-aware numpy
      # support (issue #2959).  Enabled, np.linalg.norm bridges to
      # casadi.norm_2 and returns a casadi.ArrayInterface (to_casadi() for
      # the underlying value).
      ca.GlobalOptions.setNumpyMode(1)
      try:
        f = ca.Function("f", [x], [np.linalg.norm(x).to_casadi()])  # pyright: ignore[reportArgumentType,reportCallIssue,reportAttributeAccessIssue]
        self.checkarray(f([3.0, 4.0, 0.0]), ca.DM(5))
        # A function with no casadi equivalent: dispatch returns
        # NotImplemented and numpy raises naming the function.
        with self.assertInException("argmax"):
          np.argmax(x)
      finally:
        ca.GlobalOptions.setNumpyMode(0)


  def test_quadratic(self):
    for X in [ca.SX,ca.MX]:
      x = X.sym("x")
      p = X.sym("p")
      y = X.sym("y")

      self.assertFalse(ca.is_quadratic(ca.sin(x),x))
      self.assertFalse(ca.is_quadratic(x**3,x))
      self.assertTrue(ca.is_quadratic(x**2,x))
      self.assertTrue(ca.is_quadratic(4*x,x))
      self.assertTrue(ca.is_quadratic(5,x))

      self.assertFalse(ca.is_quadratic(ca.sin(x)*p**4,x))
      self.assertFalse(ca.is_quadratic(x**3*p**4,x))
      self.assertTrue(ca.is_quadratic(x**2*p**4,x))
      self.assertTrue(ca.is_quadratic(x*p**4,x))
      self.assertTrue(ca.is_quadratic(5*p**4,x))

      self.assertFalse(ca.is_linear(ca.sin(x),x))
      self.assertFalse(ca.is_linear(x**3,x))
      self.assertFalse(ca.is_linear(x**2,x))
      self.assertTrue(ca.is_linear(3*x,x))
      self.assertTrue(ca.is_linear(5,x))

      self.assertFalse(ca.is_linear(ca.sin(x)*p**4,x))
      self.assertFalse(ca.is_linear(x**3*p**4,x))
      self.assertFalse(ca.is_linear(x**2*p**4,x))
      self.assertTrue(ca.is_linear(x*p**4,x))
      self.assertTrue(ca.is_linear(5*p**4,x))



      z = x**2+3*y**2 + 0.5*x*y + 7*x + 6*y+7
      [A,b,c] = ca.quadratic_coeff(z,ca.vertcat(x,y))

      with self.assertInException("non-quadratic"):
        [A,b,c] = ca.quadratic_coeff(x**2+3*y**2 + 0.5*x*y + 7*x + 6*y+7+ca.sin(x),ca.vertcat(x,y))

      with self.assertInException("scalar"):
        [A,b,c] = ca.quadratic_coeff(ca.vertcat(x,y),x)

      z = x**2+3*y**2 + 0.5*x*y -p*y + 7*x + 6*y+7
      [A,b,c] = ca.quadratic_coeff(z,ca.vertcat(x,y))

      xy = ca.vertcat(x,y)

      e = 0.5*ca.bilin(A,xy,xy)+ca.dot(b,xy)+c

      f = ca.Function('f',[xy,p],[z])
      f2 = ca.Function('f',[xy,p],[e])
      self.checkfunction(f,f2,inputs=[1.1,1.3])


      with self.assertInException("non-linear"):
        [A,b] = ca.linear_coeff(x**2+3*y**2 + 0.5*x*y + 7*x + 6*y+7,ca.vertcat(x,y))

      with self.assertInException("vector"):
        [A,b] = ca.linear_coeff(ca.blockcat([[x,y],[y,x]]),x)

      z = ca.vertcat(7*x + 6*y+7 ,5 -p*y )
      [A,b] = ca.linear_coeff(z,xy)

      e = A @ xy+b

      f = ca.Function('f',[xy,p],[z])
      f2 = ca.Function('f',[xy,p],[e])
      self.checkfunction(f,f2,inputs=[1.1,1.3])

  def test_evalf(self):
    x = ca.SX.sym("x")

    y = ca.SX(5)

    self.checkarray(ca.evalf(y),5)
    with self.assertInException("since variables [x] are free"):
      ca.evalf(x)
  

  def test_output_sx(self):
  
    x = ca.MX.sym("x")
    y = ca.MX.sym("y")
    z = ca.MX.sym("z")
    w = ca.sqrt(z) @ ca.sin(x*y)

    f = ca.Function('f',[x,y,z],[x*y*z,w])
    print(f)

    x = ca.SX.sym("x")
    y = ca.SX.sym("y")
    z = ca.SX.sym("z")
    
    args = [x/y,10*z,y*z]

    v = f.call(args,False,True)

    output_0 = v[0]
    output_1 = v[1]    


    self.assertEqual(output_0.dep(0).element_hash(),output_1.dep(0).element_hash())
    call_node = output_0.dep(0)
    
    print("here")
    
    v = None
    
    # Check that get_output is cached
    output_1b = call_node.get_output(1)
    self.assertEqual(output_1.element_hash(),output_1b.element_hash())
    
    h1 = output_1.element_hash()
    
    output_1 = None
    output_1b = None
    
    # Without output_0 gone, output_1b seems to get consistently reconstructed into the freshly deleted memory slot
    output_0 = None

    output_1b = call_node.get_output(1)
    
    # This test is too fragile
    #self.assertNotEqual(h1,output_1b.element_hash())


  @memory_heavy()
  def test_call_fun(self):

    A = ca.sparsify(ca.DM([[1,0,1],[0,0,6],[0,8,9]]))

    x = ca.MX.sym("x",3)
    y = ca.MX.sym("y")
    z = ca.MX.sym("z")
    w = A*ca.sqrt(z) @ ca.sin(x*y)

    f = ca.Function('f',[x,y,z],[(x*y*z)[:2],w])
    print(f)

    x = ca.SX.sym("x",3)
    y = ca.SX.sym("y")
    z = ca.SX.sym("z")
    
    args = [x/y,10*z,y*z]

    v = f(*args)
    v2 = f.call(args,False,True)
    print(v2)
    print(ca.vcat(v2))
    

    self.assertTrue(v2[0][0].is_output())
    self.assertFalse(y.is_output())
    k = 0
    for i in range(f.n_out()):
        for j in range(f.nnz_out(i)):
            self.assertEqual(v2[i][j].which_output(),k)
            k+=1
    callnode = v2[0][0].dep(0)
    self.assertTrue(callnode.is_call())
    self.assertFalse(y.is_call())
    self.assertTrue(callnode.has_output())
    self.assertFalse(y.has_output())
    self.assertTrue(callnode.which_function().__hash__()==f.__hash__())
    
    self.assertTrue(callnode.n_dep()==5)
    callnode.dep(4)
    
    e = callnode.get_output(3)
    print(e,hash(e))
    
    ee = callnode.get_output(3)
    print(ee,hash(e))
    
    
    self.assertTrue("{3}" in str(callnode.get_output(3)))
    
    print(ca.cos(v2[1] @ v2[0].T/y).shape)
    F2 = ca.Function('F',[x,y,z],[ca.sin(v2[0]*y-args[-1]),ca.cos(v2[1] @ v2[0].T/y)])
    
    with self.assertOutput("f:(i0[3],i1,i2)->(o0[2],o1[3])",[]):
        print(f)
    with self.assertOutput("[[@2,@3],[@4,@5,@7]] = f([@2,@3,@4],@5,@6);",[]):
        F2.disp(True)

    F1 = ca.Function('F',[x,y,z],[ca.sin(v[0]*y-args[-1]),ca.cos(v[1] @ v[0].T/y)])    
    for ad_weight in [True,False]:
        for ad_weight_sp in [True,False]:
            F2 = ca.Function('F',[x,y,z],[ca.sin(v2[0]*y-args[-1]),ca.cos(v2[1] @ v2[0].T/y)],{"ad_weight_sp":ad_weight_sp,"ad_weight":ad_weight})
                    
            self.checkfunction(F1,F2,inputs=[[0.1,1.7,2.3],1.13,0.11])

            self.check_serialize(F2,inputs=[[0.1,1.7,2.3],1.13,0.11])
            for avoid_stack in [True,False]:
                self.check_codegen(F2,inputs=[[0.1,1.7,2.3],1.13,0.11],opts={"avoid_stack":avoid_stack})
    
    instr = F2.instructions_sx()
    check = False
    for k in range(F2.n_instructions()):
        op = F2.instruction_id(k)
        if op==ca.OP_CALL:  # pyright: ignore[reportUndefinedVariable]
            check = True
            self.assertTrue(instr[k].is_call())
            self.assertEqual(F2.instruction_output(k),[2, 3, 4, 5, 7])
            self.assertEqual(F2.instruction_input(k),[2, 3, 4, 5, 6])

    self.assertTrue(check)

    F2 = ca.Function('F',[x,y,z],[ca.sin(v2[0]*y-args[-1]),ca.cos(v2[1] @ v2[0].T/y)],{"cse":True})
    self.checkfunction(F1,F2,inputs=[[0.1,1.7,2.3],1.13,0.11])
    
    res = F2.find_functions()
    self.assertTrue(len(res)==1)
    self.assertTrue(res[0].__hash__()==f.__hash__())

    F2 = ca.Function('F',[x,y,z],list(F2(x,y,z)))

    self.checkfunction(F1,F2,inputs=[[0.1,1.7,2.3],1.13,0.11])

    F1 = ca.Function('F',[x,y,z],list(F1(2*x,3*y,4*z)))
    F2 = ca.Function('F',[x,y,z],list(F2(2*x,3*y,4*z)))

    self.checkfunction(F1,F2,inputs=[[0.1,1.7,2.3],1.13,0.11])

    F1 = ca.Function('F',[x,y,z],list(F1(x,3*y,4*z)))
    F2 = ca.Function('F',[x,y,z],list(F2(x,3*y,4*z)))

    self.checkfunction(F1,F2,inputs=[[0.1,1.7,2.3],1.13,0.11])
    
    
    F2 = ca.Function('F',[x,y,z],[ca.sin(v2[0]*y-args[-1]),ca.cos((v2[1] @ v2[0].T)[1]/y)])
    with self.assertOutput("[[@2,@3],[NULL,@4,NULL]] = f([@2,@3,@4],@5,@6);",[]):
        F2.disp(True)

    F1 = ca.Function('F',[x,y,z],[ca.sin(v[0]*y-args[-1]),ca.cos((v[1] @ v[0].T)[1]/y)])    
    for ad_weight in [True,False]:
        for ad_weight_sp in [True,False]:
            F2 = ca.Function('F',[x,y,z],[ca.sin(v2[0]*y-args[-1]),ca.cos((v2[1] @ v2[0].T)[1]/y)],{"ad_weight_sp":ad_weight_sp,"ad_weight":ad_weight})
                    
            self.checkfunction(F1,F2,inputs=[[0.1,1.7,2.3],1.13,0.11])

            self.check_serialize(F2,inputs=[[0.1,1.7,2.3],1.13,0.11])
            for avoid_stack in [True,False]:
                self.check_codegen(F2,inputs=[[0.1,1.7,2.3],1.13,0.11],opts={"avoid_stack":avoid_stack})
    
    # multiple instances in one graph
 
  def test_call_fun2(self):
  
    print("test_call_fun2")

    for X in [ca.SX,ca.MX]:
    
        for never_inline in [True, False]:
            print(X, never_inline)
            A = ca.sparsify(ca.DM([[1,0,1],[0,0,6],[0,8,9]]))

            x = X.sym("x",3)
            y = X.sym("y")
            z = X.sym("z")
            w = A*ca.sqrt(z) @ ca.sin(x*y)
            
            f = ca.Function('fun',[x,y,z],[(x*y*z)[:2],w],{"never_inline":never_inline})

            x = ca.SX.sym("x",3)
            y = ca.SX.sym("y")
            z = ca.SX.sym("z")
            
            args = [x/y,10*z,y*z]
            res = f(*args)
            print(res)
            self.assertEqual('fun' in str(res), never_inline)
        
  def test_call_fun_nominal_out_deriv(self):
    x = ca.SX.sym("x",2)
    p = ca.SX.sym("p",2)
    rf = ca.rootfinder('rf',"newton",{'x':x,"g":ca.vertcat(ca.sin(x[0])-p[0],ca.sin(x[0]+x[1])-p[1]*p[0]),"p":p})
    
    f = ca.Function("f",[p],[ca.exp(rf(p=p**2)["x"])])
    P = ca.MX.sym("p",2)
    fref = ca.Function("f",[P],[ca.exp(rf(p=P**2)["x"])])
    
    self.checkfunction(f,fref,inputs=[[0.1,0.2]])
  
  @memory_heavy()
  def test_call_fun_copy_elision(self):
    old = ca.GlobalOptions.getCopyElisionMinSize()
    ca.GlobalOptions.setCopyElisionMinSize(0)
    
    mX = ca.MX.sym("A",3,6)
    mY = ca.MX.sym("Y")
    
    X = ca.SX.sym("A",3,6)
    Y = ca.SX.sym("Y")
    
    ca.DM.rng(1)
    
    g = ca.Function('g',[X,Y],[X,Y])
    gref = ca.Function('g',[mX,mY],[mX,mY])
    inputs = [ca.DM.rand(3,3),ca.DM.rand(1)]
    self.checkfunction(g,gref,inputs=inputs)
    for avoid_stack in [True,False]:
        self.check_codegen(g,inputs=inputs,opts={"avoid_stack":avoid_stack})
    g_roundtrip = ca.Function.deserialize(g.serialize())
    self.checkfunction(g,g_roundtrip,inputs=inputs)
    for avoid_stack in [True,False]:
        self.check_codegen(g_roundtrip,inputs=inputs,opts={"avoid_stack":avoid_stack})

    X = ca.MX.sym("A",3,3)
    Y = ca.MX.sym("Y")

    f = ca.Function('f',[X,Y],[ca.sumsqr(X)*Y],{"never_inline":True})
    fref = ca.Function('f',[X,Y],[ca.sumsqr(X)*Y])
    X = ca.SX.sym("A",3,6)
    Y = ca.SX.sym("Y")
    
    ca.DM.rng(1)
    
    g = ca.Function('g',[X,Y],[f(X[:,3:],ca.sin(Y)),X[1,3],3*X[2,4]])
    gref = ca.Function('g',[X,Y],[fref(X[:,3:],ca.sin(Y)),X[1,3],3*X[2,4]])
    inputs = [ca.DM.rand(3,3),ca.DM.rand(1)]
    self.checkfunction(g,gref,inputs=inputs)
    for avoid_stack in [True,False]:
        self.check_codegen(g,inputs=inputs,opts={"avoid_stack":avoid_stack})
    g_roundtrip = ca.Function.deserialize(g.serialize())
    self.checkfunction(g,g_roundtrip,inputs=inputs)
    for avoid_stack in [True,False]:
        self.check_codegen(g_roundtrip,inputs=inputs,opts={"avoid_stack":avoid_stack})

    X = ca.MX.sym("A",4)
    Y = ca.MX.sym("Y")

    f = ca.Function('f',[X,Y],[ca.sumsqr(X)*Y],{"never_inline":True})
    fref = ca.Function('f',[X,Y],[ca.sumsqr(X)*Y])
    X = ca.SX.sym("A",6)
    Y = ca.SX.sym("Y")
    
    ca.DM.rng(1)
    
    g = ca.Function('g',[X,Y],[f(ca.vertcat(2,X[:3]),ca.sin(Y))])
    gref = ca.Function('g',[X,Y],[fref(ca.vertcat(2,X[:3]),ca.sin(Y))])
    inputs = [ca.DM.rand(6),ca.DM.rand(1)]
    self.checkfunction(g,gref,inputs=inputs)
    for avoid_stack in [True,False]:
        self.check_codegen(g,inputs=inputs,opts={"avoid_stack":avoid_stack})
    g_roundtrip = ca.Function.deserialize(g.serialize())
    self.checkfunction(g,g_roundtrip,inputs=inputs)
    for avoid_stack in [True,False]:
        self.check_codegen(g_roundtrip,inputs=inputs,opts={"avoid_stack":avoid_stack})
        
    g = ca.Function('g',[X,Y],[f(ca.vertcat(2,X[:3]),ca.sin(Y)),3*X])
    gref = ca.Function('g',[X,Y],[fref(ca.vertcat(2,X[:3]),ca.sin(Y)),3*X])
    inputs = [ca.DM.rand(6),ca.DM.rand(1)]
    self.checkfunction(g,gref,inputs=inputs)
    for avoid_stack in [True,False]:
        self.check_codegen(g,inputs=inputs,opts={"avoid_stack":avoid_stack})
    g_roundtrip = ca.Function.deserialize(g.serialize())
    self.checkfunction(g,g_roundtrip,inputs=inputs)
    for avoid_stack in [True,False]:
        self.check_codegen(g_roundtrip,inputs=inputs,opts={"avoid_stack":avoid_stack})
        
    ca.GlobalOptions.setCopyElisionMinSize(old)

  def test_ufunc(self):
    # numpy.foo on a symbolic value is gated behind the casadi-aware numpy
    # support (issue #2959); enabled, it returns an casadi.ArrayInterface.
    ca.GlobalOptions.setNumpyMode(1)
    try:
      y = np.sin(ca.SX.sym('x'))
    finally:
      ca.GlobalOptions.setNumpyMode(0)

  def test_mmin(self):
      x = ca.SX.sym("X",2)
      f0 = ca.Function("f",[x],[(x[0]+x[1])/2])
      f1 = ca.Function("f",[x],[ca.mmin(x)])
      f2 = ca.Function("f",[x],[ca.fmin(x[0],x[1])])
      self.checkfunction(f1,f2,inputs=[[0.2,0.3]])
      self.checkfunction(f1,f2,inputs=[[2,2]])
      self.checkfunction(f1,f0,inputs=[[2,2]])
      f1 = ca.Function("f",[x],[ca.mmax(x)])
      f2 = ca.Function("f",[x],[ca.fmax(x[0],x[1])])
      self.checkfunction(f1,f2,inputs=[[0.2,0.3]])
      self.checkfunction(f1,f2,inputs=[[2,2]])
      self.checkfunction(f1,f0,inputs=[[2,2]])

      x = ca.SX.sym("X")
      f0 = ca.Function("f",[x],[x])
      f1 = ca.Function("f",[x],[ca.fmin(x,x)])
      self.checkfunction(f1,f0,inputs=[[1]])
      f1 = ca.Function("f",[x],[ca.fmax(x,x)])
      self.checkfunction(f1,f0,inputs=[[1]])

  def test_empty_broadcast(self):
    for nc in [0,2]:
      res = ca.atan2(ca.SX.sym("c",nc,1),ca.SX.sym("t",nc,3))

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

  def test_logsumexp(self):
    x = ca.SX.sym("x",3)

    f_ref = ca.Function("f_ref",[x],[ca.log(ca.exp(x[0])+ca.exp(x[1])+ca.exp(x[2]))])
    f = ca.Function("f",[x],[ca.logsumexp(x)])

    self.checkfunction(f,f_ref,inputs=[ca.vertcat(1.1,1.3,1.7)])
    self.check_codegen(f,inputs=[ca.vertcat(1.1,1.3,1.7)])
    self.checkfunction(f,f_ref,inputs=[ca.vertcat(1.1,1.3,1.3)])
    self.checkfunction(f,f_ref,inputs=[ca.vertcat(1.3,1.3,1.3)])

    self.checkarray(ca.logsumexp(ca.vertcat(1.3,1.3,1.3)),f_ref(ca.vertcat(1.3,1.3,1.3)))
    self.checkarray(ca.logsumexp(ca.vertcat(1.1,1.3,1.3)),f_ref(ca.vertcat(1.1,1.3,1.3)))
    self.checkarray(ca.logsumexp(ca.vertcat(1.1,1.3,1.7)),f_ref(ca.vertcat(1.1,1.3,1.7)))

    # Avoid overflow
    res = f(ca.vertcat(100,1000,10000))
    self.checkarray(res,10000)

    self.checkarray(ca.logsumexp(ca.vertcat(100,1000,10000)),f(ca.vertcat(100,1000,10000)))

  def test_extract_parametric_call_sx(self):
    x = ca.MX.sym("x")
    y = ca.MX.sym("y")
    z = ca.MX.sym("z")
    
    f = ca.Function("f",[x,y,z],[x*y,y*z,x*z],{"never_inline":True})
    
    x = ca.SX.sym("x")
    p = ca.SX.sym("p")
    
    [expr1,expr2,expr3] = f(1,ca.sin(p),ca.sqrt(p))
    
    expr = x*(expr1+expr3)
    
    expr_ret,symbols,parametric = ca.extract_parametric(expr,p)
    
    print("expr_ret",expr_ret)
    print("symbols",symbols)
    print("parametric",parametric)
    
    self.assertTrue("x*e_0" in str(expr_ret))
    
    print(expr_ret,symbols,parametric)
    
    self.assertFalse(ca.depends_on(expr_ret,p))
        
    expr_recreated = ca.substitute([expr_ret],symbols,parametric)[0]
    
    print(expr_recreated-expr)
    
    self.assertTrue(ca.cse(expr_recreated-expr).is_zero())


    [expr1,expr2,expr3] = f(1,ca.sin(p),x**2)
    
    expr = expr1+expr3
    

    
    expr_ret,symbols,parametric = ca.extract_parametric(expr,p)
    
    self.assertFalse(ca.depends_on(expr_ret,p))
        
    expr_recreated = ca.substitute([expr_ret],symbols,parametric)[0]
    
    f = ca.Function('f',[x,p],[expr_recreated])
    f.generate('f1.c')
    f = ca.Function('f',[x,p],[expr])
    f.generate('f2.c')
    ca.cse(expr_recreated)
    self.assertTrue(ca.cse(expr_recreated-expr).is_zero())
    
  def test_check_recursion(self):
    x = ca.SX.sym("x")


    f = ca.Function("f",[x],[x**2],{"always_inline":True})

    y = ca.MX.sym("y")
    
    f(y)

  def test_sx_eval_mx(self):
    n = 2
    y = ca.MX.sym("y",ca.Sparsity.lower(n))

    g = ca.Function('g',[y],[ca.sin(y)],{"never_inline":True})

    x = ca.SX.sym("x",ca.Sparsity.lower(n))

    f = ca.Function("f",[x],[g(x**2)-2])


    X = ca.MX.sym("X",ca.Sparsity.upper(n))

    F1 = ca.Function("F1",[X],f.call([X],True))

    ca.DM.rng(1)
    A = ca.DM.rand(n,n)

    F2 = ca.Function("F2",[X],f.call([X]))

    self.checkfunction_light(F1,F2,inputs=[A])

  def test_contains(self):
    x = ca.SX.sym("x")
    y = ca.SX.sym("y")
    z = ca.SX.sym("z")
    
    e = y*z
    
    self.assertTrue(ca.contains([x,y,z],x))
    self.assertFalse(ca.contains([x,y],z))
    self.assertTrue(ca.contains_any([x,y],[y,z]))
    self.assertFalse(ca.contains_all([x,y],[y,z]))
    self.assertTrue(ca.contains_any([x,y],[x,y]))
    self.assertTrue(ca.contains_all([x,y],[x,y]))
    self.assertTrue(ca.contains([e,x],e))
    
    with self.assertInException("Can only convert 1-by-1 matrices to scalars"):
        ca.contains([ca.vertcat(x,y)],x)
    
  def test_pow(self):
    x = ca.SX.sym("x")
    y = ca.SX(4,1)
    self.assertEqual((y**0).nnz(),4)
    y[1] = x
    self.assertEqual((y**0).nnz(),4)
    y[1] = 0
    self.assertEqual((y**0).nnz(),4)
    y[1] = 1
    self.assertEqual((y**0).nnz(),4)

  def test_set_precision(self):
    # Issue #4326: SX::set_precision was being ignored for SX printing.
    # DM and SX have their own precision settings; MX numeric printing
    # piggybacks on DM's setting (the value lives in a DM-like node).
    pi_val = 3.141592653589793239
    expected_default = "3.14159"
    expected_prec12 = "3.14159265359"

    for cls in (ca.DM, ca.SX):
      x = cls(pi_val)
      try:
        self.assertEqual(str(x), expected_default)
        cls.set_precision(12)
        self.assertEqual(str(x), expected_prec12)
        cls.set_precision(4)
        cls.set_scientific(True)
        self.assertTrue("e" in str(x).lower())
      finally:
        cls.set_precision(6)
        cls.set_scientific(False)

    # MX scalar constants follow DM precision
    mx = ca.MX(pi_val)
    try:
      self.assertEqual(str(mx), expected_default)
      ca.DM.set_precision(12)
      self.assertEqual(str(mx), expected_prec12)
    finally:
      ca.DM.set_precision(6)

    # set_precision applies to vectors/matrices too
    for cls in (ca.DM, ca.SX):
      v = cls([1.234567890123, 9.876543210987])
      try:
        self.assertEqual(str(v), "[1.23457, 9.87654]")
        cls.set_precision(12)
        self.assertEqual(str(v), "[1.23456789012, 9.87654321099]")
      finally:
        cls.set_precision(6)

  def test_linearize(self):
    x = ca.SX.sym("x")
    y = ca.sin(x)
    x0 = 0.2
    print(ca.linearize(ca.sin(x),x,x0))
    print(ca.sin(x0)+ca.cos(x0)*(x-x0))
    F = ca.Function('F',[x],[y])
    Ff = F.forward(1)
    n = F(x0)
    print(n + Ff(x0,n,x-x0))
    print(ca.taylor(y,x,x0))

if __name__ == '__main__':
    unittest.main()
