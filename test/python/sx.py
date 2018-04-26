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

import casadi as c
import numpy
from numpy import random, array, linalg, matrix, zeros, ones
import unittest
from types import *
from helpers import *
from casadi import *

scipy_available = True
try:
	from scipy.sparse import csr_matrix
except:
	scipy_available = False

class SXtests(casadiTestCase):

  def setUp(self):
    self.pool=FunctionPool()
    self.pool.append(lambda x: sqrt(x[0]),sqrt,"sqrt")
    self.pool.append(lambda x: sin(x[0]),sin,"sin")
    self.pool.append(lambda x: cos(x[0]),cos,"cos")
    self.pool.append(lambda x: tan(x[0]),tan,"tan")
    self.pool.append(lambda x: fabs(x[0]),fabs,"fabs")
    self.pool.append(lambda x: sign(x[0]),sign,"sign")
    self.pool.append(lambda x: arctan(x[0]),arctan,"arctan")
    self.pool.append(lambda x: arcsin(x[0]),arcsin,"arcsin")
    self.pool.append(lambda x: arccos(x[0]),arccos,"arccos")
    self.pool.append(lambda x: exp(x[0]),exp,"exp")
    self.pool.append(lambda x: log(x[0]),log,"log")
    self.pool.append(lambda x: x[0]**0,lambda x : x**0,"x^0",flags={'nozero'})
    self.pool.append(lambda x: x[0]**1,lambda x : x**1,"^1")
    self.pool.append(lambda x: x[0]**(-2),lambda x : x**(-2),"^-2",flags={'nozero'})
    self.pool.append(lambda x: x[0]**(0.3),lambda x : x**(0.3),"^0.3")
    self.pool.append(lambda x: floor(x[0]),floor,"floor")
    self.pool.append(lambda x: ceil(x[0]),ceil,"ceil")
    self.Jpool=FunctionPool()
    self.Jpool.append(lambda x: sqrt(x[0]),lambda x:diag(1/(2.0*sqrt(x))),"sqrt")
    self.Jpool.append(lambda x: sin(x[0]),lambda x:diag(cos(x)),"sin")
    self.Jpool.append(lambda x: fabs(x[0]),lambda x:diag(sign(x)),"fabs")
    self.Jpool.append(lambda x: sign(x[0]),lambda x:diag(x*0),"fabs")
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
    self.matrixpool.append(lambda x: norm_2(x[0]),linalg.norm,"norm_2")
    self.matrixbinarypool=FunctionPool()
    self.matrixbinarypool.append(lambda a: a[0]+a[1],lambda a: a[0]+a[1],"Matrix+Matrix")
    self.matrixbinarypool.append(lambda a: a[0]-a[1],lambda a: a[0]-a[1],"Matrix-Matrix")
    self.matrixbinarypool.append(lambda a: a[0]*a[1],lambda a: a[0]*a[1],"Matrix*Matrix")
    self.matrixbinarypool.append(lambda a: fmax(a[0],a[1]),lambda a: fmax(a[0],a[1]),"fmin")
    self.matrixbinarypool.append(lambda a: fmin(a[0],a[1]),lambda a: fmin(a[0],a[1]),"fmax")
    #self.matrixbinarypool.append(lambda a: dot(a[0],trans(a[1])),lambda a: dot(a[0].T,a[1]),name="dot(Matrix,Matrix)")
    self.matrixbinarypool.append(lambda a: mtimes(a[0],a[1].T),lambda a: np.dot(a[0],a[1].T),"dot(Matrix,Matrix.T)")

    #self.pool.append(lambda x: erf(x[0]),erf,"erf") # numpy has no erf

  def test_scalarSX(self):
      x=SX.sym("x")
      x0=0.738

      self.numpyEvaluationCheckPool(self.pool,[x],x0,name="scalarSX")

  def test_gradient(self):
      self.message("jacobian of SX**number")
      x=SX.sym("x");
      x0=1;
      p=3 # increase to 20 to showcase ticket #56
      y=x**p;
      dx=jacobian(y,x);
      dxr=p;
      self.evaluationCheck([dx],dxr,[x],x0,name="jacobian");
      dxr=1
      for i in list(range(p)):
        y=jacobian(y,x)
        dxr=dxr*(p-i)


      self.evaluationCheck([y],dxr,[x],x0,name="recursive jacobian");

  def test_gradient2(self):
      self.message("jacobian of SX**SX")
      x=SX.sym("x");
      p=SX.sym("p");
      x0=1;
      p0=3 # increase to 20 to showcase ticket #56
      y=x**p;
      dx=jacobian(y,x);
      #print dx
      dxr=p0;
      self.evaluationCheck([dx],dxr,[x,p],[x0,p0],name="jacobian");

      dxr=1
      for i in list(range(p0)):
        y=jacobian(y,x)
        dxr=dxr*(p0-i)

      self.evaluationCheck([y],dxr,[x,p],[x0,p0],name="jacobian");

  def test_SXJacobian(self):
      self.message("SX(1,1) unary operation, jacobian")
      x=SX.sym("x")
      x0=array([[0.738]])

      def fmod(f,x):
        J=f.jacobian_old(0, 0)
        return J

      self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="SX unary operations, jacobian",fmod=fmod)

  def test_SXJac(self):
      self.message("SX(1,1) unary operation, jac")
      x=SX.sym("x")
      x0=array([[0.738]])

      def fmod(f,x):
          y = f.call(x)
          J = Function('J', x, [jacobian(y[0],x[0])])
          return J

      self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="SX unary operations, jac",fmod=fmod)

  def test_SXJacobians(self):
      self.message("SX(3,1) unary operation, jacobian")
      x=SX.sym("x",3)
      x0=array([0.738,0.9,0.3])

      def fmod(f,x):
        J=f.jacobian_old(0, 0)
        return J

      self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="SX unary operations, jacobian",fmod=fmod)

  def test_SXJacobians2(self):
      self.message("SX(1,3) unary operation, jacobian")
      x=SX.sym("x",1,3)

      x0=array([0.738,0.9,0.3])

      def fmod(f,x):
        J=f.jacobian_old(0, 0)
        return J

      self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="SX unary operations, jacobian",fmod=fmod)

  def test_SX(self):
      self.message("SX unary operations")
      x=SX.sym("x",3,2)
      x0=array([[0.738,0.2],[ 0.1,0.39 ],[0.99,0.999999]])

      self.numpyEvaluationCheckPool(self.pool,[x],x0,name="SX")

      x=SX.sym("x",3,3)
      x0=array([[0.738,0.2,0.3],[ 0.1,0.39,-6 ],[0.99,0.999999,-12]])
      #self.numpyEvaluationCheck(lambda x: c.det(x[0]), lambda   x: linalg.det(x),[x],x0,name="det(SX)")
      self.numpyEvaluationCheck(lambda x: SX(c.det(x[0])), lambda   x: linalg.det(x),[x],x0,name="det(SX)")
      self.numpyEvaluationCheck(lambda x: c.inv(x[0]), lambda   x: linalg.inv(x),[x],x0,name="inv(SX)")

  def test_SXSparse(self):
      self.message("SX unary operations, sparse")
      x=SX.sym("x")
      y=SX.sym("y")
      z=SX.sym("z")
      x=SX(Sparsity(4,3,[0,2,2,3],[1,2,1]),vertcat(*[x,y,z]))
      if scipy_available:
        x0=DM(Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).sparse()

        self.numpyEvaluationCheckPool(self.pool,[x],array(x0.todense()),name="SX",setx0=x0,excludeflags={'nozero'})
      else:
        x0=DM(Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).full()

        self.numpyEvaluationCheckPool(self.pool,[x],x0,name="SX",setx0=x0,excludeflags={'nozero'})

  def test_SXbinary(self):
      self.message("SX binary operations")
      x=SX.sym("x",3,2)
      y=SX.sym("x",3,2)
      x0=array([[0.738,0.2],[ 0.1,0.39 ],[0.99,0.999999]])
      y0=array([[1.738,0.6],[ 0.7,12 ],[0,-6]])
      self.numpyEvaluationCheckPool(self.matrixbinarypool,[x,y],[x0,y0],name="SX")
      self.assertRaises(RuntimeError, lambda : mtimes(x,y))


  def test_DMbinary(self):
      self.message("SX binary operations")
      x=SX.sym("x",3,2)
      y=SX.sym("x",3,2)
      x0=array([[0.738,0.2],[ 0.1,0.39 ],[0.99,0.999999]])
      y0=array([[1.738,0.6],[ 0.7,12 ],[0,-6]])
      for f,fr,label,flags in self.matrixbinarypool.zip():
        self.checkarray(f(vertcat(*[x0,y0])),fr(vertcat(*[x0,y0])),label)

  def test_SXbinarySparse(self):
      self.message("SX binary operations")
      x=SX.sym("x")
      y=SX.sym("y")
      z=SX.sym("z")
      x2=SX.sym("x2")
      y2=SX.sym("y2")
      z2=SX.sym("z2")
      xx=SX(Sparsity(4,3,[0,2,2,3],[1,2,1]),vertcat(*[x,y,z]))
      yy=SX(Sparsity(4,3,[0,2,2,3],[0,2,3]),vertcat(*[x2,z2,y2]))

      if scipy_available:
        x0=DM(Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).sparse()
        y0=DM(Sparsity(4,3,[0,2,2,3],[0,2,3]),[1.738,0.7,-6]).sparse()

        self.numpyEvaluationCheckPool(self.matrixbinarypool,[xx,yy],[array(x0.todense()),array(y0.todense())],name="SX",setx0=[x0,y0])
      else:
        x0=DM(Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).full()
        y0=DM(Sparsity(4,3,[0,2,2,3],[0,2,3]),[1.738,0.7,-6]).full()

        self.numpyEvaluationCheckPool(self.matrixbinarypool,[xx,yy],[x0,y0],name="SX",setx0=[x0,y0])
      self.assertRaises(RuntimeError, lambda : mtimes(xx,yy))


  @known_bug()  # Test refactoring, cf. #1436
  def test_SXslicing(self):
      self.message("SX slicing/indexing")
      x=SX.sym("x",3,2)
      x0=array([[0.738,0.2],[ 0.1,0.39 ],[0.99,0.999999]])

      self.message(":dense")
      self.numpyEvaluationCheck(lambda x: SX(x[0][0,0]), lambda x: matrix(x)[0,0],[x],x0,name="x[0,0]")
      self.numpyEvaluationCheck(lambda x: SX(x[0][1,0]), lambda x: matrix(x)[1,0],[x],x0,name="x[1,0]")
      self.numpyEvaluationCheck(lambda x: SX(x[0][0,1]), lambda x: matrix(x)[0,1],[x],x0,name="x[1,0]")
      self.numpyEvaluationCheck(lambda x: SX(x[0][0,-1]), lambda x: matrix(x)[0,-1],[x],x0,name="x[0,-1]")
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

      x=SX(Sparsity(4,3,[0,2,2,3],[1,2,1]),vertcat(*[SX.sym("x"),SX.sym("y"),SX.sym("z")]))
      sx0=[0.738,0.39,0.99]
      x0=DM(Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.39,0.99]).full()
      self.numpyEvaluationCheck(lambda x: SX(x[0][0,0]), lambda x: matrix(x)[0,0],[x],x0,name="x[0,0]",setx0=[sx0])
      self.numpyEvaluationCheck(lambda x: SX(x[0][0,0]), lambda x: matrix(x)[0,0],[x],x0,name="x[0,0]",setx0=[sx0])
      self.numpyEvaluationCheck(lambda x: SX(x[0][1,0]), lambda x: matrix(x)[1,0],[x],x0,name="x[1,0]",setx0=[sx0])
      self.numpyEvaluationCheck(lambda x: SX(x[0][0,1]), lambda x: matrix(x)[0,1],[x],x0,name="x[1,0]",setx0=[sx0])
      self.numpyEvaluationCheck(lambda x: SX(x[0][0,-1]), lambda x: matrix(x)[0,-1],[x],x0,name="x[0,-1]",setx0=[sx0])
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
    x=SX.sym("x")
    y=SX.sym("y")
    f=Function("f", [vertcat(*[x,y])],[vertcat(*fun(x,y))])
    L=[2,3]
    f_in = [0]*f.n_in();f_in[0]=L
    f_out = f.call(f_in)
    z=f_out[0].full()
    zr=fun(*L)
    for i in range(3):
      self.assertAlmostEqual(z[i], zr[i],10,'SXfunction output in correct')
    self.message("SXFunction jacobian evaluation")
    J=f.jacobian_old(0, 0)
    J_in = [0]*J.n_in();J_in[0]=L
    J_out = J.call(J_in)
    Jr=matrix([[1,1],[3,2],[4,27]])
    self.checkarray(J_out[0],Jr,"SXfunction jacobian evaluates incorrectly")

  def test_SX2(self):
    self.message("SXFunction evalution 2")
    fun = lambda x,y: [3-sin(x*x)-y, sqrt(y)*x]
    # variables
    x = SX.sym("x")
    y = SX.sym("y")

    # Create function
    f = fun(x,y)
    if GlobalOptions.getSimplificationOnTheFly():
      self.assertEqual(str(f),'[SX(((3-sin(sq(x)))-y)), SX((sqrt(y)*x))]','SX representation is wrong')
    else:
      self.assertEqual(str(f),'[SX(((3-sin((x*x)))-y)), SX((sqrt(y)*x))]','SX representation is wrong'+str(f))
    fcn = Function("fcn", [vertcat(*[x,y])],[vertcat(*f)])

    self.assertEqual(repr(fcn),'Function(fcn:(i0[2])->(o0[2]) SXFunction)','SX representation is wrong')

    # Pass inputs
    L=[2,3]
    fcn_in = [0]*fcn.n_in();fcn_in[0]=L

    # Evaluate numerically
    fcn_out = fcn.call(fcn_in)

    # Get the results
    res = tuple(fcn_out[0].nonzeros())
    self.assertAlmostEqual(res[0], fun(*L)[0],10,'SXfunction evaluation wrong')
    self.assertAlmostEqual(res[1], fun(*L)[1],10,'SXfunction evaluation wrong')

  def test_SX_func(self):
    self.message("Function constructors")
    x0=SX.sym("x")
    x1=SX.sym("x")
    x2=SX.sym("x")
    x3=SX.sym("x")
    x4=SX.sym("x")
    x5=SX.sym("x")
    x6=SX.sym("x")
    y=SX.sym("y",2,3)

    f=Function("f", [y],[y])
    self.checkarray(f.size_in(0),(2,3),"Function constructors")
    self.checkarray(f.size_out(0),(2,3),"Function constructors")

    self.assertRaises(NotImplementedError,lambda: Function("f", y,[y,y]))
    self.assertRaises(NotImplementedError,lambda: Function("f", x0,[x0,x1]))

  def test_evalfail(self):
    self.message("eval fail test")
    x = SX.sym("x",2,2)
    f = Function("f", [x], [x])
    self.assertRaises(NotImplementedError,lambda: f.call(x))

  def test_SXconversion(self):
    self.message("Conversions from and to SX")
    y=SX.sym("y")
    x=SX.sym("x",3,3)
    SX(y)
    SX(x)
    c.det(x)
    y=array(DM(x))
    c.det(y)

  def test_SXbool(self):
    self.message("bool")

    x = SX.sym("x")
    y = SX.sym("y")

    f = Function("f", [vertcat(*[x,y])],[vertcat(*[logic_and(x,y),logic_or(x,y),logic_not(x)])])


    for t1 in [0,1]:
      for t2 in [0,1]:
        T1 = t1!=0
        T2 = t2!=0
        f_in = [0]*f.n_in();f_in[0]=[t1,t2]
        f_out = f.call(f_in)
        self.checkarray(f_out[0],DM([T1 and T2,T1 or T2,not T1]),"bool(%d,%d): %s" % (t1,t2,str(f_out[0])))

  def test_SXineq(self):
    self.message("SX ineq")

    x = SX.sym("x")
    y = SX.sym("y")

    f = Function("f", [vertcat(*[x,y])],[vertcat(*[x<y,x<=y,x>=y,x==y,x!=y])])


    for t1 in [-10,0.1,0,1,10]:
      for t2 in [-10,0.1,0,1,10]:
        T1 = t1
        T2 = t2
        f_in = [0]*f.n_in();f_in[0]=[t1,t2]
        f_out = f.call(f_in)
        self.checkarray(f_out[0],DM([T1 < T2,T1 <= T2, T1 >= T2, T1 == T2, T1 != T2]),"ineq(%d,%d)" % (t1,t2))



  def test_SX_func2(self):
    self.message("SXmatrix typemaps constructors")
    simplify(SX.sym("x"))
    list = [    ("number",2.3, (1,1)),
                ("SX", SX.sym("x"), (1,1))
    ];
    for name, arg,shape in list:
      self.message(":" + name)
      i=c.transpose(c.transpose(arg))
      self.assertEqual(i.shape[0],shape[0],"shape mismatch")
      self.assertEqual(i.shape[1],shape[1],"shape mismatch")
      SX(arg).is_empty()

  def test_SX_func3(self):
    self.message("vector(SXmatrix) typemaps constructors")
    y=SX.sym("y")
    x=SX.sym("x",3,1)
    vertcat(*[x,x])
    vertcat(*[y,y])
    vertcat(*[x,[]])

  def test_eval(self):
    self.message("Function eval")
    x=SX.sym("x",2,2)
    y=SX.sym("y",2,2)
    f  = Function("f", [x,y], [x*y])
    f(x,y)

  def test_symbolcheck(self):
    self.message("Check if non-symbolic inputs are caught")
    self.assertRaises(RuntimeError, lambda : Function("f", [SX(0)],[SX.sym("x")]))

  def test_sparseconstr(self):
    self.message("Check sparsity constructors")
    self.checkarray(DM.ones(Sparsity.lower(3)).full(),matrix([[1,0,0],[1,1,0],[1,1,1]]),"tril")
    self.checkarray(DM.ones(Sparsity.diag(3)).full(),matrix([[1,0,0],[0,1,0],[0,0,1]]),"diag")

  def test_subsassignment(self):
    self.message("Check subscripted assignment")

    import numpy
    numpy.random.seed(42)
    xn = numpy.random.random((3,4))

    x=DM(xn)

    y=DM(7,8)
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
    y.nz[kl]=DM(kl)

    cnt=0
    for k in kl:
      z[s.row()[k],s.get_col()[k]]=kl[cnt]
      cnt+=1
    self.checkarray(y,z,"nonzero range assignment")

  @skip(not GlobalOptions.getSimplificationOnTheFly())
  def test_substitute(self):
    self.message("Basic symbolic algebra: substitute")
    x=SX.sym("x")
    y=SX.sym("y")
    z = cos(x)*y
    self.assertTrue(depends_on(z,y))
    self.assertTrue(depends_on(z,x))
    w = substitute(z,x,0)
    self.assertTrue(w.is_symbolic())
    self.assertTrue(depends_on(w,y))
    self.assertFalse(depends_on(w,x))
    self.assertTrue(is_equal(w,y))
    r=w-y
    self.assertFalse(r.is_symbolic())
    self.assertTrue(r.is_zero())
    self.assertEqual(float(r),0)
    self.assertEqual(float(r),0)
    y = SX.sym("y",2)
    y = substitute(y+6,y,0)
    self.assertEqual(int(y[0]),6)
    self.assertEqual(int(y[1]),6)

  def test_primitivefunctions(self):
    self.message("Primitive functions")
    x=SX.sym("x")

    nums = [-2,-1.5,-1,-0.5,-0.25,0,0.25,0.5,1,1.5,2]

    def test(fun,comment,nums,reference):
      self.message(":"+comment)
      f = Function("f", [x],[fun(x)])
      for n,r in zip(nums,reference):
        f_in = [0]*f.n_in();f_in[0]=n
        f_out = f.call(f_in)
        self.assertEqual(f_out[0][0],r)

    test(casadi.sign,"sign",nums,[-1,-1,-1,-1,-1,0,1,1,1,1,1])
    test(casadi.heaviside,"heaviside",nums,[0,0,0,0,0,0.5,1,1,1,1,1])
    test(casadi.ramp,"ramp",nums,[0,0,0,0,0,0,0.25,0.50,1,1.5,2])
    test(casadi.rectangle,"rectangle",nums,[0,0,0,0.5,1,1,1,0.5,0,0,0])
    test(casadi.triangle,"triangle",nums,[0,0,0,0.5,0.75,1,0.75,0.5,0,0,0])


  def test_taylor(self):
    self.message("univariate taylor expansion")
    x=SX.sym("x")

    if GlobalOptions.getSimplificationOnTheFly():
      self.assertTrue(is_equal(taylor(sin(x),x),x))

    a_=0.13
    x_=0.15

    a = SX.sym("a")

    def test(e,r):
      f = Function("f", [x,a],[e])
      f_in = [0]*f.n_in();f_in[0]=x_
      f_in[1]=a_
      f_out = f.call(f_in)
      self.assertAlmostEqual(f_out[0][0],r,10)

    test(taylor(sin(x),x,a,0),sin(a_))
    test(taylor(sin(x),x,a,1),sin(a_)+cos(a_)*(x_-a_))
    test(taylor(sin(x),x,a,2),sin(a_)+cos(a_)*(x_-a_)-(sin(a_)*(x_-a_)**2)/2.0)
    test(taylor(sin(x),x,a,3),sin(a_)+cos(a_)*(x_-a_)-(sin(a_)*(x_-a_)**2)/2.0-(cos(a_)*(x_-a_)**3)/6.0)

    M=blockcat([[a*sin(x),a*cos(x)],[exp(a*x),a*x**2],[cos(x),0]])
    f = Function("f", [x,a],[taylor(M,x)])
    f_in = [0]*f.n_in();f_in[0]=x_
    f_in[1]=a_
    f_out = f.call(f_in)
    self.checkarray(f_out[0],matrix([[x_*a_,a_],[1+a_*x_,0],[1,0]]),"taylor on dense matrices")

  def test_null(self):
    self.message("Function null")
    x = SX.sym("x")

    f = Function("f", [x],[x**2,[]])
    f_out = f.call([0])
    self.assertTrue(f_out[1].is_empty())

    f = Function("f", [x,[]],[x**2,[]])
    f_out = f.call([0,0])
    self.assertTrue(f_out[1].is_empty())
    f_out = f.call([0,0])

    r = f.call([x,[]])
    self.assertTrue(r[1].is_empty())

    r = f.call([x,[]])
    self.assertTrue(r[1].is_empty())

    r = f.call([x,SX(0,1)])
    self.assertTrue(r[1].is_empty())

    r = f.call([x,SX(1,0)])
    self.assertTrue(r[1].is_empty())

    #self.assertRaises(Exception,lambda : f([x,x]))
    #self.assertRaises(Exception,lambda : f([[],[]]))

  def test_mtaylor(self):
    self.message("multivariate taylor expansions")
    x=SX.sym("x")
    y=SX.sym("y")
    a=SX.sym("a")
    b=SX.sym("b")

    a_=0.13
    x_=0.15

    b_=0.73
    y_=0.75

    def test(e,r):
      f = Function("f", [vertcat(*[x,y]),vertcat(*[a,b])],[e])
      f_in = [0]*f.n_in();f_in[0]=[x_,y_]
      f_in[1]=[a_,b_]
      f_out = f.call(f_in)
      self.assertAlmostEqual(f_out[0][0],r,10)

    test(mtaylor(sin(x+y),vertcat(*[x,y]),vertcat(*[a,b]),0),sin(a_+b_))
    test(mtaylor(sin(x+y),vertcat(*[x,y]),vertcat(*[a,b]),1),sin(a_+b_)+(cos(b_+a_)*(x_-a_)+cos(b_+a_)*(y_-b_)))

    def sol(x,y,a,b):
      return sin(b+a)+(cos(b+a)*(x-a)+cos(b+a)*(y-b))-(sin(b+a)*(x-a)**2+2*sin(b+a)*(y-b)*(x-a)+sin(b+a)*(y-b)**2)/2

    test(mtaylor(sin(x+y),vertcat(*[x,y]),vertcat(*[a,b]),2),sol(x_,y_,a_,b_))

    def sol(x,y,a,b):
      return sin(b+a)+(cos(b+a)*(x-a)+cos(b+a)*(y-b))-(sin(b+a)*(x-a)**2+2*sin(b+a)*(y-b)*(x-a)+sin(b+a)*(y-b)**2)/2-(cos(b+a)*(x-a)**3+3*cos(b+a)*(y-b)*(x-a)**2+3*cos(b+a)*(y-b)**2*(x-a)+cos(b+a)*(y-b)**3)/6

    test(mtaylor(sin(x+y),vertcat(*[x,y]),vertcat(*[a,b]),3),sol(x_,y_,a_,b_))

    def sol(x,y,a,b):
      return (-2*sin(b+a)*(x-a)*(y-b)-sin(b+a)*(x-a)**2)/2+cos(b+a)*(y-b)-(cos(b+a)*(x-a)**3)/6+cos(b+a)*(x-a)+sin(b+a)

    test(mtaylor(sin(x+y),vertcat(*[x,y]),vertcat(*[a,b]),3,[1,2]),sol(x_,y_,a_,b_))

    test(mtaylor(sin(x+y),vertcat(*[x,y]),vertcat(*[0,0]),4,[1,2]),(-3*x_**2*y_-x_**3)/6+y_+x_)

  def test_issue107(self):
    self.message("Regression test for issue 107: +=")
    x=SX.sym("x")
    y=SX.sym("y")

    z=x
    z+=y

    self.assertTrue(x.is_symbolic())
    self.assertFalse(z.is_symbolic())

    x=SX.sym("x")
    y=SX.sym("y")

    z=x
    z+=y

    self.assertTrue(x.is_symbolic())
    self.assertFalse(z.is_symbolic())

  def test_evalchecking(self):
    x = SX.sym("x",1,5)

    y = SX.sym("y",1,3)
    z = SX.sym("z",5,1)
    q = SX.sym("z",1,6)

    f = Function("f", [x],[x**2])

    self.assertRaises(RuntimeError, lambda : f(y))
    self.assertRaises(RuntimeError, lambda : f(q))
    f(z)

  def test_indexinglimits(self):
    self.message("Limits of indexing")
    y = casadi.SX.sym("y", 3)
    self.assertRaises(RuntimeError,lambda : y[[0, 5]] )
    try:
      y[[0, 5]] = SX.sym("a")
      self.assertTrue(False)
    except RuntimeError:
      pass
    y[[0, 2]]
    y[[0, 2]] = SX.sym("a")

  def test_issue181(self):
    self.message("Regression test #181")
    x = SX.sym("x")
    #self.assertRaises(TypeError,lambda : SX([x,None]))  # FIXME: this is leaking memory
    self.assertRaises(NotImplementedError,lambda: Function("f", [[x], [None]], [[2 * x]]))

  @known_bug()  # Not implemented
  def test_is_equal(self):
    self.message("equivalent")
    x = SX.sym("x")
    a = x*x
    b = x*x
    self.assertTrue(a.is_equal(b,1))

  @skip(not GlobalOptions.getSimplificationOnTheFly())
  def test_SXsimplifications(self):
    self.message("simplifications")
    x = SX.sym("x")

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
    ops.append(lambda x: (SX(4)-SX(4))+x)


    y = SX.sym("x")

    ops.append(lambda x: ((x+y)-(y+x))+x)
    ops.append(lambda x: ((x*y)-(y*x))+x)
    ops.append(lambda x: ((-x)-(-x))+x)

    for op in ops:
      y = op(x)
      f = Function("f", [x],[y])
      f_in = [0]*f.n_in();f_in[0]=0.3
      f_out = f.call(f_in)
      self.checkarray(f_out[0],array(DM(op(0.3))),"simplifications")
      self.assertEqual(str(y),"x")

      y = op(-x)
      f = Function("f", [x],[y])
      f_in = [0]*f.n_in();f_in[0]=0.3
      f_out = f.call(f_in)
      self.checkarray(f_out[0],array(DM(op(-0.3))),"simplifications")
      self.assertEqual(str(y),"(-x)")

  def test_truth(self):
    self.message("Truth values")
    self.assertRaises(Exception, lambda : bool(SX.sym("x")))
    self.assertRaises(Exception, lambda : bool(SX.sym("x")>0))
    self.assertTrue(bool(SX(1)))
    self.assertFalse(bool(SX(0)))
    self.assertTrue(bool(SX(0.2)))
    self.assertTrue(bool(SX(-0.2)))
    self.assertRaises(Exception, lambda : bool(SX.sym("x")))
    self.assertRaises(Exception, lambda : bool(SX.sym("x")>0))
    self.assertTrue(bool(SX(SX(1))))
    self.assertFalse(bool(SX(SX(0))))
    self.assertTrue(bool(SX(SX(0.2))))
    self.assertTrue(bool(SX(SX(-0.2))))
    self.assertRaises(Exception, lambda : bool(SX([2.0,3])))

  def test_if_else(self):
    x = SX.sym("x")
    y = if_else(x,1,2)
    f = Function("f", [x],[y])
    f_in = [0]*f.n_in();f_in[0]=1
    f_out = f.call(f_in)
    self.assertTrue(f_out[0]==1,"if_else")
    f_in = [0]*f.n_in();f_in[0]=0
    f_out = f.call(f_in)
    self.assertTrue(f_out[0]==2,"if_else")

    x0 = 2.1
    y = if_else(x>1,x**2,x**3)
    f = Function("f", [x],[y])
    f_in = [0]*f.n_in();f_in[0]=x0
    f_out = f.call(f_in)
    self.checkarray(f_out[0],x0**2,"if_else sens")

    x0 = -2.1
    f_in = [0]*f.n_in();f_in[0]=x0
    f_out = f.call(f_in)
    self.checkarray(f_out[0],x0**3,"if_else sens")


  def test_is_regular(self):
    x = SX.sym("x")

    self.assertTrue(SX(0).is_regular())
    self.assertFalse(SX(inf).is_regular())
    with self.assertRaises(Exception):
      self.assertTrue(x.nz[0])

    self.assertTrue(SX(DM([0,1])).is_regular())
    self.assertFalse(SX(DM([0,inf])).is_regular())
    self.assertFalse(vertcat(*[x,inf]).is_regular())
    with self.assertRaises(Exception):
      self.assertFalse(vertcat(*[x,x]).is_regular())


  def test_symvar(self):
    a = SX.sym("a")
    b = SX.sym("b")
    c = SX.sym("c")
    e = cos(a*b) + c
    w = symvar(e)
    self.assertEqual(len(w),3)
    if GlobalOptions.getSimplificationOnTheFly():
      self.assertTrue(is_equal(w[0],a))
      self.assertTrue(is_equal(w[1],b))
      self.assertTrue(is_equal(w[2],c))

  def test_poly_coeff(self):
    x =SX.sym("x")
    a= SX.sym("a")
    c=SX.sym("c")
    p=poly_coeff(12*x**4+x**2+a*x+c,x)
    self.assertTrue(is_equal(p[0],12))
    self.assertTrue(is_equal(p[1],0))
    self.assertTrue(is_equal(p[2],1))
    self.assertTrue(is_equal(p[3],a))
    self.assertTrue(is_equal(p[4],c))

    p=poly_coeff((x-a)*(x+a),x)
    self.assertTrue(is_equal(p[0],1))
    self.assertTrue(is_equal(p[1],0))

  def test_poly_roots(self):

    p = SX.sym("[a,b]")
    r = poly_roots(p)

    f = Function("f", [p],[r])
    f_in = [0]*f.n_in()
    f_in[0]=DM([2,7])
    a_ = f_in[0][0]
    b_ = f_in[0][1]
    f_out = f.call(f_in)
    f_out[0]
    self.checkarray(f_out[0],vertcat(*[-b_/a_]))

    p = SX.sym("[a,b]")
    r = poly_roots(vertcat(*[p,0]))

    f = Function("f", [p],[r])
    f_in = [0]*f.n_in();f_in[0]=DM([2,7])
    a_ = f_in[0][0]
    b_ = f_in[0][1]
    f_out = f.call(f_in)
    f_out[0]
    self.checkarray(f_out[0],vertcat(*[-b_/a_,0]))

    p = SX.sym("[a,b,c]")
    r = poly_roots(p)

    f = Function("f", [p],[r])
    f_in = [0]*f.n_in();f_in[0]=DM([1.13,7,3])
    a_ = f_in[0][0]
    b_ = f_in[0][1]
    c_ = f_in[0][2]
    d = b_**2-4*a_*c_
    f_out = f.call(f_in)
    x0 = (-b_-sqrt(d))/2/a_
    x1 = (-b_+sqrt(d))/2/a_
    f_out[0]
    self.checkarray(f_out[0],vertcat(*[x0,x1]))

    p = SX.sym("[a,b,c,d]")
    r = poly_roots(p)

    f = Function("f", [p],[r])
    f_in = [0]*f.n_in();f_in[0]=DM([11,1.3,-1.7,0.1])
    f_out = f.call(f_in)
    f_out[0]
    self.checkarray(f_out[0],DM([0.298028,-0.479787,0.0635774]),digits=5)

    p = SX.sym("[a,b,c,d,e]")
    r = poly_roots(p)

    f = Function("f", [p],[r])
    f_in = [0]*f.n_in();f_in[0]=DM([3,6,-123,  -126,1080])
    f_out = f.call(f_in)
    f_out[0]
    self.checkarray(f_out[0],DM([5,3,-4,-6]),digits=5)

  def test_eig_symbolic(self):
    x = SX.sym("x",2,2)
    f = Function("f", [x],[eig_symbolic(x)])
    f_in = [0]*f.n_in();f_in[0]=DM([[2,0.1],[0.3,0.7]])
    f_out = f.call(f_in)
    self.checkarray(f_out[0],DM([0.67732,2.02268]),digits=5)


    x = SX.sym("x",2)
    f = Function("f", [x],[eig_symbolic(c.diag(x))])
    f_in = [0]*f.n_in();f_in[0]=DM([3,7])
    f_out = f.call(f_in)
    self.checkarray(f_out[0],f_in[0])


    x = SX.sym("x",5)
    f = Function("f", [x],[eig_symbolic(c.diag(x))])
    f_in = [0]*f.n_in();f_in[0]=DM([3,7,2,1,6])
    f_out = f.call(f_in)
    self.checkarray(f_out[0],f_in[0])

    x = SX.sym("x",2,2)
    y = SX.sym("y",2)
    f = Function("f", [x,y],[eig_symbolic(diagcat(*[x,c.diag(y)]))])
    f_in = [0]*f.n_in();f_in[0]=DM([[2,0.1],[0.3,0.7]])
    f_in[1]=[3,7]
    f_out = f.call(f_in)
    self.checkarray(f_out[0],DM([0.67732,2.02268,3,7]),digits=5)

    x = SX.sym("x",3,3)
    x[2,0] = 0
    x[1,0] = 0

    x = sparsify(x)

    e = eig_symbolic(x)

    f = Function("f", [x],[e])
    f_in = [0]*f.n_in();f_in[0]=DM(f.sparsity_in(0),list(range(1,8)))
    f_in[0].print_dense()
    f_out = f.call(f_in)
    self.checkarray(f_out[0],DM([1,-0.29150,10.29150]),digits=5)


    x = SX.sym("x",3,3)
    x[2,0] = 0
    x[1,0] = 0
    x[2,1] = 0

    x = sparsify(x)

    e = eig_symbolic(x)

    f = Function("f", [x],[e])
    f_in = [0]*f.n_in();f_in[0]=DM(f.sparsity_in(0),list(range(1,7)))
    f_in[0].print_dense()
    f_out = f.call(f_in)
    self.checkarray(f_out[0],DM([1,3,6]),digits=5)

    x = SX.sym("x",Sparsity.upper(5))

    f = Function("f", [x],[eig_symbolic(x)])
    fin = DM(x.sparsity(),0)
    fin[Sparsity.diag(5)] = c.diag(list(range(5)))
    self.checkarray(f(fin), DM(list(range(5))))

  def test_jacobian_empty(self):
    x = SX.sym("x",3)

    s = jacobian(DM(0,0),x).shape
    self.assertEqual(s[0],0)
    self.assertEqual(s[1],3)

    s = jacobian(x,SX.sym("x",0,4)).shape
    self.assertEqual(s[0],3)
    self.assertEqual(s[1],0)

  def test_empty_SX(self):
    s = SX([]).shape
    self.assertEqual(s[0],0)
    self.assertEqual(s[1],1)
    x = vertcat(*(SX.sym("x"),SX([])))

  def test_mul_sparsity(self):

    N = 10
    x = SX.sym("x",N,N)
    y = SX.sym("y",N,N)

    x_ = self.randDM(N,N)
    y_ = self.randDM(N,N)

    filt = Sparsity.diag(N)+Sparsity.triplet(N,N,[1],[3])

    f = Function("f", [x,y],[mtimes(x,y)])
    f_in = [0]*f.n_in();f_in[0]=x_
    f_in[1]=y_
    g = Function("g", [x,y],[mac(x,y,SX.zeros(filt))])
    g_in = [0]*g.n_in();g_in[0]=x_
    g_in[1]=y_

    f_out = f.call(f_in)
    g_out = g.call(g_in)

    self.checkarray(IM.ones(filt),IM.ones(g.sparsity_out(0)))

    self.checkarray(f_out[0][filt],g_out[0])

  @skip(platform_arch==32)
  @memory_heavy()
  @unittest.skipIf(sys.version_info >= (3, 0),"pickle is not compatible")
  def test_large_hessian(self):
    import pickle

    A = pickle.load(open("../data/apoa1-2.pkl","r"))

    H = DM(A,list(range(A.nnz())))
    H = H + H.T

    H = H[:20000,:20000]

    x = SX.sym("x",H.size1())

    f = Function("f", [x],[mtimes([x.T,H,x])], {'verbose':True})
    H *= 2

    h = f.hessian_old(0, 0)
    h_out = h.call([0])

    self.assertTrue(h.sparsity_out(0)==H.sparsity())

    self.checkarray(h_out[0].nonzeros(),H.nonzeros())

  def test_mxnulloutput(self):
     a = SX(5,0)
     b = SX.sym("x",2)
     bm = MX.sym("x",2)

     f = Function("f", [b],[a])
     c = f(bm)

     self.assertEqual(c.size1(),5)
     self.assertEqual(c.size2(),0)

     c = f(b)

     self.assertEqual(c.size1(),5)
     self.assertEqual(c.size2(),0)

     a = SX(0,0)

     f = Function("f", [b],[a])

     c = f(bm)

     self.assertEqual(c.size1(),0)
     self.assertEqual(c.size2(),0)

     c = f(b)

     self.assertEqual(c.size1(),0)
     self.assertEqual(c.size2(),0)

  def test_mxnull(self):
     a = SX(5,0)
     b = SX(0,3)

     c = mtimes(a,b)

     self.assertEqual(c.nnz(),0)

     a = SX(5,3)
     b = SX(3,4)

     c = mtimes(a,b)

     self.assertEqual(c.nnz(),0)

  def  test_mxnullop(self):
    c = SX(0,0)
    x = SX.sym("x",2,3)

    with self.assertRaises(RuntimeError):
      d = x + c

    with self.assertRaises(RuntimeError):
      d = x / c

  def test_copysign(self):
    x = SX.sym("x")
    y = SX.sym("y")
    z = copysign(x,y)

    f = Function("f", [x,y],[z])

    f_in = [0]*f.n_in();f_in[0]=2
    f_in[1]=0.5
    f_out = f.call(f_in)
    self.checkarray(f_out[0],DM([2]))

    f_in = [0]*f.n_in();f_in[0]=2
    f_in[1]=-0.5
    f_out = f.call(f_in)
    self.checkarray(f_out[0],DM([-2]))

    f_in = [0]*f.n_in();f_in[0]=-2
    f_in[1]=0.5
    f_out = f.call(f_in)
    self.checkarray(f_out[0],DM([2]))

    f_in = [0]*f.n_in();f_in[0]=-2
    f_in[1]=-0.5
    f_out = f.call(f_in)
    self.checkarray(f_out[0],DM([-2]))

    f_in = [0]*f.n_in();f_in[0]=2
    f_in[1]=0
    f_out = f.call(f_in)
    self.checkarray(f_out[0],DM([2]))

    J = f.jacobian_old(0, 0)

    J_in = [0]*J.n_in();J_in[0]=2
    J_in[1]=0.5
    J_out = J.call(J_in)
    self.checkarray(J_out[0],DM([1]))

    J_in = [0]*J.n_in();J_in[0]=2
    J_in[1]=-0.5
    J_out = J.call(J_in)
    self.checkarray(J_out[0],DM([-1]))

    J_in = [0]*J.n_in();J_in[0]=-2
    J_in[1]=0.5
    J_out = J.call(J_in)
    self.checkarray(J_out[0],DM([1]))

    J_in = [0]*J.n_in();J_in[0]=-2
    J_in[1]=-0.5
    J_out = J.call(J_in)
    self.checkarray(J_out[0],DM([-1]))

    J_in = [0]*J.n_in();J_in[0]=2
    J_in[1]=0
    J_out = J.call(J_in)
    self.checkarray(J_out[0],DM([1]))

    J = f.jacobian_old(1, 0)

    J_in = [0]*J.n_in();J_in[0]=2
    J_in[1]=0.5
    J_out = J.call(J_in)
    self.checkarray(J_out[0],DM([0]))

    J_in = [0]*J.n_in();J_in[0]=2
    J_in[1]=-0.5
    J_out = J.call(J_in)
    self.checkarray(J_out[0],DM([0]))

    J_in = [0]*J.n_in();J_in[0]=-2
    J_in[1]=0.5
    J_out = J.call(J_in)
    self.checkarray(J_out[0],DM([0]))

    J_in = [0]*J.n_in();J_in[0]=-2
    J_in[1]=-0.5
    J_out = J.call(J_in)
    self.checkarray(J_out[0],DM([0]))

    J_in = [0]*J.n_in();J_in[0]=2
    J_in[1]=0
    J_out = J.call(J_in)
    self.checkarray(J_out[0],DM([0]))

  def test_depends_on(self):
    a = SX.sym("a")
    b = SX.sym("b")

    self.assertTrue(depends_on(a**2,a))
    self.assertTrue(depends_on(a,a))
    self.assertFalse(depends_on(0,a))
    self.assertTrue(depends_on(a**2,vertcat(*[a,b])))
    self.assertTrue(depends_on(a,vertcat(*[a,b])))
    self.assertFalse(depends_on(0,vertcat(*[a,b])))
    self.assertTrue(depends_on(b**2,vertcat(*[a,b])))
    self.assertTrue(depends_on(b,vertcat(*[a,b])))
    self.assertTrue(depends_on(a**2+b**2,vertcat(*[a,b])))
    self.assertTrue(depends_on(a+b,vertcat(*[a,b])))
    self.assertTrue(depends_on(vertcat(*[0,a]),a))
    self.assertTrue(depends_on(vertcat(*[a,0]),a))
    self.assertTrue(depends_on(vertcat(*[a**2,b**2]),vertcat(*[a,b])))
    self.assertTrue(depends_on(vertcat(*[a,0]),vertcat(*[a,b])))
    self.assertTrue(depends_on(vertcat(*[0,b]),vertcat(*[a,b])))
    self.assertTrue(depends_on(vertcat(*[b,0]),vertcat(*[a,b])))
    self.assertFalse(depends_on(vertcat(*[0,0]),vertcat(*[a,b])))

  @requires("is_smooth")
  def test_is_smooth(self):
    x = SX.sym("a",2,2)
    import warnings
    with warnings.catch_warnings():
      warnings.simplefilter("error",DeprecationWarning)
      with self.assertRaises(Exception):
        is_smooth(x)
      warnings.simplefilter("ignore")
      is_smooth(x)

  def test_which_depends(self):
    for X in [SX,MX]:
      x = X.sym("x")
      y = X.sym("y")

      p = X.sym("p")

      e = vertcat(0,x,y,p,2*p**3,x*y,x*p,sin(x),cos(y),sqrt(x+y),p*p*x,x*y*p)

      self.checkarray(which_depends(e, vertcat(x,y),2,True),[0, 0, 0, 0,0, 1, 0, 1, 1, 1, 0, 1])
      self.checkarray(which_depends(e, vertcat(x,y),1,True),[0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1])

      z =X.sym("z")
      e = vertcat(x*p,x+y)
      self.checkarray(which_depends(e, vertcat(x,y,p,z),2,False),[True, False, True, False])
      self.checkarray(which_depends(e, vertcat(x,y,p,z),1,False),[True, True, True, False])

      e = vertcat(x*p,x+z*y)
      self.checkarray(which_depends(e, vertcat(x,y,p),2,False),[True, False, True])
      self.checkarray(which_depends(e, vertcat(x,y,p),1,False),[True, True, True])

      e = vertcat(x*p,x+z*y)
      self.checkarray(which_depends(e, vertcat(x,y,p,z),2,False),[True, True, True, True])
      self.checkarray(which_depends(e, vertcat(x,y,p,z),1,False),[True, True, True, True])

      e = vertcat(sin(x+y)+p)
      self.checkarray(which_depends(e, vertcat(x,y,p,z),2,False),[True, True, False, False])
      self.checkarray(which_depends(e, vertcat(x,y,p,z),1,False),[True, True, True, False])

      e = vertcat(sin(x)*p**2,y**2)
      #self.checkarray(which_depends(e, vertcat(x,y,p),3,True),[True, False])
      #self.checkarray(which_depends(e, vertcat(x,y,p),3,False),[True, False, True])
      self.checkarray(which_depends(e, vertcat(x,y,p),2,True),[True, True])
      self.checkarray(which_depends(e, vertcat(x,y,p),2,False),[True, True, True])

      e = vertcat(x**2*p,y)
      #self.checkarray(which_depends(e, vertcat(x,y,p),3,True),[True, False])
      #self.checkarray(which_depends(e, vertcat(x,y,p),3,False),[True, False, False])

      self.checkarray(which_depends(e, vertcat(x,y,p),2,True),[True, False])
      self.checkarray(which_depends(e, vertcat(x,y,p),2,False),[True, False, True])

  def test_if_else_zero_sens(self):

    for X in [SX]:
      x=X.sym('x')


      a = 1+3*x+sqrt(3*x)*x+7*x
      b = 1+2*x+sin(2*x)*x +x
      z = if_else(x>0,a,b)*x

      f = Function("f",[x],[z,jacobian(z,x)])
      fa = Function("f",[x],[a*x,jacobian(a*x,x)])
      fb = Function("f",[x],[b*x,jacobian(b*x,x)])

      for i,j in zip(f([3]),fa([3])):
        self.checkarray(i,j)

      for i,j in zip(f([-3]),fb([-3])):
        self.checkarray(i,j)


      f = Function("f",[x],[z])
      fa = Function("f",[x],[a*x])
      fb = Function("f",[x],[b*x])

      self.checkfunction(f,fa,inputs=[3])
      self.checkfunction(f,fb,inputs=[-3],evals=1)

  def test_pw_const(self):
      t= SX.sym("t")

      e = pw_const(t, [0,2,3],[7,1,3,5])

      E = Function("E",[t],[e])

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
      t= SX.sym("t")

      e = pw_lin(t, [0,2,3,5], [7,1,3,2])

      E = Function("E",[t],[e])

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
      x = SX.sym("x",3)
      with self.assertInException("Use an equivalent CasADi function"):
        np.linalg.norm(x)


  def test_quadratic(self):
    for X in [SX,MX]:
      x = X.sym("x")
      p = X.sym("p")
      y = X.sym("y")

      self.assertFalse(is_quadratic(sin(x),x))
      self.assertFalse(is_quadratic(x**3,x))
      self.assertTrue(is_quadratic(x**2,x))
      self.assertTrue(is_quadratic(4*x,x))
      self.assertTrue(is_quadratic(5,x))

      self.assertFalse(is_quadratic(sin(x)*p**4,x))
      self.assertFalse(is_quadratic(x**3*p**4,x))
      self.assertTrue(is_quadratic(x**2*p**4,x))
      self.assertTrue(is_quadratic(x*p**4,x))
      self.assertTrue(is_quadratic(5*p**4,x))

      self.assertFalse(is_linear(sin(x),x))
      self.assertFalse(is_linear(x**3,x))
      self.assertFalse(is_linear(x**2,x))
      self.assertTrue(is_linear(3*x,x))
      self.assertTrue(is_linear(5,x))

      self.assertFalse(is_linear(sin(x)*p**4,x))
      self.assertFalse(is_linear(x**3*p**4,x))
      self.assertFalse(is_linear(x**2*p**4,x))
      self.assertTrue(is_linear(x*p**4,x))
      self.assertTrue(is_linear(5*p**4,x))



      z = x**2+3*y**2 + 0.5*x*y + 7*x + 6*y+7
      [A,b,c] = quadratic_coeff(z,vertcat(x,y))

      with self.assertInException("non-quadratic"):
        [A,b,c] = quadratic_coeff(x**2+3*y**2 + 0.5*x*y + 7*x + 6*y+7+sin(x),vertcat(x,y))

      with self.assertInException("scalar"):
        [A,b,c] = quadratic_coeff(vertcat(x,y),x)

      z = x**2+3*y**2 + 0.5*x*y -p*y + 7*x + 6*y+7
      [A,b,c] = quadratic_coeff(z,vertcat(x,y))

      xy = vertcat(x,y)

      e = 0.5*bilin(A,xy,xy)+dot(b,xy)+c

      f = Function('f',[xy,p],[z])
      f2 = Function('f',[xy,p],[e])
      self.checkfunction(f,f2,inputs=[1.1,1.3])


      with self.assertInException("non-linear"):
        [A,b] = linear_coeff(x**2+3*y**2 + 0.5*x*y + 7*x + 6*y+7,vertcat(x,y))

      with self.assertInException("vector"):
        [A,b] = linear_coeff(blockcat([[x,y],[y,x]]),x)

      z = vertcat(7*x + 6*y+7 ,5 -p*y )
      [A,b] = linear_coeff(z,xy)

      e = mtimes(A,xy)+b

      f = Function('f',[xy,p],[z])
      f2 = Function('f',[xy,p],[e])
      self.checkfunction(f,f2,inputs=[1.1,1.3])

  def test_evalf(self):
    x = SX.sym("x")

    y = SX(5)

    self.checkarray(evalf(y),5)
    with self.assertInException("since variables [x] are free"):
      evalf(x)


if __name__ == '__main__':
    unittest.main()
