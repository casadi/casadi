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
from numpy import *
import unittest
from types import *
from helpers import *

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
    #self.matrixbinarypool.append(lambda a: inner_prod(a[0],trans(a[1])),lambda a: dot(a[0].T,a[1]),name="inner_prod(Matrix,Matrix)") 
    self.matrixbinarypool.append(lambda a: mul(a[0],a[1].T),lambda a: dot(a[0],a[1].T),"dot(Matrix,Matrix.T)")

    #self.pool.append(lambda x: erf(x[0]),erf,"erf") # numpy has no erf
    
  def test_scalarSXElement(self):
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
        J=f.jacobian()
        J.init()
        return J
      
      self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="SX unary operations, jacobian",fmod=fmod)

      def fmod(f,x):
        #f.setOption("ad_mode","forward")
        f.init();
        J=f.jacobian()
        J.init()
        return J
        
      self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="MX unary operations, jacobian",fmod=fmod)
      
      def fmod(f,x):
        #f.setOption("ad_mode","reverse")
        f.init();
        J=f.jacobian()
        J.init()
        return J
        
      self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="MX unary operations, jacobian",fmod=fmod)
      
  def test_SXJac(self):
      self.message("SX(1,1) unary operation, jac")
      x=SX.sym("x")
      x0=array([[0.738]])

      def fmod(f,x):
        j=f.jac()
        J=SXFunction(x,[j])
        J.init()
        return J
      
      self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="SX unary operations, jac",fmod=fmod)
      
  def test_SXJacobians(self):
      self.message("SX(3,1) unary operation, jacobian")
      x=SX.sym("x",3)
      x0=array([0.738,0.9,0.3])

      def fmod(f,x):
        J=f.jacobian()
        J.init()
        return J
      
      self.numpyEvaluationCheckPool(self.Jpool,[x],x0,name="SX unary operations, jacobian",fmod=fmod)
      
  def test_SXJacobians2(self):
      self.message("SX(1,3) unary operation, jacobian")
      x=SX.sym("x",1,3)
      
      x0=array([0.738,0.9,0.3])

      def fmod(f,x):
        J=f.jacobian()
        J.init()
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
      self.numpyEvaluationCheck(lambda x: SX([c.det(x[0])]), lambda   x: linalg.det(x),[x],x0,name="det(SX)")
      self.numpyEvaluationCheck(lambda x: c.inv(x[0]), lambda   x: linalg.inv(x),[x],x0,name="inv(SX)")
        
  def test_SXSparse(self):
      self.message("SX unary operations, sparse")
      x=SXElement.sym("x")
      y=SXElement.sym("y")
      z=SXElement.sym("z")
      x=SX(Sparsity(4,3,[0,2,2,3],[1,2,1]),[x,y,z])
      if scipy_available:
        x0=DMatrix(Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).toCsc_matrix()
      
        self.numpyEvaluationCheckPool(self.pool,[x],array(x0.todense()),name="SX",setx0=x0,excludeflags={'nozero'})
      else:
        x0=DMatrix(Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).toArray()
      
        self.numpyEvaluationCheckPool(self.pool,[x],x0,name="SX",setx0=x0,excludeflags={'nozero'})
      
  def test_SXbinary(self):
      self.message("SX binary operations")
      x=SX.sym("x",3,2)
      y=SX.sym("x",3,2)
      x0=array([[0.738,0.2],[ 0.1,0.39 ],[0.99,0.999999]])
      y0=array([[1.738,0.6],[ 0.7,12 ],[0,-6]])
      self.numpyEvaluationCheckPool(self.matrixbinarypool,[x,y],[x0,y0],name="SX")
      self.assertRaises(RuntimeError, lambda : mul(x,y))

  def test_SXbinarySparse(self):
      self.message("SX binary operations")
      x=SXElement.sym("x")
      y=SXElement.sym("y")
      z=SXElement.sym("z")
      x2=SXElement.sym("x2")
      y2=SXElement.sym("y2")
      z2=SXElement.sym("z2")
      xx=SX(Sparsity(4,3,[0,2,2,3],[1,2,1]),[x,y,z])
      yy=SX(Sparsity(4,3,[0,2,2,3],[0,2,3]),[x2,z2,y2])
      
      if scipy_available:
        x0=DMatrix(Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).toCsc_matrix()
        y0=DMatrix(Sparsity(4,3,[0,2,2,3],[0,2,3]),[1.738,0.7,-6]).toCsc_matrix()
        
        self.numpyEvaluationCheckPool(self.matrixbinarypool,[xx,yy],[array(x0.todense()),array(y0.todense())],name="SX",setx0=[x0,y0])
      else:
        x0=DMatrix(Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.1,0.99]).toArray()
        y0=DMatrix(Sparsity(4,3,[0,2,2,3],[0,2,3]),[1.738,0.7,-6]).toArray()
        
        self.numpyEvaluationCheckPool(self.matrixbinarypool,[xx,yy],[x0,y0],name="SX",setx0=[x0,y0])
      self.assertRaises(RuntimeError, lambda : mul(xx,yy))


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
      
      x=SX(Sparsity(4,3,[0,2,2,3],[1,2,1]),[SXElement.sym("x"),SXElement.sym("y"),SXElement.sym("z")])
      sx0=[0.738,0.39,0.99]
      x0=DMatrix(Sparsity(4,3,[0,2,2,3],[1,2,1]),[0.738,0.39,0.99]).toArray()
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
    x=SXElement.sym("x")
    y=SXElement.sym("y")
    f=SXFunction([vertcat([x,y])],[vertcat(fun(x,y))])
    f.init()
    L=[2,3]
    f.setInput(L)
    f.evaluate()
    z=f.output(0).toArray()
    zr=fun(*L)
    for i in range(3):
      self.assertAlmostEqual(z[i], zr[i],10,'SXfunction output in correct')
    self.message("SXFunction jacobian evaluation")
    J=f.jacobian()
    J.init()
    J.setInput(L)
    J.evaluate()
    Jr=matrix([[1,1],[3,2],[4,27]])
    self.checkarray(J.getOutput(0),Jr,"SXfunction jacobian evaluates incorrectly")
          
  def test_SX2(self):
    self.message("SXFunction evalution 2")
    fun = lambda x,y: [3-sin(x*x)-y, sqrt(y)*x]
    # variables
    x = SXElement.sym("x")
    y = SXElement.sym("y")

    # Create function
    f = fun(x,y)
    if CasadiOptions.getSimplificationOnTheFly():
      self.assertEqual(str(f),'[((3-sin(sq(x)))-y), (sqrt(y)*x)]','SX representation is wrong')
    else:
      self.assertEqual(str(f),'[((3-sin((x*x)))-y), (sqrt(y)*x)]','SX representation is wrong'+str(f))
    fcn = SXFunction([vertcat([x,y])],[vertcat(f)])

    # Set some options
    fcn.setOption("name","f")

    self.assertEqual(repr(fcn),'function("f")','SX representation is wrong')
    # Initialize the function for numerical calculations
    fcn.init()

    # Pass inputs
    L=[2,3]
    fcn.setInput(L)

    # Evaluate numerically
    fcn.evaluate()

    # Get the results
    res = tuple(fcn.output().data())
    self.assertAlmostEqual(res[0], fun(*L)[0],10,'SXfunction evaluation wrong')
    self.assertAlmostEqual(res[1], fun(*L)[1],10,'SXfunction evaluation wrong')
    
  def test_SXFunctionc(self):
    self.message("SXFunction constructors")
    x0=SXElement.sym("x")
    x1=SXElement.sym("x")
    x2=SXElement.sym("x")
    x3=SXElement.sym("x")
    x4=SXElement.sym("x")
    x5=SXElement.sym("x")
    x6=SXElement.sym("x")
    y=SX.sym("y",2,3)
    
    f=SXFunction([y],[y])
    self.checkarray(f.input(0).shape,(2,3),"SXFunction constructors")
    self.checkarray(f.output(0).shape,(2,3),"SXFunction constructors")
    
    self.assertRaises(NotImplementedError,lambda: SXFunction(y,[y,y]))
    self.assertRaises(NotImplementedError,lambda: SXFunction(x0,[x0,x1]))

  def test_evalfail(self):
    self.message("eval fail test")
    x = SX.sym("x",2,2)
    f = SXFunction([x], [x])
    f.init()
    self.assertRaises(NotImplementedError,lambda: f.call(x))

  def test_SXconversion(self):
    self.message("Conversions from and to SX")
    y=SXElement.sym("y")
    x=SX.sym("x",3,3)
    SX(y)
    SX(x)
    c.det(x)
    y=array(x)
    c.det(y)
    
  def test_SXbool(self):
    self.message("bool")
    
    x = SXElement.sym("x")
    y = SXElement.sym("y")
    
    f = SXFunction([vertcat([x,y])],[vertcat([logic_and(x,y),logic_or(x,y),logic_not(x)])])
    f.init()
    
    
    for t1 in [0,1]:
      for t2 in [0,1]:
        T1 = t1!=0
        T2 = t2!=0
        f.setInput([t1,t2])
        f.evaluate()
        self.checkarray(f.getOutput(),DMatrix([T1 and T2,T1 or T2,not T1]),"bool(%d,%d): %s" % (t1,t2,str(f.getOutput())))

  def test_SXineq(self):
    self.message("SX ineq")
    
    x = SXElement.sym("x")
    y = SXElement.sym("y")
    
    f = SXFunction([vertcat([x,y])],[vertcat([x<y,x<=y,x>=y,x==y,x!=y])])
    f.init()
    
    
    for t1 in [-10,0.1,0,1,10]:
      for t2 in [-10,0.1,0,1,10]:
        T1 = t1
        T2 = t2
        f.setInput([t1,t2])
        f.evaluate()
        self.checkarray(f.getOutput(),DMatrix([T1 < T2,T1 <= T2, T1 >= T2, T1 == T2, T1 != T2]),"ineq(%d,%d)" % (t1,t2))

    
    
  def test_SXFunctionc2(self):
    self.message("SXmatrix typemaps constructors")
    with internalAPI():
      simplify(SXElement.sym("x"))                 
    SX(array([[SXElement.sym("x")]])).isEmpty()
    list = [ ("SXElement" ,SXElement.sym("x"),(1,1)),
                ("number",2.3, (1,1)),
                ("SX", SX.sym("x"), (1,1)),
                ("numpy.ndarray1D(SXElement)", array([SXElement.sym("x"),SXElement.sym("y")]), (2,1)),
                ("numpy.ndarray(SXElement)", array([[SXElement.sym("x"),SXElement.sym("y")],[SXElement.sym("w"),SXElement.sym("z")]]), (2,2)),
                ("numpy.ndarray(SX,number)", array([[SXElement.sym("x"),2.3]]), (1,2))
    ];
    for name, arg,shape in list:
      self.message(":" + name)
      i=c.transpose(c.transpose(arg))
      self.assertEqual(i.shape[0],shape[0],"shape mismatch")
      self.assertEqual(i.shape[1],shape[1],"shape mismatch")
      SX(arg).isEmpty()
      
  def test_SXconstr(self):
    self.message("SXmatrix constructors")
    list = [ ("SXElement" ,SXElement.sym("x"),(1,1)),
                ("number",2.3, (1,1)),
                ("list(SXElement)", [SXElement.sym("x"),SXElement.sym("y")], (2,1)),
                ("list(SX,number)", [SXElement.sym("x"),2.3], (2,1) ),
                ("list(list(SX,number))", [[SXElement.sym("x"),2.3],[1,SXElement.sym("y")]], (2,2) ),
                ("tuple(SXElement)", (SXElement.sym("x"),SXElement.sym("y")), (2,1)),
                ("tuple(SX,number)", (SXElement.sym("x"),2.3), (2,1)),
                ("SX", SX.sym("x"), (1,1)),
                ("numpy.ndarray1D(SXElement)", array([SXElement.sym("x"),SXElement.sym("y")]), (2,1)),
                ("numpy.ndarray(SXElement)", array([[SXElement.sym("x"),SXElement.sym("y")],[SXElement.sym("w"),SXElement.sym("z")]]), (2,2)),
                ("numpy.ndarray(SX,number)", array([[SXElement.sym("x"),2.3]]), (1,2))
    ];
    for name, arg,shape in list:
      self.message(":" + name)
      i=SX(arg)
      self.assertEqual(i.shape[0],shape[0],"shape mismatch")
      self.assertEqual(i.shape[1],shape[1],"shape mismatch")
      i.isEmpty()
    
  def test_SXFunctionc3(self):
    self.message("vector(SXmatrix) typemaps constructors")
    y=SXElement.sym("y")
    x=SX.sym("x",3,1)
    vertcat([x,x])
    vertcat([y,y])
    vertcat([x,[]])
    
  def test_eval(self):
    self.message("SXFunction eval")
    x=SX.sym("x",2,2)
    y=SX.sym("y",2,2)
    f  = SXFunction([x,y], [x*y])
    f.init()
    f.call([x,y])
    
  def test_symbolcheck(self):
    self.message("Check if non-symbolic inputs are caught")
    self.assertRaises(RuntimeError, lambda : SXFunction([SXElement(0)],[SXElement.sym("x")]))
      
  def test_sparseconstr(self):
    self.message("Check sparsity constructors")
    self.checkarray(DMatrix(Sparsity.tril(3),1).toArray(),matrix([[1,0,0],[1,1,0],[1,1,1]]),"tril")
    self.checkarray(DMatrix(Sparsity.diag(3),1).toArray(),matrix([[1,0,0],[0,1,0],[0,0,1]]),"diag")
    
  def test_subsassignment(self):
    self.message("Check subscripted assignment")

    import numpy
    numpy.random.seed(42)
    xn = numpy.random.random((3,4))

    x=DMatrix(xn)
    
    y=DMatrix.sparse(7,8)
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
      z[s.row()[k],s.getCol()[k]]=1.0
    self.checkarray(y,z,"nonzero scalar assignment")
    y.nz[kl]=DMatrix(kl)
    
    cnt=0
    for k in kl:
      z[s.row()[k],s.getCol()[k]]=kl[cnt]
      cnt+=1
    self.checkarray(y,z,"nonzero range assignment")
    
  @skip(not CasadiOptions.getSimplificationOnTheFly())
  def test_substitute(self):
    self.message("Basic symbolic algebra: substitute")
    x=SXElement.sym("x")
    y=SXElement.sym("y")
    z = cos(x)*y
    self.assertTrue(dependsOn(z,y))
    self.assertTrue(dependsOn(z,x))
    w = substitute(z,x,0)
    self.assertTrue(w.isSymbolic())
    self.assertTrue(dependsOn(w,y))
    self.assertFalse(dependsOn(w,x))
    self.assertTrue(isEqual(w,y))
    r=w-y
    self.assertFalse(r.isSymbolic())     
    self.assertTrue(r.isZero())
    self.assertEqual(r.getValue(),0)
    self.assertEqual(r.getValue(),0)
    y = SX.sym("y",2)
    y = substitute(y+6,y,0)
    self.assertEqual(int(y[0].getValue()),6)
    self.assertEqual(int(y[1].getValue()),6)
   
  def test_primitivefunctions(self):
    self.message("Primitive functions")
    x=SXElement.sym("x")
    
    nums = [-2,-1.5,-1,-0.5,-0.25,0,0.25,0.5,1,1.5,2]
    
    def test(fun,comment,nums,reference):
      self.message(":"+comment)
      f = SXFunction([x],[fun(x)])
      f.init()
      for n,r in zip(nums,reference):
        f.setInput(n)
        f.evaluate()
        self.assertEqual(f.getOutput()[0],r)
    
    test(casadi.sign,"sign",nums,[-1,-1,-1,-1,-1,0,1,1,1,1,1])
    test(casadi.heaviside,"heaviside",nums,[0,0,0,0,0,0.5,1,1,1,1,1])
    test(casadi.ramp,"ramp",nums,[0,0,0,0,0,0,0.25,0.50,1,1.5,2])
    test(casadi.rectangle,"rectangle",nums,[0,0,0,0.5,1,1,1,0.5,0,0,0])
    test(casadi.triangle,"triangle",nums,[0,0,0,0.5,0.75,1,0.75,0.5,0,0,0])
    
    
  def test_taylor(self):
    self.message("univariate taylor expansion")
    x=SXElement.sym("x")
    
    if CasadiOptions.getSimplificationOnTheFly():
      self.assertTrue(isEqual(taylor(sin(x),x),x))
      
    a_=0.13
    x_=0.15

    a = SXElement.sym("a") 
    
    def test(e,r):
      f = SXFunction([x,a],[e])
      f.init()
      f.setInput(x_,0)
      f.setInput(a_,1)
      f.evaluate()
      self.assertAlmostEqual(f.getOutput()[0],r,10)
    
    test(taylor(sin(x),x,a,0),sin(a_))
    test(taylor(sin(x),x,a,1),sin(a_)+cos(a_)*(x_-a_))
    test(taylor(sin(x),x,a,2),sin(a_)+cos(a_)*(x_-a_)-(sin(a_)*(x_-a_)**2)/2.0)
    test(taylor(sin(x),x,a,3),sin(a_)+cos(a_)*(x_-a_)-(sin(a_)*(x_-a_)**2)/2.0-(cos(a_)*(x_-a_)**3)/6.0)
    
    M=SX(matrix([[a*sin(x),a*cos(x)],[exp(a*x),a*x**2],[cos(x),0]]))
    
    f = SXFunction([x,a],[taylor(M,x)])
    f.init()
    f.setInput(x_,0)
    f.setInput(a_,1)
    f.evaluate()
    self.checkarray(f.getOutput(),matrix([[x_*a_,a_],[1+a_*x_,0],[1,0]]),"taylor on dense matrices")
    
  def test_null(self):
    self.message("SXFunction null")
    x = SX.sym("x")

    f = SXFunction([x],[x**2,[]])
    f.init()

    self.assertTrue(f.output(1).isEmpty())
    
    f = SXFunction([x,[]],[x**2,[]])
    f.init()

    self.assertTrue(f.output(1).isEmpty())
    f.evaluate()
    
    r = f.call([x,[]])
    self.assertTrue(r[1].isEmpty())

    r = f.call([x,[]])
    self.assertTrue(r[1].isEmpty())
    
    r = f.call([x,SX.sparse(0,1)])
    self.assertTrue(r[1].isEmpty())

    r = f.call([x,SX.sparse(1,0)])
    self.assertTrue(r[1].isEmpty())
    
    #self.assertRaises(Exception,lambda : f.call([x,x]))
    #self.assertRaises(Exception,lambda : f.call([[],[]]))
    
  def test_mtaylor(self):
    self.message("multivariate taylor expansions")
    x=SXElement.sym("x")
    y=SXElement.sym("y")
    a=SXElement.sym("a")
    b=SXElement.sym("b")

    a_=0.13
    x_=0.15
    
    b_=0.73
    y_=0.75
    
    def test(e,r):
      f = SXFunction([vertcat([x,y]),vertcat([a,b])],[e])
      f.init()
      f.setInput([x_,y_],0)
      f.setInput([a_,b_],1)
      f.evaluate()
      self.assertAlmostEqual(f.getOutput()[0],r,10)
    
    test(mtaylor(sin(x+y),vertcat([x,y]),vertcat([a,b]),0),sin(a_+b_))
    test(mtaylor(sin(x+y),vertcat([x,y]),vertcat([a,b]),1),sin(a_+b_)+(cos(b_+a_)*(x_-a_)+cos(b_+a_)*(y_-b_)))
    
    def sol(x,y,a,b):
      return sin(b+a)+(cos(b+a)*(x-a)+cos(b+a)*(y-b))-(sin(b+a)*(x-a)**2+2*sin(b+a)*(y-b)*(x-a)+sin(b+a)*(y-b)**2)/2
      
    test(mtaylor(sin(x+y),vertcat([x,y]),vertcat([a,b]),2),sol(x_,y_,a_,b_))
    
    def sol(x,y,a,b):
      return sin(b+a)+(cos(b+a)*(x-a)+cos(b+a)*(y-b))-(sin(b+a)*(x-a)**2+2*sin(b+a)*(y-b)*(x-a)+sin(b+a)*(y-b)**2)/2-(cos(b+a)*(x-a)**3+3*cos(b+a)*(y-b)*(x-a)**2+3*cos(b+a)*(y-b)**2*(x-a)+cos(b+a)*(y-b)**3)/6
    
    test(mtaylor(sin(x+y),vertcat([x,y]),vertcat([a,b]),3),sol(x_,y_,a_,b_))
    
    def sol(x,y,a,b):
      return (-2*sin(b+a)*(x-a)*(y-b)-sin(b+a)*(x-a)**2)/2+cos(b+a)*(y-b)-(cos(b+a)*(x-a)**3)/6+cos(b+a)*(x-a)+sin(b+a)
      
    test(mtaylor(sin(x+y),vertcat([x,y]),vertcat([a,b]),3,[1,2]),sol(x_,y_,a_,b_))
    
    test(mtaylor(sin(x+y),vertcat([x,y]),vertcat([0,0]),4,[1,2]),(-3*x_**2*y_-x_**3)/6+y_+x_)
    
  def test_issue107(self):
    self.message("Regression test for issue 107: +=")
    x=SXElement.sym("x")
    y=SXElement.sym("y")

    z=x
    z+=y
    
    self.assertTrue(x.isSymbolic())
    self.assertFalse(z.isSymbolic())
    
    x=SX.sym("x")
    y=SX.sym("y")

    z=x
    z+=y
    
    self.assertTrue(x.isSymbolic())
    self.assertFalse(z.isSymbolic())
    
  def test_evalchecking(self):
    x = SX.sym("x",1,5)
    
    y = SX.sym("y",1,3)
    z = SX.sym("z",5,1)
    q = SX.sym("z",1,6)
    
    f = SXFunction([x],[x**2])
    f.init()
    
    self.assertRaises(RuntimeError, lambda : f.call([y]))
    self.assertRaises(RuntimeError, lambda : f.call([q]))
    self.assertRaises(RuntimeError, lambda : f.call([z]))
    
  def test_indexinglimits(self):
    self.message("Limits of indexing")
    y = casadi.SX.sym("y", 3) 
    self.assertRaises(RuntimeError,lambda : y[[0, 5]] )
    try:
      y[[0, 5]] = SXElement.sym("a")
      self.assertTrue(False)
    except RuntimeError:
      pass
    y[[0, 2]]
    y[[0, 2]] = SXElement.sym("a")
    
  def test_issue181(self):
    self.message("Regression test #181")
    x = SXElement.sym("x")
    #self.assertRaises(TypeError,lambda : SX([x,None]))  # FIXME: this is leaking memory
    self.assertRaises(NotImplementedError,lambda: SXFunction([[x], [None]], [[2 * x]]))
    
  def test_printLimiting(self):
    self.message("printLimiting")

    x = SXElement.sym("x")
    for i in range(100):
      x = sin(x)*x
      
      
    self.assertTrue(len(str(x)) <  4*SX.getMaxNumCallsInPrint())
    
    SX.setMaxNumCallsInPrint(5)
    self.assertTrue(len(str(x)) <  100)
    
    SX.getMaxNumCallsInPrint()
    
  def test_isEqual(self):
    self.message("equivalent")
    x = SXElement.sym("x")
    a = x*x
    b = x*x
    self.assertTrue(a.isEqual(b,1))
    
  @skip(not CasadiOptions.getSimplificationOnTheFly())
  def test_SXsimplifications(self):
    self.message("simplifications")
    x = SXElement.sym("x")
    
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
    ops.append(lambda x: (SXElement(4)-SXElement(4))+x)
    
    
    y = SXElement.sym("x")
    
    ops.append(lambda x: ((x+y)-(y+x))+x)
    ops.append(lambda x: ((x*y)-(y*x))+x)
    ops.append(lambda x: ((-x)-(-x))+x)

    for op in ops:
      y = op(x)
      f = SXFunction([x],[y])
      f.init()
      f.setInput(0.3)
      f.evaluate()
      self.checkarray(f.getOutput(),array(op(0.3)),"simplifications")
      self.assertEqual(str(y),"x")
      
      y = op(-x)
      f = SXFunction([x],[y])
      f.init()
      f.setInput(0.3)
      f.evaluate()
      self.checkarray(f.getOutput(),array(op(-0.3)),"simplifications")
      self.assertEqual(str(y),"(-x)")

  def test_evalf(self):
    x = SXElement(3)
    y = SXElement(5)
    z = evalf(x+y)
    self.assertEqual(type(z),DMatrix)
    self.assertEqual(z,8)

  def test_evalfs(self):
    x = SXElement.sym("x")
    y = SXElement(5)
    z = evalf(x+y,x,3)
    self.assertEqual(type(z),DMatrix)
    self.assertEqual(z,8)
  
  def test_truth(self):
    self.message("Truth values")
    self.assertRaises(Exception, lambda : bool(SXElement.sym("x")))
    self.assertRaises(Exception, lambda : bool(SXElement.sym("x")>0))
    self.assertTrue(bool(SXElement(1)))
    self.assertFalse(bool(SXElement(0)))
    self.assertTrue(bool(SXElement(0.2)))
    self.assertTrue(bool(SXElement(-0.2)))
    self.assertRaises(Exception, lambda : bool(SX.sym("x")))
    self.assertRaises(Exception, lambda : bool(SX.sym("x")>0))
    self.assertTrue(bool(SX(SXElement(1))))
    self.assertFalse(bool(SX(SXElement(0))))
    self.assertTrue(bool(SX(SXElement(0.2))))
    self.assertTrue(bool(SX(SXElement(-0.2))))
    self.assertRaises(Exception, lambda : bool(SX([2.0,3])))
    
  def test_if_else(self):
    x = SXElement.sym("x")
    y = if_else(x,1,2)
    f = SXFunction([x],[y])
    f.init()
    f.setInput(1)
    f.evaluate()
    self.assertTrue(f.getOutput()==1,"if_else")
    f.setInput(0)
    f.evaluate()
    self.assertTrue(f.getOutput()==2,"if_else")
    
    x0 = 2.1
    y = if_else(x>1,x**2,x**3)
    f = SXFunction([x],[y])
    f.init()
    f.setInput(x0)
    f.evaluate()
    self.checkarray(f.getOutput(),x0**2,"if_else sens")
    
    x0 = -2.1
    f.setInput(x0)
    f.evaluate()
    self.checkarray(f.getOutput(),x0**3,"if_else sens")

    
  def test_issue548(self):
    x = SX.sym('x',100)
    f = SXFunction([x],[sum(x)**2])
    f.init()
    h = f.hessian()


  def test_isRegular(self):
    x = SX.sym("x")
    
    self.assertTrue(SXElement(0).isRegular())
    self.assertFalse(SXElement(Inf).isRegular())
    with self.assertRaises(Exception):
      self.assertTrue(x.at(0))
      
    self.assertTrue(SX(DMatrix([0,1])).isRegular())
    self.assertFalse(SX(DMatrix([0,Inf])).isRegular())
    self.assertFalse(vertcat([x,Inf]).isRegular())
    with self.assertRaises(Exception):
      self.assertFalse(vertcat([x,x]).isRegular())
      
      
  def test_getSymbols(self):
    a = SX.sym("a")
    b = SX.sym("b")
    c = SX.sym("c")
    e = cos(a*b) + c
    w = getSymbols(e)
    self.assertEqual(len(w),3)
    if CasadiOptions.getSimplificationOnTheFly():
      self.assertTrue(isEqual(w[0],a))
      self.assertTrue(isEqual(w[1],b))
      self.assertTrue(isEqual(w[2],c))
      
  def test_poly_coeff(self):
    x =SX.sym("x")
    a= SX.sym("a")
    c=SX.sym("c")
    p=poly_coeff(12*x**4+x**2+a*x+c,x)
    self.assertTrue(isEqual(p[0],12))
    self.assertTrue(isEqual(p[1],0))
    self.assertTrue(isEqual(p[2],1))
    self.assertTrue(isEqual(p[3],a))
    self.assertTrue(isEqual(p[4],c))
    
    p=poly_coeff((x-a)*(x+a),x)
    self.assertTrue(isEqual(p[0],1))
    self.assertTrue(isEqual(p[1],0))
    
  def test_poly_roots(self):
  
    p = SX.sym("[a,b]")
    r = poly_roots(p)
    
    f = SXFunction([p],[r])
    f.init()
    f.setInput([2,7])
    a_,b_ = f.input()
    f.evaluate()
    f.output()
    self.checkarray(f.output(),vertcat([-b_/a_]))

    p = SX.sym("[a,b]")
    r = poly_roots(vertcat([p,0]))
    
    f = SXFunction([p],[r])
    f.init()
    f.setInput([2,7])
    a_,b_ = f.input()
    f.evaluate()
    f.output()
    self.checkarray(f.output(),vertcat([-b_/a_,0]))
    
    p = SX.sym("[a,b,c]")
    r = poly_roots(p)
    
    f = SXFunction([p],[r])
    f.init()
    f.setInput([1.13,7,3])
    a_,b_,c_ = f.input()
    d = b_**2-4*a_*c_
    f.evaluate()
    x0 = (-b_-sqrt(d))/2/a_
    x1 = (-b_+sqrt(d))/2/a_
    f.output()
    self.checkarray(f.output(),vertcat([x0,x1]))

    p = SX.sym("[a,b,c,d]")
    r = poly_roots(p)
    
    f = SXFunction([p],[r])
    f.init()
    f.setInput([11,1.3,-1.7,0.1])
    f.evaluate()
    f.output()
    self.checkarray(f.output(),DMatrix([0.298028,-0.479787,0.0635774]),digits=5)
    
    p = SX.sym("[a,b,c,d,e]")
    r = poly_roots(p)
    
    f = SXFunction([p],[r])
    f.init()
    f.setInput([3,6,-123,  -126,1080])
    f.evaluate()
    f.output()
    self.checkarray(f.output(),DMatrix([5,3,-4,-6]),digits=5)
    
  def test_eig_symbolic(self):
    x = SX.sym("x",2,2)
    f = SXFunction([x],[eig_symbolic(x)])
    f.init()
    f.setInput(DMatrix([[2,0.1],[0.3,0.7]]))
    f.evaluate()
    self.checkarray(f.output(),DMatrix([0.67732,2.02268]),digits=5)
    
    
    x = SX.sym("x",2)
    f = SXFunction([x],[eig_symbolic(c.diag(x))])
    f.init()
    f.setInput([3,7])
    f.evaluate()
    self.checkarray(f.output(),f.input())

    
    x = SX.sym("x",5)
    f = SXFunction([x],[eig_symbolic(c.diag(x))])
    f.init()
    f.setInput([3,7,2,1,6])
    f.evaluate()
    self.checkarray(f.output(),f.input())
    
    x = SX.sym("x",2,2)
    y = SX.sym("y",2)
    f = SXFunction([x,y],[eig_symbolic(blkdiag([x,c.diag(y)]))])
    f.init()
    f.setInput(DMatrix([[2,0.1],[0.3,0.7]]),0)
    f.setInput([3,7],1)
    f.evaluate()
    self.checkarray(f.output(),DMatrix([0.67732,2.02268,3,7]),digits=5)

    x = SX.sym("x",3,3)
    x[2,0] = 0
    x[1,0] = 0

    x = sparse(x)

    e = eig_symbolic(x)
    
    f = SXFunction([x],[e])
    f.init()
    f.setInput(range(1,8))
    f.input().printDense()
    f.evaluate()
    self.checkarray(f.output(),DMatrix([1,-0.29150,10.29150]),digits=5)
    
    
    x = SX.sym("x",3,3)
    x[2,0] = 0
    x[1,0] = 0
    x[2,1] = 0
    
    x = sparse(x)

    e = eig_symbolic(x)
    
    f = SXFunction([x],[e])
    f.init()
    f.setInput(range(1,7))
    f.input().printDense()
    f.evaluate()
    self.checkarray(f.output(),DMatrix([1,3,6]),digits=5)

    x = SX.sym("x",Sparsity.triu(5))
  
    f = SXFunction([x],[eig_symbolic(x)])
    f.init()
    f.setInput(6)
    f.input()[Sparsity.diag(5)] = c.diag(range(5))
    f.evaluate()
    self.checkarray(f.output(),DMatrix(range(5)))
    
  def test_jacobian_empty(self):
    x = SX.sym("x",3)

    s = jacobian(DMatrix.sparse(0,0),x).shape
    self.assertEqual(s[0],0)
    self.assertEqual(s[1],3)

    s = jacobian(x,SX.sym("x",0,4)).shape
    self.assertEqual(s[0],3)
    self.assertEqual(s[1],0)
    
  def test_empty_SX(self):
    s = SX([]).shape
    self.assertEqual(s[0],0)
    self.assertEqual(s[1],1)
    x = SX.sym("x")
    x.append(SX([]))
    
  def test_mul_sparsity(self):

    N = 10
    x = SX.sym("x",N,N)
    y = SX.sym("y",N,N)

    x_ = self.randDMatrix(N,N)
    y_ = self.randDMatrix(N,N)

    filt = Sparsity.diag(N)+Sparsity.triplet(N,N,[1],[3])

    f = SXFunction([x,y],[mul(x,y)])
    f.init()
    f.setInput(x_,0)
    f.setInput(y_,1)
    g = SXFunction([x,y],[mul(x,y,filt)])
    g.init()
    g.setInput(x_,0)
    g.setInput(y_,1)
    
    f.evaluate()
    g.evaluate()
    
    self.checkarray(IMatrix(filt,1),IMatrix(g.output().sparsity(),1))
    
    self.checkarray(f.output()[filt],g.output())
    
  @skip(platform_arch==32)
  @memory_heavy()
  def test_large_hessian(self):
    import pickle

    A = pickle.load(file("../data/apoa1-2.pkl",'r'))

    H = DMatrix(A,range(A.size()))
    H = H + H.T
    
    H = H[:20000,:20000]
    
    x = SX.sym("x",H.size1())
    
    f = SXFunction([x],[mul([x.T,H,x])])
    H *= 2
    f.setOption("verbose",True)
    f.init()

    h = f.hessian()
    h.init()
    h.evaluate()
    
    self.assertTrue(h.output().sparsity()==H.sparsity())
    
    self.checkarray(h.output().data(),H.data())

  def test_mxnulloutput(self):
     a = SX.sparse(5,0)
     b = SX.sym("x",2)
     bm = MX.sym("x",2)
     
     f = SXFunction([b],[a])
     f.init()
     c = f.call([bm])[0]

     self.assertEqual(c.size1(),5)
     self.assertEqual(c.size2(),0)
     
     c = f.call([b])[0]

     self.assertEqual(c.size1(),5)
     self.assertEqual(c.size2(),0)
     
     a = SX.sparse(0,0)
     
     f = SXFunction([b],[a])
     f.init()
     
     c = f.call([bm])[0]

     self.assertEqual(c.size1(),0)
     self.assertEqual(c.size2(),0)
     
     c = f.call([b])[0]

     self.assertEqual(c.size1(),0)
     self.assertEqual(c.size2(),0)
     
  def test_mxnull(self):
     a = SX.sparse(5,0)
     b = SX.sparse(0,3)
     
     c = mul(a,b)
     
     self.assertEqual(c.size(),0)
     
     a = SX.sparse(5,3)
     b = SX.sparse(3,4)
     
     c = mul(a,b)
     
     self.assertEqual(c.size(),0)
     
  def  test_mxnullop(self):
    c = SX.sparse(0,0)
    x = SX.sym("x",2,3)
    
    with self.assertRaises(RuntimeError):
      d = x + c

    with self.assertRaises(RuntimeError):
      d = x / c
      
  def test_copysign(self):
    x = SX.sym("x")
    y = SX.sym("y")
    z = copysign(x,y)
    
    f = SXFunction([x,y],[z])
    f.init()
    
    f.setInput(2,0)
    f.setInput(0.5,1)
    f.evaluate()
    self.checkarray(f.output(),DMatrix([2]))

    f.setInput(2,0)
    f.setInput(-0.5,1)
    f.evaluate()
    self.checkarray(f.output(),DMatrix([-2]))
    
    f.setInput(-2,0)
    f.setInput(0.5,1)
    f.evaluate()
    self.checkarray(f.output(),DMatrix([2]))

    f.setInput(-2,0)
    f.setInput(-0.5,1)
    f.evaluate()
    self.checkarray(f.output(),DMatrix([-2]))
    
    f.setInput(2,0)
    f.setInput(0,1)
    f.evaluate()
    self.checkarray(f.output(),DMatrix([2]))
    
    J = f.jacobian()
    J.init()
    
    J.setInput(2,0)
    J.setInput(0.5,1)
    J.evaluate()
    self.checkarray(J.output(),DMatrix([1]))

    J.setInput(2,0)
    J.setInput(-0.5,1)
    J.evaluate()
    self.checkarray(J.output(),DMatrix([-1]))
    
    J.setInput(-2,0)
    J.setInput(0.5,1)
    J.evaluate()
    self.checkarray(J.output(),DMatrix([1]))

    J.setInput(-2,0)
    J.setInput(-0.5,1)
    J.evaluate()
    self.checkarray(J.output(),DMatrix([-1]))
    
    J.setInput(2,0)
    J.setInput(0,1)
    J.evaluate()
    self.checkarray(J.output(),DMatrix([1]))

    J = f.jacobian(1)
    J.init()
    
    J.setInput(2,0)
    J.setInput(0.5,1)
    J.evaluate()
    self.checkarray(J.output(),DMatrix([0]))

    J.setInput(2,0)
    J.setInput(-0.5,1)
    J.evaluate()
    self.checkarray(J.output(),DMatrix([0]))
    
    J.setInput(-2,0)
    J.setInput(0.5,1)
    J.evaluate()
    self.checkarray(J.output(),DMatrix([0]))

    J.setInput(-2,0)
    J.setInput(-0.5,1)
    J.evaluate()
    self.checkarray(J.output(),DMatrix([0]))
    
    J.setInput(2,0)
    J.setInput(0,1)
    J.evaluate()
    self.checkarray(J.output(),DMatrix([0]))
    
  def test_dependsOn(self):
    a = SX.sym("a")
    b = SX.sym("b")
    
    self.assertTrue(dependsOn(a**2,a))
    self.assertTrue(dependsOn(a,a))
    self.assertFalse(dependsOn(0,a))
    self.assertTrue(dependsOn(a**2,vertcat([a,b])))
    self.assertTrue(dependsOn(a,vertcat([a,b])))
    self.assertFalse(dependsOn(0,vertcat([a,b])))
    self.assertTrue(dependsOn(b**2,vertcat([a,b])))
    self.assertTrue(dependsOn(b,vertcat([a,b])))
    self.assertTrue(dependsOn(a**2+b**2,vertcat([a,b])))
    self.assertTrue(dependsOn(a+b,vertcat([a,b])))
    self.assertTrue(dependsOn(vertcat([0,a]),a))
    self.assertTrue(dependsOn(vertcat([a,0]),a))
    self.assertTrue(dependsOn(vertcat([a**2,b**2]),vertcat([a,b])))
    self.assertTrue(dependsOn(vertcat([a,0]),vertcat([a,b])))
    self.assertTrue(dependsOn(vertcat([0,b]),vertcat([a,b])))
    self.assertTrue(dependsOn(vertcat([b,0]),vertcat([a,b])))
    self.assertFalse(dependsOn(vertcat([0,0]),vertcat([a,b])))
    
  @requires("isSmooth")
  def test_isSmooth(self):
    x = SX.sym("a",2,2)
    import warnings
    with warnings.catch_warnings():
      warnings.simplefilter("error",DeprecationWarning)
      with self.assertRaises(Exception):
        isSmooth(x)
      warnings.simplefilter("ignore")
      isSmooth(x)
    
if __name__ == '__main__':
    unittest.main()

